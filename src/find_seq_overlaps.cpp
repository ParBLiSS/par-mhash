#include "bliss-config.hpp"

#include <string>
#include <iostream>
#include <fstream>

#include "utils/logging.h"
#include "common/alphabets.hpp"
#include "common/kmer.hpp"
#include "common/base_types.hpp"
#include "utils/kmer_utils.hpp"

#include "io/file.hpp"
#include "io/kmer_parser.hpp"
#include "io/kmer_file_helper.hpp"

#include "mxx/env.hpp"
#include "mxx/comm.hpp"
#include "mxx/utils.hpp"
#include "tclap/CmdLine.h"
#include "utils/tclap_utils.hpp"

#include "io_utils.hpp"
#include "compare_overlaps.hpp"

/// Hash Function Seed Values
#include "hash_seeds.hpp"

const static std::size_t hash_block_size = 4;
static std::size_t hash_block_count = 250;
static std::size_t hash_seeds_size;

#if (pDNA == 4)
using Alphabet = bliss::common::DNA;
#elif (pDNA == 5)
using Alphabet = bliss::common::DNA5;
#elif (pDNA == 16)
using Alphabet = bliss::common::DNA16;
#endif

// constants
#ifndef HASH_KMER_SIZE
#define HASH_KMER_SIZE 31
#endif


#if (pPARSER == FASTQ)
#define FileParser bliss::io::FASTQParser
#elif (pPARSER == FASTA)
#define FileParser bliss::io::FASTAParser
#endif


using FileReaderType = bliss::io::parallel::partitioned_file<
  bliss::io::posix_file, FileParser>;

// Kmer, FilerReader and Iterator data types
using EdgeEncoding = Alphabet;
using KmerType = bliss::common::Kmer<HASH_KMER_SIZE, Alphabet, WordType>;

template <typename Iterator, template <typename> class SeqParser>
using SeqIterType = bliss::io::SequencesIterator<Iterator, SeqParser>;

template <typename Iter>
using NonEOLIter = bliss::iterator::filter_iterator<bliss::utils::file::NotEOL, Iter>;

uint64_t load_file_data(mxx::comm& comm,
                        std::vector<std::string>& inFiles,
                        std::vector<bliss::io::file_data>& file_data){
    uint64_t total = 0;
    for (auto fn : inFiles) {
        if (comm.rank() == 0) printf("READING %s via posix\n", fn.c_str());

        FileReaderType fobj(fn, KmerType::size + 1, comm);

        file_data.push_back(fobj.read_file());
        total += file_data.back().getRange().size();
    }
    return total;
}

// Hash Value types
using HashValueType = uint64_t;
using HashBlockType = std::array<HashValueType, hash_block_size>;

//
// @brief Kmer specialization for MurmurHash.  generated hash is 128 bit.
//
template <typename KMER, bool Prefix = false>
class murmur {

protected:
    static constexpr unsigned int nBytes = (KMER::nBits + 7) / 8;

public:
    inline uint64_t operator()(const KMER & kmer, uint32_t seed_value = 42) const {
            // produces 128 bit hash.
            uint64_t h[2];
            // let compiler optimize out all except one of these.
            if (sizeof(void*) == 8)
                MurmurHash3_x64_128(kmer.getData(), nBytes, seed_value, h);
            else if (sizeof(void*) == 4)
                MurmurHash3_x86_128(kmer.getData(), nBytes, seed_value, h);
            else
                throw ::std::logic_error("ERROR: neither 32 bit nor 64 bit system");

            // use the upper 64 bits.
            if (Prefix)
                return h[1];
            else
                return h[0];
        }
};


//
// KMER : Kmer Type
// HBT : a type to store the resulting hash values, need to
//      support the following operations
//       - empty constructor
//       - copy constructor
//       - indexing operations [idx], [idx] =
//
template<typename KMER, typename HBT>
struct MinHashFunctionBlock {

    using value_type = HBT;

    HBT operator()(const KMER& tx, std::size_t sidx){
        murmur<KMER> mmur_obj;

        HBT hbx;
        // std::cout << tx.getData()[0] << " ";
        auto sdx = sidx * hash_block_size;
        for(auto i = 0; (i < hbx.size()) && (sdx < hash_seeds_size);
            i++, sdx++){
            // std::cout << hrx[i] << " ";
            hbx[i] = mmur_obj(tx, hash_seed_values[sdx]);
        }
        // std::cout << std::endl;
        return hbx;
    }
};

template<typename BT>
struct ReadMinHashBlock{
    std::size_t seq_id;
    BT hash_values;

    void reset(){
        for(auto i = 0; i < hash_values.size();i++)
            hash_values[i] = std::numeric_limits<HashValueType>::max();
    }

    ReadMinHashBlock() {
        reset();
    }

    void print(std::ostream& ofs){
        ofs << seq_id << " ";
        for(auto i = 0; i < hash_values.size();i++)
            ofs << hash_values[i] << " ";
    }

    void update_min(const BT& cx){
        for(auto i = 0;i < hash_values.size();i++) {
            if(cx[i] < hash_values[i]) hash_values[i] = cx[i];
        }
    }
};

#define WRAP_TEMPLATE(...) __VA_ARGS__

// defining your own type for structs which are non-templated
// MXX_CUSTOM_STRUCT(ReadMinHashBlock, seq_id, hash_values);
namespace mxx {
    template <typename BT>
    MXX_CUSTOM_TEMPLATE_STRUCT(WRAP_TEMPLATE(ReadMinHashBlock<BT>), \
                               seq_id, hash_values);
}

static std::size_t seed_index = 0;
//
// @tparam KmerType output value type of this parser. not necessarily the same
//                  as the map's final storage type.
//
template <typename MinHashFunctionBlockType, typename KmerType>
struct SeqMinHashGenerator {

    // type of element generated by this parser.
    using value_type = ReadMinHashBlock< typename MinHashFunctionBlockType::value_type >;
    using kmer_type = KmerType;
    static constexpr size_t window_size = kmer_type::size;

    bliss::partition::range<std::size_t> valid_range;

    SeqMinHashGenerator(::bliss::partition::range<size_t> const & _valid_range)
        : valid_range(_valid_range) {};


    //
    // @brief generate kmers from 1 sequence.  result inserted into output_iter,
    //        which may be preallocated.
    // @param read          sequence object, which has pointers to the raw byte array.
    // @param output_iter   output iterator pointing to insertion point for underlying container.
    // @return new position for output_iter
    // @tparam SeqType      type of sequence.  inferred.
    // @tparam OutputIt     output iterator type, inferred.
    template <typename SeqType, typename OutputIt>
    OutputIt operator()(SeqType & read, OutputIt output_iter) {
        MinHashFunctionBlockType hash_functions;

        static_assert(std::is_same< SeqType,
                      bliss::io::FASTQSequence<typename SeqType::IteratorType> >::value,
                      "ReadLengthParser only supports FASTQ at the moment.");

        static_assert(std::is_same<value_type,
                      typename std::iterator_traits<OutputIt>::value_type>::value,
                      "output type and output container value type are not the same");

        using Alphabet = typename KmerType::KmerAlphabet;

        /// converter from ascii to alphabet values
        using BaseCharIterator =
            bliss::iterator::transform_iterator< NonEOLIter<typename SeqType::IteratorType>,
                                                 bliss::common::ASCII2<Alphabet> >;

        /// kmer generation iterator
        using KmerIterType =
            bliss::common::KmerGenerationIterator<BaseCharIterator,
                                                  KmerType>;

        static_assert(std::is_same<typename std::iterator_traits<KmerIterType>::value_type,
                      KmerType>::value,
                      "input iterator and output iterator's value types differ");


        //== filtering iterator
        bliss::utils::file::NotEOL neol;
        NonEOLIter<typename SeqType::IteratorType> eolstart(neol, read.seq_begin,
                                                            read.seq_end);
        NonEOLIter<typename SeqType::IteratorType> eolend(neol, read.seq_end);

        //== set up the kmer generating iterators
        KmerIterType start(BaseCharIterator(eolstart, bliss::common::ASCII2<Alphabet>()),
                           true);
        KmerIterType end(BaseCharIterator(eolend, bliss::common::ASCII2<Alphabet>()),
                         false);

        // printf("First: pos %lu kmer %s\n", read.id.id,
        //       bliss::utils::KmerUtils::toASCIIString(*start).c_str());

        bliss::partition::range<size_t> seq_range(
            read.seq_global_offset(), read.seq_global_offset() + read.seq_size());
        value_type hrv;
        hrv.seq_id = read.id.get_pos();
        if (seq_range.contains(valid_range.end)) {
            // seq_range contains overlap.
            // not checking by end iterator at valid_range.end, since the NonEOLIter
            // is a filter iterator that may skip over that pos.
            int64_t valid_dist = valid_range.end - seq_range.start;
            for (auto it = start; it != end; ++it)  {
                // check tail of window -> transform iterator, get base -> non EOL
                // iterator, get base -> seq raw char iter.
                if (std::distance(read.seq_begin,
                                  it.getTrailingIterator().getBaseIterator().getBaseIterator())
                    >= valid_dist) {
                    break;
                }
                KmerType ckx = *it;
                auto rcx = ckx.reverse_complement();
                HashBlockType hx = hash_functions((*it < rcx) ? *it : rcx, seed_index);
                hrv.update_min(hx); // pick the minimum
            }
        } else {
            for (auto it = start; it != end; ++it)  {
                KmerType ckx = *it;
                auto rcx = ckx.reverse_complement();
                HashBlockType hx = hash_functions((ckx < rcx) ? ckx : rcx, seed_index);
                hrv.update_min(hx);
            }
        }
        *output_iter = hrv;
        return output_iter;
    }
};

template<typename Iterator, typename ReadIdType=uint64_t>
uint64_t generatePairs(const mxx::comm& comm,
                       Iterator start_itr,
                       Iterator end_itr,
                       std::vector<std::pair<ReadIdType, ReadIdType>>& read_pairs){
    uint64_t nsize = std::distance(start_itr, end_itr);
    // return nsize;
    for(auto outer_itr = start_itr; outer_itr != end_itr; outer_itr++){
        for(auto inner_itr = outer_itr + 1; inner_itr != end_itr; inner_itr++){
            // generate pair
            read_pairs.push_back(std::make_pair(outer_itr->seq_id,
                                                inner_itr->seq_id));
        }
    }
    return nsize;
}


template<typename BVT, typename ReadIdType=uint64_t>
uint64_t generateOverlapReadPairs(const mxx::comm& comm,
                              std::vector<BVT>&  local_rhpairs,
                              std::vector< std::pair<ReadIdType, ReadIdType> >& read_pairs){
    std::size_t start_offset, end_offset;
    std::vector<BVT> straddle_region;
    //if((*rvx_itr).hash_values != lastv.hash_values)

    // shift straddling
    shiftStraddlingRegion(comm, local_rhpairs, start_offset, end_offset,
                          straddle_region,
                          [&](const BVT& x, const BVT& y){
                              return (x.hash_values == y.hash_values);
                          });
    if(local_rhpairs.size() > 0 && start_offset >= local_rhpairs.size()){
        std::cout << "ERROR : start offset beyond size!!!" << std::endl;
        exit(1);
    }
    if(local_rhpairs.size() > 0 && end_offset > local_rhpairs.size()){
        std::cout << "ERROR : end offset beyond size!!!" << std::endl;
        exit(1);
    }
    if(start_offset > end_offset){
        std::cout << "ERROR : start offset > end offset!!!" << std::endl;
        exit(1);
    }
    // generate pairs
    auto csize = generatePairs(comm,
                               straddle_region.begin(),
                               straddle_region.end(),
                               read_pairs);
    auto max_size = csize;
    auto rbv_itr = local_rhpairs.begin() + start_offset;
    auto prev_rbv = rbv_itr;

    for(;(rbv_itr != local_rhpairs.begin() + end_offset) &&
         (rbv_itr != local_rhpairs.end());
        rbv_itr++){
        if((*rbv_itr).hash_values == (*prev_rbv).hash_values)
            continue;
        csize = generatePairs(comm, prev_rbv, rbv_itr, read_pairs);
        if(csize > max_size) max_size = csize;
        prev_rbv = rbv_itr;
    }
    if(local_rhpairs.end() != rbv_itr && 
       rbv_itr == local_rhpairs.begin() + end_offset && prev_rbv != rbv_itr)
        csize = generatePairs(comm, prev_rbv, rbv_itr, read_pairs);

    if(csize > max_size) max_size = csize;
    // return max_size;
    // remove duplicates
    comm.with_subset(read_pairs.begin() != read_pairs.end(), [&](const mxx::comm& comm){
        // sort read pairs
        mxx::sort(read_pairs.begin(), read_pairs.end(),
                  [&](const std::pair<ReadIdType, ReadIdType> x,
                      const std::pair<ReadIdType, ReadIdType> y){
                    return (x < y);
                  }, comm);

        auto prevPair = mxx::right_shift(read_pairs.back(), comm);
        if(comm.rank() == 0)
          prevPair = std::make_pair((ReadIdType) std::numeric_limits<ReadIdType>::max(),
                                    (ReadIdType) std::numeric_limits<ReadIdType>::max());
        auto cur_itr = read_pairs.begin();
        auto upd_itr = cur_itr;
        while(cur_itr != read_pairs.end()){
          if(*cur_itr != prevPair){
            *upd_itr = *cur_itr;
            upd_itr++;
          }
          prevPair = *cur_itr;
          cur_itr++;
        }
        auto psize = std::distance(read_pairs.begin(), upd_itr);
        read_pairs.resize(psize);
      });
    return read_pairs.size();
}


template<typename T>
void generateSequencePairs(mxx::comm& comm,
                           bliss::io::file_data& file_data,
                           std::vector< std::pair<T, T> >& read_pairs) {
  BL_BENCH_INIT(genpr);
  BL_BENCH_START(genpr);

  using SeqMinHashGeneratorType = SeqMinHashGenerator<
      MinHashFunctionBlock<KmerType, HashBlockType>, KmerType>;
  using RHBType = SeqMinHashGeneratorType::value_type;

  std::vector<RHBType>  local_rhpairs;

  bliss::io::KmerFileHelper::template
    parse_file_data<SeqMinHashGeneratorType, FileParser,
                    SeqIterType>(file_data, local_rhpairs, comm);

  BL_BENCH_COLLECTIVE_END(genpr, "compute_hash", local_rhpairs.size(), comm);

  BL_BENCH_START(genpr);
  auto cmp_fn = [&](const RHBType& x, const RHBType& y){
    return (x.seq_id < y.seq_id);
  };
  if(!mxx::is_sorted(local_rhpairs.begin(), local_rhpairs.end(),
                     cmp_fn, comm))
    mxx::sort(local_rhpairs.begin(), local_rhpairs.end(), cmp_fn, comm);

  auto seq_idx = mxx::scan(local_rhpairs.size(), comm);
  if(local_rhpairs.size() > 0)
    seq_idx -= local_rhpairs.size();
  for(auto fitr = local_rhpairs.begin(); fitr != local_rhpairs.end(); fitr++, seq_idx++)
    fitr->seq_id = seq_idx;

  BL_BENCH_COLLECTIVE_END(genpr, "assign_ids", local_rhpairs.size(), comm);

  BL_BENCH_START(genpr);

  mxx::sort(local_rhpairs.begin(), local_rhpairs.end(),
            [&] (const RHBType& x, const RHBType& y){
                for(auto i = 0; i < x.hash_values.size();i++){
                    if(x.hash_values[i] == y.hash_values[i]) continue;

                    if(x.hash_values[i] < y.hash_values[i]) return true;
                    else return false;
                }
                return false;
            }, comm);
  BL_BENCH_COLLECTIVE_END(genpr, "sort_records", local_rhpairs.size(), comm);
  auto total_pairs = mxx::allreduce(local_rhpairs.size());
  //if(comm.rank()  == 0)
  //    std::cout << "TOTAL RH Pairs : " << total_pairs << std::endl;

  BL_BENCH_START(genpr);
  uint64_t max_block_size = 0;
  comm.with_subset(
      local_rhpairs.begin() != local_rhpairs.end(), [&](const mxx::comm& comm){
        max_block_size =
          generateOverlapReadPairs(comm, local_rhpairs, read_pairs);
      });
  BL_BENCH_COLLECTIVE_END(genpr, "pair_gen", max_block_size, comm);
  //BL_BENCH_REPORT_MPI_NAMED(genpr, "genpair", comm);
}

template<typename T=uint64_t>
void runFSO(mxx::comm& comm,
            std::string positionFile,
            std::vector<std::string>& inFiles,
            std::string outPrefix,
            uint32_t threshold){

  // Read files
  BL_BENCH_INIT(rfso);
  std::vector<::bliss::io::file_data> file_data;
  size_t total = 0;

  BL_BENCH_START(rfso);
  load_file_data(comm, inFiles, file_data);
  BL_BENCH_COLLECTIVE_END(rfso, "read_files", total, comm);

  BL_BENCH_START(rfso);
  std::vector< std::pair<T, T> > read_pairs;
  for(auto i = 0u; i < hash_block_count; i++){
      if(seed_index >= (hash_block_count)){
        std::cout << "ERROR : seed index exceeded!!!" << std::endl;
        exit(1);
      }
      generateSequencePairs(comm, file_data.back(), read_pairs);
      seed_index++;
  }
  BL_BENCH_COLLECTIVE_END(rfso, "all_pairs", read_pairs.size(), comm);

  auto totalHashPairs = mxx::allreduce(read_pairs.size());
  if(comm.rank()  == 0)
      std::cout << "TOTAL Hash Pairs : " << totalHashPairs << std::endl;


  BL_BENCH_START(rfso);
  compareOverLaps(comm, positionFile, read_pairs, threshold);
  BL_BENCH_COLLECTIVE_END(rfso, "compare_overlaps", read_pairs.size(), comm);

  BL_BENCH_REPORT_MPI_NAMED(rfso, "rfso_app", comm);
}

void parse_args(int argc, char **argv,
                mxx::comm& comm,
                std::string& positionFile,
                std::vector<std::string>& filenames,
                std::string& outPrefix,
                uint32_t& threshold){
  try { // try-catch block for commandline

    TCLAP::CmdLine cmd("Overlap Graph Construction", ' ', "0.1");

    // MPI friendly commandline output.
    bliss::utils::tclap::MPIOutput cmd_output(comm);
    cmd.setOutput(&cmd_output);

   // position file argument
    TCLAP::ValueArg<std::string> positionArg("p", "position_file",
                                             "Position for input file (full path)",
                                             true, "", "string", cmd);



    // output file argument
    TCLAP::ValueArg<std::string> outputArg("O", "output_prefix",
                                           "Prefix for output files, including directory",
                                           false, "", "string", cmd);

    // threshold argument
    TCLAP::ValueArg<uint32_t> threshArg("d", "overlap_threshold",
                                          "Threshold for Overlap",
                                          false, 20, "int", cmd);

    // block count argument
    TCLAP::ValueArg<std::size_t> blockCountArg("B", "block_count",
                                               "Number of blocks",
                                               false, 250, "int", cmd);

    // block size argument
    //TCLAP::ValueArg<std::size_t> blockSizeArg("T", "block_size",
    //                                          "Block Size",
    //                                          false, 4, "int", cmd);


    // input files
    TCLAP::UnlabeledMultiArg<std::string> fileArg("filenames", "FASTA or FASTQ file names",
                                                  true, "string", cmd);


    // Parse the argv array.
    cmd.parse( argc, argv );

    positionFile = positionArg.getValue();
    outPrefix = outputArg.getValue();
    filenames = fileArg.getValue();
    threshold = threshArg.getValue();
    hash_block_count = blockCountArg.getValue();
    // hash_block_size = blockSizeArg.getValue();

    if(comm.rank() == 0){
      std::cout << "--------------------------------------" << std::endl;
      std::cout << "Position File   : " << positionFile << std::endl;
      std::cout << "Output File Pfx : " << outPrefix << std::endl;
      std::cout << "Input File      : " << filenames.front() << std::endl;
      std::cout << "Threshold       : " << threshold << std::endl;
      std::cout << "Block Size      : " << hash_block_size << std::endl;
      std::cout << "Block Count     : " << hash_block_count << std::endl;
      std::cout << "--------------------------------------" << std::endl;
    }

  } catch (TCLAP::ArgException &e)  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    exit(-1);
  }
}

int main(int argc, char** argv) {

  LOG_INIT(); // init logging

  mxx::env e(argc, argv); // MPI init
  mxx::comm comm;

  if (comm.rank() == 0)
    std::cout << "EXECUTING " << std::string(argv[0]) << std::endl;


  std::string positionFile;
  std::vector<std::string> filenames;
  std::string outPrefix;
  uint32_t threshold;
  outPrefix.assign("./output");

  // parse arguments
  parse_args(argc, argv, comm,
             positionFile, filenames, outPrefix, threshold);
  if(comm.rank() == 0 && filenames.size() > 0){
    for(auto fx : filenames) std::cout << fx << std::endl;
  }

  auto availableSeeeds = sizeof(hash_seed_values);
  availableSeeeds /= sizeof(hash_seed_values[0]);
  if(availableSeeeds < (hash_block_size * hash_block_count)){
      if(comm.rank() == 0)
          std::cout << "Not Enough seeds : " << std::endl;
      return 1;
  }

  // runFSO
  comm.barrier();
  auto start = std::chrono::steady_clock::now();

  if(!comm.rank())
      std::cout << "Timer started" << std::endl;

  runFSO(comm, positionFile, filenames, outPrefix, threshold);
  //std::vector< std::pair<uint64_t, uint64_t> > read_pairs;
  // compareOverLaps(comm, positionFile, read_pairs, threshold);


  comm.barrier();
  auto end = std::chrono::steady_clock::now();
  auto elapsed_time  = std::chrono::duration<double, std::milli>(end - start).count();

  if(!comm.rank())
      std::cout << "Time (ms) -> " << elapsed_time << std::endl;

 // TODO: compute elapsed time
}
