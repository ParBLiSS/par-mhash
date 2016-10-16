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

// bug in mxx left_shift in specialization for std::vector
template <typename T>
std::vector<T> mxx_left_shift(const std::vector<T>& v, const mxx::comm& comm = mxx::comm())
{
    // get datatype
    mxx::datatype dt = mxx::get_datatype<T>();

    // TODO: handle tags with MXX (get unique tag function)
    int tag = 15;
    // receive the size first
    std::vector<T> result;
    size_t right_size = mxx::left_shift(v.size(), comm);

    MPI_Request recv_req;
    // if not last processor
    // TODO: replace with comm.send/ comm.recv which automatically will resolve
    // to BIG MPI if message size is too large
    if (comm.rank() < comm.size()-1 && right_size > 0) {
        result.resize(right_size);
        MPI_Irecv(&result[0], right_size, dt.type(), comm.rank()+1, tag,
                  comm, &recv_req);
    }
    // if not first processor
    if (comm.rank() > 0 && v.size() > 0) {
        // send my most right element to the right
        MPI_Send(const_cast<T*>(&v[0]), v.size(), dt.type(), comm.rank()-1, tag, comm);
    }
    if (right_size > 0 && comm.rank() < comm.size()-1) {
        // wait for the async receive to finish
        MPI_Wait(&recv_req, MPI_STATUS_IGNORE);
    }
    return result;
}


#if (pDNA == 16)
using Alphabet = bliss::common::DNA16;
#elif (pDNA == 5)
using Alphabet = bliss::common::DNA5;
#elif (pDNA == 4)
using Alphabet = bliss::common::DNA;
#endif

#if (pPARSER == FASTQ)
#define FileParser bliss::io::FASTQParser
#elif (pPARSER == FASTA)
#define FileParser bliss::io::FASTAParser
#endif

// constants
#ifndef HASH_KMER_SIZE
#define HASH_KMER_SIZE 31
#endif

static const std::size_t hash_block_size = 4;
static const std::size_t hash_block_count = 2;
static const std::size_t hash_seeds_size = hash_block_count * hash_block_size;

/// Hash Function Seed Values
constexpr static const uint32_t
hash_seed_values[hash_seeds_size] = {
    1000000007,
    1000002043,
    1000003097,
    1000000079,
    1000005103,
    1000005203,
    1000005353,
    1000005403,
};

// Kmer, FilerReader and Iterator data types
using EdgeEncoding = Alphabet;
using KmerType = bliss::common::Kmer<HASH_KMER_SIZE, Alphabet, WordType>;

using FileReaderType = bliss::io::parallel::partitioned_file<
    bliss::io::posix_file, FileParser>;

template <typename Iterator, template <typename> class SeqParser>
using SeqIterType = bliss::io::SequencesIterator<Iterator, SeqParser>;

template <typename Iter>
using NonEOLIter = bliss::iterator::filter_iterator<bliss::utils::file::NotEOL, Iter>;

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
    std::size_t read_id;
    BT hash_values;

    void reset(){
        for(auto i = 0; i < hash_values.size();i++)
            hash_values[i] = std::numeric_limits<HashValueType>::max();
    }

    ReadMinHashBlock() {
        reset();
    }

    void print(std::ostream& ofs){
        ofs << read_id << " ";
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
// MXX_CUSTOM_STRUCT(ReadMinHashBlock, read_id, hash_values);
namespace mxx {
    template <typename BT>
    MXX_CUSTOM_TEMPLATE_STRUCT(WRAP_TEMPLATE(ReadMinHashBlock<BT>), \
                               read_id, hash_values);
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
        hrv.read_id = read.id.get_pos();
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

                HashBlockType hx = hash_functions(*it, seed_index);
                hrv.update_min(hx); // pick the minimum
            }
        } else {
            for (auto it = start; it != end; ++it)  {
                HashBlockType hx = hash_functions(*it, seed_index);
                hrv.update_min(hx);
            }
        }
        *output_iter = hrv;
        return output_iter;
    }
};

template<typename BVT>
void shiftStraddlingRegion(const mxx::comm& comm,
                           std::vector<BVT>&  local_rhpairs,
                           std::size_t& start_offset,
                           std::size_t& end_offset,
                           std::vector<BVT>& straddle_region){
    // Assumes that the local_rhpairs has at least one element
    std::vector<BVT> snd_to_left, right_region;
    // find the starting segment of local_rhpairs that straddles
    //  with the processor on the left
    auto lastv = local_rhpairs.back();
    auto prevx = mxx::right_shift(lastv, comm);

    auto fwx_itr = local_rhpairs.begin();
    if(comm.rank() > 0){
        for(;fwx_itr != local_rhpairs.end(); fwx_itr++){
            if((*fwx_itr).hash_values != prevx.hash_values)
                break;
        }
    }

    if(fwx_itr != local_rhpairs.begin()){
        auto osize = std::distance(local_rhpairs.begin(), fwx_itr);
        snd_to_left.resize(osize);
        std::copy(local_rhpairs.begin(), fwx_itr, snd_to_left.begin());
    }
    auto to_snd_size = snd_to_left.size();
    auto to_rcv_size = mxx::right_shift(to_snd_size, comm);
    // std::cout << snd_to_left.size() << std::endl;
    right_region = mxx_left_shift(snd_to_left, comm);
    start_offset = std::distance(local_rhpairs.begin(), fwx_itr);
    int soffset = (int)  (start_offset ==  local_rhpairs.size());
    auto total =  mxx::allreduce(soffset);
    if(comm.rank() == 0)
      std::cout << "Total : "<< total << std::endl;

    // find the ending segment of local_rhpairs that straddles
    //  with the processor on the right
    //  - there will be at least one value
    auto rvx_itr = local_rhpairs.rbegin();
    if(comm.rank() < comm.size() - 1) {
        for(;rvx_itr != local_rhpairs.rend();rvx_itr++){
            if((*rvx_itr).hash_values != lastv.hash_values)
                break;
        }
    }
    auto left_region_size = std::distance(local_rhpairs.rbegin(), rvx_itr);
    end_offset = local_rhpairs.size() - left_region_size;
    //std::cout << left_region_size;

    // construct straddling region from left and right region
    straddle_region.resize(left_region_size + right_region.size());
    std::copy(local_rhpairs.begin() + end_offset, local_rhpairs.end(),
              straddle_region.begin());
    std::copy(right_region.begin(), right_region.end(),
              straddle_region.begin() + left_region_size);

}

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
            read_pairs.push_back(std::make_pair(outer_itr->read_id,
                                                inner_itr->read_id));
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
    // shift straddling
    shiftStraddlingRegion(comm, local_rhpairs, start_offset, end_offset,
                          straddle_region);
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
    return max_size;
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




void generateSequencePairs(mxx::comm& comm,
                           bliss::io::file_data& file_data,
                           std::vector< std::pair<uint64_t, uint64_t> >& read_pairs) {
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
  if(comm.rank()  == 0)
      std::cout << "Total RH Pairs : " << total_pairs << std::endl;

  BL_BENCH_START(genpr);
  uint64_t max_block_size = 0;
  comm.with_subset(
      local_rhpairs.begin() != local_rhpairs.end(), [&](const mxx::comm& comm){
        max_block_size =
          generateOverlapReadPairs(comm, local_rhpairs, read_pairs);
      });
  BL_BENCH_COLLECTIVE_END(genpr, "pair_gen", max_block_size, comm);
  BL_BENCH_REPORT_MPI_NAMED(genpr, "genpair", comm);
}

void runFSO(mxx::comm& comm,
            std::vector<std::string>& inFiles, std::string outPrefix){

  // Read files
  BL_BENCH_INIT(rfso);
  std::vector<::bliss::io::file_data> file_data;
  size_t total = 0;

  BL_BENCH_START(rfso);
  load_file_data(comm, inFiles, file_data);
  BL_BENCH_COLLECTIVE_END(rfso, "read_files", total, comm);

  BL_BENCH_START(rfso);
  std::vector< std::pair<uint64_t, uint64_t> > read_pairs;
  for(auto i = 0u; i < hash_block_count; i++){
      if(seed_index >= (hash_block_count)){
        std::cout << "ERROR : seed index exceeded!!!" << std::endl;
        exit(1);
      }
      generateSequencePairs(comm, file_data.back(), read_pairs);
      seed_index++;
  }
  BL_BENCH_COLLECTIVE_END(rfso, "all_pairs", read_pairs.size(), comm);

  auto total_pairs = mxx::allreduce(read_pairs.size());
  if(comm.rank()  == 0)
      std::cout << "Total Read Pairs : " << total_pairs << std::endl;
  BL_BENCH_REPORT_MPI_NAMED(rfso, "rfso_app", comm);
}

void parse_args(int argc, char **argv,
                mxx::comm& comm,
                std::vector<std::string>& filenames,
                std::string& outPrefix){
  try { // try-catch block for commandline

    TCLAP::CmdLine cmd("Parallel de bruijn graph compaction", ' ', "0.1");

    // MPI friendly commandline output.
    bliss::utils::tclap::MPIOutput cmd_output(comm);
    cmd.setOutput(&cmd_output);

    // output file argument
    TCLAP::ValueArg<std::string> outputArg("O", "output_prefix",
                                           "Prefix for output files, including directory",
                                           false, "", "string", cmd);

    // input files
    TCLAP::UnlabeledMultiArg<std::string> fileArg("filenames", "FASTA or FASTQ file names",
                                                  true, "string", cmd);


    // Parse the argv array.
    cmd.parse( argc, argv );

    outPrefix = outputArg.getValue();
    filenames = fileArg.getValue();

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


  std::vector<std::string> filenames;
  std::string outPrefix;
  outPrefix.assign("./output");

  // parse arguments
  parse_args(argc, argv, comm, filenames, outPrefix);
  if(comm.rank() == 0 && filenames.size() > 0){
    for(auto fx : filenames) std::cout << fx << std::endl;
  }

  // runFSO
  comm.barrier();
  // auto start = std::chrono::steady_clock::now();

  // if(!comm.rank())
  //    std::cout << "Beginning computation, timer started" << std::endl;

  runFSO(comm, filenames, outPrefix);


  comm.barrier();
  // auto end = std::chrono::steady_clock::now();
  // auto elapsed_time  = std::chrono::duration<double, std::milli>(end - start).count();

  // if(!comm.rank())
  //    std::cout << "Time (ms) -> " << elapsed_time << std::endl;

 // TODO: compute elapsed time
}
