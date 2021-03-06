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

/// Run configuration
#include "run_cfg.hpp"
/// Hash Function Seed Values
#include "hash_seeds.hpp"

static struct RunArgs{

    RunArgs():min_kmer_length(8),
              max_kmer_length(21),
              min_hash_block_size(2),
              max_hash_block_size(5),
              total_blocks(0),
              processed_blocks(0),
              read_threshold(10),
              hash_block_size(2),
              hash_block_count(250),
              kmer_length(16),
              read_length(100),
              maxBucketSize(-1){}

    const int min_kmer_length;
    const int max_kmer_length;
    const int min_hash_block_size;
    const int max_hash_block_size;
    uint64_t total_blocks;
    uint64_t processed_blocks;
    unsigned read_threshold;
    int hash_block_size;
    int hash_block_count;
    int kmer_length;
    unsigned read_length;
    int maxBucketSize;

    int hash_seeds_size;
    std::string positionFile;
    std::vector<std::string> filenames;
    std::string candFile;
    std::string trueFile;
    std::string outPrefix;
    uint32_t minOverlap;
    std::string runType;

} run_params;

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
        auto sdx = sidx * run_params.hash_block_size;
        for(auto i = 0; (i < hbx.size()) && (sdx < run_params.hash_seeds_size);
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

    using hash_value_type = BT;

    std::size_t seq_id;
    BT hash_values;

    void reset(){
        seq_id = std::numeric_limits<std::size_t>::max();
        for(auto i = 0; i < hash_values.size();i++)
            hash_values[i] = std::numeric_limits<HashValueType>::max();
    }

    ReadMinHashBlock() {
        reset();
    }

    bool is_good(){
        return (hash_values[0] < std::numeric_limits<HashValueType>::max());
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
    using block_value_type = typename MinHashFunctionBlockType::value_type ;
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
        if(read.seq_size() < run_params.read_threshold){
           *output_iter = hrv;
           return output_iter;
        }
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
                block_value_type hx = hash_functions((ckx < rcx) ? ckx : rcx, seed_index);
                hrv.update_min(hx); // pick the minimum
            }
        } else {
            for (auto it = start; it != end; ++it)  {
                KmerType ckx = *it;
                auto rcx = ckx.reverse_complement();
                block_value_type hx = hash_functions((ckx < rcx) ? ckx : rcx, seed_index);
                hrv.update_min(hx);
            }
        }
        *output_iter = hrv;
        return output_iter;
    }
};
 
static double pgen_time = 0.0;
template<typename Iterator, typename ReadIdType=uint64_t>
uint64_t generatePairs(const mxx::comm& comm,
                       Iterator start_itr,
                       Iterator end_itr,
                       std::vector<std::pair<ReadIdType, ReadIdType>>& read_pairs){
    uint64_t nsize = std::distance(start_itr, end_itr);
    //auto start = std::chrono::steady_clock::now();
    // limit the maximum size of bucket
    if(run_params.maxBucketSize > 0 && nsize > run_params.maxBucketSize)
        return 0;
    for(auto outer_itr = start_itr; outer_itr != end_itr; outer_itr++){
        for(auto inner_itr = outer_itr + 1; inner_itr != end_itr; inner_itr++){
            if(outer_itr->seq_id == inner_itr->seq_id) continue;
            //nsize += 1; continue;
            // generate pair
            if(outer_itr->seq_id < inner_itr->seq_id)
                read_pairs.push_back(std::make_pair(outer_itr->seq_id,
                                                    inner_itr->seq_id));
            else
                read_pairs.push_back(std::make_pair(inner_itr->seq_id,
                                                    outer_itr->seq_id));
        }
    }
    //pgen_time += 
    return nsize;
}

template<typename T=uint64_t>
void eliminateLocalDupes(std::vector< std::pair<T, T> >& readIdPairs){
      if(readIdPairs.size () == 0) return;
      std::sort(readIdPairs.begin(), readIdPairs.end());
      auto cItr = readIdPairs.begin();
      auto vItr = cItr;
      std::pair<T,T> prevValue;
      if(cItr != readIdPairs.end()){
        prevValue = *cItr;
        cItr++;
      }
      for(;cItr != readIdPairs.end();cItr++){
        if(!(*cItr == prevValue)){
          *vItr = *cItr; vItr++;
          prevValue = *cItr;
        }
      }
      auto nPairs = std::distance(readIdPairs.begin(), vItr);
      readIdPairs.resize(nPairs);
      std::vector< std::pair<T,T> >(readIdPairs).swap(readIdPairs);
}

template<typename BVT, typename ReadIdType=uint64_t>
uint64_t generateOverlapReadPairs(const mxx::comm& comm,
                              std::vector<BVT>&  local_rhpairs,
                              std::vector< std::pair<ReadIdType, ReadIdType> >& read_pairs){
    static std::size_t size_tnow = 0;
    static std::size_t px_idx = 0;
    std::size_t start_offset, end_offset;
    std::vector<BVT> straddle_region;

    // shift straddling
    shiftStraddlingRegion(comm, local_rhpairs, start_offset, end_offset,
                          straddle_region,
                          [&](const BVT& x, const BVT& y){
                              return (x.hash_values == y.hash_values);
                          });
    bool lFlag = true;
    std::size_t isIn = 1;
    if(local_rhpairs.size() > 0 && start_offset >= local_rhpairs.size()){
        //std::cout << comm.rank() 
        //          << " ERROR : start offset beyond size!!!" << std::endl;
        isIn = 0; 
    }
    if(local_rhpairs.size() > 0 && end_offset > local_rhpairs.size()){
        //std::cout << comm.rank() 
        //          << "ERROR : end offset beyond size!!!" << std::endl;
        isIn = 0; 
    }
    if(start_offset > end_offset){
        //std::cout << comm.rank() 
        //          << "ERROR : start offset > end offset!!!" << std::endl;
        isIn = 0; 
    }
    std::size_t max_size, sub_size;
    sub_size = mxx::allreduce(isIn, comm);
    comm.with_subset((isIn == 1), [&](const mxx::comm& comm){
      auto csize = generatePairs(comm,
                                 straddle_region.begin(),
                                 straddle_region.end(),
                                 read_pairs);
      max_size = csize;
      // TODO: redistribute the rest of the region equally 
      // among all the processors
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

    // return max_size;
    // auto rmax_size = mxx::allreduce(max_size, std::greater<std::size_t>(), comm);
    // if(comm.rank() == 0)
    //   std::cout << "Maximum Size  : " << rmax_size << std::endl;
      // remove local duplicates
      eliminateLocalDupes(read_pairs);
      if(csize > max_size) max_size = csize;

      // read_pairs.clear();
      // std::vector< std::pair<ReadIdType,ReadIdType> >().swap(read_pairs);
    });
    auto rmax_size = mxx::allreduce(max_size, std::greater<std::size_t>(), comm);
    auto total_size = mxx::allreduce(read_pairs.size(), comm);
    // size_tnow += total_size;
    if(comm.rank() == 0){
        std::cout << px_idx << " : SUB COMM SIZE : " << sub_size;
       std::cout << " : MAX BUCKET : " << rmax_size ;
       std::cout << " : TOTAL PAIRS : " << total_size ;
    }
    eliminateDuplicates(comm, read_pairs);
    total_size = mxx::allreduce(read_pairs.size(), comm);
    if(comm.rank() == 0)
       std::cout << " : UNQ PAIRS : " << total_size << std::endl;
    px_idx += 1;
    return read_pairs.size();
}

template<typename T, typename KT, typename HBT>
void generateSequencePairs(mxx::comm& comm,
                           bliss::io::file_data& file_data,
                           std::vector< std::pair<T, T> >& read_pairs) {
    static std::size_t gs_idx = 0;
  BL_BENCH_INIT(genpr);
  BL_BENCH_START(genpr);

  using SeqMinHashGeneratorType = SeqMinHashGenerator<
      MinHashFunctionBlock<KT, HBT>, KT>;
  using RHBType = typename SeqMinHashGeneratorType::value_type;

  std::vector<RHBType>  local_rhpairs;

  bliss::io::KmerFileHelper::template
    parse_file_data<SeqMinHashGeneratorType, FileParser,
                    SeqIterType>(file_data, local_rhpairs, comm);

  BL_BENCH_COLLECTIVE_END(genpr, "compute_hash", local_rhpairs.size(), comm);

  auto totalPairs = mxx::allreduce(local_rhpairs.size(), comm);
  if(comm.rank() == 0)
      std::cout << gs_idx << " : NR1 : " << totalPairs;

  BL_BENCH_START(genpr);
  auto cmp_fn = [&](const RHBType& x, const RHBType& y){
    return (x.seq_id < y.seq_id);
  };

  bool isSorted = mxx::is_sorted(local_rhpairs.begin(), local_rhpairs.end(),
                                 cmp_fn, comm);
  if(comm.rank() == 0)
      std::cout << " : SRT1 : " << (isSorted ? "Y" : "N");

  if(!isSorted)
    mxx::sort(local_rhpairs.begin(), local_rhpairs.end(), cmp_fn, comm);

  isSorted = mxx::is_sorted(local_rhpairs.begin(), local_rhpairs.end(),
                            cmp_fn, comm);
  if(comm.rank() == 0)
      std::cout << " : SRT2 : " << (isSorted ? "Y" : "N");

  totalPairs = mxx::allreduce(local_rhpairs.size(), comm);
  if(comm.rank() == 0)
      std::cout << " : NR2 : " << totalPairs;


  auto seq_idx = mxx::scan(local_rhpairs.size(), comm);
  if(local_rhpairs.size() > 0)
    seq_idx -= local_rhpairs.size();
  for(auto fitr = local_rhpairs.begin(); fitr != local_rhpairs.end(); fitr++, seq_idx++)
    fitr->seq_id = seq_idx;

  std::size_t j = 0;
  for(std::size_t i = 0; i < local_rhpairs.size();i++)
      if(local_rhpairs[i].is_good())
          local_rhpairs[j++] = local_rhpairs[i];
  local_rhpairs.resize(j);

  totalPairs = mxx::allreduce(local_rhpairs.size(), comm);
  if(comm.rank() == 0)
      std::cout << " : NR3 : " << totalPairs << std::endl;

  BL_BENCH_COLLECTIVE_END(genpr, "assign_ids", local_rhpairs.size(), comm);

  BL_BENCH_START(genpr);

  comm.with_subset(
      local_rhpairs.begin() != local_rhpairs.end(), [&](const mxx::comm& comm){
        mxx::sort(local_rhpairs.begin(), local_rhpairs.end(),
            [&] (const RHBType& x, const RHBType& y){
                for(auto i = 0; i < x.hash_values.size();i++){
                    if(x.hash_values[i] == y.hash_values[i]) continue;

                    if(x.hash_values[i] < y.hash_values[i]) return true;
                    else return false;
                }
                return false;
            }, comm);
      });
  BL_BENCH_COLLECTIVE_END(genpr, "sort_records", local_rhpairs.size(), comm);
  // auto total_pairs = mxx::allreduce(local_rhpairs.size());
  // if(comm.rank()  == 0)
  //    std::cout << "FINAL PAIR CNT : " << total_pairs << std::endl;

  BL_BENCH_START(genpr);
  uint64_t max_block_size = 0;
  comm.with_subset(
      local_rhpairs.begin() != local_rhpairs.end(), [&](const mxx::comm& comm){
        max_block_size =
          generateOverlapReadPairs(comm, local_rhpairs, read_pairs);
      });
  BL_BENCH_COLLECTIVE_END(genpr, "pair_gen", max_block_size, comm);
  BL_BENCH_REPORT_MPI_NAMED(genpr, "genpair", comm);
  gs_idx++;
}

template<typename T, int KMER_LENGTH, int BLOCK_SIZE>
void runFSO(mxx::comm& comm){

    std::string positionFile = run_params.positionFile;;
    std::vector<std::string>& inFiles = run_params.filenames;
    std::string outPrefix = run_params.outPrefix;
    unsigned readLength = run_params.read_length;
    uint32_t minOverlap = run_params.minOverlap;

    using KT = bliss::common::Kmer<KMER_LENGTH, Alphabet, WordType>;
    using HBT = std::array<HashValueType, BLOCK_SIZE>;

  // Read files
  if(comm.rank() == 0){
        std::cout << " Min Overlap       : " << minOverlap
                  << " Block Size      : " << run_params.hash_block_size
                  << " Block Count     : " << run_params.hash_block_count
                  << " Kmer Length     : " << run_params.kmer_length << std::endl;
  }

  BL_BENCH_INIT(rfso);
  std::vector<::bliss::io::file_data> file_data;

  BL_BENCH_START(rfso);
  load_file_data<KT>(comm, inFiles, file_data);
  BL_BENCH_COLLECTIVE_END(rfso, "read_files", file_data.size(), comm);

  comm.barrier();
  auto start = std::chrono::steady_clock::now();
  BL_BENCH_START(rfso);
  std::vector< std::pair<T, T> > read_pairs;
  for(auto i = 0u; i < run_params.hash_block_count; i++){
      if(seed_index >= (run_params.hash_block_count)){
        std::cout << "ERROR : seed index exceeded!!!" << std::endl;
        exit(1);
      }
      generateSequencePairs<T, KT, HBT>(comm, file_data.back(), read_pairs);
      seed_index++;
  }
  BL_BENCH_COLLECTIVE_END(rfso, "all_pairs", read_pairs.size(), comm);

  auto totalHashPairs = mxx::allreduce(read_pairs.size());
  if(comm.rank()  == 0)
      std::cout << "FINAL PAIR COUNT : " << totalHashPairs << std::endl;
  std::vector<bliss::io::file_data>().swap(file_data);
  comm.barrier();
  auto end = std::chrono::steady_clock::now();
  auto elapsed_time  = std::chrono::duration<double, std::milli>(end - start).count();

  if(!comm.rank())
      std::cout << "PGEN TIME (ms)  : " << elapsed_time << std::endl;

  comm.barrier();
  start = std::chrono::steady_clock::now();
  if(outPrefix.length() > 0) {
      std::stringstream outs;
      outs << outPrefix << "_"
           << (comm.rank() < 10 ? "000" :
               (comm.rank() < 100 ? "00" :
                (comm.rank() < 1000 ? "0" : "")))
           << comm.rank() << ".txt";

      std::string outputFile = outs.str();
      std::ofstream ofs(outputFile);
      if(ofs.good())
          for(auto px : read_pairs)
              ofs << px.first << " " << px.second << std::endl;
  }
  comm.barrier();
  end = std::chrono::steady_clock::now();
  elapsed_time  = std::chrono::duration<double, std::milli>(end - start).count();
  if(!comm.rank())
      std::cout << "WRITE OUTPUT (ms)  : " << elapsed_time << std::endl;


  BL_BENCH_START(rfso);
  if(positionFile.length() > 0)
      compareOverLaps(comm, positionFile, readLength, read_pairs,  minOverlap);
  BL_BENCH_COLLECTIVE_END(rfso, "compare_overlaps", read_pairs.size(), comm);
  BL_BENCH_REPORT_MPI_NAMED(rfso, "rfso_app", comm);
}

template<typename RPM>
int parse_args(int argc, char **argv,
               mxx::comm& comm, RPM& rparams){
  try { // try-catch block for commandline

    TCLAP::CmdLine cmd("Overlap Graph Construction", ' ', "0.1");

    // MPI friendly commandline output.
    bliss::utils::tclap::MPIOutput cmd_output(comm);
    cmd.setOutput(&cmd_output);

    // position file argument
    TCLAP::ValueArg<std::string> runTypeArg("t", "run_type",
                                            "Type of run : One of 'candidate', 'true', 'eval' ",
                                            true, "", "string", cmd);

   TCLAP::ValueArg<std::string> positionArg("p", "position_file",
                                             "Position for input file (full path)",
                                             false, "", "string", cmd);



    // output file argument
    TCLAP::ValueArg<std::string> outputArg("O", "output_prefix",
                                           "Prefix for output files, including directory",
                                           false, "", "string", cmd);

    // threshold argument
    TCLAP::ValueArg<uint32_t> threshArg("d", "overlap_threshold",
                                          "Minimum Overlap",
                                          false, 20, "int", cmd);

    // block count argument
    TCLAP::ValueArg<int> blockCountArg("B", "block_count",
                                       "Number of blocks",
                                       false, 250, "int", cmd);

    // block size argument
    TCLAP::ValueArg<int> blockSizeArg("T", "block_size",
                                      "Block Size",
                                      false, 3, "int", cmd);

    // max bucket size argument
    TCLAP::ValueArg<int> maxBucketSizeArg("M", "max_bucket_size",
                                          "Max Bucket Size",
                                          false, -1, "int", cmd);

    // kmer length argument
    TCLAP::ValueArg<int> kmerLengthArg("k", "kmer_length",
                                       "Kmer Length",
                                       false, 16, "int", cmd);


    // block count argument
    TCLAP::ValueArg<unsigned> readLengthArg("r", "read_length",
                                       "read length",
                                       false, 100, "unsigned", cmd);

    // candidate files only
    TCLAP::ValueArg<std::string> candFileArg("C", "candidate_file",
                                             "File with candidate pairs",
                                             false, "", "string", cmd);

    // true parirs files only
    TCLAP::ValueArg<std::string> trueFileArg("R", "true_file",
                                             "File with true pairs",
                                             false, "", "string", cmd);

    // input files
    TCLAP::UnlabeledMultiArg<std::string> fileArg("filenames", "FASTA or FASTQ file names",
                                                  false, "string", cmd);

    // Parse the argv array.
    cmd.parse( argc, argv );

    rparams.runType = runTypeArg.getValue();
    rparams.positionFile = positionArg.getValue();
    rparams.outPrefix = outputArg.getValue();
    rparams.candFile = candFileArg.getValue();
    rparams.trueFile = trueFileArg.getValue();
    rparams.filenames = fileArg.getValue();
    rparams.minOverlap = threshArg.getValue();
    run_params.hash_block_count = blockCountArg.getValue();
    run_params.hash_block_size = blockSizeArg.getValue();
    run_params.kmer_length = kmerLengthArg.getValue();
    rparams.read_length = readLengthArg.getValue();
    rparams.maxBucketSize = maxBucketSizeArg.getValue();

    if(comm.rank() == 0){
      std::cout << "Run Type             : " << rparams.runType << std::endl;
      std::cout << "Position File        : " << rparams.positionFile << std::endl;
      std::cout << "Output File Prefix   : " << rparams.outPrefix << std::endl;
      std::cout << "Candidate Pairs File : " << rparams.candFile << std::endl;
      std::cout << "True Pairs File      : " << rparams.trueFile << std::endl;
      std::cout << "No. Files            : " << rparams.filenames.size() << std::endl;
      if(rparams.filenames.size() > 0){
          std::cout << "Input File           : ";
          for(auto& fx : rparams.filenames)
              std::cout << "["<< fx << "] ";
          std::cout << std::endl;
      }
      std::cout << "Block Size           : " << run_params.hash_block_size << std::endl;
      std::cout << "Block Count          : " << run_params.hash_block_count << std::endl;
      std::cout << "Kmer Length          : " << run_params.kmer_length  << std::endl;
      std::cout << "Read Length          : " << rparams.read_length << std::endl;
      std::cout << "Min. Overlap         : " << rparams.minOverlap << std::endl;
      std::cout << "Max. Bucket Size     : " << rparams.maxBucketSize << std::endl;
    }
  } catch (TCLAP::ArgException &e)  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }

  if(rparams.runType.length() == 0)
      return 1;
  if(rparams.runType != "true" && rparams.runType != "candidate" && rparams.runType != "eval"){
      return 1;
  }
  if(rparams.runType == "candidate"){
      if(run_params.kmer_length < run_params.min_kmer_length || run_params.kmer_length > run_params.max_kmer_length){
          if(comm.rank() == 0)
              std::cout << "Kmer length should be in range ["
                        << run_params.min_kmer_length << ", "
                        << run_params.max_kmer_length << "]" << std::endl;
          return 1;
      }
      if(run_params.hash_block_size < run_params.min_hash_block_size ||
         run_params.hash_block_size > run_params.max_hash_block_size){
          if(comm.rank() == 0)
              std::cout << "Hash block size should be in range ["
                        << run_params.min_hash_block_size << ", "
                        << run_params.max_hash_block_size << "]" << std::endl;
          return 1;
      }
  }
  // validate
  return 0;
}

using RunPGFType = std::function< void(mxx::comm& comm) >;

std::array<RunPGFType, 56> FSO_FN_ARRAY  = {
    runFSO<uint64_t, 8, 2>,
    runFSO<uint64_t, 8, 3>,
    runFSO<uint64_t, 8, 4>,
    runFSO<uint64_t, 8, 5>,
    runFSO<uint64_t, 9, 2>,
    runFSO<uint64_t, 9, 3>,
    runFSO<uint64_t, 9, 4>,
    runFSO<uint64_t, 9, 5>,
    runFSO<uint64_t, 10, 2>,
    runFSO<uint64_t, 10, 3>,
    runFSO<uint64_t, 10, 4>,
    runFSO<uint64_t, 10, 5>,
    runFSO<uint64_t, 11, 2>,
    runFSO<uint64_t, 11, 3>,
    runFSO<uint64_t, 11, 4>,
    runFSO<uint64_t, 11, 5>,
    runFSO<uint64_t, 12, 2>,
    runFSO<uint64_t, 12, 3>,
    runFSO<uint64_t, 12, 4>,
    runFSO<uint64_t, 12, 5>,
    runFSO<uint64_t, 13, 2>,
    runFSO<uint64_t, 13, 3>,
    runFSO<uint64_t, 13, 4>,
    runFSO<uint64_t, 13, 5>,
    runFSO<uint64_t, 14, 2>,
    runFSO<uint64_t, 14, 3>,
    runFSO<uint64_t, 14, 4>,
    runFSO<uint64_t, 14, 5>,
    runFSO<uint64_t, 15, 2>,
    runFSO<uint64_t, 15, 3>,
    runFSO<uint64_t, 15, 4>,
    runFSO<uint64_t, 15, 5>,
    runFSO<uint64_t, 16, 2>,
    runFSO<uint64_t, 16, 3>,
    runFSO<uint64_t, 16, 4>,
    runFSO<uint64_t, 16, 5>,
    runFSO<uint64_t, 17, 2>,
    runFSO<uint64_t, 17, 3>,
    runFSO<uint64_t, 17, 4>,
    runFSO<uint64_t, 17, 5>,
    runFSO<uint64_t, 18, 2>,
    runFSO<uint64_t, 18, 3>,
    runFSO<uint64_t, 18, 4>,
    runFSO<uint64_t, 18, 5>,
    runFSO<uint64_t, 19, 2>,
    runFSO<uint64_t, 19, 3>,
    runFSO<uint64_t, 19, 4>,
    runFSO<uint64_t, 19, 5>,
    runFSO<uint64_t, 20, 2>,
    runFSO<uint64_t, 20, 3>,
    runFSO<uint64_t, 20, 4>,
    runFSO<uint64_t, 20, 5>,
    runFSO<uint64_t, 21, 2>,
    runFSO<uint64_t, 21, 3>,
    runFSO<uint64_t, 21, 4>,
    runFSO<uint64_t, 21, 5>
};


int genOlaps(mxx::comm& comm) {

  run_params.hash_seeds_size = (run_params.hash_block_size * run_params.hash_block_count);
  auto availableSeeds = sizeof(hash_seed_values);
  availableSeeds /= sizeof(hash_seed_values[0]);
  if(comm.rank() == 0)
       std::cout << "Available seeds : " << availableSeeds <<  std::endl;
  if(availableSeeds < run_params.hash_seeds_size){
      if(comm.rank() == 0)
          std::cout << "Not Enough seeds : " 
                    << availableSeeds << " < " << run_params.hash_seeds_size
                    << std::endl;
      return 1;
  }

  // runFSO
  comm.barrier();
  auto start = std::chrono::steady_clock::now();


  int aidx = (run_params.kmer_length - run_params.min_kmer_length) *
      (1 + run_params.max_hash_block_size - run_params.min_hash_block_size);
  aidx += (run_params.hash_block_size - run_params.min_hash_block_size);
  if(comm.rank() == 0){
      std::cout << "AIDX : " << aidx << " "
                << FSO_FN_ARRAY.size() <<  " "
                << ((bool)FSO_FN_ARRAY[aidx]) << " " << std::endl;
  }

  if(aidx >= 0 && aidx < FSO_FN_ARRAY.size() && ((bool)FSO_FN_ARRAY[aidx]))
      FSO_FN_ARRAY[aidx](comm);

  // runFSO<uint64_t, 21, 3>(comm, positionFile, filenames,
  //  outPrefix, readLength, threshold);
  //std::vector< std::pair<uint64_t, uint64_t> > read_pairs;
  //compareOverLaps(comm, positionFile, read_pairs, threshold);


  comm.barrier();
  auto end = std::chrono::steady_clock::now();
  auto elapsed_time  = std::chrono::duration<double, std::milli>(end - start).count();

  auto atx = mxx::allreduce(run_params.total_blocks, comm);
  auto ptx = mxx::allreduce(run_params.processed_blocks, comm);
  if(!comm.rank()){
    std::cout << "TOTAL BLOCKS : " << atx 
              << " : PROCESSED BLOCKS : " << ptx << std::endl;
  }
  if(!comm.rank())
      std::cout << "Time (ms)  : " << elapsed_time << std::endl;

  if(comm.rank() == 0)
    std::cout << "--------------------------------------" << std::endl;
  return 0;
}

int trueOlaps(mxx::comm& comm){

    std::string positionFile = run_params.positionFile;
    std::string outPrefix =  run_params.outPrefix;
    unsigned readLength = run_params.read_length;
    uint32_t minOverlap = run_params.minOverlap;

    // generate true pairs
    std::vector<std::pair<uint64_t, uint64_t>> truePairs;
    generateTruePairs(comm, positionFile, readLength, minOverlap, truePairs);

    std::stringstream outs;
    outs << outPrefix << "_"
         << (comm.rank() < 10 ? "000" :
             (comm.rank() < 100 ? "00" :
              (comm.rank() < 1000 ? "0" : "")))
         << comm.rank() << ".txt";

    std::string outputFile = outs.str();
    std::ofstream ofs(outputFile);
    if(ofs.good())
        for(auto px : truePairs)
            ofs << px.first << " " << px.second << std::endl;

    return 0;
}

int evalOlaps(mxx::comm& comm) {
    std::string candFile = run_params.candFile;
    std::string trueFile = run_params.trueFile;
  comm.barrier();
  auto start = std::chrono::steady_clock::now();

  std::vector< std::pair<uint64_t, uint64_t> > cand_pairs, true_pairs;
  std::size_t offsetStart, offsetEnd;
  std::vector<string> bufferStore;

  compute_offsets(comm, candFile, offsetStart, offsetEnd);
  read_block(comm, candFile, offsetStart, offsetEnd, bufferStore);
  cand_pairs.resize(bufferStore.size());
  std::size_t ridx = 0;
  for(auto rcd : bufferStore){
      uint64_t p1, p2;
      std::stringstream strStream(rcd);
      strStream >> p1;
      strStream >> p2;
      if(p1 < p2)
          cand_pairs[ridx] = std::make_pair(p1, p2);
      else
          cand_pairs[ridx] = std::make_pair(p2, p1);
      ridx++;
  }

  bufferStore.clear();
  std::vector<std::string>().swap(bufferStore);

  compute_offsets(comm, trueFile, offsetStart, offsetEnd);
  read_block(comm, trueFile, offsetStart, offsetEnd, bufferStore);
  true_pairs.resize(bufferStore.size());
  ridx = 0;
  for(auto rcd : bufferStore){
      uint64_t p1, p2;
      std::stringstream strStream(rcd);
      strStream >> p1;
      strStream >> p2;
      if(p1 < p2)
          true_pairs[ridx] = std::make_pair(p1, p2);
      else
          true_pairs[ridx] = std::make_pair(p2, p1);
      ridx++;
  }
  bufferStore.clear();
  std::vector<std::string>().swap(bufferStore);

  eliminateDuplicates(comm, cand_pairs);
  eliminateDuplicates(comm, true_pairs);

  auto totalPairs = mxx::allreduce(cand_pairs.size(), comm);
  if(comm.rank() == 0)
      std::cout << " CANDIDATE TOTAL         : " << totalPairs << std::endl;
  totalPairs = mxx::allreduce(true_pairs.size(), comm);
  if(comm.rank() == 0)
      std::cout << " TRUE TOTAL              : " << totalPairs << std::endl;

  computeSetDifference(comm, cand_pairs, true_pairs);

  comm.barrier();
  auto end = std::chrono::steady_clock::now();
  auto elapsed_time  = std::chrono::duration<double, std::milli>(end - start).count();
  if(!comm.rank())
      std::cout << "Time (ms)  : " << elapsed_time << std::endl;
  return 0;
}

int main(int argc, char** argv) {

  LOG_INIT(); // init logging

  mxx::env e(argc, argv); // MPI init
  mxx::comm comm;

  if(comm.rank() == 0)
    std::cout << "--------------------------------------" << std::endl;
  if (comm.rank() == 0)
    std::cout << "EXECUTABLE  : " << std::string(argv[0]) << std::endl;


  run_params.outPrefix.assign("./output");
  run_params.total_blocks = 0;
  run_params.processed_blocks = 0;
  // parse arguments
  if(parse_args(argc, argv, comm, run_params) != 0)
      return 1;

  if(run_params.runType == "true"){
      if(comm.rank() == 0) {
          std::cout << "---- Generate True Overlaps ----" << std::endl;
      }
      return trueOlaps(comm);
  } else if(run_params.runType == "eval") {
    if(comm.rank() == 0)
      std::cout << "---- Evaluate Overlaps ----" << std::endl;
    return evalOlaps(comm);
  } else {
    if(comm.rank() == 0)
      std::cout << "---- Generate candidate Overlaps ----" << std::endl;
    return genOlaps(comm);
  }
  return 0;
}
