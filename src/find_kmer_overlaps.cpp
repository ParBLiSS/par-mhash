#include "bliss-config.hpp"

#include <cstdint>
#include <vector>
#include <string>

#include "io/file.hpp"
#include "io/kmer_parser.hpp"
#include "io/kmer_file_helper.hpp"
#include "index/kmer_index.hpp"
#include "index/kmer_hash.hpp"
#include "containers/distributed_unordered_map.hpp"

#include "mxx/env.hpp"
#include "mxx/comm.hpp"
#include "mxx/utils.hpp"
#include "tclap/CmdLine.h"
#include "utils/tclap_utils.hpp"

#include "run_cfg.hpp"
#include "io_utils.hpp"
#include "compare_overlaps.hpp"

static std::size_t readStartIndex = 0;

template<typename KmerType>
struct ReadsCountParser{
  // type of element generated by this parser.
  using value_type = char;
  using kmer_type = KmerType;
  static constexpr size_t window_size = kmer_type::size;

  bliss::partition::range<std::size_t> valid_range;

  ReadsCountParser(::bliss::partition::range<size_t> const & _valid_range)
    : valid_range(_valid_range) {};

    template <typename SeqType, typename OutputIt>
    OutputIt operator()(SeqType & read, OutputIt output_iter) {

        static_assert(std::is_same< SeqType,
                      bliss::io::FASTQSequence<typename SeqType::IteratorType> >::value,
                      "ReadLengthParser only supports FASTQ at the moment.");

        static_assert(std::is_same<value_type,
                      typename std::iterator_traits<OutputIt>::value_type>::value,
                      "output type and output container value type are not the same");

        // printf("First: pos %lu kmer %s\n", read.id.id,
        //       bliss::utils::KmerUtils::toASCIIString(*start).c_str());

        *output_iter = '1';
        return output_iter;
    }

};

// TODO: fix the maptypes

//using IdType = bliss::common::ShortSequenceKmerId;
using IdType = std::size_t;

using ValType = IdType;

template <typename KM>
using DistTrans = bliss::kmer::transform::lex_less<KM>;

template <typename KM>
using DistHash = bliss::kmer::hash::farm<KM, true>;

template <typename KM>
using StoreHash = bliss::kmer::hash::farm<KM, false>;

//template <typename Key>
//using MapParams = ::bliss::index::kmer::SingleStrandHashMapParams<Key, DistHash, StoreHash, DistTrans>;
//using SpecialKeys = ::bliss::kmer::hash::sparsehash::special_keys<KmerType, false>;
//template <typename Key>
//using MapParams = ::bliss::index::kmer::BimoleculeHashMapParams<Key, DistHash, StoreHash>;

template <typename Key>
using MapParams = ::bliss::index::kmer::CanonicalHashMapParams<Key, DistHash, StoreHash>;
using SpecialKeys = ::bliss::kmer::hash::sparsehash::special_keys<KmerType, true>;

//using MapType = dsc::densehash_multimap<KmerType, ValType, MapParams, SpecialKeys>;
using MapType = dsc::unordered_multimap<KmerType, ValType, MapParams>;


/**
 * @tparam TupleType       output value type of this parser.  not necessarily the same as the map's final storage type.
 */
template <typename TupleType>
struct KmerSeqTupleParser {

    /// type of element generated by this parser.  since kmer itself is parameterized, this is not hard coded.
  using value_type = TupleType;
  using kmer_type = typename ::std::tuple_element<0, value_type>::type;
  static constexpr size_t window_size = kmer_type::size;

  ::bliss::partition::range<size_t> valid_range;

  KmerSeqTupleParser(::bliss::partition::range<size_t> const & _valid_range) : valid_range(_valid_range) {};


  /**
   * @brief generate kmer-position pairs from 1 sequence.  result inserted into output_iter, which may be preallocated.
   * @param read          sequence object, which has pointers to the raw byte array.
   * @param output_iter   output iterator pointing to insertion point for underlying container.
   * @return new position for output_iter
   * @tparam SeqType      type of sequence.  inferred.
   * @tparam OutputIt     output iterator type, inferred.
   */
  template <typename SeqType, typename OutputIt, typename Predicate = ::fsc::TruePredicate>
  OutputIt operator()(SeqType const & read, OutputIt output_iter, Predicate const & pred = Predicate()) {

    static_assert(std::is_same<TupleType, typename ::std::iterator_traits<OutputIt>::value_type>::value,
            "output type and output container value type are not the same");
    static_assert(::std::tuple_size<TupleType>::value == 2, "kmer-pos-qual index data type should be a pair");

    using KmerType = typename std::tuple_element<0, TupleType>::type;
    using Alphabet = typename KmerType::KmerAlphabet;

    // filter out EOL characters
    using CharIter = NonEOLIter<typename SeqType::IteratorType>;
    // converter from ascii to alphabet values
    using BaseCharIterator = bliss::iterator::transform_iterator<CharIter, bliss::common::ASCII2<Alphabet> >;
    // kmer generation iterator
    using KmerIter = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;

    //== next figure out starting positions for the kmers, accounting for EOL char presenses.
    using IdType = typename std::tuple_element<1, TupleType>::type;
    // kmer position iterator type
    using IdIterType = bliss::iterator::ConstantIterator<IdType>;

    // use zip iterator to tie together the iteration of sequence raw data and id.
    using PairedIter = bliss::iterator::ZipIterator<typename SeqType::IteratorType, IdIterType>;
    using CharPosIter = NonEOLIter<PairedIter>;
    // now use 2 unzip iterators to access the values.  one of them advances.  all subsequent wrapping
    // iterators trickle back to the zip iterator above, again, one of the 2 unzip iterator will call operator++ on the underlying zip iterator
    using IdIter = bliss::iterator::AdvancingUnzipIterator<CharPosIter, 1>;

    // rezip the results
    using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIter, IdIter>;

    static_assert(std::is_same<typename std::iterator_traits<KmerIndexIterType>::value_type,
        TupleType>::value,
        "input iterator and output iterator's value types differ");


    // then compute and store into index (this will generate kmers and insert into index)
    if (::std::distance(read.seq_begin, read.seq_end) < KmerType::size) return output_iter;  // if too short...

    //== set up the kmer generating iterators.
    bliss::utils::file::NotEOL neol;
    KmerIter start(BaseCharIterator(CharIter(neol, read.seq_begin, read.seq_end), bliss::common::ASCII2<Alphabet>()), true);
    KmerIter end(BaseCharIterator(CharIter(neol, read.seq_end), bliss::common::ASCII2<Alphabet>()), false);


    //== set up the position iterators
    IdType seq_begin_id(readStartIndex);
    IdType seq_end_id(readStartIndex + 1);

    // tie chars and id together
    PairedIter pp_begin(read.seq_begin, IdIterType(seq_begin_id));
    PairedIter pp_end(read.seq_end, IdIterType(seq_end_id));

    // filter eol
    CharPosIter cp_begin(neol, pp_begin, pp_end);
    CharPosIter cp_end(neol, pp_end);

    // ==== extract new id and rezip iterators
    KmerIndexIterType index_start(start, IdIter(cp_begin));
    KmerIndexIterType index_end(end, IdIter(cp_end));


//    for (; index_start != index_end; ++index_start) {
//      auto tp = *index_start;
//
//      printf("TCP id = %lu, pos = %lu, kmer = %s\n", tp.second.get_id(), tp.second.get_pos(), bliss::utils::KmerUtils::toASCIIString(tp.first).c_str());
//
//      *output_iter = tp;
//
//    }

    ::bliss::partition::range<size_t> seq_range(read.seq_global_offset(), read.seq_global_offset() + read.seq_size());
    readStartIndex += 1;
    if (seq_range.contains(valid_range.end)) {
      // seq_range contains overlap.

      // not checking by end iterator at valid_range.end, since the NonEOLIter is a filter iterator that may skip over that pos.

      for (auto it = index_start; it != index_end; ++it, ++output_iter) {
        // check tail of window -> transform iterator, get base -> non EOL iterator, get base -> seq raw char iter.
        //int64_t valid_dist = valid_range.end - seq_range.start;
        //if ((*it).second.get_pos() >= valid_range.end) {
        //   break;
        //}
        //if (std::distance(read.seq_begin,
        //                  it.getTrailingIterator().getBaseIterator().getBaseIterator())
        //    >= valid_dist) {
        //  break;
        //}

        *output_iter = *it;
      }

      return output_iter;

    } else {

      if (::std::is_same<typename std::remove_reference<typename std::remove_cv<Predicate>::type>::type, ::fsc::TruePredicate>::value)
        return ::std::copy(index_start, index_end, output_iter);
      else
        return ::std::copy_if(index_start, index_end, output_iter, pred);
    }
  }

};

// using MapType = ::dsc::unordered_multimap_hashvec<
//  KmerType, ValType, MapParams>;
using SeqPositionIndexType = bliss::index::kmer::Index<MapType,
                                                       KmerSeqTupleParser<std::pair<typename MapType::key_type,
                                                                                    typename MapType::mapped_type> > >;


std::size_t getReadsCount(mxx::comm& comm, bliss::io::file_data& file_data){

  std::vector<char> vx;
  bliss::io::KmerFileHelper::template
    parse_file_data<ReadsCountParser<KmerType>, FileParser,
                    SeqIterType>(file_data, vx, comm);

  return vx.size();
}

void constructIndex(mxx::comm& comm, bliss::io::file_data& file_data,
                    SeqPositionIndexType& seqIdx){
  std::vector<typename SeqPositionIndexType::TupleType> tmp;
  bliss::io::KmerFileHelper::template
    parse_file_data<KmerSeqTupleParser<typename SeqPositionIndexType::TupleType>, FileParser,
                    SeqIterType>(file_data, tmp, comm);

  //std::cout << tmp.size() << std::endl;
  seqIdx.insert(tmp);
}

template<typename Iterator, typename ReadIdType=uint64_t>
uint64_t generatePairs(const mxx::comm& comm,
                       Iterator start_itr,
                       Iterator end_itr,
                       std::vector<std::pair<ReadIdType, ReadIdType>>& read_pairs){
    uint64_t nsize = std::distance(start_itr, end_itr);
    // return nsize;
    for(auto outer_itr = start_itr; outer_itr != end_itr; outer_itr++){
        auto inner_itr = outer_itr;
        for(inner_itr++; inner_itr != end_itr; inner_itr++){
            if(outer_itr->second == inner_itr->second) continue;
            // generate pair
            if(outer_itr->second < inner_itr->second)
                read_pairs.push_back(std::make_pair(outer_itr->second,
                                                    inner_itr->second));
            else
                read_pairs.push_back(std::make_pair(inner_itr->second,
                                                    outer_itr->second));
        }
    }
    return nsize;
}

template<typename T=uint64_t>
void runKSO(mxx::comm& comm,
            std::string positionFile,
            std::vector<std::string>& inFiles,
            std::string outPrefix,
            uint32_t threshold){
    BL_BENCH_INIT(kfso);
    std::vector<::bliss::io::file_data> file_data;
    size_t total = 0;

    BL_BENCH_START(kfso);
    load_file_data(comm, inFiles, file_data);
    BL_BENCH_COLLECTIVE_END(kfso, "read_files", total, comm);

    BL_BENCH_START(kfso);
    auto nLocalReads = getReadsCount(comm, file_data.back());
    readStartIndex = mxx::scan(nLocalReads, comm);
    if(nLocalReads > 0)
        readStartIndex -= readStartIndex;
    SeqPositionIndexType seqIdx(comm);
    constructIndex(comm, file_data.back(), seqIdx);
    //auto txval = std::distance(seqIdx.get_map().get_local_container().begin(), seqIdx.get_map().get_local_container().end());
    BL_BENCH_COLLECTIVE_END(kfso, "build_index", read_pairs.size(), comm);
    auto rval = mxx::allreduce(seqIdx.size(), comm);
    if(comm.rank() == 0)
        std::cout << "Index Size  : " << rval << std::endl;

    BL_BENCH_START(kfso);
    std::vector<std::pair<T, T>> read_pairs;
    auto txval = std::distance(seqIdx.get_map().get_local_container().begin(),
                               seqIdx.get_map().get_local_container().end());
    //std::size_t sum;
    //txval = generatePairs(comm, seqIdx.get_map().get_local_container().begin(),
    //                      seqIdx.get_map().get_local_container().end(),
    //                      read_pairs);

    std::size_t csize = 0, max_size = 0;
    auto rbv_itr = seqIdx.get_map().get_local_container().begin();
    auto prev_rbv = rbv_itr;

    for(;rbv_itr !=  seqIdx.get_map().get_local_container().end();rbv_itr++){
        if((*rbv_itr).first == (*prev_rbv).first)
            continue;
        csize = generatePairs(comm, prev_rbv, rbv_itr, read_pairs);
        if(csize > max_size) max_size = csize;
        prev_rbv = rbv_itr;
    }
    if(seqIdx.get_map().get_local_container().end() == rbv_itr && prev_rbv != rbv_itr)
        csize = generatePairs(comm, prev_rbv, rbv_itr, read_pairs);
    if(csize > max_size) max_size = csize;

    auto gen_pairs = mxx::allreduce(read_pairs.size(), comm);
    if(comm.rank() == 0)
        std::cout << "Total Gen Pairs : " << gen_pairs << std::endl;
    eliminateDuplicates(comm, read_pairs);
    BL_BENCH_COLLECTIVE_END(kfso, "generate_pairs", read_pairs.size(), comm);
    gen_pairs = mxx::allreduce(read_pairs.size(), comm);
    if(comm.rank() == 0)
        std::cout << "Unique Gen Pairs : " << gen_pairs << std::endl;
    auto rmax_size = mxx::allreduce(max_size, std::greater<std::size_t>(), comm);
    if(comm.rank() == 0)
        std::cout << "Maximum Size  : " << rmax_size << std::endl;

    BL_BENCH_START(kfso);
    compareOverLaps(comm, positionFile, read_pairs, threshold);
    BL_BENCH_COLLECTIVE_END(kfso, "compare_overlaps", read_pairs.size(), comm);

    BL_BENCH_REPORT_MPI_NAMED(kfso, "kfso_app", comm);
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
                                           false, "./output", "string", cmd);

    // threshold argument
    TCLAP::ValueArg<uint32_t> threshArg("d", "overlap_threshold",
                                          "Threshold for Overlap",
                                          false, 20, "int", cmd);

    // input files
    TCLAP::UnlabeledMultiArg<std::string> fileArg("filenames", "FASTA or FASTQ file names",
                                                  true, "string", cmd);


    // Parse the argv array.
    cmd.parse( argc, argv );

    positionFile = positionArg.getValue();
    outPrefix = outputArg.getValue();
    filenames = fileArg.getValue();
    threshold = threshArg.getValue();

    if(comm.rank() == 0){
      std::cout << "Position File   : " << positionFile << std::endl;
      std::cout << "Output File Pfx : " << outPrefix << std::endl;
      std::cout << "Input File      : " << filenames.front() << std::endl;
      std::cout << "Threshold       : " << threshold << std::endl;
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

  if (comm.rank() == 0) {
      std::cout << "--------------------------------------" << std::endl;
      std::cout << "EXECUTING : " << std::string(argv[0]) << std::endl;
  }


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

  comm.barrier();
  auto start = std::chrono::steady_clock::now();

  runKSO(comm, positionFile, filenames, outPrefix, threshold);

  comm.barrier();
  auto end = std::chrono::steady_clock::now();
  auto elapsed_time  = std::chrono::duration<double, std::milli>(end - start).count();

  if(!comm.rank()){
      std::cout << "Time (ms) -> : " << elapsed_time << std::endl;
      std::cout << "--------------------------------------" << std::endl;
  }
}
