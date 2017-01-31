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

//
// @tparam KmerType output value type of this parser. not necessarily the same
//                  as the map's final storage type.
//
template <typename LT, typename KmerType>
struct ReadLengthGenerator {

    // type of element generated by this parser.
    // using value_type = ReadMinHashBlock< typename MinHashFunctionBlockType::value_type >;
    using value_type = LT;
    using kmer_type = KmerType;
    static constexpr size_t window_size = kmer_type::size;

    bliss::partition::range<std::size_t> valid_range;

    ReadLengthGenerator(::bliss::partition::range<size_t> const & _valid_range)
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

        // bliss::partition::range<size_t> seq_range(
        //    read.seq_global_offset(), read.seq_global_offset() + read.seq_size());
        value_type hrv;
        hrv = read.seq_size();
        *output_iter = hrv;
        return output_iter;
    }
};


template<typename T>
void generateSequenceLengths(mxx::comm& comm,
                             bliss::io::file_data& file_data,
                             std::vector<T>& read_pairs) {
  BL_BENCH_INIT(genpr);
  BL_BENCH_START(genpr);

  using ReadLengthGeneratorType =  ReadLengthGenerator<T, KmerType>;
  std::vector<T>  local_pairs;
  auto cmp_fn = [&](const T& x, const T& y){
      return (x < y);
  };

  bliss::io::KmerFileHelper::template
      parse_file_data<ReadLengthGeneratorType, FileParser,
                    SeqIterType>(file_data, local_pairs, comm);

  BL_BENCH_COLLECTIVE_END(genpr, "compute_lengths", local_rhpairs.size(), comm);

  BL_BENCH_START(genpr);

  if(!mxx::is_sorted(local_pairs.begin(), local_pairs.end(),
                     cmp_fn, comm))
      mxx::sort(local_pairs.begin(), local_pairs.end(), cmp_fn, comm);

  BL_BENCH_COLLECTIVE_END(genpr, "sort_lengths", read_pairs.size(), comm);
  BL_BENCH_START(genpr);

  auto rsz = read_pairs.size();

  read_pairs.resize(rsz + local_pairs.size());

  std::copy(local_pairs.begin(), local_pairs.end(),
            read_pairs.begin() + rsz);
  if(!mxx::is_sorted(read_pairs.begin(), read_pairs.end(),
                     cmp_fn, comm))
      mxx::sort(read_pairs.begin(), read_pairs.end(), cmp_fn, comm);

  BL_BENCH_COLLECTIVE_END(genpr, "sort_full", local_rhpairs.size(), comm);

  // std::vector<T>(read_pairs).swap(read_pairs);
}

template<typename T=uint64_t>
void runFSO(mxx::comm& comm,
            std::vector<std::string>& inFiles,
            std::string outFile,
            uint32_t binSize){

  // Read files
  BL_BENCH_INIT(rfso);
  std::vector<::bliss::io::file_data> file_data;

  BL_BENCH_START(rfso);
  load_file_data(comm, inFiles, file_data);
  BL_BENCH_COLLECTIVE_END(rfso, "read_files", total, comm);

  comm.barrier();
  auto start = std::chrono::steady_clock::now();
  BL_BENCH_START(rfso);
  std::vector<T> read_pairs;
  for(auto i =0u; i < file_data.size(); i++)
      generateSequenceLengths(comm, file_data[i], read_pairs);
  BL_BENCH_COLLECTIVE_END(rfso, "all_pairs", read_pairs.size(), comm);

  auto totalLengths = mxx::allreduce(read_pairs.size());
  if(comm.rank()  == 0)
      std::cout << "FINAL PAIR COUNT : " << totalLengths << std::endl;
  std::vector<bliss::io::file_data>().swap(file_data);


  auto end = std::chrono::steady_clock::now();
  auto elapsed_time  = std::chrono::duration<double, std::milli>(end - start).count();
  if(!comm.rank())
      std::cout << "HGEN TIME (ms)  : " << elapsed_time << std::endl;
  BL_BENCH_START(rfso);
  if(comm.rank() == 0){
      std::cout << " Threshold       : " << binSize
                << " Kmer Length     : " << HASH_KMER_SIZE  << std::endl;
  }
  std::ofstream outf(outFile, std::ofstream::out);
  for(auto rx : read_pairs)
      outf << rx << std::endl;

  BL_BENCH_REPORT_MPI_NAMED(rfso, "rfso_app", comm);
}

void parse_args(int argc, char **argv,
                mxx::comm& comm,
                std::vector<std::string>& filenames,
                std::string& outPrefix,
                uint32_t& binSize){
  try { // try-catch block for commandline

    TCLAP::CmdLine cmd("Length Histogram", ' ', "0.1");

    // MPI friendly commandline output.
    bliss::utils::tclap::MPIOutput cmd_output(comm);
    cmd.setOutput(&cmd_output);


    // output file argument
    TCLAP::ValueArg<std::string> outputArg("O", "output_prefix",
                                           "Prefix for output files, including directory",
                                           false, "", "string", cmd);

    // bin size
    TCLAP::ValueArg<uint32_t> binSizeArg("d", "overlap_threshold",
                                         "Threshold for Overlap",
                                         false, 100, "int", cmd);

    // input files
    TCLAP::UnlabeledMultiArg<std::string> fileArg("filenames", "FASTA or FASTQ file names",
                                                  true, "string", cmd);


    // Parse the argv array.
    cmd.parse( argc, argv );

    outPrefix = outputArg.getValue();
    filenames = fileArg.getValue();
    binSize = binSizeArg.getValue();

    if(comm.rank() == 0){
      std::cout << "Output File Pfx : " << outPrefix << std::endl;
      std::cout << "Input File      : " << filenames.front() << std::endl;
      std::cout << "Threshold       : " << binSize << std::endl;
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

  if(comm.rank() == 0)
    std::cout << "--------------------------------------" << std::endl;
  if (comm.rank() == 0)
    std::cout << "EXECUTABLE  : " << std::string(argv[0]) << std::endl;


  std::vector<std::string> filenames;
  std::string outPrefix;
  uint32_t threshold;
  outPrefix.assign("./output");
  // parse arguments
  parse_args(argc, argv, comm,
             filenames, outPrefix, threshold);

  std::stringstream outs;
  outs << outPrefix << "_"
       << (comm.rank() < 10 ? "000" :
           (comm.rank() < 100 ? "00" :
            (comm.rank() < 1000 ? "0" : "")))
       << comm.rank() << ".txt";
  std::string outputFile = outs.str();

  // runFSO
  comm.barrier();
  auto start = std::chrono::steady_clock::now();

  runFSO(comm, filenames, outputFile, threshold);
  //std::vector< std::pair<uint64_t, uint64_t> > read_pairs;
  //compareOverLaps(comm, positionFile, read_pairs, threshold);

  comm.barrier();
  auto end = std::chrono::steady_clock::now();
  auto elapsed_time  = std::chrono::duration<double, std::milli>(end - start).count();

  if(!comm.rank())
      std::cout << "Time (ms)  : " << elapsed_time << std::endl;

  if(comm.rank() == 0)
    std::cout << "--------------------------------------" << std::endl;
}
