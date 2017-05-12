#ifndef RUN_CFG_HPP
#define RUN_CFG_HPP

#include "common/alphabets.hpp"
#include "common/kmer.hpp"
#include "common/base_types.hpp"
#include "utils/kmer_utils.hpp"

#include "io/file.hpp"
#include "io/sequence_iterator.hpp"
#include "io/filtered_sequence_iterator.hpp"
#include "utils/file_utils.hpp"


#if (pDNA == 4)
using Alphabet = bliss::common::DNA;
#elif (pDNA == 5)
using Alphabet = bliss::common::DNA5;
#elif (pDNA == 16)
using Alphabet = bliss::common::DNA16;
#endif

// constants
#if (pPARSER == FASTQ)
#define FileParser bliss::io::FASTQParser
#elif (pPARSER == FASTA)
#define FileParser bliss::io::FASTAParser
#endif


using FileReaderType = bliss::io::parallel::partitioned_file<
  bliss::io::posix_file, FileParser>;

// Kmer, FilerReader and Iterator data types
using EdgeEncoding = Alphabet;
using KmerType8 = bliss::common::Kmer<8, Alphabet, WordType>;
using KmerType16 = bliss::common::Kmer<16, Alphabet, WordType>;

template <typename Iterator, template <typename> class SeqParser>
using SeqIterType = bliss::io::SequencesIterator<Iterator, SeqParser>;

template <typename Iter>
using NonEOLIter = bliss::iterator::filter_iterator<bliss::utils::file::NotEOL, Iter>;

// Hash Value types
using HashValueType = uint64_t;
//using HashBlockType = std::array<HashValueType, hash_block_size>;
using HashBlockType2 = std::array<HashValueType, 2>;
using HashBlockType3 = std::array<HashValueType, 3>;
using HashBlockType4 = std::array<HashValueType, 4>;
using HashBlockType5 = std::array<HashValueType, 5>;
// using HashBlockType = std::vector<HashValueType>;

#define WRAP_TEMPLATE(...) __VA_ARGS__

#endif
