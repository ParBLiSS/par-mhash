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

const static std::size_t hash_block_size = 2;
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
#define HASH_KMER_SIZE 16
#endif

#ifndef READ_THRESHOLD
#define READ_THRESHOLD 1000
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

// Hash Value types
using HashValueType = uint64_t;
using HashBlockType = std::array<HashValueType, hash_block_size>;

#define WRAP_TEMPLATE(...) __VA_ARGS__

#endif
