# README #

This sofware implements a parallel Locality-Sensitive Hashing based heuristic algorithm to construct overlap graphs for large genomic datasets.

### Dependencies ###

* A modern, C++11 ready compiler such as `g++` version 4.7 or higher or `clang` version 3.2 or higher.
* The [cmake](www.cmake.org) build system (Version >= 2.8.11).
* A 64-bit Linux system.
* An MPI implementation. Tested with OpenMPI and MPICH only.


### Compilation ###

bruno/bliss librares are included as a submodules under the directory ext/bliss. Initialize the submodules as below, if they are not already initialized.

    git submodule init
    git submodule update

bruno/bliss also depends upon  mxx and google sparse hash, we initialize them as follows

    cd ext/bliss
    git submodule init
    git submodule update
    cd ../../

A bug in bruno is causes a error for compilation. So, we apply a patch on bruno/bliss as follows

    patch -p0 < bliss.patch

Now, we are read to build the executable. First, compile google's sparshash so that the header files are configured properly.

    cd ext/bliss/ext/sparshash
    ./configure
    make
    cd ../../../../

Next, create a build directory outside of the source directory. For example,

     mkdir build
     cd build

Finally, build the executable ***find_seq_overlaps*** as follows.

     cmake ../
     make
	 
### Usage of the executbale ###

The name of the executable is ***find_seq_overlaps***. Input files are provided in the FASTQ format. Input arguments are providede as follows:

    find_seq_overlaps  [-R <string>] [-C <string>] [-r <unsigned>] [-k
                        <int>] [-M <int>] [-T <int>] [-B <int>] [-d <int>]
                        [-O <string>] [-p <string>] -t <string> [--]
                        [--version] [-h] <file_names> ...


   Where: 

    -R <string>,  --true_file <string>
      File with true pairs

    -C <string>,  --candidate_file <string>
      File with candidate pairs

    -r <unsigned>,  --read_length <unsigned>
      read length

    -k <int>,  --kmer_length <int>
      Kmer Length

    -M <int>,  --max_bucket_size <int>
      Max Bucket Size

    -T <int>,  --block_size <int>
      Block Size

    -B <int>,  --block_count <int>
      Number of blocks

    -d <int>,  --overlap_threshold <int>
      Minimum Overlap

    -O <string>,  --output_prefix <string>
      Prefix for output files, including directory

    -p <string>,  --position_file <string>
      Position for input file (full path)

    -t <string>,  --run_type <string>
      (required)  Type of run : One of 'candidate', 'true', 'eval' 

    --,  --ignore_rest
      Ignores the rest of the labeled arguments following this flag.

    --version
      Displays version information and exits.

    -h,  --help
      Displays usage information and exits.

    <string>  (accepted multiple times)
      FASTA or FASTQ file names



