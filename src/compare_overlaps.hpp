
#ifndef COMPARE_OVERLAPS_HPP
#define COMPARE_OVERLAPS_HPP

#include "utils/logging.h"
#include "mxx/sort.hpp"
#include "mxx/comm.hpp"
#include "mxx/utils.hpp"
#include "tclap/CmdLine.h"
#include "utils/tclap_utils.hpp"

#include "io_utils.hpp"

template <typename T>
void findTruePairs (std::vector<std::pair<T,T> >& read_pos,
                    std::vector<std::pair<T,T> >& read_pairs,
                    uint64_t thresh) {
    std::size_t up=0, down=1;
    while(true) {
        if (up >= read_pos.size()) {
            break;
        }
        while (true) {
            if (down >= read_pos.size()) {
                break;
            }
            if (read_pos[down].second < (read_pos[up].second + thresh)) {
                down++;
            } else {
                break;
            }
        }
        assert(down > up);
        // gen pairs [up,down)
        for (std::size_t pos=up+1; pos < down; pos++) {
            read_pairs.push_back(std::make_pair(read_pos[up].first,
                                                read_pos[pos].first));
        }
        up++;
    }
}


template<typename T>
void loadPositionFile(const mxx::comm& comm,
                      std::string positionFile,
                      std::vector<std::pair<T, T>>& readPairs){

  // compute start and end offsets corresponding this rank
  T offsetStart, offsetEnd;
  compute_offsets(comm, positionFile, offsetStart, offsetEnd);

  // load the block
  std::vector<std::string> bufferStore;
  read_block(comm, positionFile, offsetStart, offsetEnd, bufferStore);

  // generate the pairs
  // input position file has format for paired end reads:
  // position_left_end position_right_end
  readPairs.resize(2 * bufferStore.size());
  auto startIdx = mxx::scan(bufferStore.size(), comm);
  auto nRecords = mxx::allreduce(bufferStore.size(), comm);
  if(bufferStore.size() > 0)
    startIdx = startIdx - bufferStore.size();
  auto rpItr = readPairs.begin();
  for(auto rcd : bufferStore){
    T inValue;
    std::stringstream strStream(rcd);
    strStream >> inValue;
    *rpItr = std::make_pair(startIdx, inValue);
    rpItr++; startIdx++;
    strStream >> inValue;
    *rpItr = std::make_pair(startIdx + nRecords, inValue);
    rpItr++; startIdx++;
  }

}

template<typename T>
void generateTruePairs(const mxx::comm& comm,
                       std::string positionFile,
                       uint32_t threshold,
                       std::vector<std::pair<T, T>>& truePairs){
  BL_BENCH_INIT(cmpr);
  BL_BENCH_START(cmpr);
  // load the position file
  std::vector<std::pair<T, T>> readPosPairs;
  loadPositionFile(comm, positionFile, readPosPairs);
  BL_BENCH_COLLECTIVE_END(cmpr, "load_pos", readPosPairs.size(), comm);
  if(comm.rank() == 0 && readPosPairs.size() > 0)
      std::cout << "FIRST RND : " << readPosPairs[0] << std::endl;

  // sort by positionFile
  BL_BENCH_START(cmpr);
  comm.with_subset(
      readPosPairs.begin() != readPosPairs.end(), [&](const mxx::comm& comm){
          mxx::sort(readPosPairs.begin(), readPosPairs.end(),
                    [&](const std::pair<T, T>& x,
                        const std::pair<T, T>& y){
                        return x.second < y.second;
                    }, comm);
      });
  BL_BENCH_COLLECTIVE_END(cmpr, "sort_pairs", readPosPairs.size(), comm);
  if(comm.rank() == 0 && readPosPairs.size() > 0)
      std::cout << "FIRST SRT : " << readPosPairs.front() << std::endl;
  if(comm.rank() == (comm.size() - 1) && readPosPairs.size() > 0)
      std::cout << "LAST SRT : " << readPosPairs.back() << std::endl;
  auto totalPosPairs = mxx::allreduce(readPosPairs.size(), comm);
  if(comm.rank() == 0)
     std::cout << "POS TOTAL : " << totalPosPairs << std::endl;

  BL_BENCH_START(cmpr);
  // get the straddling region
  uint64_t startOffset, endOffset;
  auto posCompare = [&](const std::pair<T,T>& x,
                        const std::pair<T,T>& y){
      return ((x.second > y.second) ?
              (x.second - y.second) : (y.second - x.second)) < threshold;
  };

  std::vector<std::pair<T,T>> straddleRegion;
  shiftStraddlingRegion(comm, readPosPairs,
                        startOffset, endOffset,
                        straddleRegion, posCompare);
  BL_BENCH_COLLECTIVE_END(cmpr, "shift_straddle", straddleRegion.size(), comm);

  BL_BENCH_START(cmpr);
  std::vector<std::pair<T,T>> localRegion;
  localRegion.resize(endOffset + straddleRegion.size()); // TODO: check if this is correct
  std::copy(readPosPairs.begin(), readPosPairs.begin() + endOffset,
            localRegion.begin());
  std::copy(straddleRegion.begin(), straddleRegion.end(),
            localRegion.begin() + endOffset);
  BL_BENCH_COLLECTIVE_END(cmpr, "local_region", localRegion.size(), comm);

  // generate pairs
  BL_BENCH_START(cmpr);
  findTruePairs(localRegion, truePairs, threshold);
  BL_BENCH_COLLECTIVE_END(cmpr, "true_pairs", truePairs.size(), comm);
  BL_BENCH_REPORT_MPI_NAMED(cmpr, "cmpr_app", comm);
  return;
}

template<typename T>
void compareOverLaps(const mxx::comm& comm,
                     std::string positionFile,
                     std::vector<std::pair<T, T>>& candidatePairs,
                     uint32_t threshold){

    // generate true pairs
    std::vector<std::pair<T, T>> truePairs;
    generateTruePairs(comm, positionFile, threshold, truePairs);

    auto totalTruePairs = mxx::allreduce(truePairs.size(), comm);
    if(comm.rank() == 0)
       std::cout << "TRUE TOTAL : " << totalTruePairs << std::endl;
}


#endif
