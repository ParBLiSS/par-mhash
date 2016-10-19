
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
    rpItr++;
    strStream >> inValue;
    *rpItr = std::make_pair(startIdx + nRecords, inValue);
    rpItr++; startIdx++;
  }

}

template<typename T>
void eliminateDuplicates(const mxx::comm& comm,
                         std::vector<std::pair<T, T>>& truePairs){
  // sort and eliminate duplicates
  comm.with_subset(truePairs.begin() != truePairs.end(), [&](const mxx::comm& comm){
      mxx::sort(truePairs.begin(), truePairs.end(), comm);
      auto prevValue = mxx::right_shift(truePairs.back(), comm);
      auto cItr = truePairs.begin();
      auto vItr = cItr;
      if(comm.rank() == 0 && cItr != truePairs.begin()){
        prevValue = *cItr;
        cItr++;
      }
      for(;cItr != truePairs.end();cItr++){
        if(!(*cItr == prevValue)){
          *vItr = *cItr; vItr++;
          prevValue = *cItr;
        }
      }
      auto nPairs = std::distance(truePairs.begin(), vItr);
      truePairs.resize(nPairs);
      std::vector< std::pair<T,T> >(truePairs).swap(truePairs);
    });
}

template<typename T>
void generateTruePairs(const mxx::comm& comm,
                       std::string positionFile,
                       uint32_t threshold,
                       std::vector< std::pair<T, T> >& truePairs){
  BL_BENCH_INIT(cmpr);
  BL_BENCH_START(cmpr);
  // load the position file
  std::vector<std::pair<T, T>> readPosPairs;
  loadPositionFile(comm, positionFile, readPosPairs);
  BL_BENCH_COLLECTIVE_END(cmpr, "load_pos", readPosPairs.size(), comm);
  if(comm.rank() == 0 && readPosPairs.size() > 2)
      std::cout << "FIRST RND : " << readPosPairs[2] << std::endl;

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
      return ((x.second > y.second) ? (x.second - y.second) :
              (y.second - x.second)) < threshold;
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

  BL_BENCH_START(cmpr);
  eliminateDuplicates(comm, truePairs);
  BL_BENCH_COLLECTIVE_END(cmpr, "unique_pairs", truePairs.size(), comm);
  BL_BENCH_REPORT_MPI_NAMED(cmpr, "cmpr_app", comm);
  return;
}

template<typename T>
struct ReadPair{
    T first, second;
    char indicator;
    bool operator==(const ReadPair& other){
        return (first == other.first) &&
            (second == other.second);
    }
};

#define CMP_WRAP_TEMPLATE(...) __VA_ARGS__
namespace mxx {
    template <typename T>
    MXX_CUSTOM_TEMPLATE_STRUCT(CMP_WRAP_TEMPLATE(ReadPair<T>), \
                               first, second, indicator);
}

template<typename T>
void computeSetDifference(const mxx::comm& comm,
                          std::vector<std::pair<T, T>>& candidatePairs,
                          std::vector<std::pair<T, T>>& truePairs){
    BL_BENCH_INIT(cmpr);
    BL_BENCH_START(cmpr);
    uint64_t setDiff[2] = {0, 0};
    uint64_t intersectPairs = 0;
    comm.with_subset(
      candidatePairs.size() + truePairs.size() > 0, [&](const mxx::comm& comm){
          std::vector<ReadPair<T>> mergedPairs(candidatePairs.size() + truePairs.size());
          auto rdItr = mergedPairs.begin();
          for(auto px : candidatePairs){
              rdItr->first = px.first; rdItr->second = px.second;
              rdItr->indicator = 0; rdItr++;
          }
          for(auto px : truePairs){
              rdItr->first = px.first; rdItr->second = px.second;
              rdItr->indicator = 1; rdItr++;
          }

          mxx::sort(mergedPairs.begin(), mergedPairs.end(),
                    [&](const ReadPair<T>& x, const ReadPair<T>& y){
                        return (x.first < y.first) ||
                            ((x.first == y.first) && (x.second < y.second));
                    }, comm);

          // assuming a pair can appear atmost once
          auto curVal = mxx::right_shift(mergedPairs.back(), comm);
          auto nextVal = mxx::left_shift(mergedPairs.back(), comm);
          auto itrx = mergedPairs.begin();
          if(comm.rank() == 0 && itrx != mergedPairs.end()){
              curVal = *itrx; itrx++;
          } else if( itrx != mergedPairs.end() ){
              if(!(curVal == *itrx)){
                  curVal = *itrx; itrx++;
              }
          }
          for(; itrx != mergedPairs.end(); itrx++){
              if(!(*itrx == curVal)){
                  setDiff[curVal.indicator] += 1;
                  curVal = *itrx;
              } else {
                  intersectPairs += 1; itrx++;
                  if(itrx != mergedPairs.end())
                      curVal = *itrx;
                  else
                      return;
              }
          }
          if(comm.rank() < (comm.size() - 1)){
              if(!(curVal == nextVal))
                  setDiff[curVal.indicator] += 1;
          } else { // TODO: verify this
              setDiff[curVal.indicator] += 1;
          }
      });
    auto inCandidates = mxx::allreduce(setDiff[0], comm);
    auto inTruth = mxx::allreduce(setDiff[1], comm);
    auto totalIntersect = mxx::allreduce(intersectPairs, comm);
    if(comm.rank() == 0){
        std::cout << "|CANDIDATES \\   TRUTH|  : " << inCandidates << std::endl;
        std::cout << "|TRUTH \\   CANDIDATES|  : " << inTruth << std::endl;
        std::cout << "|TRUTH \\cup CANDIDATES|  : " << totalIntersect << std::endl;
    }
    BL_BENCH_COLLECTIVE_END(cmpr, "count_setdiff", truePairs.size(), comm);
    BL_BENCH_REPORT_MPI_NAMED(cmpr, "cmpr_app", comm);
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

    //
    computeSetDifference(comm, candidatePairs, truePairs);
}


#endif
