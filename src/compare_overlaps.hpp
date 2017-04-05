
#ifndef COMPARE_OVERLAPS_HPP
#define COMPARE_OVERLAPS_HPP

#include "utils/logging.h"
#include "mxx/sort.hpp"
#include "mxx/comm.hpp"
#include "mxx/utils.hpp"
#include "tclap/CmdLine.h"
#include "utils/tclap_utils.hpp"

#include "io_utils.hpp"

static uint64_t RPOS_MISSING;
static uint32_t RLEN_MISSING;

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
struct ReadIPL{
    T rid, rpos;
    unsigned rlen;

    inline void set(T id, T pos, unsigned rl){
        rid = id; rpos = pos; rlen = rl;
    }
};

#define CMP_WRAP_TEMPLATE(...) __VA_ARGS__
namespace mxx {
    template <typename T>
    MXX_CUSTOM_TEMPLATE_STRUCT(CMP_WRAP_TEMPLATE(ReadIPL<T>), \
                               rid, rpos, rlen);
}

template<typename T>
struct PosReadPair{
    T first, second;
    uint64_t first_position, second_position;
    uint32_t first_length, second_length;
    bool operator==(const PosReadPair& other){
        return (first == other.first) &&
            (second == other.second);
    }
    PosReadPair(){
        first_position = second_position = RPOS_MISSING;
        first_length = second_length = RLEN_MISSING;
    }

    bool is_pos_bad() const{
        return (first_position == RPOS_MISSING) ||
            (second_position == RPOS_MISSING);
    }
};

#define CMP_WRAP_TEMPLATE(...) __VA_ARGS__
namespace mxx {
    template <typename T>
    MXX_CUSTOM_TEMPLATE_STRUCT(CMP_WRAP_TEMPLATE(PosReadPair<T>), \
                               first, second, first_position, second_position, first_length, second_length);
}


bool has_suffix(const std::string &str, const std::string &suffix)
{
    return str.size() >= suffix.size() &&
        str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

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
            if(read_pos[up].first == read_pos[pos].first) continue;
            if(read_pos[up].first < read_pos[pos].first)
                read_pairs.push_back(std::make_pair(read_pos[up].first,
                                                    read_pos[pos].first));
            else
                read_pairs.push_back(std::make_pair(read_pos[pos].first,
                                                    read_pos[up].first));
        }
        up++;
    }
}

template<typename T>
void findTruePairs(std::vector< ReadIPL<T> >& read_pos,
                   std::vector<std::pair<T,T> >& read_pairs,
                   uint64_t minOverlap) {
    long unsigned int up=0, down=1;
    assert(read_pos.size() == read_len.size());
    while(true) {
        if (up >= read_pos.size()) {
            break;
        }
        while (true) {
            if (down >= read_pos.size()) {
                break;
            }
            //if (read_pos[down].second - read_pos[up].second < thresh) {
            if (read_pos[up].rpos + read_pos[up].rlen >= read_pos[down].rpos + minOverlap) {
                down++;
            }
            else {
                break;
            }
        }
        assert(down > up);
        // gen pairs [up,down)
        for (std::size_t pos=up+1; pos < down; pos++) {
            read_pairs.emplace_back(std::make_pair(read_pos[up].rid,
                                                   read_pos[pos].rid));
        }
        up++;
    }
}


template<typename T>
void eliminateDuplicates(const mxx::comm& comm,
                         std::vector<std::pair<T, T>>& readIdPairs){
  // sort and eliminate duplicates
  comm.with_subset(readIdPairs.begin() != readIdPairs.end(), [&](const mxx::comm& comm){
      mxx::sort(readIdPairs.begin(), readIdPairs.end(), comm);
      auto prevValue = mxx::right_shift(readIdPairs.back(), comm);
      auto cItr = readIdPairs.begin();
      auto vItr = cItr;
      if(comm.rank() == 0 && cItr != readIdPairs.end()){
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
    });
}

template<typename T>
void loadSimulatedPosFile(const mxx::comm& comm,
                          std::string positionFile,
                          unsigned rlen,
                          std::vector< ReadIPL<T> >& readPairs){

  if(comm.rank() == 0)
     std::cout << "---- Start Loading SIM POS File : [" << positionFile << "] ----" << std::endl;

  // compute start and end offsets corresponding this rank
  T offsetStart, offsetEnd;
  compute_offsets(comm, positionFile, offsetStart, offsetEnd);

  // load the block
  std::vector<std::string> bufferStore ;
  read_block(comm, positionFile, offsetStart, offsetEnd, bufferStore);
  auto totalPosLines = mxx::allreduce(bufferStore.size(), comm);
  if(comm.rank() == 0)
      std::cout << "POS FILE RECORDS : " << totalPosLines << std::endl;
  if(totalPosLines == 0)
      return;
  auto vx = mxx::left_shift(bufferStore.front(), comm);

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
    rpItr->set(startIdx, inValue, rlen);
    rpItr++;
    strStream >> inValue;
    rpItr->set(startIdx + nRecords, inValue, rlen);
    rpItr++; startIdx++;
  }

  totalPosLines = mxx::allreduce(readPairs.size(), comm);

  if(comm.rank() == 0)
     std::cout << "POS DATA : " << totalPosLines << std::endl;
  if(comm.rank() == 0)
     std::cout << "---- Finish Loading SIM POS File ----" << std::endl;
}

template<typename T>
void loadSAMPosFile(const mxx::comm& comm,
                    std::string positionFile,
                    std::vector<ReadIPL<T>>& readPairs){

    if(comm.rank() == 0)
       std::cout << "---- Start Loading SAM POS File: [" << positionFile <<  "] ----" << std::endl;

    // compute start and end offsets corresponding this rank
    T offsetStart, offsetEnd;
    compute_offsets(comm, positionFile, offsetStart, offsetEnd);

    // load the block
    std::vector<std::string> bufferStore ;
    read_block(comm, positionFile, offsetStart, offsetEnd, bufferStore);
    auto totalPosLines = mxx::allreduce(bufferStore.size(), comm);
    if(comm.rank() == 0)
        std::cout << "POS FILE RECORDS : " << totalPosLines << std::endl;
    if(totalPosLines == 0)
        return;
    auto vx = mxx::left_shift(bufferStore.front(), comm);
    // generate the pairs
    // input position file has format for paired end reads:
    // position_left_end position_right_end
    readPairs.resize(bufferStore.size());
    std::size_t rix = 0;
    for(auto rcd : bufferStore){
        T readIdx, inValue;
        unsigned rlen;
        std::stringstream strStream(rcd);
        strStream >> readIdx;
        strStream >> inValue;
        strStream >> rlen;
        readPairs[rix].set(readIdx, inValue, rlen);
        rix++;
    }

    auto totalRecords = mxx::allreduce(rix, comm);
    if(comm.rank() == 0)
        std::cout << "POS DATA : " << totalRecords << std::endl;
    if(comm.rank() == 0)
       std::cout << "---- Finish Loading SAM POS File ----" << std::endl;
}

template<typename T>
void loadPositionFile(const mxx::comm& comm,
                      std::string positionFile,
                      unsigned readLength,
                      std::vector< ReadIPL<T> >& readPairs){
    if(has_suffix(positionFile, "map")){
        loadSAMPosFile(comm, positionFile, readPairs);
    } else {
        loadSimulatedPosFile(comm, positionFile, readLength, readPairs);
    }
}

template<typename T>
void generateTruePairs(const mxx::comm& comm,
                       std::string positionFile,
                       unsigned rlen,
                       uint32_t minOverlap,
                       std::vector< std::pair<T, T> >& truePairs){
  BL_BENCH_INIT(cmpr);
  BL_BENCH_START(cmpr);
  // load the position file
  std::vector<ReadIPL<T>> readPosPairs;
  loadPositionFile(comm, positionFile, rlen, readPosPairs);
  BL_BENCH_COLLECTIVE_END(cmpr, "load_pos", readPosPairs.size(), comm);
  auto totalPosLines = mxx::allreduce(readPosPairs.size(), comm);
  if(totalPosLines == 0)
      return;

  // sort by positionFile
  BL_BENCH_START(cmpr);
  comm.with_subset(
      readPosPairs.begin() != readPosPairs.end(), [&](const mxx::comm& comm){
          mxx::sort(readPosPairs.begin(), readPosPairs.end(),
                    [&](const ReadIPL<T>& x,
                        const ReadIPL<T>& y){
                        return x.rpos < y.rpos;
                    }, comm);
      });
  BL_BENCH_COLLECTIVE_END(cmpr, "sort_pairs", readPosPairs.size(), comm);
  auto totalPosPairs = mxx::allreduce(readPosPairs.size(), comm);
  if(comm.rank() == 0)
     std::cout << "POS TOTAL : " << totalPosPairs << std::endl;

  BL_BENCH_START(cmpr);

  // get the straddling region
  uint64_t startOffset, endOffset;
  auto posCompare = [&](const ReadIPL<T>& x,
                        const ReadIPL<T>& y){
      // return ((x.rpos > y.rpos) ? (x.rpos - y.rpos) :
      //         (y.rpos - x.rpos)) < threshold;
      T pos_delta = ((x.rpos > y.rpos) ? (x.rpos - y.rpos) :
                     (y.rpos - x.rpos));
      T o_delta = ((x.rlen > minOverlap) ?
                   (x.rlen - minOverlap) : (minOverlap - x.rlen));
      return pos_delta < o_delta ;
  };

  std::vector<ReadIPL<T>> straddleRegion;
  shiftStraddlingRegion(comm, readPosPairs,
                        startOffset, endOffset,
                        straddleRegion, posCompare);
  BL_BENCH_COLLECTIVE_END(cmpr, "shift_straddle", straddleRegion.size(), comm);

  BL_BENCH_START(cmpr);
  std::vector<ReadIPL<T>> localRegion;
  localRegion.resize(endOffset + straddleRegion.size()); // TODO: check if this is correct
  std::copy(readPosPairs.begin(), readPosPairs.begin() + endOffset,
            localRegion.begin());
  std::vector<ReadIPL<T>>().swap(readPosPairs);
  std::copy(straddleRegion.begin(), straddleRegion.end(),
            localRegion.begin() + endOffset);
  std::vector<ReadIPL<T>>().swap(straddleRegion);
  BL_BENCH_COLLECTIVE_END(cmpr, "local_region", localRegion.size(), comm);

  // generate pairs
  BL_BENCH_START(cmpr);
  findTruePairs(localRegion, truePairs, minOverlap);
  BL_BENCH_COLLECTIVE_END(cmpr, "true_pairs", truePairs.size(), comm);

  BL_BENCH_START(cmpr);
  eliminateDuplicates(comm, truePairs);
  BL_BENCH_COLLECTIVE_END(cmpr, "unique_pairs", truePairs.size(), comm);
  BL_BENCH_REPORT_MPI_NAMED(cmpr, "cmpr_app", comm);
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
        float fpRatio, fnRatio;
        fnRatio = (float) inTruth;
        fnRatio /= (float)(inTruth + totalIntersect);
        fpRatio = (float) inCandidates;
        fpRatio /= (float)(totalIntersect);
        std::cout << "|CANDIDATES \\   TRUTH|  : " << inCandidates;
        std::cout << " |TRUTH \\   CANDIDATES|  : " << inTruth;
        std::cout << " |TRUTH \\cup CANDIDATES|  : " << totalIntersect;
        std::cout << " FN Ratio : "  << fnRatio;
        std::cout << " FP Ratio : "  << fpRatio << std::endl;
    }
    BL_BENCH_COLLECTIVE_END(cmpr, "count_setdiff", truePairs.size(), comm);
    BL_BENCH_REPORT_MPI_NAMED(cmpr, "cmpr_app", comm);
}

template<typename T>
void compareOverLaps(const mxx::comm& comm,
                     std::string positionFile,
                     unsigned readLength,
                     std::vector<std::pair<T, T>>& candidatePairs,
                     uint32_t minOverlap){

    // generate true pairs
    std::vector<std::pair<T, T>> truePairs;
    generateTruePairs(comm, positionFile, readLength, minOverlap, truePairs);

    auto totalTruePairs = mxx::allreduce(truePairs.size(), comm);
    if(comm.rank() == 0)
       std::cout << "TRUE TOTAL : " << totalTruePairs << std::endl;

    //
    computeSetDifference(comm, candidatePairs, truePairs);
}


template<typename T>
void loadSimulatedPosFile(const mxx::comm& comm,
                          std::string positionFile,
                          std::vector<PosReadPair<T>>& readPairs){

  if(comm.rank() == 0)
     std::cout << "---- Start Loading SIM POS File : [" << positionFile << "] ----" << std::endl;

  // compute start and end offsets corresponding this rank
  T offsetStart, offsetEnd;
  compute_offsets(comm, positionFile, offsetStart, offsetEnd);

  // load the block
  std::vector<std::string> bufferStore ;
  read_block(comm, positionFile, offsetStart, offsetEnd, bufferStore);
  auto totalPosLines = mxx::allreduce(bufferStore.size(), comm);
  if(comm.rank() == 0)
      std::cout << "POS FILE RECORDS : " << totalPosLines << std::endl;
  if(totalPosLines == 0)
      return;
  auto vx = mxx::left_shift(bufferStore.front(), comm);

  // generate the pairs
  // input position file has format for paired end reads:
  // position_left_end position_right_end
  readPairs.resize(2 * bufferStore.size());
  auto startIdx = mxx::scan(bufferStore.size(), comm);
  auto nRecords = mxx::allreduce(bufferStore.size(), comm);
  if(bufferStore.size() > 0)
    startIdx = startIdx - bufferStore.size();
  auto rpItr = readPairs.begin();
  auto nReads = startIdx - startIdx;
  for(auto rcd : bufferStore){
    uint64_t inValue;
    std::stringstream strStream(rcd);
    strStream >> inValue;
    rpItr->first = rpItr->second = startIdx;
    rpItr->first_position = rpItr->second_position = inValue;
    rpItr->first_length = rpItr->second_length = 0;
    rpItr++;
    strStream >> inValue;
    rpItr->first = rpItr->second = startIdx + nRecords;
    rpItr->first_position = rpItr->second_position = inValue;
    rpItr->first_length = rpItr->second_length = 0;
    rpItr++; startIdx++;
    nReads += 2;
  }

  totalPosLines = mxx::allreduce(nReads, comm);
  if(comm.rank() == 0)
     std::cout << "POS DATA : " << totalPosLines << std::endl;
  if(comm.rank() == 0)
     std::cout << "---- Finish Loading SIM POS File ----" << std::endl;
}


template<typename T>
void loadSAMPosFile(const mxx::comm& comm,
                    std::string positionFile,
                    std::vector<PosReadPair<T>>& readPairs){

    if(comm.rank() == 0)
       std::cout << "---- Start Loading SAM POS File: [" << positionFile <<  "] ----" << std::endl;

    // compute start and end offsets corresponding this rank
    T offsetStart, offsetEnd;
    compute_offsets(comm, positionFile, offsetStart, offsetEnd);

    // load the block
    std::vector<std::string> bufferStore ;
    read_block(comm, positionFile, offsetStart, offsetEnd, bufferStore);
    auto totalPosLines = mxx::allreduce(bufferStore.size(), comm);
    if(comm.rank() == 0)
        std::cout << "POS FILE RECORDS : " << totalPosLines << std::endl;
    if(totalPosLines == 0)
        return;
    auto vx = mxx::left_shift(bufferStore.front(), comm);
    // generate the pairs
    // input position file has format for paired end reads:
    // position_left_end position_right_end
    readPairs.resize(bufferStore.size());
    std::size_t rix = 0;
    for(auto rcd : bufferStore){
        T readIdx;
        uint64_t inPos;
        uint32_t inLen;
        std::stringstream strStream(rcd);
        strStream >> readIdx;
        strStream >> inPos;
        strStream >> inLen;
        readPairs[rix].first =  readPairs[rix].second = readIdx;
        readPairs[rix].first_position =  readPairs[rix].second_position = inLen;
        readPairs[rix].first_length =  readPairs[rix].second_length = inLen;
        rix++;
    }

    auto totalRecords = mxx::allreduce(rix, comm);
    if(comm.rank() == 0)
        std::cout << "POS DATA : " << totalRecords << std::endl;
    if(comm.rank() == 0)
       std::cout << "---- Finish Loading SAM POS File ----" << std::endl;
}

template<typename T>
void loadPositionFile(const mxx::comm& comm,
                      std::string positionFile,
                      std::vector<PosReadPair<T>>& readPairs){
    if(has_suffix(positionFile, "map")){
        loadSAMPosFile(comm, positionFile, readPairs);
    } else {
        loadSimulatedPosFile(comm, positionFile, readPairs);
    }
}


template<typename T>
void updateReadPL(const mxx::comm& comm,
                  std::vector<PosReadPair<T>>& readPairs){

    auto fst_cmp_fn = [&] (const PosReadPair<T>& x, const PosReadPair<T>& y){
        return ((x.first < y.first) ? true : (x.first_position < y.first_position));
    };

    auto snd_cmp_fn = [&] (const PosReadPair<T>& x, const PosReadPair<T>& y){
        return (x.second < y.second) ? true : (x.second_position < y.second_position);
    };

    auto val_fn = [&](const PosReadPair<T>& x, const PosReadPair<T>& y){
        return y.is_pos_bad() ? x : y;
    };


    if(!mxx::is_sorted(readPairs.begin(), readPairs.end(),
                       fst_cmp_fn, comm))
        mxx::sort(readPairs.begin(), readPairs.end(), fst_cmp_fn, comm);

    PosReadPair<T> fxp;
    for(auto rpitr = readPairs.rbegin(); rpitr != readPairs.rend(); rpitr++)
        if(rpitr->first_position != RPOS_MISSING){
            fxp = *rpitr;
            break;
        }

    auto fyp = mxx::scan(fxp, val_fn, comm);;

    // update first
    for(auto rpx : readPairs){
        if(rpx.first_position == RPOS_MISSING){
            rpx.first_position = fyp.first_position;
            rpx.first_length = fyp.first_length;
        } else {
            fyp = rpx;
        }
    }

    if(!mxx::is_sorted(readPairs.begin(), readPairs.end(),
                       snd_cmp_fn, comm))
        mxx::sort(readPairs.begin(), readPairs.end(), snd_cmp_fn, comm);

    PosReadPair<T> sxp;
    for(auto rpitr = readPairs.rbegin(); rpitr != readPairs.rend(); rpitr++)
        if(rpitr->first_position != RPOS_MISSING){
            sxp = *rpitr;
            break;
        }
    auto syp = mxx::scan(sxp, val_fn, comm);;

    // update second
    for(auto rpx : readPairs){
        if(rpx.second_position == RPOS_MISSING){
            rpx.second_position = syp.second_position;
            rpx.second_length = syp.second_length;
        } else {
            syp = rpx;
        }
    }

}

template<typename T>
void computeSetDiff(const mxx::comm& comm,
                        std::vector<PosReadPair<T>>& readPairs,
                        uint32_t threshold){
    // ABANDON THIS
    // REMOVE THIS LATER
}

template<typename T>
void comparePosOverLaps(const mxx::comm& comm,
                        std::string positionFile,
                        std::vector<std::pair<T, T>>& candidatePairs,
                        uint32_t threshold) {

    std::vector<PosReadPair<T>> readPairs;
    // generate true pairs
    loadPositionFile(comm, positionFile, readPairs);

    auto nsize = readPairs.size();
    readPairs.resize(nsize + candidatePairs.size());
    for(auto idx = nsize, jdx = (nsize - nsize); idx < readPairs.size(); idx++, jdx++){
        readPairs[idx].first = candidatePairs[jdx].first;
        readPairs[idx].second = candidatePairs[jdx].second;
    }

    auto totalTruePairs = mxx::allreduce(readPairs.size(), comm);
    if(comm.rank() == 0)
        std::cout << "ALL TOTAL : " << totalTruePairs << std::endl;

    updateReadPL(comm, readPairs);

    computeSetDiff(comm, readPairs, threshold);
}


#endif
