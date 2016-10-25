#include <cstdint>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

#include "mxx/comm.hpp"

#include "io_utils.hpp"
#include "run_cfg.hpp"

// UTILITY FUNCTIONS -------------------------------------------------
// trim taken from stack overflow
// http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
static inline std::string &ltrim(std::string &s) {
  s.erase(s.begin(),
          std::find_if(s.begin(), s.end(),
                       std::not1(std::ptr_fun<int, int>(std::isspace))));
  return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       std::not1(std::ptr_fun<int, int>(std::isspace))).base(),
          s.end());
  return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
  return ltrim(rtrim(s));
}
void updateStrStore(const std::string& in_str,
                    std::vector<char>& StrStore,
                    std::vector<uint64_t>& Offset,
                    unsigned long& position){
    StrStore.resize(StrStore.size() + in_str.length() + 1);
    memcpy(&StrStore[position], in_str.c_str(), in_str.length() + 1);
    Offset.push_back(position);
    position += in_str.length() + 1;
}

void compute_offsets(const mxx::comm& comm,
                     std::string inFileName,
                     uint64_t& offsetStart,
                     uint64_t& offsetEnd){
  // offset computation for static work load
  std::ifstream fin(inFileName.c_str());
  fin.seekg(0,std::ios::end);
  auto fileSize = fin.tellg();
  //uint64_t tmp = (fileSize/((uint64_t)comm.size())) * (((uint64_t)comm.rank()) + 1);
  uint64_t tmp = block_high(comm.rank(), comm.size(), fileSize);
  fin.seekg(tmp, std::ios::beg);
  if(comm.rank() == (comm.size() - 1)){
    tmp = fileSize;
  } else {
    std::string rd_line;
    std::getline(fin, rd_line);
    tmp = fin.tellg();
    //tmp += 1;
  }
  offsetEnd = tmp;
  uint64_t tmp2 = mxx::right_shift(tmp, comm);
  if(comm.rank() == 0)
    offsetStart = 0;
  else
    offsetStart = tmp2;
  //offsetStart = block_low(comm.size(), comm.rank(), fileSize);
  //offsetEnd = block_high(comm.size(), comm.rank() + 1, fileSize);
  
}

std::size_t get_file_size(std::string inFileName){
  std::ifstream fin(inFileName.c_str());
  fin.seekg(0,std::ios::end);
  auto fileSize = fin.tellg();
  return fileSize;
}


void read_block(const mxx::comm& comm,
                std::string inFileName,
                uint64_t offsetStart,
                uint64_t offsetEnd,
                std::vector<std::string>& readStore){

  std::ifstream in_stream(inFileName.c_str());
  in_stream.seekg(offsetStart,std::ios::beg);

  if(comm.rank() > 0){
      std::string rd_line;
      std::getline(in_stream, rd_line);
  }
  
  while(in_stream.good()){
    std::string rd_line;
    std::getline(in_stream, rd_line);
    std::string trimValue = trim(rd_line);
    if(trimValue.length() > 0)
      readStore.push_back(trimValue);
    if(!in_stream.good() ||
       in_stream.tellg() > (std::streamoff) offsetEnd)
      break;
  }
  // auto totalPosLines = mxx::allreduce(readStore.size(), comm);
  // if(comm.rank() == 0)
  //    std::cout << "RDBLOCK TOTAL : " << totalPosLines << std::endl;
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
