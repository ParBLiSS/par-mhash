#include <cstdint>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

#include "mxx/comm.hpp"


#include "io_utils.hpp"

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
  offsetStart = (fileSize/comm.size()) * comm.rank();
  offsetEnd = (fileSize/comm.size()) * (comm.rank() + 1);
}


void read_block(const mxx::comm& comm,
                std::string inFileName,
                uint64_t offsetStart,
                uint64_t offsetEnd,
                std::vector<std::string>& readStore){


  std::ifstream in_stream(inFileName.c_str());
  in_stream.seekg(offsetStart,std::ios::beg);

  // ignore first line, will be read by rank+1
  if(comm.rank() > 0 && in_stream.good()){
    std::string rd_line;
    std::getline(in_stream, rd_line);
    // right on the new line character, then read another line
    if(trim(rd_line).length() == 0 && in_stream.good())
      std::getline(in_stream, rd_line);
  }

  while(in_stream.good()){
    std::string rd_line;
    std::getline(in_stream, rd_line);
    readStore.push_back(trim(rd_line));
    if(!in_stream.good() ||
       in_stream.tellg() > (std::streamoff) offsetEnd)
      break;
  }
}
