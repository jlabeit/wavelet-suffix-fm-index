// This code is part of the paper "Parallel Wavelet Tree Construction"
// Copyright (c) 2014 Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#include <fstream>
#include <iostream>
#include <algorithm>
#include "gettime.h"
#include "getmemory.h"
#include "parallel.h"
#include "IO.h"
#include "parseCommandLine.h"
#include "WT.h"
#include "sequence.h"

using namespace std;
using namespace benchIO;

//naive linear time rank
intT naiveRank(long* bitmap, long start, long finish, int bit) {
  intT rank = 0; 
  for(long i=start;i<finish;i++) {
    long w = bitmap[i/64];
    if((1 & (w >> (i % 64))) == bit) rank++;
  }
  return rank;
}

//access queries---use for testing only; for performance, need to
//implement a fast rank structure
intT query(pair<WTnode*,long*> R, long index, long sigma, long n) {
  long levels = max(1,utils::log2Up(sigma));
  WTnode* nodes = R.first;
  long* bitmap = R.second;
  WTnode curr = nodes[0]; 
  long min = 0, max = 1 << utils::log2Up(sigma);

  //loop until no children nodes left
  while(max-min > 2 && (curr.leftChild != UINT_T_MAX || curr.rightChild != UINT_T_MAX)) {
    long offset = curr.bitmapPtr;
    long w = bitmap[(index+offset)/64];
    int bit = (w >> ((index+offset) % 64)) & 1;
    if(bit) {
      //subtract number of zeros
      index -= naiveRank(bitmap, offset, offset+index, 0);
      min += (long) 1 << (utils::log2Up(max-min)-1);
      curr = nodes[curr.rightChild];
    }
    else {
      index -= naiveRank(bitmap, offset, offset+index, 1);
      max = min + ((long) 1 << (utils::log2Up(max-min)-1)); 
      curr = nodes[curr.leftChild];
    }
  }
  long offset = curr.bitmapPtr;
  long w = bitmap[(index+offset)/64];
  if(1 & (w >> ((index+offset) % 64))) return min+1;
  else return min;
} 

void timeWT(symbol* s, long n, int rounds, char* outFile, int check) {
  //sigma
  long k = 1 + sequence::reduce(s, n, utils::maxF<symbol>());
  int* A = newA(int,k+1);
  parallel_for(long i=0;i<k+1;i++) A[i] = 0;
  parallel_for(long i=0;i<n;i++) {
    if(!A[s[i]]) A[s[i]] = 1;
  }
  long sigma = sequence::plusScan(A,A,k+1);
  //cout << "n = " << n << endl;
  //cout << "sigma = " << sigma << endl;
  int* reverseMap = newA(int,sigma);
  parallel_for(long i=0;i<k;i++) 
    if(A[i] != A[i+1]) reverseMap[A[i]] = i;

  parallel_for(long i=0;i<n;i++) s[i] = A[s[i]];
  free(A);

  // Version for sdsl
  sdsl::int_vector<sizeof(intT)*8> input(n, 0, sizeof(symbol)*8);
  for (uint64_t i = 0; i < n; i++) {
	input[i] = s[i];
  }
  string input_file = "@input.sdsl"; 
  sdsl::store_to_file(input, input_file);
  sdsl::wt_int<> wt;
  sdsl::construct(wt, input_file);
  sdsl::rank_support_v<1> rs;
  for (int i=0; i < rounds; i++) {
    startTime();
    sdsl::util::init_support(rs, wt.get_m_tree());
    nextTimeN();
  }
  //cout<<"Peak-memory: " << getPeakRSS() / (1024*1024)<< endl;

  if(check) {
	  std::cout<<"Checking...";
	  //TODO 
  }
  if (outFile != NULL) {
	sdsl::store_to_file(rs, outFile);
  }
  // if(outFile != NULL) {
  //   intT* foo = newA(intT,n); 
  //   parallel_for(long i=0;i<n;i++) foo[i] = (intT) s[i];
  //   ofstream out(outFile, ofstream::out | ios::binary);
  //   out.write((char*)foo, sizeof(intT)*n);
  //   free(foo);
  //   out.close();
  // }

  //if (outFile != NULL) writeIntArrayToFile((intT*) R, (intT) n, outFile);
}

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);
  bool binary = P.getOption("-b");
  int check = P.getOptionIntValue("-c",0);
  if(binary) {
    ifstream in(iFile,ifstream::in |ios::binary);
    in.seekg(0,ios::end);
    long n = in.tellg();
    in.seekg(0);
    char* s = newA(char, n);
    in.read(s,n);
    in.close(); 
    timeWT((symbol*)s, n/sizeof(symbol), rounds, oFile, check);
    free(s);
  }
  else {
#ifdef INT
    _seq<uintT> S = readIntArrayFromFile<uintT>(iFile);
    uintT n = S.n;
    timeWT(S.A, n, rounds, oFile, check);
#else
    _seq<char> S = readStringFromFile(iFile);
    uintT n = S.n;
    timeWT((unsigned char*) S.A, n, rounds, oFile, check);
#endif
    S.del();
  }
}
