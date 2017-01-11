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

void timeWT(sdsl::bit_vector& input, int rounds, char* outFile, int check) {
  sdsl::rank_support_v<1> rs;
  sdsl::util::init_support(rs, &input);
  for (int i=0; i < rounds; i++) {
    startTime();
    sdsl::util::init_support(rs, &input);
    nextTimeN();
  }
  cout<<"Peak-memory: " << getPeakRSS() / (1024*1024)<< endl;

  if(check) {
	  std::cout<<"Checking...";
	  //TODO 
  }
  if (outFile != NULL) {
	sdsl::store_to_file(rs, outFile);
  }
}

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);
  bool binary = P.getOption("-b");
  int check = P.getOptionIntValue("-c",0);
  sdsl::bit_vector input;
  sdsl::load_from_file(input, iFile);
  timeWT(input, rounds, oFile, check);
}
