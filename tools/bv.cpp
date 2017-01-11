#include <iostream>
#include "parseCommandLine.h"
#include "sequence.h"
#include "IO.h"
#include <sdsl/suffix_trees.hpp>
#include <sdsl/wt_blcd.hpp>


using namespace benchIO;
using namespace std;

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<inFile> <outFile>");
  char* iFile = P.getArgument(1);
  char* oFile = P.getArgument(0);
  _seq<char> S = readStringFromFile(iFile);
  sdsl::int_vector<8> data(S.n);
  string input_file = "@input.sdsl";
  for (size_t i = 0; i < S.n; ++i) {
    data[i] = S.A[i];
  }
  sdsl::store_to_file(data, input_file);
  sdsl::wt_blcd<> wt;
  sdsl::construct(wt, input_file);
  sdsl::store_to_file(wt.m_bv, oFile);
  std::cout << oFile << std::endl;
}
