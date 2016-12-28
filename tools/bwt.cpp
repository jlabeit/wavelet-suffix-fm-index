#include <iostream>
#include "parseCommandLine.h"
#include "sequence.h"
#include "IO.h"

#include "divsufsort.h"

using namespace benchIO;
using namespace std;

void bwt(unsigned char* s, int64_t n, char* outFile) {
  int64_t* R;
  R = (int64_t*)malloc(sizeof(int64_t)*n);
  divsufsort(s, R, n);
  unsigned char* BWT = (unsigned char*)malloc(sizeof(unsigned char)*n);
  int64_t to_add[2] = {(int64_t)-1,(int64_t)n-1};
  parallel_for(int64_t i = 0; i < n; ++i) {
      BWT[i] = s[R[i]+to_add[R[i]==0]];
  }
  if (outFile != NULL) writeStringToFile((char*) BWT, n, outFile);
  free(R);
}

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<inFile> <outFile>");
  char* iFile = P.getArgument(1);
  char* oFile = P.getArgument(0);
  std::cout << iFile << endl;
  std::cout << oFile << endl;
  _seq<char> S = readStringFromFile(iFile);
  bwt((unsigned char*) S.A, S.n, oFile);
}
