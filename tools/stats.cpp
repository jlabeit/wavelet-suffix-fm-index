#include <iostream>
#include "parseCommandLine.h"
#include "sequence.h"
#include "IO.h"
#include <iomanip>

#include "divsufsort.h"

using namespace benchIO;
using namespace std;

void printStats(unsigned char* s, int64_t n) {
  int64_t* SA = (int64_t*)malloc(sizeof(int64_t)*n);
  divsufsort(s, SA, n);
  int64_t* ISA = (int64_t*)malloc(sizeof(int64_t)*n);
  parallel_for(int64_t i=0; i < n; i++)
      ISA[SA[i]]=i;
  uint64_t sum = 0;
  int64_t k = 0;
  for(int64_t i = 0; i < n; i++, k?k--:0) {
      if (ISA[i] == n-1) {
          k=0;
          continue;
      }
      int j = SA[ISA[i]+1];
      while( i + k < n && j + k < n && s[i+k] == s[j+k])
          k++;
      //LCP[ISA[i]] = k;
      sum += k;
  }
  cout << setprecision(1) << fixed << (double)sum / (double)n << endl;
  free(SA);
  free(ISA);
}

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<inFile>");
  char* iFile = P.getArgument(0);
  std::cout << iFile << endl;
  _seq<char> S = readStringFromFile(iFile);
  printStats((unsigned char*) S.A, S.n);
}
