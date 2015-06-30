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
#include "parallel.h"
#include "utils.h"
#include <iostream>
#include "sequence.h"
#include "WT.h"
#include "blockRadixSort.h"
#include "gettime.h"
using namespace std;

#define THRESHOLD 10000

template <class E>
struct topL {
  uint shift;
  topL(uint _s) { shift = _s; }
  E operator() (E a) {return a >> shift;} ; 
};

template <class E>
struct topLCmp {
  uint shift;
  topLCmp(uint _s) { shift = _s; }
  E operator() (E a, E b) {return (a >> shift) < (b >> shift);} ; 
};

inline void writeOr(long *a, long b) {
  volatile long newV, oldV; 
  do {oldV = *a; newV = oldV | b;}
  while ((oldV != newV) && !utils::CAS(a, oldV, newV));
}

struct notMax { bool operator() (uintT i) {return i != UINT_T_MAX;}};

pair<WTnode*,long*> WT(symbol* s, uintT n, uintT sigma) {
  int levels = max(1,utils::log2Up(sigma));
  //cout << "levels = " << levels << endl;

  uintT numNodes = (long)1 << (1+levels);
  WTnode* nodes = newA(WTnode,numNodes);
  nodes[0].length = n; nodes[0].bitmapPtr = 0;
#ifdef POINTERS
  nodes[0].parent = UINT_T_MAX;
  parallel_for(long i=1;i<numNodes;i++) nodes[i].parent = UINT_T_MAX-1; //indicates not in tree yet
#endif
  long* wt = newA(long,((long)n*levels+63)/64);
  parallel_for(long i=0;i<((long)n*levels+63)/64; i++) wt[i] = 0;
  //cout << "wt length = " << (n*levels+63)/64 << endl;;
  //cout << "sigma = " << sigma << endl;

  uintT* space = newA(uintT,n); //temp space
  uintT* space2 = newA(uintT,n); //temp space

  symbol* s2 = newA(symbol,n); //space for reordered character string
  parallel_for(long i=0;i<n;i++) s2[i] = s[i];

  for(int l = 0; l < levels; l++) {
    //cout << "l = " << l << endl;
    int mask = (long)1 << (levels - l - 1);

    if(l == 0) {
      uintT endWord = (n-1)/64;
      parallel_for(uintT k=0;k<=endWord;k++) {
	uintT b = 64*k;
	for(uintT i=b;i<min(b+64,n);i++) 
	  if(s2[i] & mask) wt[k] |= (long)1 << i % 64;
      }
    } else {
      intOffset wtOffset = l*n;

      uintT numNodesLevel = (long)1 << l;
      uint shift = levels-l;
      
      //stably sort on top l bits to put characters in correct position
      //compSort(s2,n,topLCmp<symbol>(shift));
      //JJ:why sort after top shift symbols, shouldnt sort after one bit be enough in msortWT vs. sortWT
      intSort::iSort(s2, n, numNodesLevel, topL<symbol>(shift));

      //create bitmaps for this level
      intOffset s = wtOffset, e = wtOffset+n-1, word = s/64;
      while(s % 64) {
	if(s2[s-wtOffset] & mask) wt[word] |= (long)1 << s % 64;	
	s++;
      } 
      intOffset startWord = s/64, endWord = e/64;
      parallel_for(intOffset k=startWord;k<=endWord;k++) {
	intOffset b = 64*k;
	for(intOffset i=b;i<min(wtOffset+n,b+64);i++) 
	  if(s2[i-wtOffset] & mask) wt[k] |= (long)1 << i % 64;
      }

      //determine offsets and compute number of elements per node
      //JJ:calc offsets and number of elemnts per node without using n*log(n) space
      space[0] = 0;
      parallel_for(long i=1;i<n;i++) {
	space[i] = ((s2[i] >> shift) != (s2[i-1] >> shift)) ? i : UINT_T_MAX; 
      }

      uintT actualNumNodes = sequence::filter(space,space2,n,notMax());
      //cout << "numNodesLevel = " << numNodesLevel << " actualNumNodes = " << actualNumNodes << endl;
      if(actualNumNodes < THRESHOLD) {
	for(long i=0;i<actualNumNodes;i++) {
	  uintT nodeID = numNodesLevel-1+(s2[space2[i]] >> shift);
	  uintT length = (i == actualNumNodes-1) ? (n-space2[i]) : (space2[i+1]-space2[i]);
	  nodes[nodeID].length = length;
	  nodes[nodeID].bitmapPtr = wtOffset+space2[i];
#ifdef POINTERS
	  nodes[nodeID].parent = (nodeID-1)/2; //put node in tree
#endif
	}
      } 
      else {
	parallel_for(long i=0;i<actualNumNodes;i++) {
	  uintT nodeID = numNodesLevel-1+(s2[space2[i]] >> shift);
	  uintT length = (i == actualNumNodes-1) ? (n-space2[i]) : (space2[i+1]-space2[i]);
	  nodes[nodeID].length = length;
	  nodes[nodeID].bitmapPtr = wtOffset+space2[i];
#ifdef POINTERS
	  nodes[nodeID].parent = (nodeID-1)/2; //put node in tree
#endif
	}
      }
    }
    //cout << "finishing level " << l << endl;
  }

#ifdef POINTERS
  //fix child pointers of tree
  parallel_for(long i=0;i<numNodes/2-1;i++) {
    if(nodes[2*i+1].parent == i) nodes[i].leftChild = 2*i+1; 
    else nodes[i].leftChild = UINT_T_MAX;
    if(nodes[2*i+2].parent == i) nodes[i].rightChild = 2*i+2; 
    else nodes[i].rightChild = UINT_T_MAX;
  }
#endif
  
  free(space); free(space2); 
  free(s2); 
  return make_pair(nodes,wt);
}
