// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
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

#include "config.h"
#include "divsufsort_private.h"
#include "parallel.h"
#include <sdsl/select_support_mcl.hpp>
#include <sdsl/int_vector.hpp>

#include <iostream>
#include <algorithm>
#include "gettime.h"
#include "sequence.h"
#include "blockRadixSort.h"
#include "quickSort.h"

using namespace std;


//#define printInfo

#ifdef printInfo
#define nextTimeM(_str) nextTime(_str)
#else
#define nextTimeM(_str) 
#endif


struct seg {
  saidx_t start;
  saidx_t length;
  seg(saidx_t s, saidx_t l) : start(s), length(l) {}
};

struct isSeg {bool operator() (seg s) {return s.length > 1;}};

struct cmp_offset {
	saidx_t* SA,*ISA;
	saidx_t n,offset;
	cmp_offset(saidx_t* SA_, saidx_t* ISA_, saidx_t n_, saidx_t offset_ ) : SA(SA_), ISA(ISA_), n(n_), offset(offset_) {}
	// Return rank of the suffix offset characters after the suffix at position pos
	saidx_t operator() (const saidx_t& pos)const {
		return (pos + offset >= n) ? (n-pos-offset) : ISA[pos+offset]; 
	}
	bool operator()(const saidx_t& a, const saidx_t& b)const {
		return (*this)(a) < (*this)(b);
	}
};


void splitSegment(saidx_t start, saidx_t l, saidx_t* ISA, saidx_t* ISA_buf, saidx_t *SA,
		  bool addRanks, cmp_offset F, sdsl::bit_vector& segBounds) {
  if (l < 1000) { // sequential version

    if (addRanks) {
      // if following two loops are fused performance goes way down?
      saidx_t name = start;
      ISA_buf[SA[start]] = name + 1;
      for (saidx_t i=start+1; i < start + l; i++) {
	if (F(SA[i-1]) != F(SA[i])) name = i; 
	ISA_buf[SA[i]] = name + 1;
      }
    }

    saidx_t name = start;
    for (saidx_t i=start+1; i < start + l; i++) {
      if (F(SA[i-1]) != F(SA[i])) {
	segBounds[i-1] = 1;
	name = i;
      } else segBounds[i-1] = 0;
    }
    segBounds[start+l-1] = 1;

  } else { // parallel version
    uintT *names = newA(uintT,l);

    parallel_for (uintT i = start+1;  i < l;  i++) 
	names[i] = (F(SA[i]) != F(SA[i-1])) ? (i-start) : 0;
    names[0] = 0;
    sequence::scanI(names,names,l,utils::maxF<uintT>(),(uintT)0);
    if (addRanks) 
      parallel_for (uintT i = start;  i < start + l;  i++) 
	ISA_buf[SA[i]] = names[i-start]+start+1;

    parallel_for (uintT i =1;  i < l;  i++)
      if (names[i] == i) 
	segBounds[i-1 + start] = 1;
      else segBounds[i-1 + start] = 0;
    segBounds[l-1+start] = 1;

    free(names);
  }
}  

uintT getStartSeg(uintT i, sdsl::select_support_mcl<1,1>& ss) {
	return ss.select(2*i+1);
}
uintT getEndSeg(uintT i, sdsl::select_support_mcl<1,1>& ss) {
	return ss.select(2*i+2)+1;
}

bool brokenCilk(saidx_t *SA, saidx_t offset, saidx_t n, saidx_t* ISA, saidx_t* ISA_buf, sdsl::bit_vector& segBounds, sdsl::bit_vector& segBoundsBuf) {
  // TODO do this in parallel and fast
  uintT nSegs = 0;
  for (uintT i = 0; i < n; i++) {
	if (segBounds[i] & 1) nSegs++;
  }
  nSegs /= 2;
  if(nSegs == 0) return false;
  // Init select/rank support on segBounds
  sdsl::select_support_mcl<1,1> ss(&segBounds);
  parallel_for (uintT i=0; i < nSegs; i++) {
    saidx_t start = getStartSeg(i, ss);
    saidx_t *SAi = SA + start;
    saidx_t l = getEndSeg(i, ss) - start;
    if (l >= 256) 
      intSort::iSort(SAi, l, n , cmp_offset(SA, ISA, n, offset));
    else
      quickSort(SAi,l,cmp_offset(SA, ISA, n, offset));
  }
  // Write new seg bounds to buffer
  memset(segBoundsBuf.data(), 0, segBoundsBuf.bit_size() / 8); // Fill with zeros

  parallel_for (uintT i=0; i < nSegs; i++) {
    uintT start = getStartSeg(i, ss);
    uintT length = getEndSeg(i, ss)-start;
    splitSegment(start, length, 
		 ISA, ISA_buf, SA, 1, cmp_offset(SA, ISA, n, offset), segBoundsBuf);
  }
  parallel_for(saidx_t i = 0; i < n; i++) ISA[i] = ISA_buf[i];

  swap(segBounds, segBoundsBuf);
  return true;
}

void paralleltrsort(saidx_t* ISA, saidx_t* SA, saidx_t n, saidx_t* buf, saidx_t buffer_len) {
	saidx_t offset = 1;
	saidx_t *ISA_buf = (buffer_len >= n) ? buf : newA(saidx_t,n);
	parallel_for(saidx_t i = 0; i < n; i++) ISA_buf[i] = ++ISA[i]; // Make 1 based index
	sdsl::bit_vector segBounds(n,0);
	sdsl::bit_vector segBoundsBuf(n,0); 
	// Init segBounds
	splitSegment(0, n, ISA, ISA_buf, SA, 0, cmp_offset(SA, ISA, n, 0), segBoundsBuf);
	while (true) {
		utils::myAssert(offset <= n,  "Suffix Array:  Too many rounds");
		if (!brokenCilk(SA, offset, n, ISA, ISA_buf, segBounds, segBoundsBuf)) break;
		// double offset		
	    	offset = 2 * offset;			
	}
	// Use 0 based index again
	parallel_for (saidx_t i = 0; i < n; i++) {
		ISA[i]--;
	}
	if (buffer_len < n) free(ISA_buf);
}
