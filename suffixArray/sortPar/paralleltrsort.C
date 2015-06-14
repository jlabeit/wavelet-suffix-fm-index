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


void splitSegment(seg *segOut, saidx_t start, saidx_t l, saidx_t* ISA, saidx_t* ISA_buf, saidx_t *SA,
		  bool addRanks, cmp_offset F) {
  if (l < 1000) { // sequential version

    if (addRanks) {
      // if following two loops are fused performance goes way down?
      saidx_t name = 0;
      ISA_buf[SA[0]] = name + start + 1;
      for (saidx_t i=1; i < l; i++) {
	if (F(SA[i-1]) != F(SA[i])) name = i; // TODO use comparison function 
	ISA_buf[SA[i]] = name + start + 1;
      }
    }

    saidx_t name = 0;
    for (saidx_t i=1; i < l; i++) {
      if (F(SA[i-1]) != F(SA[i])) {
	segOut[i-1] = seg(name+start,i-name);
	name = i;
      } else segOut[i-1] = seg(0,0);
    }
    segOut[l-1] = seg(name+start,l-name);

  } else { // parallel version
    uintT *names = newA(uintT,l);

    parallel_for (uintT i = 1;  i < l;  i++) 
	names[i] = (F(SA[i]) != F(SA[i-1])) ? i : 0;
    names[0] = 0;
    sequence::scanI(names,names,l,utils::maxF<uintT>(),(uintT)0);
    if (addRanks) 
      parallel_for (uintT i = 0;  i < l;  i++) 
	ISA_buf[SA[i]] = names[i]+start+1;

    parallel_for (uintT i = 1;  i < l;  i++)
      if (names[i] == i) 
	segOut[i-1] = seg(start+names[i-1],i-names[i-1]);
      else segOut[i-1] = seg(0,0);
    segOut[l-1] = seg(start+names[l-1],l-names[l-1]);

    free(names);
  }
}  



void brokenCilk(uintT nSegs, seg *segments, saidx_t *SA, saidx_t offset, saidx_t n, saidx_t* ISA,saidx_t* ISA_buf, seg *segOuts, saidx_t* offsets) {
  parallel_for (uintT i=0; i < nSegs; i++) {
    saidx_t start = segments[i].start;
    saidx_t *SAi = SA + start;
    saidx_t l = segments[i].length;
    if (l >= 256) 
      intSort::iSort(SAi, l, n , cmp_offset(SA, ISA, n, offset));
    else
      quickSort(SAi,l,cmp_offset(SA, ISA, n, offset));
  }

  parallel_for (uintT i=0; i < nSegs; i++) {
    uintT start = segments[i].start;
    splitSegment(segOuts + offsets[i], start, segments[i].length, 
		 ISA, ISA_buf, SA + start, 1, cmp_offset(SA, ISA, n, offset));
  }
  parallel_for(saidx_t i = 0; i < n; i++) ISA[i] = ISA_buf[i];
}

void paralleltrsort(saidx_t* ISA, saidx_t* SA, saidx_t n, saidx_t* buf, saidx_t buffer_len) {
	saidx_t offset = 1;
	seg *segOuts = newA(seg,n); 
  	seg *segments= newA(seg,n/2);
	saidx_t *offsets = newA(saidx_t,n/2);
	saidx_t *ISA_buf = (buffer_len >= n) ? buf : newA(saidx_t,n);
	splitSegment(segOuts, 0, n, ISA, ISA_buf, SA, 1, cmp_offset(SA, ISA, n, 0));
	parallel_for(saidx_t i = 0; i < n; i++) ISA[i] = ISA_buf[i];
	saidx_t nKeys = n;
	while (true) {
		utils::myAssert(offset <= n,  "Suffix Array:  Too many rounds");
		// Calculate segments 
		saidx_t nSegs = sequence::filter(segOuts,segments,nKeys,isSeg());
		if (nSegs == 0) break;

		parallel_for (uintT i=0; i < nSegs; i++)
			offsets[i] = segments[i].length;

		nKeys = sequence::scan(offsets,offsets,nSegs,utils::addF<uintT>(),(saidx_t)0);

		// Sort segments and update ranks
		brokenCilk(nSegs, segments, SA, offset, n, ISA, ISA_buf, segOuts, offsets);
		// double offset		
	    	offset = 2 * offset;			
	}
	// Use 0 based index again
	parallel_for (saidx_t i = 0; i < n; i++) {
		ISA[i]--;
	}
	free(segOuts);
	free(segments);
	free(offsets);
	if (buffer_len < n) free(ISA_buf);
}
