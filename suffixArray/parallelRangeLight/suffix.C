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

// SDSL includes
#include "../../select/select_support_mcl_par.hpp"
//#include <sdsl/select_support_mcl.hpp>
#include "../../rank/rank_support_v_par.hpp"
//#include <sdsl/rank_support_v.hpp>
#include <sdsl/int_vector.hpp>

#include <iostream>
#include <algorithm>
#include "gettime.h"
#include "sequence.h"
#include "blockRadixSort.h"
#include "quickSort.h"
#include "parallel.h"
#include "SA.h"


using namespace std;

//#define printInfo

#ifdef printInfo
#define nextTimeM(_str) nextTime(_str)
#else
#define nextTimeM(_str) 
#endif

typedef unsigned int uint;
typedef unsigned char uchar;


inline uintT grabChars(uchar *s, uint bits, uintT nChars) {
  uintT r = s[0];
  for (uintT i=1; i < nChars; i++) r = r<<bits | s[i];
  return r; 
}

inline uintT grabCharsEnd(uchar *s, uint bits, uintT nChars, uintT end) {
  uintT r = s[0];
  for (uintT i=1; i < nChars; i++) 
    r = r<<bits | ((i < end) ? s[i] : 0);
  return r; 
}

struct cmp_offset {
	uint* SA,*ISA;
	uint n,offset;
	cmp_offset(uint* SA_, uint* ISA_, uint n_, uint offset_ ) : SA(SA_), ISA(ISA_), n(n_), offset(offset_) {}
	// Return rank of the suffix offset characters after the suffix at position pos
	uint operator() (const uint& pos)const {
		return (pos + offset >= n) ? (n-pos-offset) : ISA[pos+offset]; 
	}
	bool operator()(const uint& a, const uint& b)const {
		return (*this)(a) < (*this)(b);
	}
};


inline void setStartSeg(uint i, sdsl::bit_vector& segBounds) {
	segBounds[i] = 1;
}
inline void setEndSeg(uint i, sdsl::bit_vector& segBounds) {
	segBounds[i] = segBounds[i] & 1 ? 0 : 1;	
}
inline uint getStartSeg(uintT i, sdsl::select_support_mcl<1,1>& ss) {
	return ss.select(2*i+1);
}
inline uint getEndSeg(uintT i, sdsl::select_support_mcl<1,1>& ss) {
	return ss.select(2*i+2)+1;
}

const uint BLOCK_SIZE = 128*1024;


void splitSegmentsParallel(sdsl::bit_vector& segBounds, sdsl::select_support_mcl<1,1>& ss, sdsl::rank_support_v<1,1>& rs, uint* ISA, uint* SA, uint n, uint offset, cmp_offset F) {
	uintT nSegs = rs.rank(n) / 2;
	sdsl::bit_vector segBoundsBuf(n,0);
	uint num_blocks = n / BLOCK_SIZE +1;
	uint* names = new uint[num_blocks]; // Keep track of last name of segments spanning blocks
	bool* first_comp = new bool[num_blocks];

	// First round do all comparisons
	parallel_for (uint b = 0; b < num_blocks; b++) {
		uint bs = b * BLOCK_SIZE;
		uint be = std::min(n, (b+1)*BLOCK_SIZE);
		// Get first/last segment in block
		uint sb = rs.rank(bs)/2;
		uint se = (rs.rank(be)+1) / 2;
		names[b] = -1;// -1 means not set
		for (uint seg = sb; seg < se; seg++) { // Go over all segments
			uint start = getStartSeg(seg,ss);
			uint end = std::min(be, getEndSeg(seg,ss));
			if (start >= bs) {
				// By definition first value is 1
				segBoundsBuf[start] = 1;
				names[b] = start;
			} else {  // segments starts in previous block
				start = bs-1;
				first_comp[b] = F(SA[start]) != F(SA[start+1]); // Save comparison with element from previous block
			}
			for (uint i = start + 1; i < end; i++) {
				if (F(SA[i-1]) != F(SA[i])) { segBoundsBuf[i] = 1; names[b] = i;}
			}
		}
	}
	// Pass names across multiple blocks TODO : in paralle
	for (uint b = 1; b < num_blocks; b++) {
		if (names[b] == -1)  names[b] = names[b-1];
	}
	
	// Second round mark beginning end end of new segments, update ISA 
	parallel_for (uint b = 0; b < num_blocks; b++) {
		uint bs = b * BLOCK_SIZE;
		uint be = std::min(n, (b+1)*BLOCK_SIZE);
		// Get first/last segment in block
		uint sb = rs.rank(bs)/2;
		uint se = (rs.rank(be)+1) / 2;
		bool eval;
		uint name, start, end;
		uint tmp_seg_end;
		for (uint seg = sb; seg < se; seg++) {
			start = std::max(bs, getStartSeg(seg,ss));
			tmp_seg_end = getEndSeg(seg,ss);
			end = std::min(be, tmp_seg_end);
			name = b == 0 ? start : names[b-1];
			for (uint i = start; i < end-1; i++) {
				eval = segBoundsBuf[i] & 1;
				// Update name and isa
				if (eval) name = i;
				ISA[SA[i]] = name;
				// Is beginning of zero run
				if (eval && !(segBoundsBuf[i+1]&1)) segBoundsBuf[i] = 1;
				// End of zero run
				else if( !eval && (segBoundsBuf[i+1]&1)) segBoundsBuf[i] = 1;
				else if (eval) segBoundsBuf[i] = 0;
			}
			eval = segBoundsBuf[end-1];
			if (eval) name = end-1;
			ISA[SA[end-1]] = name;
			// Extra case for end
			if (end == tmp_seg_end) { // If really end of segment
				segBoundsBuf[end-1] = eval ? 0 : 1;
			} else { // if segment continues in next block
				if (eval && !first_comp[b+1]) segBoundsBuf[end-1] = 1;
				else if (!eval && first_comp[b+1]) segBoundsBuf[end-1] = 1;
				else if (eval) segBoundsBuf[end-1] = 0;
			}
		}
	}
	swap(segBounds, segBoundsBuf);
	delete [] names;
	delete [] first_comp;
}

void sortSegmentsParallel(uintT nSegs, sdsl::select_support_mcl<1,1>& ss, uint* ISA, uint* SA, uint n, uint offset) {
	parallel_for (uintT i=0; i < nSegs; i++) {
		uint start = getStartSeg(i, ss);
		uint *SAi = SA + start;
		uint l = getEndSeg(i, ss) - start;
		if (l >= 256) 
			intSort::iSort(SAi, l, n , cmp_offset(SA, ISA, n, offset));
		else
			quickSort(SAi,l,cmp_offset(SA, ISA, n, offset));
	}
}

bool brokenCilk(uint *SA, uint offset, uint n, uint* ISA, sdsl::bit_vector& segBounds, sdsl::select_support_mcl<1,1>& ss, sdsl::rank_support_v<1,1>& rs) {
  uint nSegs = rs.rank(n) / 2;
  // Init select/rank support on segBounds
  sortSegmentsParallel(nSegs, ss, ISA, SA, n, offset);
  // Write new seg bounds to buffer
  splitSegmentsParallel(segBounds, ss, rs, ISA, SA, n, offset, cmp_offset(SA, ISA, n, offset));
  // Rebuild rank and select support
  rs = sdsl::rank_support_v<1,1>(&segBounds);
  if (rs.rank(n) == 0) return false;
  ss = sdsl::select_support_mcl<1,1>(&segBounds);
  return true;
}

uintT* suffixArrayInternal(unsigned char* ss, long n) { 
  // following line is used to fool icpc into starting the scheduler
  if (n < 0) cilk_spawn printf("ouch");
  uintT *SA = newA(uintT,n);
  uintT *ISA = newA(uintT,n);
  uchar* s = newA(uchar, n);

  uintT flags[256];
  for (uintT i=0; i < 256; i++) flags[i] = 0;
  parallel_for (uintT i=0; i < n; i++) 
    if (!flags[ss[i]]) flags[ss[i]] = 1;

  // renumber characters densely
  // start at 1 so that end-of-string is 0
  uintT m = sequence::scan(flags,flags,256,utils::addF<uintT>(),(uintT)1);
  parallel_for (uintT i=0; i < n; i++) 
    s[i] = flags[ss[i]];
  #ifdef printInfo
  cout << "m = " << m << endl;
  #endif

  uintT bits = max(1,utils::log2Up(m));
  uintT nchars = 31/bits;

  // pack characters into word in chunks of "bits"
  startTime();
  if(n+1 > nchars) {
    parallel_for (uintT i=0; i < n-nchars+1; i++) {
      ISA[i] = grabChars(s+i,bits,nchars); 
      SA[i]  = i;
    }

    for (uintT i=n-nchars+1; i < n; i++) {
      ISA[i] = grabCharsEnd(s+i,bits,nchars,n-i); 
      SA[i] = i;
    }
  } else {
    for (uintT i=0; i < n; i++) {
      ISA[i] = grabCharsEnd(s+i,bits,nchars,n-i); 
      SA[i] = i;
    }
  }
  free(s);

  //intSort::iSort(C,n,(uintT)1 << bits*nchars,utils::firstF<uintT,uintT>());
  intSort::iSort(SA, n, (uintT)1 << bits*nchars, cmp_offset(SA, ISA, n, 0));
  sdsl::bit_vector segBounds(n,0);
  segBounds[0] = 1; segBounds[n-1] = 1;
  sdsl::select_support_mcl<1,1> selectsupport(&segBounds);
  sdsl::rank_support_v<1,1> ranksupport(&segBounds);
  // Init segBounds
  splitSegmentsParallel(segBounds, selectsupport, ranksupport, ISA, SA, n, 0, cmp_offset(SA, ISA, n, 0)); 
  selectsupport = sdsl::select_support_mcl<1,1>(&segBounds);
  ranksupport = sdsl::rank_support_v<1,1>(&segBounds);

  uintT offset = nchars;
  
  uint round =0;
  while (1) {
    utils::myAssert(round++ < 40, "Suffix Array:  Too many rounds");
    if (!brokenCilk(SA, offset, n, ISA, segBounds, selectsupport, ranksupport)) break;
    offset = 2 * offset;
  }
  free(ISA);
  return SA;
}

uintT* suffixArray(unsigned char* ss, long n) { 
  return (uintT*)suffixArrayInternal(ss, n);
}
