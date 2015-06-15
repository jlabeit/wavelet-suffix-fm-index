/*
 * divsufsort.c for libdivsufsort
 * Copyright (c) 2003-2008 Yuta Mori All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "config.h"
#include "gettime.h"
#include "divsufsort_private.h"
#include "parallel.h"
#include "sequence.h"
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support_v.hpp>
#ifdef _OPENMP
# include <omp.h>
#endif
/*- Private Functions -*/
void calculateBucketOffsets(saidx_t * bucket_A, saidx_t* bucket_B) {
/*
note:
  A type B* suffix is lexicographically smaller than a type B suffix that
  begins with the same first two characters.
*/
  saidx_t i,j,t;
  saint_t c0,c1;
  /* Calculate the index of start/end point of each bucket. */
  for(c0 = 0, i = 0, j = 0; c0 < ALPHABET_SIZE; ++c0) {
    t = i + BUCKET_A(c0);
    BUCKET_A(c0) = i + j; /* start point */
    i = t + BUCKET_B(c0, c0);
    for(c1 = c0 + 1; c1 < ALPHABET_SIZE; ++c1) {
      j += BUCKET_BSTAR(c0, c1);
      BUCKET_BSTAR(c0, c1) = j; /* end point */
      i += BUCKET_B(c0, c1);
    }
  }
}


void initBStarBuckets(const sauchar_t *T, saidx_t *SA, saidx_t* bucket_B, saidx_t n, saidx_t m, saidx_t* PAb) {
    saidx_t num_blocks = getWorkers();
    saidx_t block_size = (m-1) / num_blocks + 1;    
    saidx_t* block_bucket_cnt = new saidx_t[num_blocks*BUCKET_B_SIZE];
    memset(block_bucket_cnt, 0, sizeof(saidx_t)*num_blocks*BUCKET_B_SIZE); // TODO is this faster than parallel memset?
    // First pass count buckets for each block
    parallel_for (saidx_t b = num_blocks-1; 0 <= b; b--) {
	saidx_t start = std::min((m-1), (b+1)*block_size);
	saidx_t end = b*block_size;
	saidx_t *bucket_B = block_bucket_cnt + b * BUCKET_B_SIZE; 
	saidx_t t;
	sauchar_t c0,c1;
	for (saidx_t i = start-1; end <= i; --i) {
		t = PAb[i], c0 = T[t], c1 = T[t+1];
		BUCKET_BSTAR(c0,c1)--;	
	}
    }
    // Prefix sum
    parallel_for (saidx_t i = 0; i < BUCKET_B_SIZE; i++) {
	saidx_t sum = bucket_B[i];	
	for (saidx_t b = num_blocks-1; 0 <= b; b--) {
		sum += block_bucket_cnt[b*BUCKET_B_SIZE + i];
		block_bucket_cnt[b*BUCKET_B_SIZE + i] = sum - block_bucket_cnt[b*BUCKET_B_SIZE + i];
	}
    }
    // Second pass to fill the actual buckets
    parallel_for (saidx_t b = num_blocks-1; 0 <= b; b--) {
	saidx_t start = std::min(m-1, (b+1)*block_size);
	saidx_t end = b*block_size;
	saidx_t *bucket_B = block_bucket_cnt + b*BUCKET_B_SIZE;
	saidx_t t;
	sauchar_t c0,c1;
	for (saidx_t i = start -1; end <= i; --i) {
		t = PAb[i], c0 = T[t], c1 = T[t+1];
		SA[--BUCKET_BSTAR(c0,c1)] = i;
	}
    }
    // Init again correctly
    parallel_for (saidx_t i = 0; i < BUCKET_B_SIZE; i++)  bucket_B[i] = block_bucket_cnt[i];
    delete [] block_bucket_cnt;
    saidx_t t;
    sauchar_t c0,c1;
    t = PAb[m - 1], c0 = T[t], c1 = T[t + 1];
    SA[--BUCKET_BSTAR(c0, c1)] = m - 1;
}


void initBuckets(const sauchar_t *T, saidx_t *SA,
               saidx_t *bucket_A, saidx_t *bucket_B,
               saidx_t n, saidx_t& m,
	       saidx_t num_blocks, saidx_t block_size, saidx_t* bstar_count,
	       sdsl::bit_vector& bstar_flags) {
  
	saidx_t* tempBA = new saidx_t[num_blocks*BUCKET_A_SIZE];
	memset(tempBA, 0, sizeof(saidx_t)*num_blocks*BUCKET_A_SIZE);
	saidx_t* tempBB = new saidx_t[num_blocks*BUCKET_B_SIZE];
	memset(tempBB, 0, sizeof(saidx_t)*num_blocks*BUCKET_B_SIZE);
	parallel_for (saidx_t b = 0; b < num_blocks; b++) {
		// Init values with 0
		saidx_t reducer_m = 0;
		saidx_t* reducer_bucket_A = tempBA + b * BUCKET_A_SIZE;
		saidx_t* reducer_bucket_B = tempBB + b * BUCKET_B_SIZE;

		saidx_t s = std::min(n, block_size * (b+1)) - 1;
		saidx_t e = block_size * b;
		sauchar_t c0,c1;
		// Go until definitly finding A type suffix 
		if (s < n-1) // If not the last block
			while (e < s && T[s] <= T[s+1]) s--; 
		// it is ensured that s is set to a position of a A-type suffix
		// loop goes past e until finding A-type suffix after a B-type suffix
		for(saidx_t i = s, c0 = T[s]; e <= i; ) {
			do { 
				if (i < e && c0 > c1) { // If next block can be sure it's a A-type suffix
					i = -1;
					break;
				}
				++RED_BUCKET_A(c1 = c0); 
			} while((0 <= --i) && ((c0 = T[i]) >= c1));
			if(0 <= i) {
				/* type B* suffix. */
				++RED_BUCKET_BSTAR(c0, c1);
				bstar_flags[i] = 1; // TODO is blocksize multiple of 64!!
				reducer_m++;
				/* type B suffix. */
				for(--i, c1 = c0; (0 <= i) && ((c0 = T[i]) <= c1); --i, c1 = c0) {
					++RED_BUCKET_B(c0, c1);
				}
			}	
		}
		*(bstar_count + b) = reducer_m;
	}
	m = 0; // inclusive prefix sum
	for (int b = 0; b < num_blocks; b++) {
		m += bstar_count[b];
		bstar_count[b] = m;
	} 
	
	saidx_t* SAb = SA + n - m;
	parallel_for(saidx_t i_ = 0; i_ < BUCKET_A_SIZE; ++i_) { bucket_A[i_] = 0; for (int b = 0; b < num_blocks; b++) bucket_A[i_] += tempBA[i_ + b*BUCKET_A_SIZE]; } 
	parallel_for(saidx_t i_ = 0; i_ < BUCKET_B_SIZE; ++i_) { bucket_B[i_] = 0; for (int b = 0; b < num_blocks; b++) bucket_B[i_] += tempBB[i_ + b*BUCKET_B_SIZE]; } 
	cilk_spawn calculateBucketOffsets(bucket_A, bucket_B);	// Buckets offsets can be calculated during the second pass
	delete [] tempBA;
	delete [] tempBB;
	// Write position of BSTAR suffixes to the end of SA array
	// Pack all elements i from [0,n-1] to SA+n-m such that i is a BSTAR suffix	
	parallel_for (saidx_t b = 0; b < num_blocks; b++) {
		saidx_t s = std::min(n, block_size * (b+1)) - 1;
		saidx_t e = block_size * b;
		sauchar_t c0,c1;
		saidx_t m_ = bstar_count[b];
		// Go until definitly finding A type suffix 
		if (s < n-1) // If not the last block
			while (e < s && T[s] <= T[s+1]) s--; 
		// it is ensured that s is set to a position of a A-type suffix
		// loop goes past e until finding A-type suffix after a B-type suffix
		for(saidx_t i = s, c0 = T[s]; e <= i; ) {
			do { 
				if (i < e && c0 > c1) { // If next block can be sure it's a A-type suffix
					i = -1;
					break;
				}
				c1 = c0;
			} while((0 <= --i) && ((c0 = T[i]) >= c1));
			if(0 <= i) {
				/* type B* suffix. */
				SAb[--m_] = i;
				/* type B suffix. */
				for(--i, c1 = c0; (0 <= i) && ((c0 = T[i]) <= c1); --i, c1 = c0) { }
			}	
		}
	}
}


/* Sorts suffixes of type B*. */
static
saidx_t
sort_typeBstar(const sauchar_t *T, saidx_t *SA,
               saidx_t *bucket_A, saidx_t *bucket_B,
               saidx_t n, sdsl::bit_vector& bstar_flags) {
  saidx_t *PAb, *ISAb, *buf;
  saidx_t i, j, k, t, m, bufsize;
  saint_t c0, c1;

  /* Initialize bucket arrays. */
  //parallel_for(saidx_t i_ = 0; i_ < BUCKET_A_SIZE; ++i_) { bucket_A[i_] = 0; }
  //parallel_for(saidx_t i_ = 0; i_ < BUCKET_B_SIZE; ++i_) { bucket_B[i_] = 0; }
  memset(bucket_A, 0, sizeof(saidx_t)*BUCKET_A_SIZE);
  memset(bucket_B, 0, sizeof(saidx_t)*BUCKET_B_SIZE);

  /* Count the number of occurrences of the first one or two characters of each
     type A, B and B* suffix. Moreover, store the beginning position of all
     type B* suffixes into the array SA. */
  saidx_t num_blocks = getWorkers();
  saidx_t block_size = n / num_blocks + 1;    
  saidx_t* bstar_count = new saidx_t[num_blocks];
  memset(bstar_count, 0, sizeof(saidx_t)*num_blocks);
  initBuckets(T, SA, bucket_A, bucket_B, n, m, num_blocks, block_size, bstar_count, bstar_flags);

    PAb = SA + n - m; ISAb = SA + m;
  if(0 < m) {
    initBStarBuckets(T, SA, bucket_B, n, m, PAb);
    nextTime("BSTARSORT, Init buck\t\t");

    /* Sort the type B* substrings using sssort. */
    buf = SA + m, bufsize = n - (2 * m);
    bufsize = 0; // Dont use buffer when multithreadding
    parallel_for(saint_t c0_ = 0; c0_ < ALPHABET_SIZE-1; ++c0_) {
      parallel_for(saint_t c1_ = c0_+1; c1_ < ALPHABET_SIZE; ++c1_) {
        i = BUCKET_BSTAR(c0_, c1_);
	if (c1_ < ALPHABET_SIZE -1) {
		j = BUCKET_BSTAR(c0_, c1_+1);
	} else if (c0_ < ALPHABET_SIZE - 2) {
		j = BUCKET_BSTAR(c0_+1, c0_+2);
	} else {
		j = m;
	}
        if(1 < (j - i)) {
         sssort(T, PAb, SA + i, SA + j, buf, bufsize, 2, n, *(SA + i) == (m - 1));
        }
      }
    }
  nextTime("BSTARSORT, sssort\t\t");
    /* Compute ranks of type B* substrings. */
    saidx_t block_size = m / num_blocks + 1;
    saidx_t* block_start_rank = new saidx_t[num_blocks];
    block_start_rank[0] = 0;
    // First pass calculate block start rank
    parallel_for (saidx_t b = 1; b < num_blocks; b++) {
		saidx_t s = std::min(m, block_size * (b+1)) - 1;
		saidx_t e = block_size * b;
		if (SA[e] < 0) {
			while (e <= s && SA[e] < 0) e++;
			if (e > s && e < m && SA[e] < 0) block_start_rank[b] = -1; // block starts in previous chunks 
			else block_start_rank[b] = e-1;	 // block starts in this chunk at position e-1 
		} else {
			block_start_rank[b] = 0; // there is no block starting in this chunk which ends in later chunks
		}
    }
    for (saidx_t b = num_blocks-2; b > 0; b--) {
	if (block_start_rank[b] == -1)
		block_start_rank[b] = block_start_rank[b+1];
    }
    // Second pass for actuall rank calculation
    parallel_for (saidx_t b = 0; b < num_blocks; b++) {
	saidx_t s = std::min(m, block_size * (b+1)) - 1;
	saidx_t e = block_size * b;
	saidx_t i,j;
	// Deal with block crossing the boarder
	if (b < num_blocks-1 && block_start_rank[b+1] !=  0) {
		j = block_start_rank[b+1];
		while (s >= e && SA[s] < 0) { ISAb[SA[s] = ~SA[s]] = j; s--; }
		if (e <= s) {
			ISAb[SA[s]] = j;
			s--;
		}
	}
	for(i = s; e <= i; --i) {
		while (e <= i && 0 <= SA[i]) { ISAb[SA[i]] = i; i--;}
		j = i;
		while (e <= i && SA[i] < 0) { ISAb[SA[i] = ~SA[i]] = j; i--; }
		if (e <= i)
			ISAb[SA[i]] = j; // End of the bucket with equal suffixes
	}
    }
    delete []block_start_rank;
  }

    nextTime("BSTARSORT, ranks\t\t");
    buf = SA + (2*m);
    bufsize = n - (2*m);
    paralleltrsort(ISAb, SA, m, buf, bufsize);
    nextTime("BSTARSORT, trsort\t\t");
  return m;
}





/* Constructs the suffix array by using the sorted order of type B* suffixes. */
static
void
construct_SA(const sauchar_t *T, saidx_t *SA,
             saidx_t *bucket_A, saidx_t *bucket_B,
             saidx_t n, saidx_t m, sdsl::bit_vector& bstar_flags) {
	saidx_t num_blocks = getWorkers();
	// Make copy of ranks of B* suffixes
	saidx_t* ISAb = new saidx_t[m];
	memcpy(ISAb, SA + m, sizeof(saidx_t) *m);
	// Init with bucket sort
	saidx_t buckets[ALPHABET_SIZE*ALPHABET_SIZE+1];	
	memset(buckets, 0, sizeof(buckets));
	// Count
	buckets[T[n-1] * ALPHABET_SIZE]++;	
	for (saidx_t i = 0; i < n-1; i++) {
		buckets[T[i]*ALPHABET_SIZE + T[i+1]]++;
	}
	// Exclusive prefix sum
	saidx_t sum = 0;
	for (saidx_t i = 0; i <= ALPHABET_SIZE * ALPHABET_SIZE; i++) {
		sum += buckets[i];		
		buckets[i] = sum - buckets[i];
	}	
	// Filling buckets
	SA[buckets[T[n-1]*ALPHABET_SIZE]++] = n-1;
	for (saidx_t i = 0; i < n-1; i++) {
		SA[buckets[T[i]*ALPHABET_SIZE + T[i+1]]++] = i;
	}
	// Move buckets positions to start again
	for (saidx_t i = ALPHABET_SIZE*ALPHABET_SIZE; i > 0; i--) {
		buckets[i] = buckets[i-1];
	}
	buckets[0] = 0;
	printf("%d %d\n", n, buckets[ALPHABET_SIZE * ALPHABET_SIZE]);
	// Init rank support
	sdsl::rank_support_v<1> rs;
	sdsl::util::init_support(rs, &bstar_flags);
	// Sort all of them
	parallel_for (saidx_t b = 0; b < ALPHABET_SIZE*ALPHABET_SIZE; b++) {
		std::sort(SA + buckets[b], SA + buckets[b+1], 
		[&] (saidx_t a, saidx_t b) -> bool {
			a += 2; b += 2;
			while (a < n && b < n) {
				if (T[a] != T[b]) return T[a] < T[b];
				// Look if bstar rank can break the tie
				if (1 & (bstar_flags[a] | bstar_flags[b]) && a < n-1 && b < n-1) {
					a++; b++;
					if (T[a] != T[b]) return T[a] < T[b];
					if (1 & bstar_flags[a-1] && 1 & bstar_flags[b-1]) {
						// Compare ranks of both
						return ISAb[rs.rank(a-1)] < ISAb[rs.rank(b-1)];
					} // Only one flag is set
					if (bstar_flags[a-1] & 1) return true;
					return false;
				}					
				++a;++b;
			}
			if (a == n && b < n) return true;
			return false;
		});
	}
	delete [] ISAb;
}

/* Constructs the burrows-wheeler transformed string directly
   by using the sorted order of type B* suffixes. */
static
saidx_t
construct_BWT(const sauchar_t *T, saidx_t *SA,
              saidx_t *bucket_A, saidx_t *bucket_B,
              saidx_t n, saidx_t m) {
  saidx_t *i, *j, *k, *orig;
  saidx_t s;
  saint_t c0, c1, c2;

  if(0 < m) {
    /* Construct the sorted order of type B suffixes by using
       the sorted order of type B* suffixes. */
    for(c1 = ALPHABET_SIZE - 2; 0 <= c1; --c1) {
      /* Scan the suffix array from right to left. */
      for(i = SA + BUCKET_BSTAR(c1, c1 + 1),
          j = SA + BUCKET_A(c1 + 1) - 1, k = NULL, c2 = -1;
          i <= j;
          --j) {
        if(0 < (s = *j)) {
          assert(T[s] == c1);
          assert(((s + 1) < n) && (T[s] <= T[s + 1]));
          assert(T[s - 1] <= T[s]);
          c0 = T[--s];
          *j = ~((saidx_t)c0);
          if((0 < s) && (T[s - 1] > c0)) { s = ~s; }
          if(c0 != c2) {
            if(0 <= c2) { BUCKET_B(c2, c1) = k - SA; }
            k = SA + BUCKET_B(c2 = c0, c1);
          }
          assert(k < j);
          *k-- = s;
        } else if(s != 0) {
          *j = ~s;
#ifndef NDEBUG
        } else {
          assert(T[s] == c1);
#endif
        }
      }
    }
  }

  /* Construct the BWTed string by using
     the sorted order of type B suffixes. */
  k = SA + BUCKET_A(c2 = T[n - 1]);
  *k++ = (T[n - 2] < c2) ? ~((saidx_t)T[n - 2]) : (n - 1);
  /* Scan the suffix array from left to right. */
  for(i = SA, j = SA + n, orig = SA; i < j; ++i) {
    if(0 < (s = *i)) {
      assert(T[s - 1] >= T[s]);
      c0 = T[--s];
      *i = c0;
      if((0 < s) && (T[s - 1] < c0)) { s = ~((saidx_t)T[s - 1]); }
      if(c0 != c2) {
        BUCKET_A(c2) = k - SA;
        k = SA + BUCKET_A(c2 = c0);
      }
      assert(i < k);
      *k++ = s;
    } else if(s != 0) {
      *i = ~s;
    } else {
      orig = i;
    }
  }

  return orig - SA;
}


/*---------------------------------------------------------------------------*/

/*- Function -*/

saint_t
divsufsort(const sauchar_t *T, saidx_t *SA, saidx_t n) {
  saidx_t *bucket_A, *bucket_B;
  saidx_t m;
  saint_t err = 0;


  /* Check arguments. */
  if((T == NULL) || (SA == NULL) || (n < 0)) { return -1; }
  else if(n == 0) { return 0; }
  else if(n == 1) { SA[0] = 0; return 0; }
  else if(n == 2) { m = (T[0] < T[1]); SA[m ^ 1] = 0, SA[m] = 1; return 0; }

  bucket_A = (saidx_t *)malloc(BUCKET_A_SIZE * sizeof(saidx_t));
  bucket_B = (saidx_t *)malloc(BUCKET_B_SIZE * sizeof(saidx_t));
  sdsl::bit_vector bstar_flags(n,0);
  
  /* Suffixsort. */
  if((bucket_A != NULL) && (bucket_B != NULL)) {
    startTime();    
    m = sort_typeBstar(T, SA, bucket_A, bucket_B, n, bstar_flags);
    construct_SA(T, SA, bucket_A, bucket_B, n, m, bstar_flags);
    nextTime();
  } else {
    err = -2;
  }

  free(bucket_B);
  free(bucket_A);

  return err;
}

saidx_t
divbwt(const sauchar_t *T, sauchar_t *U, saidx_t *A, saidx_t n) {
  saidx_t *B;
  saidx_t *bucket_A, *bucket_B;
  saidx_t m, pidx, i;

  /* Check arguments. */
  if((T == NULL) || (U == NULL) || (n < 0)) { return -1; }
  else if(n <= 1) { if(n == 1) { U[0] = T[0]; } return n; }

  if((B = A) == NULL) { B = (saidx_t *)malloc((size_t)(n + 1) * sizeof(saidx_t)); }
  bucket_A = (saidx_t *)malloc(BUCKET_A_SIZE * sizeof(saidx_t));
  bucket_B = (saidx_t *)malloc(BUCKET_B_SIZE * sizeof(saidx_t));
  sdsl::bit_vector bstar_flags(n, 0);

  /* Burrows-Wheeler Transform. */
  if((B != NULL) && (bucket_A != NULL) && (bucket_B != NULL)) {
    m = sort_typeBstar(T, B, bucket_A, bucket_B, n, bstar_flags);
    pidx = construct_BWT(T, B, bucket_A, bucket_B, n, m);

    /* Copy to output string. */
    U[0] = T[n - 1];
    for(i = 0; i < pidx; ++i) { U[i + 1] = (sauchar_t)B[i]; }
    for(i += 1; i < n; ++i) { U[i] = (sauchar_t)B[i]; }
    pidx += 1;
  } else {
    pidx = -2;
  }

  free(bucket_B);
  free(bucket_A);
  if(A == NULL) { free(B); }

  return pidx;
}

const char *
divsufsort_version(void) {
  return PROJECT_VERSION_FULL;
}
