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
#include "divsufsort_private.h"
#include "parallel.h"
#include "sequence.h"
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
//#include <cilk/reducer_opadd.h>
#ifdef _OPENMP
# include <omp.h>
#endif
#include <atomic>

/*- Private Functions -*/

/* Sorts suffixes of type B*. */
static
saidx_t
sort_typeBstar(const sauchar_t *T, saidx_t *SA,
               saidx_t *bucket_A, saidx_t *bucket_B,
               saidx_t n) {
  saidx_t *PAb, *ISAb, *buf;
  saidx_t i, j, k, t, m, bufsize;
  saint_t c0, c1;

  /* Initialize bucket arrays. */
  parallel_for(saidx_t i_ = 0; i_ < BUCKET_A_SIZE; ++i_) { bucket_A[i_] = 0; }
  parallel_for(saidx_t i_ = 0; i_ < BUCKET_B_SIZE; ++i_) { bucket_B[i_] = 0; }

  /* Count the number of occurrences of the first one or two characters of each
     type A, B and B* suffix. Moreover, store the beginning position of all
     type B* suffixes into the array SA. */
  saidx_t  block_size = 1024*1024*16;
  saidx_t num_blocks = n / block_size + 1;	
  saidx_t* bstar_count = new saidx_t[num_blocks];
  memset(bstar_count, 0, sizeof(saidx_t)*num_blocks);
  { // Count A,B,BSTAR types
	saidx_t* tempBA = new saidx_t[num_blocks*BUCKET_A_SIZE];
	memset(tempBA, 0, sizeof(saidx_t)*num_blocks*BUCKET_A_SIZE);
	saidx_t* tempBB = new saidx_t[num_blocks*BUCKET_B_SIZE];
	memset(tempBB, 0, sizeof(saidx_t)*num_blocks*BUCKET_B_SIZE);
	parallel_for (saidx_t b = 0; b < num_blocks; b++) {
		// Init values with 0
		saidx_t *reducer_m = (bstar_count + b);
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
				(*reducer_m)++;
				/* type B suffix. */
				for(--i, c1 = c0; (0 <= i) && ((c0 = T[i]) <= c1); --i, c1 = c0) {
					++RED_BUCKET_B(c0, c1);
				}
			}	
		}
	}
	m = 0; // inclusive prefix sum
	for (int b = 0; b < num_blocks; b++) {
		m += bstar_count[b];
		bstar_count[b] = m;
	} 
	saidx_t* SAb = SA + n - m;
	parallel_for(saidx_t i_ = 0; i_ < BUCKET_A_SIZE; ++i_) { bucket_A[i_] = 0; for (int b = 0; b < num_blocks; b++) bucket_A[i_] += tempBA[i_ + b*BUCKET_A_SIZE]; } 
	parallel_for(saidx_t i_ = 0; i_ < BUCKET_B_SIZE; ++i_) { bucket_B[i_] = 0; for (int b = 0; b < num_blocks; b++) bucket_B[i_] += tempBB[i_ + b*BUCKET_B_SIZE]; } 
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
/*
note:
  A type B* suffix is lexicographically smaller than a type B suffix that
  begins with the same first two characters.
*/

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

  if(0 < m) {
    /* Sort the type B* suffixes by their first two characters. */
    PAb = SA + n - m; ISAb = SA + m;
    for(i = m - 2; 0 <= i; --i) {
      t = PAb[i], c0 = T[t], c1 = T[t + 1];
      SA[--BUCKET_BSTAR(c0, c1)] = i;
    }
    t = PAb[m - 1], c0 = T[t], c1 = T[t + 1];
    SA[--BUCKET_BSTAR(c0, c1)] = m - 1;

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
    //for (int i = 0; i < n; i++) {
	//printf("%d ", SA[i]);
    //}printf("\n");
    
    /* Compute ranks of type B* substrings. */
    for(i = m - 1; 0 <= i; --i) {
      if(0 <= SA[i]) {
        j = i;
        do { ISAb[SA[i]] = i; } while((0 <= --i) && (0 <= SA[i]));
        //SA[i + 1] = i - j; // skip values for already sorted sequence (negative!)
        if(i <= 0) { break; }
      }
      j = i;
      do { ISAb[SA[i] = ~SA[i]] = j; } while(SA[--i] < 0);
      ISAb[SA[i]] = j; // End of the bucket with equal suffixes
    }
    // Construct the inverse suffix array of type B* suffixes using trsort. 
    //trsort(ISAb, SA, m, 1);
    buf = SA + (2*m);
    bufsize = n - (2*m);
    paralleltrsort(ISAb, SA, m, buf, bufsize);

    /* Set the sorted order of type B* suffixes. */
	/*    
    for(i = n - 1, j = m, c0 = T[n - 1]; 0 <= i;) {
      for(--i, c1 = c0; (0 <= i) && ((c0 = T[i]) >= c1); --i, c1 = c0) { } // A suffix
      if(0 <= i) {
        t = i; // B* suffix
        for(--i, c1 = c0; (0 <= i) && ((c0 = T[i]) <= c1); --i, c1 = c0) { } // B suffix
        SA[ISAb[--j]] = ((t == 0) || (1 < (t - i))) ? t : ~t;
      }
    }*/
    parallel_for (saidx_t b = 0; b < num_blocks; b++) {
		saidx_t s = std::min(n, block_size * (b+1)) - 1;
		saidx_t e = block_size * b;
		sauchar_t c0,c1;
		saidx_t t,j,i;
		// Go until definitly finding A type suffix 
		if (s < n-1) // If not the last block
			while (e < s && T[s] <= T[s+1]) s--; 
		for(i = s, j = bstar_count[b], c0 = T[s]; e <= i;) {
			for(--i, c1 = c0; (0 <= i) && ((c0 = T[i]) >= c1); --i, c1 = c0) { 
				if (i < e && c0 > c1) { // If next block can be sure it's a A-type suffix
					i = -1;
					break;
				}
			} // A suffix
			if(0 <= i) {
				t = i; // B* suffix
				for(--i, c1 = c0; (0 <= i) && ((c0 = T[i]) <= c1); --i, c1 = c0) { } // B suffix
				SA[ISAb[--j]] = ((t == 0) || (1 < (t - i))) ? t : ~t;
			}
		}
    }
    delete [] bstar_count;

    /* Calculate the index of start/end point of each bucket. */
    BUCKET_B(ALPHABET_SIZE - 1, ALPHABET_SIZE - 1) = n; /* end point */
    for(c0 = ALPHABET_SIZE - 2, k = m - 1; 0 <= c0; --c0) {
      i = BUCKET_A(c0 + 1) - 1;
      for(c1 = ALPHABET_SIZE - 1; c0 < c1; --c1) {
        t = i - BUCKET_B(c0, c1);
        BUCKET_B(c0, c1) = i; /* end point */

        /* Move all type B* suffixes to the correct position. */
        for(i = t, j = BUCKET_BSTAR(c0, c1);
            j <= k;
            --i, --k) { SA[i] = SA[k]; }
      }
      BUCKET_BSTAR(c0, c0 + 1) = i - BUCKET_B(c0, c0) + 1; /* start point */
      BUCKET_B(c0, c0) = i; /* end point */
    }
  }

  return m;
}

/* Constructs the suffix array by using the sorted order of type B* suffixes. */
static
void
construct_SA(const sauchar_t *T, saidx_t *SA,
             saidx_t *bucket_A, saidx_t *bucket_B,
             saidx_t n, saidx_t m) {
  saidx_t *i, *j, *k;
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
          *j = ~s;
          c0 = T[--s];
          if((0 < s) && (T[s - 1] > c0)) { s = ~s; }
          if(c0 != c2) {
            if(0 <= c2) { BUCKET_B(c2, c1) = k - SA; }
            k = SA + BUCKET_B(c2 = c0, c1);
          }
          assert(k < j);
          *k-- = s;
        } else {
          assert(((s == 0) && (T[s] == c1)) || (s < 0));
          *j = ~s;
        }
      }
    }
  }
  
  /* Construct the suffix array by using
     the sorted order of type B suffixes. */
  k = SA + BUCKET_A(c2 = T[n - 1]);
  *k++ = (T[n - 2] < c2) ? ~(n - 1) : (n - 1);
  /* Scan the suffix array from left to right. */
  for(i = SA, j = SA + n; i < j; ++i) {
    if(0 < (s = *i)) {
      assert(T[s - 1] >= T[s]);
      c0 = T[--s];
      if((s == 0) || (T[s - 1] < c0)) { s = ~s; }
      if(c0 != c2) {
        BUCKET_A(c2) = k - SA;
        k = SA + BUCKET_A(c2 = c0);
      }
      assert(i < k);
      *k++ = s;
    } else {
      assert(s < 0);
      *i = ~s;
    }
  }
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

  /* Suffixsort. */
  if((bucket_A != NULL) && (bucket_B != NULL)) {
    m = sort_typeBstar(T, SA, bucket_A, bucket_B, n);
    construct_SA(T, SA, bucket_A, bucket_B, n, m);
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

  /* Burrows-Wheeler Transform. */
  if((B != NULL) && (bucket_A != NULL) && (bucket_B != NULL)) {
    m = sort_typeBstar(T, B, bucket_A, bucket_B, n);
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
