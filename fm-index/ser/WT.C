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
#include <divsufsort.h>
#include "WT.h"
#include "parallel.h"
#include "utils.h"
#include <iostream>
#include "sequence.h"


#define THRESHOLD 10000


void FMIndex(symbol* s, long n, sdsl::int_vector<sizeof(symbol)*8>& bwt) {
	// Calculate SA
    int32_t *data = newA(int32_t, n); 
    divsufsort(s, (int32_t*)data, n);
	// Calculate BWT
	int32_t to_add[2] = {(int32_t)-1,(int32_t)n-1};
    for (long i=0; i < n; ++i) {
		bwt[i] = s[data[i]+to_add[data[i]==0]];
	}
    free(data);
	// Construct WT
	sdsl::wt_huff<> wt(bwt,n);
}
