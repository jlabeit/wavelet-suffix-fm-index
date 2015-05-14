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
#include "select_support_mcl_par.hpp"
#include "rank_support_v_par.hpp"
//#include "wt_int_par.hpp"
#include "wt_huff_par.hpp"
#include <sdsl/construct.hpp>
//#include <sdsl/wt_huff.hpp>

using namespace std;

#define THRESHOLD 10000




pair<WTnode*,long*> WT(symbol* s, uintT n, uintT sigma) {
  sdsl::int_vector<> input;
  input.resize(n);
  symbol x = 0;
  for (uintT i = 0; i < n; i++) 
	  input[i] = (uint64_t)s[i];
  sdsl::util::bit_compress(input);  
  sdsl::store_to_file(input, "input_vector");
  sdsl::wt_huff<> wt;
  sdsl::construct(wt, "input_vector");
  // output node
  sdsl::store_to_file(wt, "output_wt");


  long* empty_result = (long*)malloc(sizeof(long));
  WTnode* emptyNode = (WTnode*)malloc(sizeof(WTnode)); 
  return make_pair((WTnode*)emptyNode,(long*)empty_result);
}
