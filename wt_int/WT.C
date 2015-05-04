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
#include "wt_int_par.hpp"

using namespace std;

#define THRESHOLD 10000




pair<WTnode*,long*> WT(symbol* s, uintT n, uintT sigma) {
  sdsl::int_vector<sizeof(symbol)> input(n);
  for (uintT i = 0; i < n; i++)
	  input[i] = s[i];
  sdsl::wt_int<> wt(input);
  long* result_tree = (long*)malloc(sizeof(symbol)*n);
  memcpy((void*)result_tree, (void*)wt.m_tree.data(), sizeof(symbol)*n);
  WTnode* emptyNode = (WTnode*)malloc(sizeof(WTnode)); 
  return make_pair((WTnode*)emptyNode,(long*)result_tree);
}
