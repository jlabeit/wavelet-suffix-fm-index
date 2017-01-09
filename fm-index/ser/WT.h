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

#define POINTERS 1
//#define INT 1

#include "gettime.h"
#include "wt_pc.hpp"
#include <sdsl/csa_wt.hpp>
//#include <sdsl/construct.hpp>
//#include <sdsl/util.hpp>

//#include <sdsl/wt_int.hpp>
//typedef ulong intOffset;

#ifdef INT
typedef uintT symbol;
#else
typedef unsigned char symbol;
#endif

//pair<WTnode*,long*> WT(symbol* s, uintT n, uintT sigma);
void FMIndex(symbol* s, long n, sdsl::int_vector<sizeof(symbol)*8>& bwt);
