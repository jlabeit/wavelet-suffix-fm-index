#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <ctime>
#include <unistd.h>
#include <omp.h>
#include "gettime.h"
#include "getmemory.h"

#include "inmem_psascan.h"
#include "utils.h"
#include "uint40.h"


template<typename saidx_t>
void test(std::uint8_t *text, std::uint64_t text_length, std::uint64_t n_threads, std::uint64_t extra_bytes_per_symbol) {
  unsigned char *tab = (unsigned char *)malloc(text_length * (sizeof(saidx_t) + 1));
  psascan_private::inmem_psascan_private::inmem_psascan<saidx_t>(text, text_length,
      tab, n_threads, false, false, NULL, n_threads, extra_bytes_per_symbol);
  free(tab);
}

template<typename saidx_t>
void test_file(const char *filename, std::uint64_t n_threads, std::uint64_t extra_bytes_per_symbol, std::uint64_t runs = 3) {
  // Allocate and read text.
  std::uint64_t text_length = utils::file_size(std::string(filename));
  std::uint8_t *text = new std::uint8_t[text_length];
  utils::read_from_file(text, text_length, std::string(filename));

  // Run tests.
  test<saidx_t>(text, text_length, n_threads, extra_bytes_per_symbol);
  for (std::uint64_t i = 0; i < runs; ++i) {
    startTime();
    test<saidx_t>(text, text_length, n_threads, extra_bytes_per_symbol);
    nextTimeN();
  }
  std::cout<<"Peak-memory: " <<getPeakRSS() / (1024*1024) << std::endl;

  // Free text.
  delete[] text;
}

int main(int argc, char **argv) {
  if (argc == 1) {
    fprintf(stderr, "Usage: %s FILE1 FILE2 ...\nBenchmark inmem_psascan on provided text files.\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  std::uint64_t n_threads = omp_get_max_threads();
    //test_file<std::int32_t>(argv[i], n_threads, 1);    // most space efficient
  test_file<std::int32_t>(argv[1], n_threads, 5);    // fastest
    //test_file<uint40>(argv[i], n_threads, 3);          // default

    // Remaining combinations (not tested):
    // test_file<std::int32_t>(argv[i], n_threads, 2);
    // test_file<std::int32_t>(argv[i], n_threads, 3);
    // test_file<std::int32_t>(argv[i], n_threads, 4);
    // test_file<uint40>(argv[i], n_threads, 1);
    // test_file<uint40>(argv[i], n_threads, 2);
    // test_file<uint40>(argv[i], n_threads, 4);
    // test_file<uint40>(argv[i], n_threads, 5);
}
