#include "miscellaneous/sieve.hpp"
#include <benchmark/benchmark.h>
#include <cstdint> // integer types

using namespace mathjoy;

static void BM_MathJoySift_bit(benchmark::State &state) {
  for (auto _ : state) {
    long long length = 1e6;
    std::vector<bool> masks(length);
    sift(masks.begin(), length);
  }
}

static void BM_MathJoySift_uint8(benchmark::State &state) {
  for (auto _ : state) {
    long long length = 1e6;
    std::vector<uint8_t> masks(length);
    sift(masks.begin(), length);
  }
}

static void BM_MathJoySift_uint16(benchmark::State &state) {
  for (auto _ : state) {
    long long length = 1e6;
    std::vector<uint16_t> masks(length);
    sift(masks.begin(), length);
  }
}

static void BM_MathJoySift_uint32(benchmark::State &state) {
  for (auto _ : state) {
    long long length = 1e6;
    std::vector<uint32_t> masks(length);
    sift(masks.begin(), length);
  }
}

static void BM_MathJoySift_uint64(benchmark::State &state) {
  for (auto _ : state) {
    long long length = 1e6;
    std::vector<uint64_t> masks(length);
    sift(masks.begin(), length);
  }
}

BENCHMARK(BM_MathJoySift_bit);
BENCHMARK(BM_MathJoySift_uint8);
BENCHMARK(BM_MathJoySift_uint16);
BENCHMARK(BM_MathJoySift_uint32);
BENCHMARK(BM_MathJoySift_uint64);

BENCHMARK_MAIN();
