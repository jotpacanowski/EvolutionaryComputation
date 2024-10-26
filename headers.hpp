#pragma once

// NOLINTBEGIN
#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

using namespace std;
namespace fs = std::filesystem;  // NOLINT

// NOLINTEND

static const int LARGE_SCORE = 10'000'000;

class Stopwatch {
   private:
    chrono::high_resolution_clock::time_point start_tp;

   public:
    Stopwatch() { this->reset(); }
    void reset() { this->start_tp = chrono::high_resolution_clock::now(); }

    [[nodiscard]] long long count_nanos() const
    {
        const auto end = chrono::high_resolution_clock::now();
        // chrono::duration<long long, std::milli> dur = end - this->start_tp;
        return chrono::nanoseconds(end - start_tp).count();
    }

    [[nodiscard]] std::string pretty_print() const
    {
        const long long t = this->count_nanos();
        if (t <= 90'000'000) {
            // 90 ms = 0.090 s
            return std::format("{:.3} ms", (double)t * 1e-6);
        }
        const double s = (double)t * 1e-9;
        // TODO: if s > 60: "1 min 0 s"
        return std::format("{:.3} s", s);
    }
};
