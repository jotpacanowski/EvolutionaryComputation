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

static inline std::string format_with_spaces(long long number)
{
    auto str = std::to_string(number);
    int n = str.length() - 3;
    while (n > 0) {
        str.insert(n, " ");
        n -= 3;
    }
    return str;
}

static inline std::string format_with_spaces(double number)
{
    auto [integral, fractional] = std::div(static_cast<long long>(number * 1000), 1000ll);
    auto int_part = format_with_spaces(integral);
    return std::format("{}.{:03}", int_part, std::abs(fractional));
}

static inline ofstream openOutFile(string_view filename, string_view prefix = "")
{
    fs::path path;
    if (prefix.empty()) {
        path = fs::path(filename);
    }
    else {
        path = fs::path(prefix) / filename;
    }
    std::ofstream out(path);
    if (!out.is_open()) {
        // perror("Error");
        std::error_code ec(errno, std::generic_category());
        cerr << "Error: " << ec.message() << " (code: " << ec.value() << ")\n";

        fs::path parent_dir = fs::path(path).parent_path();
        if (!fs::exists(parent_dir)) {
            cerr << "Creating missing directory: " << parent_dir << endl;
            std::error_code ec;
            fs::create_directories(parent_dir, ec);
            // Check error without exceptions
            if (ec.value() != 0) {
                cerr << "Failed: " << ec.message() << endl;
                exit(1);
            }
            // Try again
            out.open(path);
            if (!out.is_open()) {
                perror("second open");
                exit(1);
            }
        }
        else {
            exit(1);
        }
    }
    return out;
}

static inline void printResults(const vector<int>& sol, bool toFile = false,
                                string_view filename = "solution.txt",
                                string_view prefix = "data/results")
{
    if (toFile) {
        ofstream out = openOutFile(filename, prefix);
        std::copy(sol.begin(), sol.end(), std::ostream_iterator<int>(out, "\n"));
    }
    else {
        std::copy(sol.begin(), sol.end(), std::ostream_iterator<int>(std::cout, "\n"));
    }
}

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

struct SolutionStats {
    vector<int> best_sol;
    int best_sol_value = LARGE_SCORE;

    vector<int> worst_sol;
    int worst_sol_value = 0;

    long long average_numerator = 0;
    long long count = 0;

    [[nodiscard]] double average() const { return (double)(average_numerator) / count; }

    void track(const vector<int>& sol, int value)
    {
        if (value < best_sol_value) {
            best_sol_value = value;
            best_sol = sol;
        }
        if (value > worst_sol_value) {
            worst_sol_value = value;
            worst_sol = sol;
        }
        average_numerator += value;
        count += 1;
    }

    vector<long long> timings;

    void add_time(long long nanosec) { timings.push_back(nanosec); }

    [[nodiscard]] string format_latex_3() const
    {
        return std::format("{:>7} & {:>12} & {:>7}",
                           format_with_spaces((long long)this->best_sol_value),
                           format_with_spaces(this->average()),
                           format_with_spaces((long long)this->worst_sol_value));
    }
    [[nodiscard]] string format_latex_one_field() const
    {
        return std::format("{} ({} - {})", format_with_spaces((long long)this->average()),
                           format_with_spaces((long long)this->best_sol_value),
                           format_with_spaces((long long)this->worst_sol_value));
    }
};
