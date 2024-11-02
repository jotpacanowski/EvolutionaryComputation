#pragma once

#include "headers.hpp"

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

template <typename T>
struct VectorSummary {
    size_t count;
    double mean;
    double std_dev;
    T min;
    T percentile_25;
    T median;
    T percentile_75;
    T max;
    double sum;
};

// VectorSummary<T> describe_vec(const vector<T>& data)
template <typename T>
VectorSummary<T> describe_vec(vector<T>& data)
{
    if (data.empty()) {
        cerr << "Data vector is empty.\n";
        exit(5);
    }

    auto count = data.size();

    VectorSummary<T> summary = {
        // initialize fields
        .count = data.size()
        // ...
    };

    summary.mean = accumulate(data.begin(), data.end(), 0.0) / count;

    double sum_of_squares =
        accumulate(data.begin(), data.end(), 0.0, [summary](double sum, double val) {
            return sum + (val - summary.mean) * (val - summary.mean);
        });
    // Note: Bessel's correction: divide by N-1 for sample standard deviation
    // https://en.wikipedia.org/wiki/Standard_deviation
    // Quote: dividing by n would underestimate the variability
    summary.std_dev = std::sqrt(sum_of_squares / (count - 1));

    // vector<T> sorted_data = data;
    // std::sort(sorted_data.begin(), sorted_data.end());
    std::sort(data.begin(), data.end());
    summary.min = data.front();
    summary.max = data.back();

    auto percentile = [&](double p) {
        double idx = p * (count - 1);
        size_t idx_below = floor(idx);
        size_t idx_above = ceil(idx);
        if (idx_below == idx_above) return data[idx_below];
        // for median: interpolate between middle elements if count is even
        T r = (double)data[idx_below];
        r += (double)(data[idx_above] - data[idx_below]) / 2.0;
        return r;
    };

    summary.percentile_25 = percentile(0.25);
    summary.median = percentile(0.50);
    summary.percentile_75 = percentile(0.75);

    summary.sum = accumulate(data.begin(), data.end(), 0.0);

    return summary;
}

template <typename T>
void print_summary(const VectorSummary<T>& summary)
{
    // cout << format("Count: {}\n", summary.count);
    // cout << format("Mean: {:.2f}\n", summary.mean);
    // cout << format("Std Dev: {:.2f}\n", summary.std_dev);
    // cout << format("Min: {}\n", summary.min);
    // cout << format("25%: {}\n", summary.percentile_25);
    // cout << format("50% (Median): {}\n", summary.median);
    // cout << format("75%: {}\n", summary.percentile_75);
    // cout << format("Max: {}\n", summary.max);
    // cout << format("Sum: {:.2f}\n", summary.sum);

    cout << format("{:>6} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n", "Count",
                   "Mean", "Std. Dev.", "Min", "25%", "50%", "75%", "Max");
    cout << format("{:>6} {:>10.2f} {:>10.2f} {:>10} {:>10} {:>10} {:>10} {:>10}\n",
                   summary.count, summary.mean, summary.std_dev, summary.min,
                   summary.percentile_25, summary.median, summary.percentile_75,
                   summary.max);
}
