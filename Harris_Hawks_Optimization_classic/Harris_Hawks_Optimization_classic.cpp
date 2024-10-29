#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <cmath>
#include "HHO.h"
#include "benchmark_function.h"

using namespace std;
using namespace std::chrono;

#define QUANTITY_RUNS 10

long double calculate_deviation(const vector<long double>& values, const long double mean) {
    size_t size = values.size();
    long double deviation = 0.0l;

    for (size_t i = 0; i < size; i++) {
        deviation += pow(values[i] - mean, 2.0l);
    }
    return sqrt(deviation / size);
}

long double calculate_median(vector<long double>& values) {
    size_t size = values.size();
    sort(values.begin(), values.end());
    if (size % 2) return values[size / 2];
    return (values[size / 2 - 1] + values[size / 2]) / 2.0l;
}

int main() {
    ofstream file("data.txt");
    if (!file.is_open()) {
        cerr << "err opening file" << endl;
        return 1;
    }

    vector<int> dimension = { 5, 10, 30 };
    vector<benchmark_function> test_function = get_all_functions();

    cout << "QUANTITY_RUNS: " << QUANTITY_RUNS << endl;
    file << "QUANTITY_RUNS: " << QUANTITY_RUNS << endl;
    for (int i = 0; i < dimension.size(); i++) {
        cout << endl << "Dimension: " << dimension[i] << endl;
        file << endl << "Dimension: " << dimension[i] << endl;
        vector<long double> borders = get_borders(dimension[i]);
        for (int j = 0; j < QUANTITY_BENCHMARK_FUNCTIONS; j++) {
            long double best_value = INFINITY;
            long double average = 0.0l;
            vector<long double> answers;
            answers.resize(QUANTITY_RUNS);
            double total_duration = 0.0;
            for (int k = 0; k < QUANTITY_RUNS; k++) {
                int T = 1000, size = 500;
                high_resolution_clock::time_point start = high_resolution_clock::now();
                vector<long double> answer = harris_hawks_optimazation(T, size, borders[j], -borders[j], dimension[i], test_function[j]);
                high_resolution_clock::time_point end = high_resolution_clock::now();
                duration<double> time_span = duration_cast<duration<double>>(end - start);
                total_duration += time_span.count();

                long double temp = test_function[j](answer);
                average += temp;
                best_value = min(best_value, temp);
                answers[k] = temp;
            }
            average /= QUANTITY_RUNS;
            long double deviation = calculate_deviation(answers, average);
            long double median = calculate_median(answers);
            double average_time_s = total_duration / QUANTITY_RUNS;
            cout << "Function: " << j + 1
                << ", Average: " << average
                << ", Best: " << best_value
                << ", Deviation: " << deviation
                << ", Median: " << median
                << ", Average Time: " << average_time_s << " s" << endl;

            file << "Function: " << j + 1
                << ", Average: " << average
                << ", Best: " << best_value
                << ", Deviation: " << deviation
                << ", Median: " << median
                << ", Average Time: " << average_time_s << " s" << endl;
        }
    }
    cout << endl << "End of the test for: " << QUANTITY_RUNS << " runs" << endl;

    return 0;
}
