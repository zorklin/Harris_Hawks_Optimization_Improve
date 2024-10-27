#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include "HHO.h"
#include "benchmark_function.h"

using namespace std;
using namespace std::chrono;

#define QUANTITY_RUNS 5

long double calculate_deviation(vector<long double> values, long double mean) {
    size_t size = values.size();
    long double deviation = 0.0l;

    for (size_t i = 0; i < size; i++) {
        deviation += pow(values[i] - mean, 2.0l);
    }
    return sqrt(deviation / size);
}

int main() {

    ofstream file("data.txt");
    if (!file.is_open()) {
        cerr << "err open file" << endl;
        return 1;
    }

    vector<int> dimension = { 5, 10, 30 };
    vector<benchmark_function> test_function = get_all_functions();

    cout << "QUANTITY_RUNS: " << QUANTITY_RUNS << endl;
    file << "QUANTITY_RUNS: " << QUANTITY_RUNS << endl;
    long double total_rank = 0.0l;
    for (int i = 0; i < dimension.size(); i++) {
        cout << endl << "Dimension: " << dimension[i] << endl;
        file << endl << "Dimension: " << dimension[i] << endl;
        vector<long double> borders = get_borders(dimension[i]);
        long double rank_dimension = 0.0l;
        for (int j = 0; j < QUANTITY_BENCHMARK_FUNCTIONS; j++) {
            long double best_value = INFINITY;
            long double averege = 0.0l;
            vector<long double> answers;
            answers.resize(QUANTITY_RUNS);
            double total_duration = 0.0;
            for (int k = 0; k < QUANTITY_RUNS; k++) {
                int T = dimension[i] * 100, size = dimension[i] * 20;
                high_resolution_clock::time_point start = high_resolution_clock::now();
                vector<long double> answer = harris_hawks_optimazation(T, size, borders[j], -borders[j], dimension[i], test_function[j]);
                high_resolution_clock::time_point end = high_resolution_clock::now();
                duration<double> time_span = duration_cast<duration<double>>(end - start);
                total_duration += time_span.count();

                long double temp = test_function[j](answer);
                averege += temp;
                best_value = min(best_value, temp);
                answers[k] = temp;
            }
            averege /= QUANTITY_RUNS;
            long double deviation = calculate_deviation(answers, averege);
            double average_time_s = total_duration / QUANTITY_RUNS;
            rank_dimension += averege;
            cout << "Function: " << j + 1
                << ", Averege: " << averege
                << ", Best: " << best_value
                << ", Deviation: " << deviation
                << ", Average Time: " << average_time_s << " s" << endl;

            file << "Function: " << j + 1
                << ", Averege: " << averege
                << ", Best: " << best_value
                << ", Deviation: " << deviation
                << ", Average Time: " << average_time_s << " s" << endl;
        }
        total_rank += rank_dimension / QUANTITY_RUNS;
        cout << "Rand dimension:" << rank_dimension;
    }
    cout << endl << "End of the test for: " << QUANTITY_RUNS << " runs" << endl << "Total rank: " << total_rank / dimension.size() << endl;

    return 0;
}