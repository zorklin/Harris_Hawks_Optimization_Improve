#include <iostream>
#include <vector>
#include <fstream>
#include "HHO.h"
#include "benchmark_function.h"

using namespace std;

#define QUANTITY_RUNS 10

long double calculate_deviation(vector<long double> values, long double mean) {
    size_t size = values.size();
    long double deviation = 0.0;

    for (size_t i = 0; i < size; i++) {
        deviation += pow(values[i] - mean, 2);
    }
    return sqrt(deviation / size);
}

int main(){

    ofstream file("data.txt");
    if (!file.is_open()) {
        cerr << "err open file" << endl;
        return 1;
    }

    vector<int> dimension = { 5, 10, 30, 50 };
    vector<benchmark_function> test_function = get_all_functions();

    cout << "QUANTITY_RUNS: " << QUANTITY_RUNS << endl;
    file << "QUANTITY_RUNS: " << QUANTITY_RUNS << endl;
    for (int i = 0; i < dimension.size(); i++) {
        cout << endl << "Dimension: " << dimension[i] << endl;
        file << endl << "Dimension: " << dimension[i] << endl;
        vector<long double> borders = get_borders(dimension[i]);
        for (int j = 0; j < QUANTITY_BENCHMARK_FUNCTIONS; j++) {
            long double best_value = INFINITY;
            long double averege = 0.0l;
            vector<long double> answers;
            answers.resize(QUANTITY_RUNS);
            for (int k = 0; k < QUANTITY_RUNS; k++) {
                int T = dimension[i] * 100, size = dimension[i] * 20;
                vector<long double> answer = harris_hawks_optimazation(T, size, borders[j], -borders[j], dimension[i], test_function[j]);
                long double temp = test_function[j](answer);
                averege += temp;
                best_value = min(best_value, temp);
                answers[k] = temp;
            }
            averege /= dimension[i];
            long double deviation = calculate_deviation(answers, averege);
            cout << "Function: " << j + 1 << ", Averege: " << averege << ", Best: " << best_value << ", Deviation: " << deviation << endl;
            file << "Function: " << j + 1 << ", Averege: " << averege << ", Best: " << best_value << ", Deviation: " << deviation << endl;
        }
    }

    return 0;
}