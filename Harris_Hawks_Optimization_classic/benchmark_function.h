#pragma once

#include <cmath>
#include <vector>
#include <random>

#define SIZE_BENCHMARK_FUNCTION 5

using benchmark_function = double(*)(std::vector<double>);

double benchmark_f1(std::vector<double> point) {
	double sum = 0.0;
	for (int i = 0; i < point.size(); i++) {
		sum += pow(point[i], 2.0);
	}
	return sum;
}

double benchmark_f2(std::vector<double> point) {
	double sum = 0.0;
	for (int i = 0; i < point.size(); i++) {
		sum += i * pow(point[i], 2.0);
	}
	return sum;
}

double benchmark_f3(std::vector<double> point) {
	double sum = 0.0;
	for (int i = 0; i < point.size(); i++) {
		sum += pow(fabs(point[i]), i + 2);
	}
	return sum;
}

double benchmark_f4(std::vector<double> point) {
	double sum = 0.0;
	for (int i = 0; i < point.size(); i++) {
		sum += fabs(pow(point[i], 5.0) - 3.0 * pow(point[i], 4.0) + 4 * pow(point[i], 3.0) + 2 * pow(point[i], 2.0) - 10 * point[i] - 4.0);
	}
	return sum;
}

double benchmark_f5(std::vector<double> point) {
	double sum = 0.0;
	for (int i = 0; i < point.size(); i++) {
		sum += pow(point[i], 2.0);
	}
	return 1.0 - ((1.0 + cos(12.0 * sqrt(sum))) / (0.5 * sum + 2.0));
}

std::vector<benchmark_function> initialize_benchmark_functions() {
	return { benchmark_f1, benchmark_f2, benchmark_f3, benchmark_f4, benchmark_f5 };
}