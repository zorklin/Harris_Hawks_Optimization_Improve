#pragma once

#include <vector>

#define QUANTITY_BENCHMARK_FUNCTIONS 30

using benchmark_function = long double(*)(std::vector<long double>);

std::vector<benchmark_function> get_all_functions();
std::vector<long double> get_borders(int dimension);

long double benchmark_f1(std::vector<long double> values);
long double benchmark_f2(std::vector<long double> values);
long double benchmark_f3(std::vector<long double> values);
long double benchmark_f4(std::vector<long double> values);
long double benchmark_f5(std::vector<long double> values);
long double benchmark_f6(std::vector<long double> values);
long double benchmark_f7(std::vector<long double> values);
long double benchmark_f8(std::vector<long double> values);
long double benchmark_f9(std::vector<long double> values);
long double benchmark_f10(std::vector<long double> values);

long double benchmark_f11(std::vector<long double> values);
long double benchmark_f12(std::vector<long double> values);
long double benchmark_f13(std::vector<long double> values);
long double benchmark_f14(std::vector<long double> values);
long double benchmark_f15(std::vector<long double> values);
long double benchmark_f16(std::vector<long double> values);
long double benchmark_f17(std::vector<long double> values);
long double benchmark_f18(std::vector<long double> values);
long double benchmark_f19(std::vector<long double> values);
long double benchmark_f20(std::vector<long double> values);

long double benchmark_f21(std::vector<long double> values);
long double benchmark_f22(std::vector<long double> values);
long double benchmark_f23(std::vector<long double> values);
long double benchmark_f24(std::vector<long double> values);
long double benchmark_f25(std::vector<long double> values);
long double benchmark_f26(std::vector<long double> values);
long double benchmark_f27(std::vector<long double> values);
long double benchmark_f28(std::vector<long double> values);
long double benchmark_f29(std::vector<long double> values);
long double benchmark_f30(std::vector<long double> values);