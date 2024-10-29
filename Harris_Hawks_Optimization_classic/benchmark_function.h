#pragma once

#include <vector>

#define QUANTITY_BENCHMARK_FUNCTIONS 30

using benchmark_function = long double(*)(const std::vector<long double>&);

std::vector<benchmark_function> get_all_functions();
std::vector<long double> get_borders(int dimension);

long double benchmark_f1(const std::vector<long double>& values);
long double benchmark_f2(const std::vector<long double>& values);
long double benchmark_f3(const std::vector<long double>& values);
long double benchmark_f4(const std::vector<long double>& values);
long double benchmark_f5(const std::vector<long double>& values);
long double benchmark_f6(const std::vector<long double>& values);
long double benchmark_f7(const std::vector<long double>& values);
long double benchmark_f8(const std::vector<long double>& values);
long double benchmark_f9(const std::vector<long double>& values);
long double benchmark_f10(const std::vector<long double>& values);

long double benchmark_f11(const std::vector<long double>& values);
long double benchmark_f12(const std::vector<long double>& values);
long double benchmark_f13(const std::vector<long double>& values);
long double benchmark_f14(const std::vector<long double>& values);
long double benchmark_f15(const std::vector<long double>& values);
long double benchmark_f16(const std::vector<long double>& values);
long double benchmark_f17(const std::vector<long double>& values);
long double benchmark_f18(const std::vector<long double>& values);
long double benchmark_f19(const std::vector<long double>& values);
long double benchmark_f20(const std::vector<long double>& values);

long double benchmark_f21(const std::vector<long double>& values);
long double benchmark_f22(const std::vector<long double>& values);
long double benchmark_f23(const std::vector<long double>& values);
long double benchmark_f24(const std::vector<long double>& values);
long double benchmark_f25(const std::vector<long double>& values);
long double benchmark_f26(const std::vector<long double>& values);
long double benchmark_f27(const std::vector<long double>& values);
long double benchmark_f28(const std::vector<long double>& values);
long double benchmark_f29(const std::vector<long double>& values);
long double benchmark_f30(const std::vector<long double>& values);