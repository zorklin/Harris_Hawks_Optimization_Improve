#pragma once

#include "benchmark_function.h"

constexpr long double PI = 3.141592653589793238462643383279502884l;
constexpr long double E = 2.71828182845904523536028747135266249775724709369995L;

long double benchmark_f1(const std::vector<long double>& values) {
	size_t size = values.size();
	long double sum = 0.0l;
	for (int i = 0; i < size; i++) {
		sum += pow(values[i], 2.0l);
	}
	return sum;
}

long double benchmark_f2(const std::vector<long double>& values) {
	size_t size = values.size();
	long double sum = 0.0l;
	for (int i = 0; i < size; i++) {
		sum += (1.0 + i) * pow(values[i], 2.0l);
	}
	return sum;
}

long double benchmark_f3(const std::vector<long double>& values) {
	size_t size = values.size();
	long double sum = 0.0l;
	for (int i = 0; i < size; i++) {
		sum += fabs(pow(values[i], 2.0l + i));
	}
	return sum;
}

long double benchmark_f4(const std::vector<long double>& values) {
	size_t size = values.size();
	long double sum = 0.0l;
	for (int i = 0; i < size; i++) {
		sum += fabs(pow(values[i], 5.0l) - 3.0l * pow(values[i], 4.0l) + 4.0l * pow(values[i], 3.0l) + 2.0l * pow(values[i], 2.0l) - 10l * values[i] - 4.0l);
	}
	return sum;
}

long double benchmark_f5(const std::vector<long double>& values) {
	size_t size = values.size();
	long double sum = 0.0l;
	for (int i = 0; i < size; i++) {
		sum += pow(values[i], 2.0l);
	}
	return 1.0l - (1.0l + cos(12.0l * sqrt(sum))) / (0.5l * sum + 2.0l);
}

long double benchmark_f6(const std::vector<long double>& values) {
	const long double a = 0.5l;
	const long double b = 3.0l;
	const int k_max = 20;

	long double f_sum = 0.0l;
	size_t size = values.size();

	std::vector<long double> a_powers(k_max);
	std::vector<long double> b_powers(k_max);
	for (int k = 0; k < k_max; k++) {
		a_powers[k] = pow(a, k);
		b_powers[k] = pow(b, k);
	}

	for (size_t i = 0; i < size; i++) {
		long double temp = 0.0l;
		for (int k = 0; k < k_max; k++) {
			temp += a_powers[k] * cos(2.0l * PI * b_powers[k] * (values[i] + 0.5l));
		}
		f_sum += temp;
	}

	long double s_sum = 0.0l;
	for (int k = 0; k < k_max; k++) {
		s_sum += a_powers[k] * cos(PI * b_powers[k]);
	}

	return f_sum - size * s_sum;
}

long double benchmark_f7(const std::vector<long double>& values) {
	long double sum = 0.0;
	size_t size = values.size();
	for (int i = 0; i < size; i++) {
		sum += fabs(values[i] * sin(values[i]) + 0.1l * values[i]);
	}
	return sum;
}

long double benchmark_f8(const std::vector<long double>& values) {
	long double f_sum = 0.0, s_sum = 0.0l;
	size_t size = values.size();
	for (int i = 0; i < size; i++) {
		f_sum += pow(values[i], 2.0l);
		s_sum += cos(2.0l * PI * values[i]);
	}
	return -20.0 * pow(E, -0.2l * sqrt(1.0l / size * f_sum)) - pow(E, 1.0l / size * s_sum) + E + 20.0l;
}

long double benchmark_f9(const std::vector<long double>& values) {
	long double sum = 0.0, product = 1.0l;
	size_t size = values.size();
	for (int i = 0; i < size; i++) {
		sum += pow(values[i], 2.0l);
		product *= cos(values[i] / sqrt(1.0l + i));
	}
	return 1.0l / 4000.0l * sum - product + 1.0l;
}

long double benchmark_f10(const std::vector<long double>& values) {
	long double sum = 0.0;
	size_t size = values.size();
	for (int i = 0; i < size; i++) {
		sum += pow(values[i], 2.0l) - 10.0l * cos(2.0l * PI * values[i]);
	}
	return 10.0l * size + sum;
}

long double benchmark_f11(const std::vector<long double>& values) {
	const long double a = 0.125l;
	long double f_sum = 0.0l, s_sum = 0.0l;
	size_t size = values.size();
	for (int i = 0; i < size; i++) {
		f_sum += pow(values[i], 2.0l);
		s_sum += values[i];
	}
	return pow(pow(f_sum - size, 2.0l), a) + ((0.5l * f_sum + s_sum) / size) + 0.5l;
}

long double benchmark_f12(const std::vector<long double>& values) {
	long double f_sum = 0.0l, s_sum = 0.0l;
	size_t size = values.size();
	for (int i = 0; i < size; i++) {
		f_sum += pow(values[i], 2.0l);
		s_sum += values[i];
	}
	return sqrt(pow(f_sum, 2.0l) - pow(s_sum, 2.0l)) + ((0.5l * f_sum + s_sum) / size) + 0.5l;
}

long double benchmark_f13(const std::vector<long double>& values) {
	long double sum = 0.0l;
	size_t size = values.size() - 1;
	for (int i = 0; i < size; i++) {
		sum += 100.0l * pow(values[i + 1] - pow(values[i], 2.0l), 2.0l) + pow(values[i] - 1.0l, 2.0l);
	}
	return sum;
}

long double benchmark_f14(const std::vector<long double>& values) {
	long double sum = 0.0l;
	size_t size = values.size();
	for (int i = 0; i < size; i++) {
		sum += pow(1000000.0l, i / (size - 1)) * pow(values[i], 2.0l);
	}
	return sum;
}

long double benchmark_f15(const std::vector<long double>& values) {
	long double sum = 0.0l;
	size_t size = values.size();
	for (int i = 1; i < size; i++) {
		sum += pow(values[i], 2.0l);
	}
	return 1000000.0l * pow(values[0], 2.0l) + sum;
}

long double benchmark_f16(const std::vector<long double>& values) {
	long double sum = 0.0l;
	size_t size = values.size();
	for (int i = 1; i < size; i++) {
		sum += 1000000.0l * pow(values[i], 2.0l);
	}
	return pow(values[0], 2.0l) + sum;
}

long double benchmark_f17(const std::vector<long double>& values) {
	const long double b = 0.5l;
	long double sum = 0.0l;
	size_t size = values.size();

	for (size_t i = 0; i < size; ++i) {
		long double temp_sum = 0.0;
		long double pow_i = pow(i + 1, i + 1);

		for (size_t j = 0; j < size; ++j) {
			temp_sum += (pow(j + 1, i + 1) + b) * (pow(values[j] / (j + 1), i + 1) - 1.0l);
		}

		sum += pow(temp_sum, 2.0l);
	}

	return sum;
}

/*long double benchmark_f17(const std::vector<long double>& values) {
	const long double b = 0.5l;
	long double sum = 0.0l;
	size_t size = values.size();

	for (int i = 0; i < size; i++) {
		long double temp_sum = 0.0;
		for (int j = 0; j < size; j++) {
			temp_sum += (pow(j + 1, i + 1) + b) * (pow(values[j] / (j + 1), i + 1) - 1.0l);
		}
		sum += pow(temp_sum, 2.0l);
	}
	return sum;
}*/

long double benchmark_f18(const std::vector<long double>& values) {
	long double sum = 0.0l, s = 0.0l;
	size_t size = values.size() - 1;
	for (int i = 0; i < size; i++) {
		s = sqrt(pow(values[i], 2.0l) + pow(values[i + 1], 2.0l));
		sum += sqrt(s) + sqrt(s) * pow(sin(50.0l * pow(s, 0.2l)), 2.0l);
	}
	return pow((1.0l / size) * sum, 2.0l);
}

long double benchmark_f19(const std::vector<long double>& values) {
	long double sum = 0.0l, g = 0.0l, temp = 0.0;
	size_t size = values.size() - 1;
	for (int i = 0; i < size; i++) {
		temp = pow(values[i], 2.0l) + pow(values[i + 1], 2.0l);
		g = (pow(sin(sqrt(temp)), 2.0l) - 0.5l) / pow(1.0l + 0.001 * temp, 2.0l) + 0.5l;
		sum += g;
	}
	temp = pow(values[size - 1], 2.0l) + pow(values[0], 2.0l);
	g = (pow(sin(sqrt(temp)), 2.0l) - 0.5l) / pow(1.0l + 0.001 * temp, 2.0l) + 0.5l;
	return sum + g;
}

long double benchmark_f20(const std::vector<long double>& values) {
	long double sum = 0.0l;
	size_t size = values.size();
	for (int i = 0; i < size; i++) {
		sum += pow(values[i], 2.0l) * (size - i);
	}
	return sum;
}

long double benchmark_f21(const std::vector<long double>& values) {
	long double sum = 0.0l;
	size_t size = values.size();
	for (int i = 0; i < size; i++) {
		sum += values[i] * sin(sqrt(fabs(values[i])));
	}
	return -sum + 418.9828872724337l * size;
}

long double benchmark_f22(const std::vector<long double>& values) {
	long double sum = 0.0l;
	size_t size = values.size();
	for (int i = 0; i < size; i++) {
		sum += pow(fabs(values[i]), 2.0l + (4.0l * (i / (size - 1))));
	}
	return sqrt(sum);
}

long double benchmark_f23(const std::vector<long double>& values) {
	long double f_sum = 0.0l, s_sum = 0.0l;
	size_t size = values.size();
	for (int i = 0; i < size; i++) {
		f_sum += fabs(values[i]);
		s_sum += sin(pow(values[i], 2.0l));
	}
	return f_sum * pow(E, -s_sum);
}

long double benchmark_f24(const std::vector<long double>& values) {
	long double max = -INFINITY;
	size_t size = values.size();
	for (int i = 0; i < size; i++) {
		long double temp = fabs(values[i]);
		if (temp > max) max = temp;
	}
	return max;
}

long double benchmark_f25(const std::vector<long double>& values) {
	long double sum = 0.0l, product = 1.0l;
	size_t size = values.size();
	for (int i = 0; i < size; i++) {
		sum += fabs(values[i]);
		product *= fabs(values[i]);
	}
	return sum + product;
}

long double benchmark_f26(const std::vector<long double>& values) {
	long double sum = 0.0l;
	size_t size = values.size();
	for (int i = 0; i < size; i++) {
		sum += pow(values[i], 2.0l);
	}
	return 1.0l - cos(2.0l * PI * sqrt(sum)) + 0.1l * sqrt(sum);
}

long double benchmark_f27(const std::vector<long double>& values) {
	const long double d = 2.0l, a = 0.1l;
	long double sum = 0.0l;
	size_t size = values.size();
	for (int i = 1; i < size; i++) {
		sum += fabs(values[i]);
	}
	return fabs(values[0]) + d * pow(sum, a);
}

long double benchmark_f28(const std::vector<long double>& values) {
	long double f_sum = 0.0l, s_sum = 0.0l;
	size_t size = values.size();
	for (int i = 0; i < size; i++) {
		f_sum += pow(values[i], 2.0l);
		s_sum += 0.5l * (1.0l + i) * values[i];
	}
	return f_sum + pow(s_sum, 2.0) + pow(s_sum, 4.0);
}

long double benchmark_f29(const std::vector<long double>& values) {
	long double f_sum = 0.0l, s_sum = 0.0l, product = 1.0l;
	size_t size = values.size();
	for (int i = 0; i < size; i++) {
		f_sum += pow(values[i] / 15.0l, 10.0l);
		s_sum += pow(values[i], 2.0l);
		product *= pow(cos(values[i]), 2.0l);
	}
	return 10000.0l * (1.0l + (pow(E, -f_sum) - 2.0l * pow(E, -s_sum)) * product);
}

long double benchmark_f30(const std::vector<long double>& values) {
	long double f_sum = 0.0l, s_sum = 0.0l, t_sum = 0.0l;
	size_t size = values.size();
	for (int i = 0; i < size; i++) {
		f_sum += pow(sin(values[i]), 2.0l);
		s_sum += pow(values[i], 2.0l);
		t_sum += pow(sin(sqrt(fabs(values[i]))), 2.0l);
	}
	return 10000.0l * (1.0l + (f_sum - pow(E, -s_sum)) * pow(E, -t_sum));
}

std::vector<benchmark_function> get_all_functions() {
	return {
		benchmark_f1,
		benchmark_f2,
		benchmark_f3,
		benchmark_f4,
		benchmark_f5,
		benchmark_f6,
		benchmark_f7,
		benchmark_f8,
		benchmark_f9,
		benchmark_f10,
		benchmark_f11,
		benchmark_f12,
		benchmark_f13,
		benchmark_f14,
		benchmark_f15,
		benchmark_f16,
		benchmark_f17,
		benchmark_f18,
		benchmark_f19,
		benchmark_f20,
		benchmark_f21,
		benchmark_f22,
		benchmark_f23,
		benchmark_f24,
		benchmark_f25,
		benchmark_f26,
		benchmark_f27,
		benchmark_f28,
		benchmark_f29,
		benchmark_f30,
	};
}

std::vector<long double> get_borders(int dimension) {
	return {
		100.0l,
		100.0l,
		10.0l,
		20.0l,
		5.12l,
		0.5l,
		10.0l,
		32.768l,
		100.0l,
		5.12l,

		20.0l,
		15.0l,
		10.0l,
		100.0l,
		100.0l,
		100.0l,
		static_cast<long double>(dimension),
		100.0l,
		100.0l,
		100.0l,

		500.0l,
		10.0l,
		2.0l * PI,
		100.0l,
		100.0l,
		20.0l,
		100.0l,
		10.0l,
		20.0l,
		100.0l
	};
}