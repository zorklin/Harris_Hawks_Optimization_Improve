#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>

struct hawk {
    std::vector<long double> X;
    long double fitness;

    hawk() : fitness(0.0) {}
    hawk(long double value, size_t size) : X(size, value), fitness(0) {}
    hawk(const std::vector<long double>& values) : X(values), fitness(0) {}

    hawk operator+(const hawk& other) const {
        hawk new_hawk = *this;

        std::transform(X.begin(), X.end(), other.X.begin(), new_hawk.X.begin(), std::plus<>());
        return new_hawk;
    }

    hawk operator-(const hawk& other) const {
        hawk new_hawk = *this;
        std::transform(X.begin(), X.end(), other.X.begin(), new_hawk.X.begin(), std::minus<>());
        return new_hawk;
    }

    hawk operator+(long double scalar) const {
        hawk new_hawk = *this;
        transform(X.begin(), X.end(), new_hawk.X.begin(), [scalar](long double value) { return value + scalar; });
        return new_hawk;
    }

    hawk operator-(long double scalar) const {
        hawk new_hawk = *this;
        transform(X.begin(), X.end(), new_hawk.X.begin(), [scalar](double value) { return value - scalar; });
        return new_hawk;
    }

    hawk operator*(long double scalar) const {
        hawk new_hawk = *this;
        transform(X.begin(), X.end(), new_hawk.X.begin(), [scalar](long double value) { return value * scalar; });
        return new_hawk;
    }

    hawk operator/(long double scalar) const {
        hawk new_hawk = *this;
        transform(X.begin(), X.end(), new_hawk.X.begin(), [scalar](long double value) { return value / scalar; });
        return new_hawk;
    }

    hawk module_number() const {
        hawk new_hawk = *this;
        transform(X.begin(), X.end(), new_hawk.X.begin(), [](long double value) { return std::fabs(value); });
        return new_hawk;
    }
};

struct HHO {
    int T, size, dimension;
    long double max, min;
    long double (*fitness)(const std::vector<long double>&);
    std::vector<hawk> hawks;
};

std::mt19937 initialization_rand();
long double rand_real(std::mt19937& generator, const long double low, const long double high);
int rand_int(std::mt19937& generator, const int low, const int high);
HHO initialization_hawks(const int T, int size, const long double max, const long double min, const int dimension, long double (*function)(const std::vector<long double>&));
void border_correction(HHO& hho, const int index);
hawk best_hawk(const HHO& hho);
long double calculate_energy(const HHO& hho, const int t, std::mt19937& generator);
hawk calculate_Xm(const HHO& hho);
void very_high_energy(HHO& hho, const int index, hawk leader_hawk, hawk averege_hawk, std::mt19937& generator);
void high_energy_high_chance(HHO& hho, const int index, const long double E, hawk leader_hawk, std::mt19937& generator);
void low_energy_high_chance(HHO& hho, const int index, const long double E, hawk leader_hawk);
hawk better(const HHO& hho, const hawk f_hawk, const hawk s_hawk);
long double levy_flight(std::mt19937& generator);
void high_energy_low_chance(HHO& hho, const int index, const long double E, hawk leader_hawk, std::mt19937& generator);
void low_energy_low_chance(HHO& hho, const int index, const long double E, hawk leader_hawk, hawk averege_hawk, std::mt19937& generator);
void elite_opposition_based_learning(HHO& hho, std::mt19937& generator);
long double gaissian_rand(const long double expectations, const long double deviation, std::mt19937& generator);
void gaussian_walk_learning(HHO& hho, const int t, hawk leader_hawk, std::mt19937& generator);
std::vector<long double> harris_hawks_optimazation(const int T, const int size, const long double x_max, const long double x_min, const int dimension, long double (*function)(const std::vector<long double>&));