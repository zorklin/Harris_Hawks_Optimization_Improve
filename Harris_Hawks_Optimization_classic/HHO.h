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
    long double (*fitness)(std::vector<long double>);
    std::vector<hawk> hawks;
};

std::mt19937 initialization_rand();
long double rand_real(std::mt19937& generator, long double low, long double high);
int rand_int(std::mt19937& generator, int low, int high);
HHO* initialization_hawks(int T, int size, long double max, long double min, int dimension, long double (*function)(std::vector<long double>));
void border_correction(HHO* hho, int index);
hawk best_hawk(HHO* hho);
long double calculate_energy(HHO* hho, int t, std::mt19937& generator);
hawk calculate_Xm(HHO* hho);
void very_high_energy(HHO* hho, int index, hawk leader_hawk, hawk averege_hawk, std::mt19937& generator);
void high_energy_high_chance(HHO* hho, int index, long double E, hawk leader_hawk, std::mt19937& generator);
void low_energy_high_chance(HHO* hho, int index, long double E, hawk leader_hawk);
hawk better(HHO* hho, hawk f_hawk, hawk s_hawk);
long double levy_flight(HHO* hho, std::mt19937& generator);
void high_energy_low_chance(HHO* hho, int t, int index, long double E, hawk leader_hawk, std::mt19937& generator);
void low_energy_low_chance(HHO* hho, int t, int index, long double E, hawk leader_hawk, hawk averege_hawk, std::mt19937& generator);
void elite_opposition_based_learning(HHO* hho, std::mt19937& generator);
long double gaissian_rand(long double expectations, long double deviation, std::mt19937& generator);
void gaussian_walk_learning(HHO* hho, int t, hawk leader_hawk, std::mt19937& generator);
void delete_data(HHO* hho);
std::vector<long double> harris_hawks_optimazation(int T, int size, long double x_max, long double x_min, int dimension, long double (*function)(std::vector<long double>));