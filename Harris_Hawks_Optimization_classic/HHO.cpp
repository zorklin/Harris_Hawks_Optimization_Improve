#pragma once

#include "halton_sequence.h"
#include "HHO.h"
#include <algorithm>

constexpr long double beta = 1.5l;
constexpr long double PI = 3.141592653589793238462643383279502884l;
constexpr long double E = 2.71828182845904523536028747135266249775724709369995l;
constexpr long double POWER_STAGNATION = 15.0l;

std::mt19937 initialization_rand() {
    std::random_device rand_number;
    return std::mt19937(rand_number());
}

long double rand_real(std::mt19937& generator, const long double low, const long double high) {
    std::uniform_real_distribution<long double> range(low, high);
    return range(generator);
}

int rand_int(std::mt19937& generator, const int low, const int high) {
    std::uniform_int_distribution<int> range(low, high);
    return range(generator);
}

HHO initialization_hawks(const int T, const int size, const long double max, const long double min, const int dimension, long double (*function)(const std::vector<long double>&)) {
    HHO hho;
    hho.T = T;
    hho.size = size;
    hho.max = max;
    hho.min = min;
    hho.dimension = dimension;
    hho.fitness = function;
    hho.hawks.resize(size);

    const std::vector<int> bases = sieve_eratosthenes(dimension);
    for (int i = 0; i < size; i++) {
        hawk temp = halton_sequence(i, bases);
        hho.hawks[i] = temp * (max - min) + min;
        hho.hawks[i].fitness = hho.fitness(hho.hawks[i].X);
    }

    return hho;
}
 
void border_correction(HHO& hho, const int index) {
    for (int i = 0; i < hho.dimension; i++) {
        hho.hawks[index].X[i] = std::max(std::min(hho.hawks[index].X[i], hho.max), hho.min);
    }
}

hawk best_hawk(const HHO& hho) {
    hawk best_hawk = hho.hawks[0];
    for (int i = 1; i < hho.size; i++) {
        if (hho.hawks[i].fitness < best_hawk.fitness) {
            best_hawk = hho.hawks[i];
        }
    }
    return best_hawk;
}

hawk calculate_Xm(const HHO& hho) {
    hawk average_hawk(0.0l, hho.dimension);
    for (int i = 0; i < hho.size; i++) {
        average_hawk = average_hawk + hho.hawks[i];
    }
    return average_hawk / hho.size;
}

long double calculate_energy(const HHO& hho, const int t, std::mt19937& generator) {
    long double E = pow(cos(PI * pow(static_cast<long double>(t) / hho.T, 2.0l / 3.0l)), 5.0l) + 1.0l;
    return E * rand_real(generator, -1.0, 1.0);
}

void very_high_energy(HHO& hho, int index, hawk leader_hawk, hawk averege_hawk, std::mt19937& generator) {
    long double q = rand_real(generator, 0.0l, 1.0l);
    if (q >= 0.5) {
        int rand_index = rand_int(generator, 0, hho.size - 1);
        hawk rand_hawk = hho.hawks[rand_index];
        hawk hawk_deviation = (rand_hawk - hho.hawks[index] * 2.0l * rand_real(generator, 0.0l, 1.0l)).module_number();
        hho.hawks[index] = rand_hawk - hawk_deviation * rand_real(generator, 0.0l, 1.0l);
    }
    else {
        long double correction = hho.min + (hho.max - hho.min) * rand_real(generator, 0.0l, 1.0l);
        hho.hawks[index] = leader_hawk - averege_hawk - correction * rand_real(generator, 0.0l, 1.0l);
    }
}

void high_energy_high_chance(HHO& hho, const int index, const long double E, hawk leader_hawk, std::mt19937& generator) {
    long double J = rand_real(generator, 0.0l, 1.0l);
    hawk diff_hawk = leader_hawk - hho.hawks[index];
    hho.hawks[index] = diff_hawk - (leader_hawk * J - hho.hawks[index]).module_number() * E;
}

void low_energy_high_chance(HHO& hho, const int index, const long double E, hawk leader_hawk) {
    hawk diff_hawk = leader_hawk - hho.hawks[index];
    hho.hawks[index] = leader_hawk - diff_hawk.module_number() * E;
}

hawk better(const HHO& hho, const hawk f_hawk, const hawk s_hawk) {
    if (hho.fitness(f_hawk.X) < hho.fitness(s_hawk.X)) return f_hawk;
    return s_hawk;
}

long double levy_flight(std::mt19937& generator) {
    const long double u = rand_real(generator, 0.0l, 1.0l), v = rand_real(generator, 0.0l, 1.0l);
    const long double numerator = tgamma(1.0l + beta) * sin((PI * beta) / 2.0l);
    const long double denominator = tgamma((1.0l + beta) / 2.0l) * beta * pow(2.0l, (beta - 1.0l) / 2.0l);
    const long double sigma = pow(numerator / denominator, 1.0 / beta);
    return 0.01l * ((u * sigma) / pow(fabs(v), 1.0l / beta));
}

void high_energy_low_chance(HHO& hho, const int index, const long double E, hawk leader_hawk, std::mt19937& generator) {
    long double J = rand_real(generator, 0.0l, 1.0l);
    hawk y_hawk = leader_hawk - (leader_hawk * J - hho.hawks[index]).module_number() * E;
    hawk rand_hawk(0.0l, hho.dimension);
    for (int i = 0; i < hho.dimension; i++) {
        rand_hawk.X[i] = rand_real(generator, hho.min, hho.max);
    }
    hawk z_hawk = y_hawk + rand_hawk * levy_flight(generator);
    hho.hawks[index] = better(hho, better(hho, y_hawk, hho.hawks[index]), z_hawk);
}

void low_energy_low_chance(HHO& hho, const int index, const long double E, hawk leader_hawk, hawk averege_hawk, std::mt19937& generator) {
    long double J = rand_real(generator, 0.0l, 1.0l);
    hawk rand_hawk(0.0l, hho.dimension);
    for (int i = 0; i < hho.dimension; i++) {
        rand_hawk.X[i] = rand_real(generator, hho.min, hho.max);
    }
    hawk y_hawk = leader_hawk - (leader_hawk * J - averege_hawk).module_number() * E;
    hawk z_hawk = y_hawk + rand_hawk * levy_flight(generator);
    hho.hawks[index] = better(hho, better(hho, y_hawk, hho.hawks[index]), z_hawk);
}

void elite_opposition_based_learning(HHO& hho, std::mt19937& generator) {
    hawk min_hawk(INFINITY, hho.dimension);
    hawk max_hawk(-INFINITY, hho.dimension);

    for (int i = 0; i < hho.dimension; i++) {
        for (int j = 0; j < hho.size; j++) {
            min_hawk.X[i] = std::min(min_hawk.X[i], hho.hawks[j].X[i]);
            max_hawk.X[i] = std::max(max_hawk.X[i], hho.hawks[j].X[i]);
        }
    }

    for (int i = 0; i < hho.size; i++) {
        hawk temp_hawk = hho.hawks[i];
        for (int j = 0; j < hho.dimension; j++) {
            temp_hawk.X[j] = rand_real(generator, min_hawk.X[j], max_hawk.X[j]);
            temp_hawk.fitness = hho.fitness(temp_hawk.X);
            if (temp_hawk.fitness < hho.hawks[i].fitness) hho.hawks[i] = temp_hawk;
        }
    }
}

long double gaissian_rand(const long double expectations, const long double deviation, std::mt19937& generator) {
    std::normal_distribution<> dist(expectations, deviation);
    return dist(generator);
}

void gaussian_walk_learning(HHO& hho, const int t, hawk leader_hawk, std::mt19937& generator) {
    hho.hawks[0] = leader_hawk;
    for (int i = 1; i < hho.size; i++) {
        int index = rand_int(generator, 0, hho.size - 1);
        hawk rand_hawk = hho.hawks[index];
        hawk tau_hawk = ((hho.hawks[i] - rand_hawk) * cos(PI / 2.0l * pow(static_cast<long double>(t) / hho.T, 2.0l))).module_number();
        for (int j = 0; j < hho.dimension; j++) {
            long double deviation = std::max(tau_hawk.X[j], std::numeric_limits<long double>::min());
            hho.hawks[i].X[j] = gaissian_rand(hho.hawks[i].X[j], deviation, generator);
        }
    }
}

std::vector<long double> harris_hawks_optimazation(const int T, const int size, const long double max, const long double min, const int dimension, long double (*function)(const std::vector<long double>&)) {
    std::mt19937 generator = initialization_rand();
    HHO hho = initialization_hawks(T, size, max, min, dimension, function);
    hawk best_solution(INFINITY, dimension);
    best_solution.fitness = INFINITY;
    long double stagnation = 0;
    for (int t = 0; t < T; t++) {
        elite_opposition_based_learning(hho, generator);
        hawk leader_hawk = best_hawk(hho);
        hawk average_hawk = calculate_Xm(hho);
        if (stagnation > POWER_STAGNATION) {
            gaussian_walk_learning(hho, t, leader_hawk, generator);
            for (int i = 0; i < size; i++) {
                very_high_energy(hho, i, leader_hawk, average_hawk, generator);
                border_correction(hho, i);
                hho.hawks[i].fitness = hho.fitness(hho.hawks[i].X);
            }
            stagnation = 0;
        }
        else {
            for (int i = 0; i < size; i++) {
                long double E = calculate_energy(hho, t, generator);
                if (fabs(E) >= 1.0l) very_high_energy(hho, i, leader_hawk, average_hawk, generator);
                else {
                    long double q = rand_real(generator, 0.0l, 1.0l);
                    if (fabs(E) >= 0.5l && q >= 0.5l) high_energy_high_chance(hho, i, E, leader_hawk, generator);
                    else if (fabs(E) < 0.5l && q >= 0.5l) low_energy_high_chance(hho, i, E, leader_hawk);
                    else if (fabs(E) >= 0.5l && q < 0.5l) high_energy_low_chance(hho, i, E, leader_hawk, generator);
                    else if (fabs(E) < 0.5l && q < 0.5l) low_energy_low_chance(hho, i, E, leader_hawk, average_hawk, generator);
                }
                border_correction(hho, i);
                hho.hawks[i].fitness = hho.fitness(hho.hawks[i].X);
            }
        }
        hawk iter_solution = leader_hawk;
        stagnation += cos(PI / 3.0l * cbrt(static_cast<long double>(t) / T));
        if (iter_solution.fitness < best_solution.fitness) {
            stagnation = 0;
            best_solution = iter_solution;
        }
    }
    return best_solution.X;
}