#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include "halton_sequence.h"
#include "benchmark_function.h"

#define beta 1.5l

using namespace std;

struct hawk {
    vector<long double> X;
    long double fitness;

    hawk() : fitness(0.0) {}
    hawk(long double value, size_t size) : X(size, value), fitness(0) {}
    hawk(const vector<long double>& values) : X(values), fitness(0) {}

    hawk operator+(const hawk& other) const {
        hawk new_hawk = *this;
        transform(X.begin(), X.end(), other.X.begin(), new_hawk.X.begin(), plus<>());
        return new_hawk;
    }

    hawk operator-(const hawk& other) const {
        hawk new_hawk = *this;
        transform(X.begin(), X.end(), other.X.begin(), new_hawk.X.begin(), minus<>());
        return new_hawk;
    }

    hawk operator+(long double scalar) const {
        hawk new_hawk = *this;
        transform(X.begin(), X.end(), new_hawk.X.begin(), [scalar](long double value) { return value + scalar; });
        return new_hawk;
    }

    hawk operator-(long double scalar) const {
        hawk new_hawk = *this;
        transform(X.begin(), X.end(), new_hawk.X.begin(), [scalar](long double value) { return value - scalar; });
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

    hawk fabs() const {
        hawk new_hawk = *this;
        transform(X.begin(), X.end(), new_hawk.X.begin(), [](long double value) { return std::fabs(value); });
        return new_hawk;
    }

    friend ostream& operator<<(ostream& os, const hawk& h) {
        os << "Fitness: " << h.fitness << ", X: (";
        for (size_t i = 0; i < h.X.size(); ++i) {
            os << h.X[i];
            if (i < h.X.size() - 1) {
                os << ", ";
            }
        }
        os << ")";
        return os;
    }
};

struct HHO {
    int T, size, dimension;
    long double x_max, x_min;
    vector<hawk> hawks;
};

mt19937 initialization_rand() {
    random_device rand_number;
    return mt19937(rand_number());
}

long double rand_real(mt19937& generator, long double low, long double high) {
    uniform_real_distribution<long double> range(low, high);
    return range(generator);
}

int rand_int(mt19937& generator, int low, int high) {
    uniform_int_distribution<int> range(low, high);
    return range(generator);
}

long double fitness(HHO* hho, hawk hawk) {
    long double sum = 0.0, prodact = 1.0;
    for (int i = 0; i < hho->dimension; i++) {
        sum += pow(hawk.X[i], 2.0);
        prodact *= cos(hawk.X[i] / sqrt(i + 1));
        //sum += pow(hawk.X[i], 2.0);
        //sum += fabs(pow(hawk.X[i], 5.0) - 3.0 * pow(hawk.X[i], 4.0) + 4.0 * pow(hawk.X[i], 3.0) + 2.0 * pow(hawk.X[i], 2.0) - 10 * pow(hawk.X[i], 5.0) - 4.0);
    }
    //return sum;
    return (1.0 / 4000.0) * sum - prodact + 1.0;
}

HHO* initialization_hawks(int T, int size, long double x_max, long double x_min, int dimension) {
    HHO* hho = new HHO();
    if (hho == nullptr) exit(EXIT_FAILURE);
    hho->T = T;
    hho->size = size;
    hho->x_max = x_max;
    hho->x_min = x_min;
    hho->dimension = dimension;
    hho->hawks.resize(size);

    vector<int> bases = sieve_eratosthenes(dimension);
    for (int i = 0; i < size; i++) {
        hawk temp = halton_sequence(i, bases);
        hho->hawks[i] = temp * (x_max - x_min) + x_min;
        hho->hawks[i].fitness = fitness(hho, hho->hawks[i]);
    }

    return hho;
}

hawk best_hawk(HHO* hho) {
    hawk best_hawk = hho->hawks[0];
    for (int i = 0; i < hho->size; i++) {
        if (hho->hawks[i].fitness < best_hawk.fitness) { 
            best_hawk = hho->hawks[i];
        }
    }
    return best_hawk;
}

long double calculate_energy(HHO* hho, int t, mt19937& generator) {
    long double E;
    if (t <= hho->T / 2.0) E = cos(_Pi * (static_cast<long double>(t) / hho->T + 0.5) + 2.0);
    else E = cos(_Pi * pow((static_cast<long double>(t) / hho->T - 0.5), 1.0 / 3.0));
    return E * (2.0 * rand_real(generator, 0.0, 1.0) - 1.0);
}

hawk calculate_Xm(HHO* hho) {
    hawk average_hawk(0.0, hho->dimension);
    for (int i = 0; i < hho->size; i++) {
        average_hawk = average_hawk + hho->hawks[i];
    }
    return average_hawk / hho->size;
}

void very_high_energy(HHO* hho, int index, mt19937& generator) {
    long double q = rand_real(generator, 0.0, 1.0);
    if (q >= 0.5) {
        int rand_index = rand_int(generator, 0, hho->size - 1);
        hawk rand_hawk = hho->hawks[rand_index];
        hawk hawk_deviation = (rand_hawk - hho->hawks[index] * 2.0 * rand_real(generator, 0.0, 1.0)).fabs();
        hho->hawks[index] = rand_hawk - hawk_deviation * rand_real(generator, 0.0, 1.0);
    }
    else {
        hawk averege_hawk = calculate_Xm(hho), leader_hawk = best_hawk(hho);
        long double correction = hho->x_min + (hho->x_max - hho->x_min) * rand_real(generator, 0.0, 1.0);
        hho->hawks[index] = leader_hawk - averege_hawk - correction * rand_real(generator, 0.0, 1.0);
    }
}

void high_energy_high_chance(HHO* hho, int index, long double E, mt19937& generator) {
    long double J = rand_real(generator, 0.0, 1.0);
    hawk leader_hawk = best_hawk(hho);
    hawk diff_hawk = leader_hawk - hho->hawks[index];
    hho->hawks[index] = diff_hawk - (leader_hawk * J - hho->hawks[index]).fabs() * E;
}

void low_energy_high_chance(HHO* hho, int index, long double E) {
    hawk leader_hawk = best_hawk(hho);
    hawk diff_hawk = leader_hawk - hho->hawks[index];
    hho->hawks[index] = leader_hawk - diff_hawk.fabs() * E;
}

hawk better(HHO* hho, hawk f_hawk, hawk s_hawk) {
    if (fitness(hho, f_hawk) < fitness(hho, s_hawk)) return f_hawk;
    return s_hawk;
}

long double levy_flight(HHO* hho, mt19937& generator) {
    long double u = rand_real(generator, 0.0, 1.0), v = rand_real(generator, 0.0, 1.0);
    long double numerator = tgamma(1.0 + beta) * sin((_Pi * beta) / 2.0);
    long double denominator = tgamma((1.0 + beta) / 2.0) * beta * pow(2.0, (beta - 1.0) / 2.0);
    long double sigma = pow(numerator / denominator, 1.0 / beta);
    return 0.01 * ((u * sigma) / pow(fabs(v), 1.0 / beta));
}

void high_energy_low_chance(HHO* hho, int index, long double E, mt19937& generator) {
    long double J = rand_real(generator, 0.0, 1.0);
    hawk leader_hawk = best_hawk(hho);
    hawk y_hawk = leader_hawk - (leader_hawk * J - hho->hawks[index]).fabs() * E;
    hawk rand_hawk(0.0, hho->dimension);
    for (int i = 0; i < hho->dimension; i++) {
        rand_hawk.X[i] = rand_real(generator, hho->x_min, hho->x_max);
    }
    hawk z_hawk = y_hawk + rand_hawk * levy_flight(hho, generator);
    hho->hawks[index] = better(hho, y_hawk, z_hawk);
}

void low_energy_low_chance(HHO* hho, int index, long double E, mt19937& generator) {
    long double J = rand_real(generator, 0.0, 1.0);
    hawk leader_hawk = best_hawk(hho);
    hawk averege_hawk = calculate_Xm(hho);
    hawk rand_hawk(0.0, hho->dimension);
    for (int i = 0; i < hho->dimension; i++) {
        rand_hawk.X[i] = rand_real(generator, hho->x_min, hho->x_max);
    }
    hawk y_hawk = leader_hawk - (leader_hawk * J - averege_hawk).fabs() * E;
    hawk z_hawk = y_hawk + rand_hawk * levy_flight(hho, generator);
    hho->hawks[index] = better(hho, y_hawk, z_hawk);
}

long double gaissian_rand(long double expectations, long double deviation, mt19937& generator) {
    normal_distribution<> dist(expectations, deviation);
    return dist(generator);
}

void gaussian_walk_learning(HHO* hho, int t, mt19937& generator) {
    for (int i = 0; i < hho->size; i++) {
        int index = rand_int(generator, 0, hho->size - 1);
        hawk rand_hawk = hho->hawks[index];
        hawk tau_hawk = ((hho->hawks[i] - rand_hawk) * cos(_Pi / 2.0 * pow(t / hho->T, 2.0))).fabs();
        for (int j = 0; j < hho->dimension; j++) {
            long double deviation = max(tau_hawk.X[j], numeric_limits<long double>::min());
            hho->hawks[i].X[j] = gaissian_rand(hho->hawks[i].X[j], deviation, generator);
        }
    }
}

void delete_data(HHO* hho) {
    if (hho != nullptr) {
        hho->hawks.clear();
        delete hho;
    }
}

void harris_hawks_optimazation(int T, int size, long double x_max, long double x_min, int dimension) {
    mt19937 generator = initialization_rand();
    HHO* hho = initialization_hawks(T, size, x_max, x_min, dimension);
    hawk best_solution(INFINITY, dimension);
    best_solution.fitness = INFINITY;
    int stagnation = 0;
    for (int t = 0; t < T; t++) {
        //_Pi / 2.0 * pow(t / T, 1.0 / 2.0))
        if (stagnation >= 20.0) {
            gaussian_walk_learning(hho, t, generator);
            long double E = calculate_energy(hho, t, generator);
            for (int i = 0; i < size; i++) {
                very_high_energy(hho, i, generator);
                high_energy_low_chance(hho, i, pow(E, 2.0), generator);
                hho->hawks[i].fitness = fitness(hho, hho->hawks[i]);
            }
            stagnation = 0;
        }
        else {
            long double E = calculate_energy(hho, t, generator);
            for (int i = 0; i < size; i++) {
                if (fabs(E) >= 1) very_high_energy(hho, i, generator);
                else {
                    long double q = rand_real(generator, 0.0, 1.0);
                    if (fabs(E) >= 0.5 && q >= 0.5) high_energy_high_chance(hho, i, E, generator);
                    else if (fabs(E) < 0.5 && q >= 0.5) low_energy_high_chance(hho, i, E);
                    else if (fabs(E) >= 0.5 && q < 0.5) high_energy_low_chance(hho, i, E, generator);
                    else if (fabs(E) < 0.5 && q < 0.5) low_energy_low_chance(hho, i, E, generator);
                }
                hho->hawks[i].fitness = fitness(hho, hho->hawks[i]);
            }
        }
        hawk iter_solution = best_hawk(hho);
        stagnation += cos(_Pi / 2.0 * pow(t / T, 3.0 / 2.0));
        if (iter_solution.fitness < best_solution.fitness) { 
            stagnation = 0;
            best_solution = iter_solution;
        }
        cout << "iteration: " << t << ", solution: " << best_solution << endl;
    }
    cout << "best solution: " << best_solution << endl;
    delete_data(hho);
}

int main(){
    int T = 500, size = 100, dimension = 10;
    long double x_max = 5.0, x_min = -5.0;
    harris_hawks_optimazation(T, size, x_max, x_min, dimension);

    return 0;
}