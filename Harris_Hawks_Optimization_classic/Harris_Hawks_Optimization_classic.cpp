#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include "halton_sequence.h"

//make J logic
//#define J 1.0
#define beta 1.5

using namespace std;

struct hawk {
    vector<double> X;
    double fitness;

    hawk() : fitness(0.0) {}
    hawk(double value, size_t size) : X(size, value), fitness(0) {}
    hawk(const vector<double>& values) : X(values), fitness(0) {}

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

    hawk operator+(double scalar) const {
        hawk new_hawk = *this;
        transform(X.begin(), X.end(), new_hawk.X.begin(), [scalar](double value) { return value + scalar; });
        return new_hawk;
    }

    hawk operator-(double scalar) const {
        hawk new_hawk = *this;
        transform(X.begin(), X.end(), new_hawk.X.begin(), [scalar](double value) { return value - scalar; });
        return new_hawk;
    }

    hawk operator*(double scalar) const {
        hawk new_hawk = *this;
        transform(X.begin(), X.end(), new_hawk.X.begin(), [scalar](double value) { return value * scalar; });
        return new_hawk;
    }

    hawk operator/(double scalar) const {
        hawk new_hawk = *this;
        transform(X.begin(), X.end(), new_hawk.X.begin(), [scalar](double value) { return value / scalar; });
        return new_hawk;
    }


    hawk fabs() const {
        hawk new_hawk = *this;
        transform(X.begin(), X.end(), new_hawk.X.begin(), [](double value) { return std::fabs(value); });
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
    double x_max, x_min;
    vector<hawk> hawks;
};

mt19937 initialization_rand() {
    random_device rand_number;
    return mt19937(rand_number());
}

double rand_real(mt19937& generator, double low, double high) {
    uniform_real_distribution<double> range(low, high);
    return range(generator);
}

int rand_int(mt19937& generator, int low, int high) {
    uniform_int_distribution<int> range(low, high);
    return range(generator);
}

double fitness(HHO* hho, hawk hawk) {
    double sum = 0.0, prodact = 1.0;
    for (int i = 0; i < hho->dimension; i++) {
        //sum += pow(hawk.X[i], 2.0);
        prodact *= cos(hawk.X[i] / sqrt(i + 1));
        //sum += cos(pow(hawk.X[i], 1.0 / pow(hawk.X[i], hawk.X[i])));
        sum += pow(hawk.X[i], 2.0);
        //sum += fabs(pow(hawk.X[i], 5.0) - 3.0 * pow(hawk.X[i], 4.0) + 4.0 * pow(hawk.X[i], 3.0) + 2.0 * pow(hawk.X[i], 2.0) - 10 * pow(hawk.X[i], 5.0) - 4.0);
    }
    //return sum;
    return (1.0 / 4000.0) * sum - prodact + 1.0;
}

HHO* initialization_hawks(int T, int size, double x_max, double x_min, int dimension) {
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
        hho->hawks[i].X = halton_sequence(i, bases);
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

double calculate_energy(HHO* hho, int t, mt19937& generator) {
    double E;
    if (t <= hho->T / 2.0) E = cos(_Pi * (static_cast<double>(t) / hho->T + 0.5) + 2.0);
    else E = cos(_Pi * pow((static_cast<double>(t) / hho->T - 0.5), 1.0 / 3));
    return fabs(E * (2.0 * rand_real(generator, 0.0, 1.0) - 1.0));
}

hawk calculate_Xm(HHO* hho) {
    hawk average_hawk(0.0, hho->dimension);
    for (int i = 0; i < hho->size; i++) {
        average_hawk = average_hawk + hho->hawks[i];
    }
    return average_hawk / hho->size;
}

void very_high_energy(HHO* hho, int index, mt19937& generator) {
    double q = rand_real(generator, 0.0, 1.0);
    if (q >= 0.5) {
        int rand_index = rand_int(generator, 0, hho->size - 1);
        hawk rand_hawk = hho->hawks[rand_index];
        hawk hawk_deviation = (rand_hawk - hho->hawks[index] * 2.0 * rand_real(generator, 0.0, 1.0)).fabs();
        hho->hawks[index] = rand_hawk - hawk_deviation * rand_real(generator, 0.0, 1.0);
    }
    else {
        hawk averege_hawk = calculate_Xm(hho), leader_hawk = best_hawk(hho);
        double correction = rand_real(generator, 0.0, 1.0) * rand_real(generator, 0.0, 1.0) * (hho->x_max - hho->x_min);
        hho->hawks[index] = leader_hawk - averege_hawk - correction;
    }
}

void high_energy_high_chance(HHO* hho, int index, double E, mt19937& generator) {
    double J = rand_real(generator, 0.0, 1.0);
    hawk leader_hawk = best_hawk(hho);
    hawk diff_hawk = leader_hawk - hho->hawks[index];
    hho->hawks[index] = diff_hawk - (leader_hawk * J - hho->hawks[index]).fabs() * E;
}

void low_energy_high_chance(HHO* hho, int index, double E) {
    hawk leader_hawk = best_hawk(hho);
    hawk diff_hawk = leader_hawk - hho->hawks[index];
    hho->hawks[index] = leader_hawk - diff_hawk.fabs() * E;
}

hawk better(HHO* hho, hawk f_hawk, hawk s_hawk) {
    if (fitness(hho, f_hawk) < fitness(hho, s_hawk)) return f_hawk;
    return s_hawk;
}

double levy_flight(HHO* hho, int dimension, mt19937& generator) {
    double u = rand_real(generator, hho->x_min, hho->x_max), v = rand_real(generator, hho->x_min, hho->x_max);
    double numerator = tgamma(1.0 + beta) * sin((_Pi * beta) / 2.0);
    double denominator = tgamma((1.0 + beta) / 2.0) * beta * pow(2.0, (beta - 1.0) / 2.0);
    double sigma = pow(numerator / denominator, 1.0 / beta);
    return 0.01 * ((u * sigma) / pow(fabs(v), 1.0 / beta));
}

void high_energy_low_chance(HHO* hho, int index, double E, mt19937& generator) {
    double J = rand_real(generator, 0.0, 1.0);
    hawk leader_hawk = best_hawk(hho);
    hawk y_hawk = leader_hawk - (leader_hawk * J - hho->hawks[index]).fabs() * E;
    hawk rand_hawk(0.0, hho->dimension);
    for (int i = 0; i < hho->dimension; i++) {
        double rand_num = rand_real(generator, hho->x_min, hho->x_max);
        rand_hawk.X[i] = rand_num;
    }
    hawk z_hawk = y_hawk + rand_hawk * levy_flight(hho, hho->dimension, generator);
    hho->hawks[index] = better(hho, y_hawk, z_hawk);
}

void low_energy_low_chance(HHO* hho, int index, double E, mt19937& generator) {
    double J = rand_real(generator, 0.0, 1.0);
    hawk leader_hawk = best_hawk(hho);
    hawk averege_hawk = calculate_Xm(hho);
    hawk rand_hawk(0.0, hho->dimension);
    for (int i = 0; i < hho->dimension; i++) {
        rand_hawk.X[i] = rand_real(generator, hho->x_min, hho->x_max);
    }
    hawk y_hawk = leader_hawk - (leader_hawk * J - averege_hawk).fabs() * E;
    hawk z_hawk = y_hawk + rand_hawk * levy_flight(hho, hho->dimension, generator);
    hho->hawks[index] = better(hho, y_hawk, z_hawk);
}

void opposition_based_learning(HHO* hho) {
    for (int i = 0; i < hho->size; i++) {
        hawk revers_hawk(0.0, hho->dimension);
        for (int j = 0; j < hho->dimension; j++) {
            revers_hawk.X[j] = hho->x_min + hho->x_max - hho->hawks[i].X[j];
        }
        revers_hawk.fitness = fitness(hho, revers_hawk);
        hho->hawks[i] = better(hho, hho->hawks[i], revers_hawk);
    }
}

double gaissian_function(double expectations, double deviation, mt19937& generator) {
    normal_distribution<> dist(expectations, deviation);
    return dist(generator);
}

void gaussian_walk_learning(HHO* hho, int t, mt19937& generator) {
    for (int i = 0; i < hho->size; i++) {
        int index = rand_int(generator, 0, hho->size - 1);
        hawk rand_hawk = hho->hawks[index];
        hawk tau_hawk = ((hho->hawks[i] - rand_hawk) * cos(_Pi / 2.0 * pow(t / hho->T, 2.0))).fabs();
        for (int j = 0; j < hho->dimension; j++) {
            double deviation = max(tau_hawk.X[j], numeric_limits<double>::min());
            hho->hawks[i].X[j] = gaissian_function(hho->hawks[i].X[j], deviation, generator);
        }
    }
}

void delete_data(HHO* hho) {
    if (hho != nullptr) {
        hho->hawks.clear();
        delete hho;
    }
}

void harris_hawks_optimazation(int T, int size, double x_max, double x_min, int dimension) {
    mt19937 generator = initialization_rand();
    HHO* hho = initialization_hawks(T, size, x_max, x_min, dimension);
    hawk best_solution(INFINITY, dimension);
    best_solution.fitness = INFINITY;
    int stagnation = 1000;
    for (int t = 0; t < T; t++) {
        opposition_based_learning(hho);
        if (stagnation >= pow(T, 0.5)) { 
            gaussian_walk_learning(hho, t, generator);
            stagnation = 0;
        }
        else {
            for (int i = 0; i < size; i++) {
                double r = rand_real(generator, 0.0, 1.0);
                double E = calculate_energy(hho, t, generator);
                if (E >= 1) very_high_energy(hho, i, generator);
                else {
                    double q = rand_real(generator, 0.0, 1.0);
                    if (E >= 0.5 && q >= 0.5) high_energy_high_chance(hho, i, E, generator);
                    else if (E < 0.5 && q >= 0.5) low_energy_high_chance(hho, i, E);
                    else if (E >= 0.5 && q < 0.5) high_energy_low_chance(hho, i, E, generator);
                    else if (E < 0.5 && q < 0.5) low_energy_low_chance(hho, i, E, generator);
                }
                hho->hawks[i].fitness = fitness(hho, hho->hawks[i]);
            }
        }
        hawk iter_solution = best_hawk(hho);
        stagnation++;
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
    int T = 1000, size = 200, dimension = 3;
    double x_max = 1.0, x_min = -1.0;
    harris_hawks_optimazation(T, size, x_max, x_min, dimension);

    return 0;
}