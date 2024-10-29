#pragma once

#include "halton_sequence.h"

std::vector<int> sieve_eratosthenes(const int size) {
    int limit;
    if (size > 6) limit = static_cast<int>(size * (log(size) + log(log(size))));
    else limit = 15;
    std::vector<bool> number(limit, true);
    number[0] = false;
    number[1] = false;
    std::vector<int> primes;

    for (int i = 2; i < limit; i++) {
        if (number[i]) primes.push_back(i);
        if (primes.size() == size) break;
        for (int j = i * i; j < limit; j += i) {
            number[j] = false;
        }
    }
    return primes;
}

std::vector<long double> halton_sequence(const int index, const std::vector<int>& bases) {
    std::vector<long double> point(bases.size());

    for (size_t i = 0; i < bases.size(); i++) {
        int base = bases[i];
        long double position = 1.0;
        long double result = 0.0;
        int current = index + 1;

        while (current > 0) {
            position /= base;
            result += position * (current % base);
            current /= base;
        }

        point[i] = result;
    }
    return point;
}