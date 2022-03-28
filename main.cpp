#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

using namespace std;

float count_distrib(float * hists, unsigned num_hists, float const * values, unsigned num_values, float dev_from) {
    auto min = numeric_limits<float>::max();
    auto max = numeric_limits<float>::min();
    for (int j = 0; j < num_values; j++) {
        if (values[j] < min) {
            min = values[j];
        }
        if (values[j] > max) {
            max = values[j];
        }
    }
    auto dv = (max - min) / (float) (num_hists - 1);
    for (int j = 0; j < num_values; j++) {
        hists[(int) round((values[j] - min) / dv)] += 1;
    }
    return (dev_from - min) / dv;
}

float pdf(float * hists, unsigned num_hists, float const * values, unsigned num_values, float dev_from) {
    auto value_idx = count_distrib(hists, num_hists, values, num_values, dev_from);
    for (unsigned j = 0; j < num_hists; j++) {
        hists[j] /= (float) num_values;
    }
    return value_idx;
}

void print_int_distrib(float const * hists, unsigned num_hists) {
    for (int i = 0; i < num_hists; i++) {
        for (int j = 0; j < (int) hists[i]; j++) {
            cout << "|";
        }
        cout << endl;
    }
}

float dispersion(float const * values, unsigned num_values, float dev_from) {
    float sum = 0;
    for (unsigned i = 0; i < num_values; i++) {
        auto delta = values[i] - dev_from;
        sum = fma(delta, delta, sum);
    }
    return sum / (float) num_values;
}

unsigned read_data(float * values) {
    ifstream file("../raw.csv");
    int num_values = 0;
    if (file.is_open()) {
        while (file >> values[num_values]) {
            num_values++;
        }
    }
    file.close();
    return num_values;
}

float maxwell_pdf(float * hists, unsigned num_hists, float scale, float T) {
    auto max_idx = (float) (num_hists - 1);
    auto max_dev_idx = max_idx / 2;
    auto dv = scale / max_idx;
    for (int i = 0; i < num_hists; i++) {
        auto v = ((float) (2 * i) - max_idx) * dv;
        hists[i] = (exp(-(v * v) / T) / sqrt(T * (float) numbers::pi)) * dv;
    }
    return max_idx / 2;
}

float mean_direct(float const * psi, float const * pdf, unsigned size) {
    float sum = 0;
    for (auto i = 0; i < size; i++) {
        sum += pdf[i] * psi[i];
    }
    return sum;
}

float mean_fma(float const * psi, float const * pdf, unsigned size) {
    float sum = 0;
    for (auto i = 0; i < size; i++) {
        sum = fma(pdf[i], psi[i], sum);
    }
    return sum;
}

float mean_recurse(float const * psi, float const * pdf, unsigned s) {
    if (s == 1) return psi[0] * pdf[0];
    auto h = s / 2;
    return mean_recurse(psi, pdf, h) + mean_recurse(&psi[h], &pdf[h], s - h);
}

float mean_kahan(float const * psi, float const * pdf, unsigned s) {
    float sum = 0;
    float err = 0;
    for (auto i = 0; i < s; i++) {
        auto x = fma(pdf[i], psi[i], -err);
        auto sum_new = sum + x;
        err = (sum_new - sum) - x;
        sum = sum_new;
    }
    return sum;
}

void test_mean(float const * pdf, unsigned size, float v0, float T) {
    cout << "mean " << v0 << " from 0 to " << size - 1 << endl;
    auto psi = new float[size];
    for (int i = 0; i < size; i++) {
        psi[i] = 1;
//        psi[i] = abs(v0 - (float) i);
    }
    cout << "T: " << T << endl;
    cout << "expected mean: " << sqrt(T / numbers::pi) << endl;
    cout << "mean directly: " << mean_direct(psi, pdf, size) << endl;
    cout << "mean by fma: " << mean_fma(psi, pdf, size) << endl;
    cout << "mean by recursion: " << mean_recurse(psi, pdf, size) << endl;
    cout << "mean by Kahan: " << mean_kahan(psi, pdf, size) << endl;
    delete[] psi;
}

float test_data(unsigned num_hists) {
    cout << " --- testing raw.csv --- " << endl;

    float values[273];
    auto num_values = read_data(values);
    float expected = 500;
    auto T = 2 * dispersion(values, num_values, expected);

    auto hists = new float[num_hists] {0};
    auto expected_idx = pdf(hists, num_hists, values, num_values, expected);
    test_mean(hists, num_hists, expected_idx, T);
    delete[] hists;

    cout << endl;
    return T;
}

void test_maxwell(unsigned num_hists, float T) {
    cout << " --- testing maxwell distribution --- " << endl;

    auto hists = new float[num_hists];
    auto value_idx = maxwell_pdf(hists, num_hists, 100, T);
    test_mean(hists, num_hists, value_idx, T);

    cout << endl;
}

int main() {
    auto num_hists = 20;
    auto T = test_data(num_hists);
    test_maxwell(num_hists, T);
}