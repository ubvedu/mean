#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

using namespace std;

double count_distrib(double * hists, unsigned num_hists, double const * values, unsigned num_values, double dev_from) {
    auto min = numeric_limits<double>::max();
    auto max = numeric_limits<double>::min();
    for (int j = 0; j < num_values; j++) {
        if (values[j] < min) {
            min = values[j];
        }
        if (values[j] > max) {
            max = values[j];
        }
    }
    auto dv = (max - min) / (double) (num_hists - 1);
    for (int j = 0; j < num_values; j++) {
        hists[(int) round((values[j] - min) / dv)] += 1;
    }
    return (dev_from - min) / dv;
}

double pdf(double * hists, unsigned num_hists, double const * values, unsigned num_values, double dev_from) {
    auto value_idx = count_distrib(hists, num_hists, values, num_values, dev_from);
    for (unsigned j = 0; j < num_hists; j++) {
        hists[j] /= (double) num_values;
    }
    return value_idx;
}

void print_int_distrib(double const * hists, unsigned num_hists) {
    for (int i = 0; i < num_hists; i++) {
        for (int j = 0; j < (int) hists[i]; j++) {
            cout << "|";
        }
        cout << endl;
    }
}

double dispersion(double const * values, unsigned num_values, double dev_from) {
    double sum = 0;
    for (unsigned i = 0; i < num_values; i++) {
        auto delta = values[i] - dev_from;
        sum = fma(delta, delta, sum);
    }
    return sum / (double) num_values;
}

unsigned read_data(double * values) {
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

double maxwell_pdf(double * hists, unsigned num_hists, double T) {
    auto dev_from_idx = 0.5 * (double) (num_hists - 1);
    for (int i = 0; i < num_hists; i++) {
        auto v = dev_from_idx - (double) i;
        hists[i] = exp(-v * v / T) / sqrt(T * numbers::pi);
    }
    return dev_from_idx;
}

double mean_direct(double const * psi, double const * pdf, unsigned size) {
    double sum = 0.;
    for (auto i = 0; i < size; i++) {
        sum += pdf[i] * psi[i];
    }
    return sum;
}

double mean_fma(double const * psi, double const * pdf, unsigned size) {
    double sum = 0.;
    for (auto i = 0; i < size; i++) {
        sum = fma(pdf[i], psi[i], sum);
    }
    return sum;
}

double mean_recurse(double const * psi, double const * pdf, unsigned s) {
    if (s == 1) return psi[0] * pdf[0];
    auto h = s / 2;
    return mean_recurse(psi, pdf, h) + mean_recurse(&psi[h], &pdf[h], s - h);
}

double mean_kahan(double const * psi, double const * pdf, unsigned s) {
    auto sum = 0.;
    auto err = 0.;
    for (auto i = 0; i < s; i++) {
        auto x = fma(pdf[i], psi[i], -err);
        auto sum_new = sum + x;
        err = sum_new - sum - x;
        sum = sum_new;
    }
    return sum;
}

void test_mean(double const * pdf, unsigned size, double v0, double T) {
    cout << "mean " << v0 << " from 0 to " << size - 1 << endl;
    auto psi = new double[size];
    for (int i = 0; i < size; i++) {
        psi[i] = abs(v0 - (double) i);
    }
    cout << "T: " << T << endl;
    cout << "expected mean: " << sqrt(T / numbers::pi) << endl;
    cout << "mean directly: " << mean_direct(psi, pdf, size) << endl;
    cout << "mean by fma: " << mean_fma(psi, pdf, size) << endl;
    cout << "mean by recursion: " << mean_recurse(psi, pdf, size) << endl;
    cout << "mean by Kahan: " << mean_kahan(psi, pdf, size) << endl;
    delete[] psi;
}

double test_data(unsigned num_hists) {
    cout << " --- testing raw.csv --- " << endl;

    double values[273];
    auto num_values = read_data(values);
    double expected = 500;
    auto T = 2 * dispersion(values, num_values, expected);

    auto hists = new double[num_hists] {0};
    auto expected_idx = pdf(hists, num_hists, values, num_values, expected);
    test_mean(hists, num_hists, expected_idx, T);
    delete[] hists;

    cout << endl;
    return T;
}

void test_maxwell(unsigned num_hists, double T) {
    cout << " --- testing maxwell distribution --- " << endl;

    auto hists = new double[num_hists];
    auto value_idx = maxwell_pdf(hists, num_hists, T);
    test_mean(hists, num_hists, value_idx, T);

    cout << endl;
}

int main() {
    auto s = 100;
    auto T = test_data(s);
    test_maxwell(s, T);
}