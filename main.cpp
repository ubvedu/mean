#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

using namespace std;

void opd(double * hists, int num_hists, double const * values, int num_values) {
    for (int j = 0; j < num_hists; j++) {
        hists[j] = 0;
    }
    double min = numeric_limits<double>::max();
    double max = numeric_limits<double>::min();
    for (int i = 0; i < num_values; i++) {
        if (values[i] < min) {
            min = values[i];
        }
        if (values[i] > max) {
            max = values[i];
        }
    }
    auto dv = (max - min) / (double) (num_hists - 1);
    for (int i = 0; i < num_values; i++) {
        hists[(int) roundevenf((values[i] - min) / dv)] += 1;
    }
}

void pdf(double * hists, int num_hists, double const * values, int num_values) {
    opd(hists, num_hists, values, num_values);
    for (int j = 0; j < num_hists; j++) {
        hists[j] /= (double) num_values;
    }
}

void print_opd(double const * hists, int num_hists) {
    for (int i = 0; i < num_hists; i++) {
        for (int j = 0; j < (int) hists[i]; j++) {
            cout << "|";
        }
        cout << endl;
    }
}

double std_dev(double const * values, int num_values, double mean) {
    double std_dev = 0;
    for (int i = 0; i < num_values; i++) {
        auto delta = values[i] - mean;
        std_dev = fma(delta, delta, std_dev);
    }
    return sqrt(std_dev / (double) num_values);
}

double read_pdf(double * hists, int num_hists) {
    ifstream file("../data.csv");
    auto values = new double[300];
    int num_values = 0;
    if (file.is_open()) {
        while (file >> values[num_values]) {
            num_values++;
        }
    }
    file.close();
    pdf(hists, num_hists, values, num_values);
    auto sigma = std_dev(values, num_values, 500);
    return 2 * sigma * sigma;
}

double maxwell_pdf(double * hists, int num_hists, double T) {
    for (int i = 0; i < num_hists; i++) {
        auto v = (double) i - 0.5 * (double) (num_hists - 1);
        hists[i] = exp(-v * v / T) / sqrt(T * numbers::pi);
    }
    return T;
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

double mean_rec(double const * psi, double const * pdf, unsigned s) {
    if (s == 1) return psi[0] * pdf[0];
    auto h = s / 2;
    return mean_rec(psi, pdf, h) + mean_rec(&psi[h], &pdf[h], s - h);
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

int main() {
    auto s = 100;
    auto distrib = new double[s];
    auto T = maxwell_pdf(distrib, s, 5);
    cout << "T: " << T << endl;
    cout << "expectation: " << sqrt(T / numbers::pi) << endl;
    auto psi = new double[s];
    for (int i = 0; i < s; i++) {
        auto v = (double) i - 0.5 * (double) (s - 1);
        psi[i] = abs(v);
//        psi[i] = 1;
    }
    cout << "mean direct: " << mean_direct(psi, distrib, s) << endl;
    cout << "mean by fma: " << mean_fma(psi, distrib, s) << endl;
    cout << "mean by recursion: " << mean_rec(psi, distrib, s) << endl;
    cout << "mean by Kahan: " << mean_kahan(psi, distrib, s) << endl;
    delete[] psi;
    delete[] distrib;
}