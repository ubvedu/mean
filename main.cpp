#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

using namespace std;

void opd(float * hists, int num_hists, float const * values, int num_values) {
    for (int j = 0; j < num_hists; j++) {
        hists[j] = 0;
    }
    float min = numeric_limits<float>::max();
    float max = numeric_limits<float>::min();
    for (int i = 0; i < num_values; i++) {
        if (values[i] < min) {
            min = values[i];
        }
        if (values[i] > max) {
            max = values[i];
        }
    }
    auto dv = (max - min) / (float) (num_hists - 1);
    for (int i = 0; i < num_values; i++) {
        hists[(int) roundevenf((values[i] - min) / dv)] += 1;
    }
}

void pdf(float * hists, int num_hists, float const * values, int num_values) {
    opd(hists, num_hists, values, num_values);
    for (int j = 0; j < num_hists; j++) {
        hists[j] /= (float) num_values;
    }
}

void print_opd(float const * hists, int num_hists) {
    for (int i = 0; i < num_hists; i++) {
        for (int j = 0; j < (int) hists[i]; j++) {
            cout << "|";
        }
        cout << endl;
    }
}

float std_dev(float const * values, int num_values, float mean) {
    float std_dev = 0;
    for (int i = 0; i < num_values; i++) {
        auto delta = values[i] - mean;
        std_dev = fma(delta, delta, std_dev);
    }
    return sqrtf(std_dev / (float) num_values);
}

float read_pdf(float * hists, int num_hists) {
    ifstream file("../data.csv");
    auto values = new float[300];
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

float maxwell_pdf(float * hists, int num_hists, float T) {
    for (int i = 0; i < num_hists; i++) {
        auto v = 0.5f + (float) i - 0.5f * (float) num_hists;
        hists[i] = exp(-v * v / T) / sqrt((float) (T * numbers::pi));
    }
    return T;
}

float bias(float v, float v0) { return abs(v - v0); }

float mean_direct(float const * psi, float const * pdf, unsigned size) {
    float sum = 0.f;
    for (auto i = 0; i < size; i++) {
        sum += pdf[i] * psi[i];
    }
    return sum;
}

float mean_fma(float const * psi, float const * pdf, unsigned size) {
    float sum = 0.f;
    for (auto i = 0; i < size; i++) {
        sum = fma(pdf[i], psi[i], sum);
    }
    return sum;
}

float mean_recursive(float const * psi, float const * pdf, unsigned size) {
    float sum = 0.f;
    if (size == 1) {
        return psi[0] * pdf[0];
    }
    return mean_recursive(psi, pdf, size / 2) + mean_recursive(&psi[size / 2], &pdf[size / 2], size - size / 2);
}

int main() {
    auto s = 100;
    auto distrib = new float[s];
    auto T = maxwell_pdf(distrib, s, 2);
    cout << "T: " << T << endl;
    cout << "expectation: " << sqrt(T / numbers::pi) << endl;
    auto psi = new float[s];
    for (int i = 0; i < s; i++) {
        psi[i] = bias(distrib[i], 0);
//        psi[i] = 1;
    }
    cout << "mean direct: " << mean_direct(psi, distrib, s) << endl;
    cout << "mean by fma: " << mean_fma(psi, distrib, s) << endl;
    cout << "mean by recursion: " << mean_recursive(psi, distrib, s) << endl;
    delete[] psi;
    delete[] distrib;
}