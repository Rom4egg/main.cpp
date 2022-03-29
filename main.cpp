#include <iostream>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

double IE(double x, double L) {
    return sin(4 * M_PI * x / L);
}

double Left_Triangle(double Umn, double Um_1n, double tau, double h) {
    return Umn - tau / h * (Umn - Um_1n);
}

int main() {
    double L = 10;
    double T = 18;
    double Co = 0.5;
    double h = 0.03;
    double tau = Co * h;
    int Nt = T / tau + 1;
    int Nx = L / h + 1;
    vector<double> U(Nt * Nx);
    for (int i = 0; i < Nx; i++) {
        U[i] = IE(i * h, L);
    }
    for (int j = 1; j < Nt; j++) {
        for (int i = 1; i < Nx; i++) {
            U[j * Nx + i] = Left_Triangle(U[(j - 1) * Nx + i], U[(j - 1) * Nx + i - 1], tau, h);
        }
        U[j * Nx] = U[j * Nx + Nx - 1];
    }
    for (int j = 530; j < 600; j++) {
        fstream fout("a" + std::to_string(j) + ".csv", ios::out);
        fout << "time, x, U" << endl;
        for (int i = 0; i < Nx; i++) {
            fout << j * tau << ", " << i * h << ", " << U[j * Nx + i] << endl;
        }


    }
    return 0;
}
