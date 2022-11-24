#include "ArraySupport.hpp"
#include "EnergyLandscape.hpp"
#include "EnergyLandscape.cpp"
#include <iostream>

using namespace std;

int SIZE = 3;
double C = 1;
double BETA = 1;
double E[3] = {1, 0, 0};
double LOWERBARRIER = 1;
double ATTEMPTFREQ = 1;
double EPSILON = 0.00000001;

EnergyLandscape* energyLandscape = new EnergyLandscape(SIZE, BETA);
vector<vector<double>> jumpRates = energyLandscape->getJumpRates(ATTEMPTFREQ, LOWERBARRIER);
vector<vector<vector<double>>> eqOccNum = energyLandscape->getEqOccNum(C, EPSILON);

void K1(int i[3]) {
    double term1 = 0;
    for (auto k : getNNidx(i)) {
        term1 += -(acc(jumpRates, i, k));
    }
}

int main() {
    setLength(SIZE);
    int i[3] = {1,1,1};
    K1(i);
    return 0;
};