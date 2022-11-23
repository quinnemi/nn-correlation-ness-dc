#include "EnergyLandscape.hpp"
#include "ArraySupport.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

int L;
double beta;
double chemPot;
vector<vector<vector<double>>> energies;

EnergyLandscape::EnergyLandscape(int L, double beta) {
    this->L = L;
    this->beta = beta;
    this->energies = vector<vector<vector<double>>>(L, vector<vector<double>>(L, vector<double>(L)));
    createCheckerPattern(1);
    //if (pattern == "checker") {
    //}
}

void EnergyLandscape::createCheckerPattern(double scale) {
    if (this->L % 2 == 0) {
        throw invalid_argument("Length must be odd for checker pattern");
    }
    
    for(int i=0; i<(L*L*L); i++) {
        set(energies, indexToTuple(i), (double)(i%2));
    }
}

vector<vector<double>> EnergyLandscape::getJumpRates(double attemptFreq, double lowerBarrier) {
    double nu = attemptFreq * exp(-(this->beta * lowerBarrier));

    auto jumpRate = [this, nu, attemptFreq](int* i, int* j) {
        vector<double> term1 = {1.0, exp(-(this->beta * (acc(this->energies, j) - acc(this->energies, i))))};
        return nu * *min_element(term1.begin(), term1.end());
    };

    vector<vector<double>> jumpRates = vector<vector<double>>(pow(L, 3), vector<double>(pow(L, 3)));
    for (int i=0; i<pow(L, 3); i++) {
        for (int j=0; j<pow(L, 3); j++) {
            if (distance(indexToTuple(i), indexToTuple(j)) <= 1) {
                jumpRates[i][j] = jumpRate(indexToTuple(i), indexToTuple(j));
            } else {
                jumpRates[i][j] = 0;
            }
        }
    }

    return jumpRates;
}

vector<vector<vector<double>>> EnergyLandscape::getEqOccNum(double concentration, double epsilon) {
    double chemPot = getChemPot(concentration, epsilon);

    auto eqOccNum = [this, chemPot](int* i) {
        return 1 / (exp(-this->beta * (acc(this->energies, i) - chemPot)) + 1);  
    };

    vector<vector<vector<double>>> eqOccNums = vector<vector<vector<double>>>(L, vector<vector<double>>(L, vector<double>(L)));
    for (int i=0; i<pow(L, 3); i++) {
        set(eqOccNums, indexToTuple(i), eqOccNum(indexToTuple(i)));
    }

    return eqOccNums;
}

double EnergyLandscape::getChemPot(double concentration, double epsilon) {
    if (this->chemPot) {
        return this->chemPot;
    }

    // calculate site energy distribution g(E), along with sorted values of E
    vector<vector<double>> energyDist = getEnergyDistribution(epsilon);

    double a = -9.0;
    double b = 10.0;
    double tol = 1e-8;

    this->chemPot = bisection(energyDist, concentration, 10*a, 10*b, tol, 1000);

    return this->chemPot;
}

vector<vector<double>> EnergyLandscape::getEnergyDistribution(double epsilon) {
    vector<double> sortedEnergies = sortEnergies();
    vector<double> probDist = {1};
    
    // combine epsilon-equal energy values in sortedEnergies and count up in probDist
    for (int i=0; i<sortedEnergies.size(); i++) {
        if (compare(sortedEnergies[i], sortedEnergies[i+1], epsilon) && i < sortedEnergies.size()-1) {
            sortedEnergies.erase(next(sortedEnergies.begin(), i));
            probDist[i] = probDist[i] + 1.0;
            i--;
        } else if (i < sortedEnergies.size()-1) {
            probDist.push_back(1.0);
        }
    }
    for (int i=0; i<probDist.size(); i++) {
        probDist[i] = probDist[i] / (double)pow(L, 3);
    }

    double sum = 0;
    for (double a: probDist) {
        sum += a;
    }

    vector<vector<double>> ret;
    ret.push_back(sortedEnergies);
    ret.push_back(probDist);
    return ret; 
}

vector<double> EnergyLandscape::sortEnergies() {
    vector<double> sortedEnergies = vector<double>(pow(L, 3), 0);
    int* indexTuple;
    for (int i=0; i<pow(L, 3); i++) {
        indexTuple = indexToTuple(i);
        sortedEnergies[i] = this->energies[indexTuple[0]][indexTuple[1]][indexTuple[2]];
    }
    sort(sortedEnergies.begin(), sortedEnergies.end());

    return sortedEnergies;
}

double EnergyLandscape::approx(double mu, vector<vector<double>> energyDist, double concentration) {
    double sum = 0;
    for(int i=0; i<energyDist[0].size(); i++) {
        sum += energyDist[1][i] / (exp(this->beta * (energyDist[0][i] - mu)) + 1);
    }
    return sum-concentration;
}

double EnergyLandscape::bisection(vector<vector<double>> energyDist, double concentration, double a, double b, double tol, int nMax) {
    int n = 1;
    double c;
    double fc;
    double fa = approx(a, energyDist, concentration);
    double fb = approx(b, energyDist, concentration);

    if (!(((fa < 0) && (fb > 0)) || ((fa > 0) && (fb < 0)))) {
        cerr << "Bisection: f(a) < 0 and f(b) > 0 or f(a) > 0 and f(b) < 0 not fulfilled" << endl;
    }

    while (n <= nMax) {
        c = (a + b) / 2.0;
        fc = approx(c, energyDist, concentration);
        if ((fc == 0) || (abs((b - a)) / 2 < tol)) {
            fa = approx(a, energyDist, concentration);
            fb = approx(b, energyDist, concentration);
            cout << "bisection ended after n=" << n << " with mu=" << c << " f(mu)=" << fc << " tol=" << tol << endl;
            return c;
        }
        n++;
        if (approx(a, energyDist, concentration) * approx(c, energyDist, concentration) >= 0) {
            a = c;
        } else {
            b = c;
        }
    }
    cerr << "Chemical potential could not be calculated within nMax=" << nMax << " with tol=" << tol << ", a=" << a << ", b=" << b << endl; 
    throw std::overflow_error("");
}

int EnergyLandscape::distance(int a[3], int b[3]) {
        return cyclDist(abs(a[0]-b[0])) + cyclDist(abs(a[1]-b[1])) + cyclDist(abs(a[2]-b[2]));
}

int EnergyLandscape::cyclDist(int d) {
    if (d > L / 2.0) {
        d = L - d;
    }
    return d;
}

int* EnergyLandscape::indexToTuple(int i) {
    int* a = new int[3];
    a[2] = i % L;
    a[1] = i / L % L;
    a[0] = i / (L*L) % L;

    return a;
}