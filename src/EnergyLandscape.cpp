#include "EnergyLandscape.hpp"
#include "ArraySupport.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <boost/multiprecision/cpp_bin_float.hpp>
using namespace boost::multiprecision;
using namespace std;

int L;
cpp_bin_float_100 beta;
cpp_bin_float_100 chemPot;
vector<vector<vector<cpp_bin_float_100>>> energies;

EnergyLandscape::EnergyLandscape(int L, cpp_bin_float_100 beta) {
    this->L = L;
    this->beta = beta;
    this->energies = vector<vector<vector<cpp_bin_float_100>>>(L, vector<vector<cpp_bin_float_100>>(L, vector<cpp_bin_float_100>(L)));
    createCheckerPattern(0.001);
}

void EnergyLandscape::createCheckerPattern(cpp_bin_float_100 scale) {
    if (this->L % 2 == 0) {
        throw invalid_argument("Length must be odd for checker pattern");
    }
    
    for(int i=0; i<(L*L*L); i++) {
        setVal(energies, indexToTuple(i), (cpp_bin_float_100)(i%2));
    }
}

vector<vector<cpp_bin_float_100>> EnergyLandscape::getJumpRates(cpp_bin_float_100 attemptFreq, cpp_bin_float_100 lowerBarrier) {
    cpp_bin_float_100 nu = attemptFreq * exp(-(this->beta * lowerBarrier));
    auto jumpRate = [this, nu, attemptFreq](int* i, int* j) {
        vector<cpp_bin_float_100> term1 = {1.0, exp(-(this->beta * (acc(this->energies, j) - acc(this->energies, i))))};
        return nu * *min_element(term1.begin(), term1.end());
    };

    vector<vector<cpp_bin_float_100>> jumpRates = vector<vector<cpp_bin_float_100>>(pow(L, 3), vector<cpp_bin_float_100>(pow(L, 3)));
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

vector<vector<vector<cpp_bin_float_100>>> EnergyLandscape::getEqOccNum(cpp_bin_float_100 concentration, cpp_bin_float_100 epsilon) {
    cpp_bin_float_100 chemPot = getChemPot(concentration, epsilon);
    cpp_bin_float_100 one = pow(2,-112);
    auto eqOccNum = [this, chemPot, one](int* i) {
        //cout << (exp(this->beta * (acc(this->energies, i) - chemPot))) << " " << one << endl;
        return 1 / (exp(this->beta * (acc(this->energies, i) - chemPot)) + one); //! 
    };

    vector<vector<vector<cpp_bin_float_100>>> eqOccNums = vector<vector<vector<cpp_bin_float_100>>>(L, vector<vector<cpp_bin_float_100>>(L, vector<cpp_bin_float_100>(L)));
    for (int i=0; i<pow(L, 3); i++) {
        setVal(eqOccNums, indexToTuple(i), eqOccNum(indexToTuple(i)));
    }

    return eqOccNums;
}

cpp_bin_float_100 EnergyLandscape::getChemPot(cpp_bin_float_100 concentration, cpp_bin_float_100 epsilon) {
    if (this->chemPot) {
        return this->chemPot;
    }

    // calculate site energy distribution g(E), along with sorted values of E
    vector<vector<cpp_bin_float_100>> energyDist = getEnergyDistribution(epsilon);

    cpp_bin_float_100 a = -100.0;
    cpp_bin_float_100 b = 200.0;
    cpp_bin_float_100 tol = 1e-8;

    this->chemPot = bisection(energyDist, concentration, a, b, tol, 1000);

    return this->chemPot;
}

vector<vector<cpp_bin_float_100>> EnergyLandscape::getEnergyDistribution(cpp_bin_float_100 epsilon) {
    vector<cpp_bin_float_100> sortedEnergies = sortEnergies();
    vector<cpp_bin_float_100> probDist = {1};
    
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
        probDist[i] = probDist[i] / (cpp_bin_float_100)pow(L, 3);
    }

    cpp_bin_float_100 sum = 0;
    for (cpp_bin_float_100 a: probDist) {
        sum += a;
    }

    vector<vector<cpp_bin_float_100>> ret;
    ret.push_back(sortedEnergies);
    ret.push_back(probDist);
    return ret; 
}

vector<cpp_bin_float_100> EnergyLandscape::sortEnergies() {
    vector<cpp_bin_float_100> sortedEnergies = vector<cpp_bin_float_100>(pow(L, 3), 0);
    int* indexTuple;
    for (int i=0; i<pow(L, 3); i++) {
        indexTuple = indexToTuple(i);
        sortedEnergies[i] = this->energies[indexTuple[0]][indexTuple[1]][indexTuple[2]];
    }
    sort(sortedEnergies.begin(), sortedEnergies.end());

    return sortedEnergies;
}

cpp_bin_float_100 EnergyLandscape::approx(cpp_bin_float_100 mu, vector<vector<cpp_bin_float_100>> energyDist, cpp_bin_float_100 concentration) {
    cpp_bin_float_100 sum = 0;
    for(int i=0; i<energyDist[0].size(); i++) {
        sum += energyDist[1][i] / (exp(this->beta * (energyDist[0][i] - mu)) + 1);
    }
    return sum-concentration;
}

cpp_bin_float_100 EnergyLandscape::bisection(vector<vector<cpp_bin_float_100>> energyDist, cpp_bin_float_100 concentration, cpp_bin_float_100 a, cpp_bin_float_100 b, cpp_bin_float_100 tol, int nMax) {
    int n = 1;
    cpp_bin_float_100 c;
    cpp_bin_float_100 fc;
    cpp_bin_float_100 fa = approx(a, energyDist, concentration);
    cpp_bin_float_100 fb = approx(b, energyDist, concentration);

    if (!(((fa < 0) && (fb > 0)) || ((fa > 0) && (fb < 0)))) {
        cerr << "Bisection: f(a) < 0 and f(b) > 0 or f(a) > 0 and f(b) < 0 not fulfilled" << endl;
    }

    while (n <= nMax) {
        c = (a + b) / 2.0;
        fc = approx(c, energyDist, concentration);
        if ((fc == 0) || (abs((b - a)) / 2 < tol)) {
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

