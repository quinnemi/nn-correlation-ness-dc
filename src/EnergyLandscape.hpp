#ifndef ENERGYLANDSCAPE_HPP
#define ENERGYLANDSCAPE_HPP
#include "ArraySupport.hpp"
#include <vector>
using namespace std;

class EnergyLandscape {
    public:
        int L;
        double beta;
        double chemPot;
        vector<vector<vector<double>>> energies;

        EnergyLandscape(int L, double beta);
        void createCheckerPattern(double scale);
        vector<vector<double>> getJumpRates(double attemptFreq, double lowerBarrier);
        vector<vector<vector<double>>> getEqOccNum(double concentration, double epsilon);
        double getChemPot(double concentration, double epsilon);
        vector<vector<double>> getEnergyDistribution(double epsilon);
    
    private:
        vector<double> sortEnergies();
        double approx(double mu, vector<vector<double>> energyDist, double concentration);
        double bisection(vector<vector<double>> energyDist, double concentration, double a, double b, double tol, int nMax);
        int distance(int a[3], int b[3]);
        int cyclDist(int d);
        int* indexToTuple(int i);
};
#endif