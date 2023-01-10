#ifndef ENERGYLANDSCAPE_HPP
#define ENERGYLANDSCAPE_HPP
#include "ArraySupport.hpp"
using namespace std;

class EnergyLandscape {
    public:
        int L;
        cpp_bin_float_100 beta;
        cpp_bin_float_100 chemPot;
        vector<vector<vector<cpp_bin_float_100>>> energies;

        EnergyLandscape(int L, cpp_bin_float_100 beta);
        void createCheckerPattern(cpp_bin_float_100 scale);
        vector<vector<cpp_bin_float_100>> getJumpRates(cpp_bin_float_100 attemptFreq, cpp_bin_float_100 lowerBarrier);
        vector<vector<vector<cpp_bin_float_100>>> getEqOccNum(cpp_bin_float_100 concentration, cpp_bin_float_100 epsilon);
        cpp_bin_float_100 getChemPot(cpp_bin_float_100 concentration, cpp_bin_float_100 epsilon);
        vector<vector<cpp_bin_float_100>> getEnergyDistribution(cpp_bin_float_100 epsilon);
        cpp_bin_float_100 approx(cpp_bin_float_100 mu, vector<vector<cpp_bin_float_100>> energyDist, cpp_bin_float_100 concentration);
    
    private:
        vector<cpp_bin_float_100> sortEnergies();
        cpp_bin_float_100 bisection(vector<vector<cpp_bin_float_100>> energyDist, cpp_bin_float_100 concentration, cpp_bin_float_100 a, cpp_bin_float_100 b, cpp_bin_float_100 tol, int nMax);
        int distance(int a[3], int b[3]);
        int cyclDist(int d);
        int* indexToTuple(int i);
};
#endif