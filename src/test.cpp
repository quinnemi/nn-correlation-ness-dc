#include "ArraySupport.hpp"
#include "EnergyLandscape.hpp"
#include "EnergyLandscape.cpp"
#include <iostream>
#include <fstream>
using namespace std;

int main(){
    ofstream file;
    file.open("../../data/BisectionTest/bisectionData.txt");
    double mu, mu0;
    EnergyLandscape* el = new EnergyLandscape(5, 1);
    for (double mu=-200; mu<=300; mu+=1) {
        mu0 = el->approx(mu, el->getEnergyDistribution(0.00000001), 0.05);
        file << mu << "\t" << mu0 << endl; 
    }
    file.close();
    return 0;
}