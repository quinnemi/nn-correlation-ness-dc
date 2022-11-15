#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

class EnergyLandscape {
    public:
        int L;
        double chemPot;
        vector<vector<vector<double>>> energies;

        EnergyLandscape(int L, string pattern) {
            this->L = L;
            this->energies = vector<vector<vector<double>>>(L, vector<vector<double>>(L, vector<double>(L)));
            if (pattern == "checker") {
                createCheckerPattern();
            }
        }

        void createCheckerPattern(double scale = 1) {
            if (this->L % 2 == 0) {
                throw invalid_argument("Length must be odd for checker pattern");
            }
            
            for(int i=0; i<(L*L*L); i++) {
                set(energies, indexToTuple(i), (double)(i%2));
            }
        }
    
        int* indexToTuple(int i) {
            static int a[3];
            a[2] = i % L;
            a[1] = i / L % L;
            a[0] = i / (L*L) % L;

            return a;
        }

        double acc(double** a[], int i[3]) {
            return a[i[0]][i[1]][i[2]] ;
        }

        double acc(vector<vector<vector<double>>> a, int i[3]) {
            return a[i[0]][i[1]][i[2]] ;
        }

        void set(double** a[], int i[3], double value) {
            a[i[0]][i[1]][i[2]] = value;
        }

        void set(vector<vector<vector<double>>> &a, int i[3], double value) {
            a[i[0]][i[1]][i[2]] = value;
        }

    private:

        void displayArray(int a[], int size) {
            for (int i=0; i<size; i++) {
                cout << a[i] << " ";
            }
            cout << endl;
        }

        void displayArray(vector<double> a, int size) {
            for (int i=0; i<size; i++) {
                cout << (double)a[i] << " ";
            }
            cout << endl;
        }

};

int main() {
    EnergyLandscape e = EnergyLandscape(3, "checker");

    return 0;
}