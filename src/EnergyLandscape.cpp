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

        double approx(double c) {
            return c-1;
        }

        void testBisection() {
            double result = bisection(-10, 10, 0.001, 1000);
            cout << endl << result << endl;
            cout << approx(result) << endl;
        }

        double bisection(double a, double b, double tol, int nMax) {
            int n = 1;
            double c;
            double fc;
            double fa = approx(a);
            double fb = approx(b);
            if (!(((fa < 0) && (fb > 0)) || ((fa > 0) && (fb < 0)))) {
                cerr << "Bisection: f(a) < 0 and f(b) > 0 or f(a) > 0 and f(b) < 0 not fulfilled" << endl;
            }

            while (n <= nMax) {
                c = (a + b) / 2;
                fc = approx(c);
                if ((fc == 0) || ((b - a) / 2 < tol)) {
                    return c;
                }
                n++;
                cout << a << " " << b << " " << c << endl;
                if (approx(a) * approx(c) > 0) {
                    a = c;
                } else {
                    b = c;
                }
            }
            cerr << "Chemical potential could not be calculated within nMax=" << nMax << " with tol=" << tol << ", a=" << a << ", b=" << b << endl; 
            throw std::overflow_error("");
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
    e.testBisection();
    return 0;
}