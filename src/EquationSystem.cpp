#include "ArraySupport.hpp"
#include "EnergyLandscape.hpp"
#include "EnergyLandscape.cpp"
#include <iostream>
#include <cfloat>

using namespace std;

int SIZE = 5;
double C = 0.95;
double BETA = 100;
double E[3] = {1, 0, 0};
double LOWERBARRIER = 1;
double ATTEMPTFREQ = 1;
double EPSILON = 0.00000001;

int NSITES = pow(SIZE, 3);

EnergyLandscape* energyLandscape = new EnergyLandscape(SIZE, BETA);
vector<vector<double>> jumpRates = energyLandscape->getJumpRates(ATTEMPTFREQ, LOWERBARRIER);
vector<vector<vector<double>>> eqOccNum = energyLandscape->getEqOccNum(C, EPSILON);

vector<int> P;

double D(int i[3], int j[3]) {
    double term1 = 0;
    for (vector<int> k: getNNidx(i)) {
        if (!eqIdx(k, j)) {
            term1 += acc(jumpRates, k, i) * acc(eqOccNum, k);
        }
    }
    term1 /= acc(eqOccNum, i);

    double term2 = 0;
    for (vector<int> l: getNNidx(j)) {
        if (!eqIdx(l, i)) {
            term2 += acc(jumpRates, l, j) * acc(eqOccNum, l);
        }
    }
    term2 /= acc(eqOccNum, j);

    return term1 + term2;
}

double D(int i[3], vector<int> j) {
    int k[3] = {j[0], j[1], j[2]};
    return D(i,k);
}

double b(int i[3], int j[3]) {
    int diff[3] = {j[0]-i[0], j[1]-i[1], j[2]-i[2]};
    double dotProduct = diff[0]*E[0] + diff[1]*E[1] + diff[2]*E[2];
    
    return dotProduct;
}

double b(vector<int> i, int j[3]) {
    int k[3] = {i[0], i[1], i[2]};
    return b(k, j);
}

double b(vector<int> i, vector<int> j) {
    int k[3] = {i[0], i[1], i[2]};
    int l[3] = {j[0], j[1], j[2]};
    return b(k, l);
}

double K1(int i[3]) {
    
    double term1 = 0;
    for (vector<int> k : getNNidx(i)) {
        term1 += -(acc(jumpRates, i, k));
    }

    double term2 = 0;
    for (vector<int> k: getNNidx(i)) {
        double term2_1 = (acc(jumpRates, i, k) - acc(jumpRates, k, i)) / D(i, k);
        double term2_2 = 0; 
        for (vector<int> l: getNNidx(k)) {
            if (!eqIdx(i, l)) {
                term2_2 += acc(jumpRates, l, k) * acc(eqOccNum, l);
            }
        }
        term2 += term2_1 * term2_2;
    }

    return term1 + term2;
}

double K2(int i[3], vector<int> j) {
    double term1 = acc(jumpRates, j, i);
    
    double term2_1 = (acc(jumpRates, i, j) - acc(jumpRates, j, i)) / D(i, j);
    double term2_2 = 0;
    for (vector<int> k: getNNidx(i)) {
        if (!eqIdx(k, j)) {
            term2_2 = acc(jumpRates, k, i) * acc(eqOccNum, k);
        }
    }
    double term2 = term2_1 * term2_2;

    double term3 = 0;
    double term3_1 = acc(jumpRates, j, i) * (1 - acc(eqOccNum, i)) / (1 - acc(eqOccNum, j));
    double term3_2 = 0;
    for (vector<int> k: getNNidx(i)) {
        term3_2 += (acc(jumpRates, i, k) - acc(jumpRates, k, i)) / D(i, k) * acc(eqOccNum, k);
    }
    term3 = term3_1 * term3_2;

    return term1 + term2 + term3;
}

double K3(int i[3], vector<int> j) {
    double term1_1 = acc(eqOccNum, i) / (1 - acc(eqOccNum, j));
    double term1_2 = 0;
    for (vector<int> k: getNNidx(i)) {
        term1_2 += (acc(jumpRates, i ,k) - acc(jumpRates, k, i)) / D(i, k) * acc(jumpRates, j, k) * (1 - acc(eqOccNum, k));
    }

    return term1_1 * term1_2;
}

double hi(int i[3]) {
    double term1 = 0;
    for (vector<int> j: getNNidx(i)) {
        term1 += b(j, i) * acc(jumpRates, j, i) * acc(eqOccNum, j);
    }
    term1 *= -(1 - acc(eqOccNum, i));

    double term2 = 0;
    for (vector<int> j: getNNidx(i)) {
        double term2_1 = (acc(jumpRates, j, i) - acc(jumpRates, i, j)) / D(i, j);
        double term2_2_1 = 0;
        for (vector<int> k: getNNidx(i)) {
            if (!eqIdx(k, j)) {
                term2_2_1 += b(k, i) * acc(jumpRates, k, i) * acc(eqOccNum, k);
            }
        }
        term2_2_1 *= (1 - acc(eqOccNum, i)) * acc(eqOccNum, j);

        double term2_2_2 = 0;
        for (vector<int> l: getNNidx(j)) {
            if (!eqIdx(l, i)) {
                term2_2_2 = b(l, j) * acc(jumpRates, l, j) * acc(eqOccNum, l);
            }
        }
        term2_2_2 *= acc(eqOccNum, i) * (1 - acc(eqOccNum, j));

        term2 += term2_1 * (term2_2_1 + term2_2_2);
    }

    return term1 + term2;
}

vector<double> h() {
    vector<double> ret = vector<double>(NSITES);
    for (int i=0; i<NSITES; i++) {
        ret[i] = hi(indexToTuple(i));
    }

    return ret;
}

vector<vector<double>> K() {
    vector<vector<double>> Kmat = vector<vector<double>>(NSITES, vector<double>(NSITES));

    for (vector<double> row: Kmat) {
        row.end();
        fill(row.begin(), row.end(), 0);
    }

    for(int idxI=0; idxI<NSITES; idxI++) {
        int* i = indexToTuple(idxI);
        vector<vector<int>> nn = getNNidx(i);
        vector<vector<int>> nnn = getNNNidx(i);
        
        Kmat[idxI][idxI] = K1(i);
        
        for (vector<int> j: nn) {
            Kmat[idxI][tupleToIndex(j)] = K2(i, j);
        }
        for (vector<int> j: nnn) {
            Kmat[idxI][tupleToIndex(j)] = K3(i, j);
        }
    }

    return Kmat;
}

vector<vector<double>> LUPDecomposition(vector<vector<double>> A, bool debug) {
    // https://en.wikipedia.org/wiki/LU_decomposition
    int N = A.size();
    for (int i=0; i<=N; i++) {
        P.push_back(i);
    }

    int imax, i, j, k;
    double maxA, absA;
    vector<double> row;
    for (i=0; i<N; i++) {
        maxA = 0;
        imax = i;
        for (k=i; k<N; k++) {
            if ((absA = abs(A[k][i])) > maxA) {
                maxA = absA;
                imax = k;
            }
        }

        if (debug) {
            cout << "LUPDecomposition pivot element: " << maxA << endl;
        }

        if (imax != i) {
            j = P[i];
            P[i] = imax;
            P[imax] = j;

            row = A[i];
            A[i] = A[imax];
            A[imax] = row;

            P[N]++;
        }

        for (j=i+1; j<N; j++) {
            A[j][i] /= A[i][i];
            for (k=i+1; k<N; k++) {
                A[j][k] -= A[j][i] * A[i][k];
            }
        }
    }

    return A;
}

vector<double> LUPSolve(vector<vector<double>> A, vector<double> b) {
    int i, k;
    int N = A.size();
    vector<double> x = vector<double>(N);
    
    for (i=0; i<N; i++) {
        x[i] = b[P[i]];
        for (k=0; k<i; k++) {
            x[i] -= A[i][k] * x[k];
        }
    }

    for (i=N-1; i>=0; i--) {
        for (k=i+1; k<N; k++) {
            x[i] -= A[i][k] * x[k];
        }
        x[i] /= A[i][i];
    }

    return x;

}

int main() {
    setLength(SIZE);
    cout << "setting up K" << endl;
    vector<vector<double>> A = K();
    displayArray(A);
    cout << "setting up h" << endl;
    vector<double> b = h();
    cout << "solving equation system" << endl;
    vector<double> res = LUPSolve(LUPDecomposition(A, false), b);
    displayArray(res, res.size());

    return 0;
};