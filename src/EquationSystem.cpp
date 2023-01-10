#include "ArraySupport.hpp"
#include "EnergyLandscape.hpp"
#include "EnergyLandscape.cpp"
#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>
using namespace boost::multiprecision;
using namespace std;

cpp_bin_float_100 one(1);

int SIZE = 11;
cpp_bin_float_100 C = 0.95;
cpp_bin_float_100 BETA = 100;
cpp_bin_float_100 E[3] = {1, 0, 0};
cpp_bin_float_100 LOWERBARRIER = 1;
cpp_bin_float_100 ATTEMPTFREQ = 1;
cpp_bin_float_100 EPSILON = 0.00000001;
cpp_bin_float_100 B = 1;

int NSITES = pow(SIZE, 3);

EnergyLandscape* energyLandscape = new EnergyLandscape(SIZE, BETA);
vector<vector<cpp_bin_float_100>> jumpRates = energyLandscape->getJumpRates(ATTEMPTFREQ, LOWERBARRIER);
vector<vector<vector<cpp_bin_float_100>>> eqOccNumTup = energyLandscape->getEqOccNum(C, EPSILON);
vector<cpp_bin_float_100> eqOccNum;

vector<int> P;

cpp_bin_float_100 D(int itup[3], int jtup[3]) {
    cpp_bin_float_100 term1 = 0;
    int i = tupleToIndex(itup);
    int j = tupleToIndex(jtup);
    int k;
    for (vector<int> ktup: getNNidx(itup)) {
        k = tupleToIndex(ktup);
        if (k != j) {
            term1 += jumpRates[k][i] * eqOccNum[k];
        }
    }
    term1 /= eqOccNum[i];

    cpp_bin_float_100 term2 = 0;
    int l;
    for (vector<int> ltup: getNNidx(jtup)) {
        l = tupleToIndex(ltup);
        if (l != i) {
            term2 += jumpRates[l][j] * eqOccNum[l];
        }
    }
    term2 /= eqOccNum[j];

    return term1 + term2;
}

cpp_bin_float_100 D(int i[3], vector<int> j) {
    int k[3] = {j[0], j[1], j[2]};
    return D(i,k);
}

cpp_bin_float_100 b(int i[3], int j[3]) {
    int diff[3] = {j[0]-i[0], j[1]-i[1], j[2]-i[2]};
    cpp_bin_float_100 dotProduct = diff[0]*E[0] + diff[1]*E[1] + diff[2]*E[2];
    
    return dotProduct;
}

cpp_bin_float_100 b(vector<int> i, int j[3]) {
    int k[3] = {i[0], i[1], i[2]};
    return b(k, j);
}

cpp_bin_float_100 b(vector<int> i, vector<int> j) {
    int k[3] = {i[0], i[1], i[2]};
    int l[3] = {j[0], j[1], j[2]};
    return b(k, l);
}

cpp_bin_float_100 K1(int itup[3]) {
    //cout << "K1: " << itup[0] << " " << itup[1] << " " << itup[2] << ", " << endl;
    int i = tupleToIndex(itup);
    cpp_bin_float_100 term1 = 0;
    int k;
    for (vector<int> ktup : getNNidx(itup)) {
        k = tupleToIndex(ktup);
        term1 += -jumpRates[i][k];
    }

    cpp_bin_float_100 term2 = 0;
    for (vector<int> ktup: getNNidx(itup)) {
        k = tupleToIndex(ktup);
        cpp_bin_float_100 term2_1 = (jumpRates[i][k] - jumpRates[k][i]) / D(itup, ktup);
        cpp_bin_float_100 term2_2 = 0;
        int l;
        for (vector<int> ltup: getNNidx(ktup)) {
            l = tupleToIndex(ltup);
            if (i != l) {
                term2_2 += jumpRates[l][k] * eqOccNum[l];
            }
        }
        term2 += term2_1 * term2_2;
    }

    return term1 + term2;
}

cpp_bin_float_100 term1;
cpp_bin_float_100 term2_1;
cpp_bin_float_100 term2_2;
cpp_bin_float_100 term2;
cpp_bin_float_100 term3;
cpp_bin_float_100 term3_1_a;
cpp_bin_float_100 term3_1_b;
cpp_bin_float_100 term3_1_c;
cpp_bin_float_100 term3_1;
cpp_bin_float_100 term3_2;
cpp_bin_float_100 ret;
cpp_bin_float_100 K2(int itup[3], vector<int> jtup) {
    //cout << "K2: " << itup[0] << " " << itup[1] << " " << itup[2] << ", " << jtup[0] << " " << jtup[1] << " " << jtup[2] << endl;
    int i = tupleToIndex(itup);
    int j = tupleToIndex(jtup);
    term1 = jumpRates[j][i];
    
    term2_1 = (jumpRates[i][j] - jumpRates[j][i]) / D(itup, jtup);
    term2_2 = 0;
    int k;
    for (vector<int> ktup: getNNidx(itup)) {
        k = tupleToIndex(ktup);
        if (k != j) {
            term2_2 = jumpRates[k][i] * eqOccNum[k];
        }
    }
    term2 = term2_1 * term2_2;

    term3 = 0;
    term3_1_a = jumpRates[j][i];
    term3_1_b = (one - eqOccNum[i]);
    term3_1_c = (one - eqOccNum[j]);
    term3_1 = term3_1_a * term3_1_b / term3_1_c;
    term3_2 = 0;
    for (vector<int> ktup: getNNidx(itup)) {
        k = tupleToIndex(ktup);
        term3_2 += (jumpRates[i][k] - jumpRates[k][i]) / D(itup, ktup) * eqOccNum[k];
    }
    term3 = term3_1 * term3_2;
    cpp_bin_float_100 ret = term1 + term2 + term3;;
    if (isinf(abs(ret)) or isnan(abs(ret))) {
        displayArray(eqOccNum);
        cout << ret << endl;
    }
    return ret;
}
cpp_bin_float_100 term1_1;
cpp_bin_float_100 term1_2;
cpp_bin_float_100 K3(int itup[3], vector<int> jtup) {
    //cout << "K3: " << itup[0] << " " << itup[1] << " " << itup[2] << ", " << jtup[0] << " " << jtup[1] << " " << jtup[2] << endl;
    int i = tupleToIndex(itup);
    int j = tupleToIndex(jtup);
    term1_1 = eqOccNum[i] / (one - eqOccNum[j]);
    term1_2 = 0;
    int k;
    for (vector<int> ktup: getNNidx(itup)) {
        k = tupleToIndex(ktup);
        term1_2 += (jumpRates[i][k] - jumpRates[k][i]) / D(itup, ktup) * jumpRates[j][k] * (one - eqOccNum[k]);
    }

    return term1_1 * term1_2;
}

cpp_bin_float_100 term2_2_1;
cpp_bin_float_100 term2_2_2;
cpp_bin_float_100 hi(int itup[3]) {
    term1 = 0;
    int i = tupleToIndex(itup);
    int j;
    for (vector<int> jtup: getNNidx(itup)) {
        j = tupleToIndex(jtup);
        term1 += b(jtup, itup) * jumpRates[j][i] * eqOccNum[j];
    }
    term1 *= -(one - eqOccNum[i]);

    term2 = 0;
    int k;
    int l;
    for (vector<int> jtup: getNNidx(itup)) {
        j = tupleToIndex(jtup);
        term2_1 = (jumpRates[j][i] - jumpRates[i][j]) / D(itup, jtup);
        term2_2_1 = 0;
        for (vector<int> ktup: getNNidx(itup)) {
            k = tupleToIndex(ktup);
            if (k != j) {
                term2_2_1 += b(ktup, itup) * jumpRates[k][i] * eqOccNum[k];
            }
        }
        term2_2_1 *= (one - eqOccNum[i]) * eqOccNum[j];

        term2_2_2 = 0;
        for (vector<int> ltup: getNNidx(jtup)) {
            l = tupleToIndex(ltup);
            if (l != i) {
                term2_2_2 = b(ltup, jtup) * jumpRates[l][j] * eqOccNum[l];
            }
        }
        term2_2_2 *= eqOccNum[i] * (one - eqOccNum[j]);

        term2 += term2_1 * (term2_2_1 + term2_2_2);
    }

    return term1 + term2;
}

vector<cpp_bin_float_100> h() {
    vector<cpp_bin_float_100> ret = vector<cpp_bin_float_100>(NSITES);
    for (int i=0; i<NSITES; i++) {
        ret[i] = hi(indexToTuple(i));
    }

    return ret;
}

vector<vector<cpp_bin_float_100>> K() {
    vector<vector<cpp_bin_float_100>> Kmat = vector<vector<cpp_bin_float_100>>(NSITES, vector<cpp_bin_float_100>(NSITES));

    for (vector<cpp_bin_float_100> row: Kmat) {
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

vector<vector<cpp_bin_float_100>> LUPDecomposition(vector<vector<cpp_bin_float_100>> A, bool debug) {
    // https://en.wikipedia.org/wiki/LU_decomposition
    cout << "LUP decomposition" << endl;
    int N = A.size();
    for (int i=0; i<=N; i++) {
        P.push_back(i);
    }

    int imax, i, j, k;
    cpp_bin_float_100 maxA, absA;
    vector<cpp_bin_float_100> row;
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

vector<cpp_bin_float_100> LUPSolve(vector<vector<cpp_bin_float_100>> A, vector<cpp_bin_float_100> b) {
    cout << "Solving LUP decomposed system" << endl;
    int i, k;
    int N = A.size();
    vector<cpp_bin_float_100> x = vector<cpp_bin_float_100>(N);
    
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

cpp_bin_float_100 term2_1_1;
cpp_bin_float_100 term2_1_2;
cpp_bin_float_100 term3_1_1;
cpp_bin_float_100 term3_1_2;
cpp_bin_float_100 I(int itup[3], int jtup[3], vector<cpp_bin_float_100> kappa) {
    int i = tupleToIndex(itup);
    int j = tupleToIndex(jtup);
    term1 = jumpRates[i][j] * b(itup, jtup) * eqOccNum[i] * (1 - eqOccNum[j]) + jumpRates[i][j] * kappa[i] + jumpRates[j][i] * kappa[j];
    term2 = 0;
    
    term2_1 = 0;
    term2_1_1 = 0;
    int k;
    for (vector<int> ktup : getNNidx(itup)) {
        k = tupleToIndex(ktup);       
        if (k != j) {
            term2_1_1 += b(ktup, itup) * jumpRates[k][i] * eqOccNum[k] + jumpRates[k][i] / (1-eqOccNum[k]) * kappa[k];
        }
    }
    term2_1_1 *= eqOccNum[j] * (1-eqOccNum[i]);

    term2_1_2 = 0;
    for (vector<int> ktup : getNNidx(itup)) {
        k = tupleToIndex(ktup);
        if (k != j) {
            term2_1_2 += jumpRates[k][i] * eqOccNum[k];
        }
    }
    term2_1_2 *= kappa[j];

    term2_1 = (jumpRates[i][j] - jumpRates[j][i])/D(itup, jtup) * (term2_1_1 + term2_1_2);

    term3_1 = 0;
    term3_1_1 = 0;
    int l;
    for (vector<int> ltup : getNNidx(jtup)) {
        l = tupleToIndex(ltup);
        if (i != l) {
            term3_1_1 += b(ltup, jtup) * jumpRates[l][j] * eqOccNum[l] + jumpRates[l][j] / (1-eqOccNum[l]) * kappa[j];
        }
    }
    term3_1_1 *= eqOccNum[i] * (1-eqOccNum[j]);

    term3_1_2 = 0;
    for (vector<int> ltup : getNNidx(jtup)) {
        l = tupleToIndex(ltup);
        if (i != l) {
            term3_1_2 += jumpRates[l][j] * eqOccNum[l];
        }
    }
    term3_1_2 *= kappa[i];

    term3_1 = (jumpRates[j][i] - jumpRates[i][j]) / D(jtup, itup) * (term3_1_1 + term3_1_2);

    return B * (term1 + term2_1 - term3_1);
}

cpp_bin_float_100 sigmaDCequation(vector<cpp_bin_float_100> kappa) {
    cpp_bin_float_100 term1 = 0;
    cpp_bin_float_100 nu = ATTEMPTFREQ * exp(-(BETA * LOWERBARRIER));
    int* itup;
    int* jtup;
    for (int i = 0; i < NSITES; i++) {
        itup = indexToTuple(i);
        for (int j = 0; j < NSITES; j++) {
            jtup = indexToTuple(j);
            term1 += nu * B * b(itup, jtup) * I(itup, jtup, kappa);
        }
    }
    term1 /= NSITES;
    return term1;
}

int main() {
    setLength(SIZE);
    for (vector<vector<cpp_bin_float_100>> vv : eqOccNumTup) {
        for (vector<cpp_bin_float_100> v : vv) {
            for (cpp_bin_float_100 n : v) {
                eqOccNum.push_back(n);
            }
        }
    }
    eqOccNumTup.clear();
    cout << "setting up K" << endl;
    vector<vector<cpp_bin_float_100>> A = K();
    displayArray(A);
    cout << "setting up h" << endl;
    vector<cpp_bin_float_100> b = h();
    cout << "solving equation system" << endl;
    vector<cpp_bin_float_100> res = LUPSolve(LUPDecomposition(A, false), b);
    displayArray(res);
    int j[3] = {1,1,1};
    int k[3] = {1,1,0};
    cpp_bin_float_100 iij = I(j, k, res);
    cout << endl << sigmaDCequation(res) << endl;
    return 0;
};