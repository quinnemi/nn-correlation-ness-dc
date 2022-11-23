#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

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

bool compare(double a, double b, double epsilon) {
    double diff = a - b;
    return (diff < epsilon) && (-diff < epsilon);
}

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

void displayArray(vector<vector<double>> a, int size) {
    for (int i=0; i<size; i++) {
        displayArray(a[i], a[i].size());
        cout << endl;
    }
}

void displayArray(vector<int> a, int size) {
    for (int i=0; i<size; i++) {
        cout << (int)a[i] << " ";
    }
    cout << endl;
}