#ifndef ARRAYSUPPORT_HPP
#define ARRAYSUPPORT_HPP
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

int length;

void setLength(int l) {
    length = l;
}

int tupleToIndex(int i[3]) {
    return i[0] * i[1]*length + i[2]*length*length;
}

int tupleToIndex(vector<int> i) {
    return i[0] * i[1]*length + i[2]*length*length;
}

int pIdx(int i) {
    return i % length;
}
double acc(double** a[], int i[3]) {
    return a[i[0]][i[1]][i[2]] ;
}

double acc(vector<vector<vector<double>>> a, int i[3]) {
    return a[i[0]][i[1]][i[2]] ;
}

double acc(vector<vector<double>> a, int i[3], int j[3]) {
    return a[tupleToIndex(i)][tupleToIndex(j)];
}

double acc(vector<vector<double>> a, int i[3], vector<int> j) {
    return a[tupleToIndex(i)][tupleToIndex(j)];
}

double acc(vector<vector<double>> a, vector<int> i, int j[3]) {
    return a[tupleToIndex(i)][tupleToIndex(j)];
}

double acc(vector<vector<double>> a, vector<int> i, vector<int> j) {
    return a[tupleToIndex(i)][tupleToIndex(j)];
}

void set(double** a[], int i[3], double value) {
    a[i[0]][i[1]][i[2]] = value;
}

void set(vector<vector<vector<double>>> &a, int i[3], double value) {
    a[i[0]][i[1]][i[2]] = value;
}

vector<vector<int>> getNNidx(int i[3]) {
    vector<vector<int>> ret;
    ret.push_back(vector<int>({pIdx(i[0]-1), i[1], i[2]}));
    ret.push_back(vector<int>({pIdx(i[0]+1), i[1], i[2]}));
    ret.push_back(vector<int>({i[0], pIdx(i[1]-1), i[2]}));
    ret.push_back(vector<int>({i[0], pIdx(i[1]+1), i[2]}));
    ret.push_back(vector<int>({i[0], i[1], pIdx(i[2]-1)}));
    ret.push_back(vector<int>({i[0], i[1], pIdx(i[2]+1)}));

    return ret;
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

void displayArray(vector<vector<double>> a) {
    for (int i=0; i<a.size(); i++) {
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
#endif