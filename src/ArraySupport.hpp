#ifndef ARRAYSUPPORT_HPP
#define ARRAYSUPPORT_HPP
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <boost/multiprecision/cpp_bin_float.hpp>
using namespace boost::multiprecision;
using namespace std;

int length;
int length2;

void setLength(int l) {
    length = l;
    length2 = length * length;
}

int tupleToIndex(int i[3]) {
    return i[0] + i[1]*length + i[2]*length2;
}

int tupleToIndex(vector<int> i) {
    return i[0] + i[1]*length + i[2]*length2;
}

int* indexToTuple(int i) {
    int* a = new int[3];
    a[2] = i % length;
    a[1] = i / length % length;
    a[0] = i / length2 % length;

    return a;
}

int mod(int i) {
    int r = i % length;
    return r < 0 ? r + length : r;
}
    
cpp_bin_float_100 acc(cpp_bin_float_100** a[], int i[3]) {
    return a[i[0]][i[1]][i[2]] ;
}

cpp_bin_float_100 acc(vector<vector<vector<cpp_bin_float_100>>> a, int i[3]) {
    return a[i[0]][i[1]][i[2]] ;
}

cpp_bin_float_100 acc(vector<vector<vector<cpp_bin_float_100>>> a, vector<int> i) {
    return a[i[0]][i[1]][i[2]] ;
}

cpp_bin_float_100 acc(vector<vector<cpp_bin_float_100>> a, int i[3], int j[3]) {
    return a[i[0] + i[1]*length + i[2]*length2][j[0] + j[1]*length + j[2]*length2];
}

cpp_bin_float_100 acc(vector<cpp_bin_float_100> a, int i[3]) {
    return a[i[0] + i[1]*length + i[2]*length2];
}

cpp_bin_float_100 acc(vector<cpp_bin_float_100> a, vector<int> i) {
    return a[i[0] + i[1]*length + i[2]*length2];
}

cpp_bin_float_100 acc(vector<vector<cpp_bin_float_100>> a, int i[3], vector<int> j) {
    return a[i[0] + i[1]*length + i[2]*length2][j[0] + j[1]*length + j[2]*length2];
}

cpp_bin_float_100 acc(vector<vector<cpp_bin_float_100>> a, vector<int> i, int j[3]) {
    return a[i[0] + i[1]*length + i[2]*length2][j[0] + j[1]*length + j[2]*length2];
}

cpp_bin_float_100 acc(vector<vector<cpp_bin_float_100>> a, vector<int> i, vector<int> j) {
    return a[i[0] + i[1]*length + i[2]*length2][j[0] + j[1]*length + j[2]*length2];
}

void setVal(cpp_bin_float_100** a[], int i[3], cpp_bin_float_100 value) {
    a[i[0]][i[1]][i[2]] = value;
}

void setVal(vector<vector<vector<cpp_bin_float_100>>> &a, int i[3], cpp_bin_float_100 value) {
    a[i[0]][i[1]][i[2]] = value;
}

vector<vector<int>> getNNidx(int i[3]) {
    vector<vector<int>> ret;
    ret.push_back(vector<int>({mod(i[0]-1), i[1], i[2]}));
    ret.push_back(vector<int>({mod(i[0]+1), i[1], i[2]}));
    ret.push_back(vector<int>({i[0], mod(i[1]-1), i[2]}));
    ret.push_back(vector<int>({i[0], mod(i[1]+1), i[2]}));
    ret.push_back(vector<int>({i[0], i[1], mod(i[2]-1)}));
    ret.push_back(vector<int>({i[0], i[1], mod(i[2]+1)}));

    return ret;
}

vector<vector<int>> getNNidx(vector<int> i) {
    int j[3] = {i[0], i[1], i[2]};
    return getNNidx(j);
}

vector<vector<int>> getNNNidx(int i[3]) {
    vector<vector<int>> ret;
    
    ret.push_back(vector<int>({mod(i[0]-2), i[1], i[2]}));
    ret.push_back(vector<int>({mod(i[0]-1), mod(i[1]-1), i[2]}));
    ret.push_back(vector<int>({mod(i[0]-1), mod(i[1]+1), i[2]}));
    ret.push_back(vector<int>({mod(i[0]-1), i[1], mod(i[2]-1)}));
    ret.push_back(vector<int>({mod(i[0]-1), i[1], mod(i[2]+1)}));
    ret.push_back(vector<int>({i[0], mod(i[1]-2), i[2]}));
    ret.push_back(vector<int>({mod(i[0]+1), mod(i[1]-1), i[2]}));
    ret.push_back(vector<int>({i[0], mod(i[1]-1), mod(i[2]-1)}));
    ret.push_back(vector<int>({i[0], mod(i[1]-1), mod(i[2]+1)}));
    ret.push_back(vector<int>({mod(i[0]+2), i[1], i[2]}));
    ret.push_back(vector<int>({mod(i[0]+1), mod(i[1]+1), i[2]}));
    ret.push_back(vector<int>({mod(i[0]+1), i[1], mod(i[2]-1)}));
    ret.push_back(vector<int>({mod(i[0]+1), i[1], mod(i[2]+1)}));
    ret.push_back(vector<int>({i[0], mod(i[1]+2), i[2]}));
    ret.push_back(vector<int>({i[0], mod(i[1]+1), mod(i[2]-1)}));
    ret.push_back(vector<int>({i[0], mod(i[1]+1), mod(i[2]+1)}));
    ret.push_back(vector<int>({i[0], i[1], mod(i[2]-2)}));
    ret.push_back(vector<int>({i[0], i[1], mod(i[2]+2)}));

    return ret;
}

bool eqIdx(vector<int> a, vector<int> b) {
    return (a[0] == b[0]) && (a[1] == b[1]) && (a[2] == b[2]);
}

bool eqIdx(int a[3], vector<int> b) {
    return (a[0] == b[0]) && (a[1] == b[1]) && (a[2] == b[2]);
}

bool eqIdx(vector<int> a, int b[3]) {
    return (a[0] == b[0]) && (a[1] == b[1]) && (a[2] == b[2]);
}

bool eqIdx(int a[3], int b[3]) {
    return (a[0] == b[0]) && (a[1] == b[1]) && (a[2] == b[2]);
}

bool compare(cpp_bin_float_100 a, cpp_bin_float_100 b, cpp_bin_float_100 epsilon) {
    cpp_bin_float_100 diff = a - b;
    return (diff < epsilon) && (-diff < epsilon);
}

void displayArray(int a[], int size) {
    for (int i=0; i<size; i++) {
        cout << a[i] << " ";
    }
    cout << endl;
}

void displayArray(vector<cpp_bin_float_100> a) {
    for (int i=0; i<a.size(); i++) {
        cout << (cpp_bin_float_100)a[i] << " ";
    }
    cout << endl;
}

void displayArray(vector<vector<cpp_bin_float_100>> a) {
    for (int i=0; i<a.size(); i++) {
        displayArray(a[i]);
        cout << endl;
    }
}

void displayArray(vector<vector<vector<cpp_bin_float_100>>> a) {
    for (int i=0; i<a.size(); i++) {
        displayArray(a[i]);
        cout << endl << endl;
    }
}

void displayArray(vector<int> a, int size) {
    for (int i=0; i<size; i++) {
        cout << (int)a[i] << " ";
    }
    cout << endl;
}
#endif