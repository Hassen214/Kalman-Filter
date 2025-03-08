#pragma once
#include <vector>
#include <iostream>

using namespace std;

class Matrix {
private:
    vector<vector<double>> mat;
    int rows, cols;

public:
    // Constructors
    Matrix(int r, int c);
    Matrix() ;
    Matrix(int r, int c, const vector<vector<double>>& values);
    Matrix(const Matrix& other);  // Copy constructor
    Matrix(int row) ;
    Matrix(Matrix&& other);       // Move constructor
    // Display function
    void display() const;
    int row() const ;
    int columns() const  ;
    double& operator()(int i,int j);

    // Operator overloads
    bool operator==(const Matrix& other) const;
    Matrix operator=(const Matrix& other);
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(const Matrix& other) const;
    Matrix transpose() const ;
    // Matrix inverse function
    Matrix inverse() const ;
    double element(int i,int j)  const ;
    Matrix identity(int row) const;
    void assign(int i,int j,double t) ;
    Matrix identity(int row) ;

};
