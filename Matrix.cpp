#include <iostream>
#include <vector>
#include <iomanip>
#include "Matrix.h"
#include <stdexcept>
#include <cmath>
#include <utility>  // for std::move
#include <stdexcept>

using namespace std;
Matrix::Matrix() : rows(0), cols(0), mat({}) {} ;
// Default Constructor
Matrix::Matrix(int r, int c) : rows(r), cols(c), mat(r, vector<double>(c, 0)) {}

// Parameterized Constructor
Matrix::Matrix(int r, int c, const vector<vector<double>>& values) : rows(r), cols(c), mat(values) {}

// Copy Constructor
Matrix::Matrix(const Matrix& other) : rows(other.rows), cols(other.cols), mat(other.mat) {}

// Move Constructor
Matrix::Matrix(Matrix&& other) : rows(other.rows), cols(other.cols), mat(move(other.mat)) {
    other.rows = 0;
    other.cols = 0;
}
int Matrix::row()const {
    return rows;
}
int Matrix::columns() const {
    return cols;
}
// Display matrix
void Matrix::display() const {
    for (const auto &row : mat) {
        for (double val : row)
            cout << setw(8) << fixed << setprecision(2) << val << " ";
        cout << endl;
    }
}

// Matrix addition
Matrix Matrix::operator+(const Matrix &other) const {
    if (rows != other.rows || cols != other.cols) {
        throw invalid_argument("Matrix dimensions do not match for addition.");
    }
    Matrix result(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result.mat[i][j] = mat[i][j] + other.mat[i][j];
    return result;
}

// Matrix subtraction
Matrix Matrix::operator-(const Matrix &other) const {
    if (rows != other.rows || cols != other.cols) {
        throw invalid_argument("Matrix dimensions do not match for subtraction.");
    }
    Matrix result(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result.mat[i][j] = mat[i][j] - other.mat[i][j];
    return result;
}


// Matrix equality
bool Matrix::operator==(const Matrix& other) const {
    if (rows != other.rows || cols != other.cols) {
        return false;
    }
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            if (mat[i][j] != other.mat[i][j])
                return false;
        }
    }
    return true;
}

// Matrix multiplication
Matrix Matrix::operator*(const Matrix &other) const {
    if (cols != other.rows) {
        throw invalid_argument("Matrix dimensions do not match for multiplication.");
    }
    Matrix result(rows, other.cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < other.cols; j++) {
            for (int k = 0; k < cols; k++) {
                result.mat[i][j] += mat[i][k] * other.mat[k][j];
            }
        }
    }
    return result;
}

// Compute inverse using Gaussian elimination
Matrix Matrix::inverse() const {
    if (rows != cols) {
        throw invalid_argument("Only square matrices can be inverted.");
    }

    int n = rows;
    Matrix aug(n, 2 * n);

    // Create augmented matrix [A | I]
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            aug.mat[i][j] = mat[i][j]; // Copy original matrix
        for (int j = n; j < 2 * n; j++)
            aug.mat[i][j] = (j - n == i) ? 1 : 0; // Identity matrix
    }

    // Perform Gaussian elimination
    for (int i = 0; i < n; i++) {
        // Make diagonal element 1
        double diag = aug.mat[i][i];
        if (fabs(diag) < 1e-9)
            throw runtime_error("Matrix is singular and cannot be inverted.");
        for (int j = 0; j < 2 * n; j++)
            aug.mat[i][j] /= diag;

        // Make other rows' column i zero
        for (int k = 0; k < n; k++) {
            if (k == i) continue;
            double factor = aug.mat[k][i];
            for (int j = 0; j < 2 * n; j++)
                aug.mat[k][j] -= factor * aug.mat[i][j];
        }
    }

    // Extract inverse matrix
    Matrix inv(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv.mat[i][j] = aug.mat[i][j + n];

    return inv;
}
Matrix Matrix::transpose()  const{
    Matrix A = Matrix(cols,rows) ;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            A.mat[j][i] = mat[i][j];
        }
    }
    return A;
}
double Matrix::element(int i,int j) const {
    if (i>=rows or j>=cols) {
        throw std::out_of_range("Matrix element index out of range");

    }
    return mat[i][j] ;
}
void Matrix::assign(int i,int j,double t){
    mat[i][j] = t ;
}

Matrix Matrix::operator=(const Matrix& other) {
    int row= other.row();
    int col = other.columns();
    Matrix A(row, col);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            A.mat[i][j] = other.mat[i][j];
        }
    }
    return A;
}

Matrix::Matrix(int row) {

    rows = row;
    cols = row;
    vector<vector<double>> values(row, vector<double>(row, 0));
    for (int i = 0; i < row; i++) {
        values[i][i] = 1;
    }
    mat = values ;
} ;
double& Matrix::operator()(int i,int j){
    return mat[i][j] ;
}