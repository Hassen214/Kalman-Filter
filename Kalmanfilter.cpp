#include "Matrix.h"
#include <vector>
#include "Kalmanfilter.h"
using namespace std;

void Kalmanfilter::predict(const Matrix& u) {
    // Predict state
    Matrix F_prime = F.transpose() ;
    x_hat = F * x_hat + B * u;
    // Predict error covariance
    P = F * P * F_prime + Q;
}


void Kalmanfilter::update(const Matrix y, double dt, const Matrix A) {
    // Compute residual

    Matrix A_prime = A.transpose() ;


    A_prime = P*A_prime;
    Matrix residual = y - A * x_hat;

    Matrix S = A * A_prime + R;


    // Compute Kalman gain
    Matrix K = P * A.transpose() * S.inverse();
    // Update state estimate
    x_hat = x_hat + K * residual;
    // Update error covariance
    Matrix I(n);
    P = (I - K * A) * P;
    // Update time
    t += dt;
}


