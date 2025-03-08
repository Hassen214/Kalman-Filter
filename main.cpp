#include <iostream>
#include "Matrix.h"
#include <random>
#include "Kalmanfilter.h"
using namespace std;

// Helper function to generate Gaussian noise
double gaussianNoise(double mean, double stddev) {
    static default_random_engine generator;
    normal_distribution<double> dist(mean, stddev);
    return dist(generator);
}

int main() {
    double dt = 0.1;  // Time step (0.1 seconds)

    // System matrices ----------------------------------------------------------
    // State transition matrix (position & velocity)
    Matrix F = Matrix(2);
    F.assign(0,1,dt);

    // Control matrix (no control input)
    Matrix B(2, 1);  // Zero matrix

    // Measurement matrix (only measures position)
    Matrix H(1, 2);
    H(0, 0) = 1.0;

    // Process noise covariance (small uncertainties)
    Matrix Q(2, 2);
    Q(0, 0) = 0.01;  // Position process noise
    Q(1, 1) = 0.01;  // Velocity process noise

    // Measurement noise covariance
    Matrix R(1, 1);
    R(0, 0) = 1.0;   // Position measurement noise

    // Initial state [position; velocity]
    Matrix x0(2, 1);
    x0(0, 0) = 0.0;  // Initial position
    x0(1, 0) = 1.0;  // Initial velocity (1 m/s)

    // Initial error covariance
    Matrix P(2, 2);
    P(0, 0) = 1.0;    // Large initial position uncertainty
    P(1, 1) = 1.0;    // Large initial velocity uncertainty

    // Create Kalman filter -----------------------------------------------------
    Kalmanfilter kf(F, B, H, Q, R, P, x0);

    // Simulation parameters ----------------------------------------------------
    Matrix true_state = x0;      // True object state (hidden)
    Matrix u(1, 1);              // Control input (none)
    const int steps = 10;        // Number of simulation steps

    // Store true position and measurements for comparison
    vector<double> true_positions, measurements, estimates;

    // Simulation loop ----------------------------------------------------------
    for(int i = 0; i < steps; ++i) {

        // True state evolution (constant velocity with slight process noise)
        true_state = F * true_state;
        true_state(0, 0) += gaussianNoise(0, 0.05);  // Add position noise
        true_state(1, 0) += gaussianNoise(0, 0.01);  // Add velocity noise

        // Generate noisy measurement
        Matrix y(1, 1);
        y(0, 0) = (H * true_state)(0,0) + gaussianNoise(0, 1.0);

        // Kalman filter steps
        kf.predict(u);
        // Prediction
        kf.update(y, dt, H);
        // Update with measurement

        // Store results
        true_positions.push_back(true_state(0, 0));
        measurements.push_back(y(0, 0));
        estimates.push_back(kf.value()(0, 0));

        // Print results
        cout << "Step " << i + 1 << ":\n"
             << "True position: " << true_state(0, 0) << "\n"
             << "Measurement: " << y(0, 0) << "\n"
             << "Estimate: " << kf.value()(0, 0) << "\n"
             << "Velocity estimate: " << kf.value()(1, 0) << "\n\n";
    }

    return 0;
}