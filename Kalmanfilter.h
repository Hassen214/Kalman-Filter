#include "Matrix.h"
using namespace std;

class Kalmanfilter {
private:
    Matrix F, B, H, Q, R, P, x_hat, x0;
    int n, m;
    double t, dt, t0;

public:
    Kalmanfilter(const Matrix& F, const Matrix& B, const Matrix& H,
                const Matrix& Q, const Matrix& R, const Matrix& P,
                const Matrix& x_hat)
        : F(F), B(B), H(H), Q(Q), R(R), P(P), x_hat(x_hat), x0(x_hat) {
        n = x_hat.row();
        m = H.row();
        t0 = 0.0;
        t = t0;
        dt = 0.0;
    }
	void predict(const Matrix& u) ;
    void update(const Matrix y, double dt, const Matrix A);
    Matrix value() { return x_hat; }
    double time() { return t; }
};