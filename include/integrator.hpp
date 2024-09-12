#ifndef INTEGRATOR

#define INTEGRATOR

#include <Eigen/Core>

template <typename T>
class ODEIntegrator {
   public:
    double timestep_init = 0.5;
    double timestep = timestep_init;
    double timestep_min = 1e-6;

    /**
     * @brief Runge-Kutta 4/5 explicit method
     *
     * @param x0 - the integration point
     * @param t - the current time
     * @param gradient - the function which computes gradients, of signature (x,
     * t, prealloc (size x))
     */
    void RKF45(Eigen::ArrayBase<T>& x0, const double t,
               const std::function<T(Eigen::Ref<Eigen::ArrayBase<T>>, double,
                                     Eigen::Ref<Eigen::ArrayBase<T>> prealloc)>
                   gradient);

    T FixedTimeRK5(const T&, const double t,
                   const std::function<T(T, double)> gradient);

    ODEIntegrator(int r, int c) {
        this->prealloc1 = Eigen::ArrayXd(r, c);
        this->prealloc2 = Eigen::ArrayXd(r, c);
        this->prealloc3 = Eigen::ArrayXd(r, c);
        this->prealloc4 = Eigen::ArrayXd(r, c);
        this->prealloc5 = Eigen::ArrayXd(r, c);
        this->prealloc6 = Eigen::ArrayXd(r, c);
    };

    void set_timestep(double t) { timestep = t; }

   private:
    T prealloc1, prealloc2, prealloc3, prealloc4, prealloc5, prealloc6;
    double errtol = 1e-6;
};

template <typename T>
void ODEIntegrator<T>::RKF45(
    Eigen::ArrayBase<T>& x, const double t,
    const std::function<T(Eigen::Ref<Eigen::ArrayBase<T>>, double,
                          Eigen::Ref<Eigen::ArrayBase<T>> prealloc)>
        gradient) {
    Eigen::ArrayBase<T> k0, k1, k2, k3, k4, k5, k6, grad;

    double h = 5 * this->timestep;
    double err_estimate = 0.0;
    double htemp = 0.0;
    while (true) {
        k0 = x;
        k1 = h * gradient(k0, t, prealloc1);

        k2 = x + 1.0 / 5.0 * k1;
        k2 = h * gradient(k2, t + 1.0 / 5.0 * h, prealloc2);

        k3 = x + 3.0 / 40.0 * k1 + 9.0 / 40.0 * k2;
        k3 = h * gradient(k3, t + 3.0 / 10.0 * h, prealloc3);

        k4 = x + 3.0 / 10.0 * k1 - 9.0 / 10.0 * k2 + 6.0 / 5.0 * k3;
        k4 = h * gradient(k4, t + 3.0 / 5.0 * h, prealloc4);

        k5 = x - 11.0 / 54.0 * k1 + 5.0 / 2.0 * k2 - 70.0 / 27.0 * k3 +
             35.0 / 27.0 * k4;
        k5 = h * gradient(k5, t + 1.0 * h, prealloc5);

        k6 = x + 1631.0 / 55296.0 * k1 + 175.0 / 512.0 * k2 +
             575.0 / 13824.0 * k3 + 44275.0 / 110592.0 * k4 +
             253.0 / 4096.0 * k5;
        k6 = h * gradient(k6, t + 7.0 / 8.0 * h, prealloc6);

        err_estimate =
            ((37.0 / 378.0 - 2825.0 / 27468.0) * k1 + 0.0 * k2 +
             (250.0 / 621.0 - 18575.0 / 48384.0) * k3 +
             (125.0 / 594.0 - 13525.0 / 55296.0) * k4 +
             (0.0 - 277.0 / 14336.0) * k5 + (512.0 / 1771.0 - 1.0 / 4.0) * k6)
                .cwiseAbs()
                .maxCoeff();

        err_estimate /= errtol;
        if (err_estimate <= 1.0) {
            break;
            h = std::min(h, 60.0);
        }

        htemp = 0.9 * h * pow(err_estimate, -0.20);
        h = (h >= 0.0 ? std::max(htemp, 0.1 * h) : std::min(htemp, 0.1 * h));
    }

    this->timestep = h;
    grad = 37.0 / 378.0 * k1 + 0.0 * k2 + 250.0 / 621.0 * k3 +
           125.0 / 594.0 * k4 + 0.0 * k5 + 512.0 / 1771.0 * k6;

    x += grad;
};

template <typename T>
T ODEIntegrator<T>::FixedTimeRK5(const T& x, const double t,
                                 std::function<T(T, double)> gradient) {
    double h = this->timestep;
    T k0, k1, k2, k3, k4, k5, k6, grad;

    k0 = x;
    k1 = h * gradient(k0, t, prealloc1);

    k2 = x + 1.0 / 5.0 * k1;
    k2 = h * gradient(k2, t + 1.0 / 5.0 * h, prealloc2);

    k3 = x + 3.0 / 40.0 * k1 + 9.0 / 40.0 * k2;
    k3 = h * gradient(k3, t + 3.0 / 10.0 * h, prealloc3);

    k4 = x + 3.0 / 10.0 * k1 - 9.0 / 10.0 * k2 + 6.0 / 5.0 * k3;
    k4 = h * gradient(k4, t + 3.0 / 5.0 * h, prealloc4);

    k5 = x - 11.0 / 54.0 * k1 + 5.0 / 2.0 * k2 - 70.0 / 27.0 * k3 +
         35.0 / 27.0 * k4;
    k5 = h * gradient(k5, t + 1.0 * h, prealloc5);

    k6 = x + 1631.0 / 55296.0 * k1 + 175.0 / 512.0 * k2 + 575.0 / 13824.0 * k3 +
         44275.0 / 110592.0 * k4 + 253.0 / 4096.0 * k5;
    k6 = h * gradient(k6, t + 7.0 / 8.0 * h, prealloc6);

    grad = 37.0 / 378.0 * k1 + 0.0 * k2 + 250.0 / 621.0 * k3 +
           125.0 / 594.0 * k4 + 0.0 * k5 + 512.0 / 1771.0 * k6;

    return x + grad;
}

#endif
