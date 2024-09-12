#include "onedimintegrator.hpp"

#include <cmath>

double OneDimQuadrature::AdaptiveTrapezoidal(
    const std::function<double(double)> &function, const double &low,
    const double &high) {
    double prev, curr;
    int j = 4;
    do {
        j++;
        prev = this->trapzd(function, low, high, j);
        curr = this->trapzd(function, low, high, j + 1);
    } while (std::abs(prev - curr) <= 1e-6);
    return curr;
}

double OneDimQuadrature::GaussLegendre32Point(
    const std::function<double(double)> &function, const double &low,
    const double &high) {
    double width = (high - low);
    Eigen::Vector<double, 32> fvals =
        this->GL32_abscissae
            .unaryExpr(
                [&high, &width](const double x) { return high + width * 0.5 * (x-1) ; })
            .unaryExpr([&function](double y) { return function(y); });
    return width * 0.5 * this->GL32_weights.cwiseProduct(fvals).sum();
}

double OneDimQuadrature::trapzd(
    const std::function<double(const double &)> &function, const double &low,
    const double &high, const int &N) {
    if (N == 1) {
        return (this->s =
                    0.5 * (high - low) * (function(low) + function(high)));
    } else {
        int it, j;
        for (it = 1, j = 1; j < N - 1; j++) it <<= 1;
        double tnm = it;
        double del = (high - low) / tnm;
        double x = low + 0.5 * del;

        double sum;
        for (sum = 0.0, j = 0; j < it; j++, x += del) sum += function(x);
        return (this->s = 0.5 * (this->s + (high - low) * sum / tnm));
    }
}
