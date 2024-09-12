#include "twodimintegrator.hpp"

double TwoDimQuadrature::RectangleIntegrator(
    const std::function<double(const double &, const double &)> function,
    const double &x, const double &y, const double &L_x, const double &L_y) {
    const double x_width = L_x;
    const double x_hi = x + x_width;

    const double y_width = L_y;
    const double y_hi = y + y_width;

    /*Since this is a rectangle, computing the range of x and y values is very
     * easy*/
    Eigen::Vector<double, 32> xvals = this->GL32_abscissae.unaryExpr(
        [&x_hi, &x_width](double x) { return x_hi + x_width * 0.5 * (x - 1); });

    Eigen::Vector<double, 32> yvals = this->GL32_abscissae.unaryExpr(
        [&y_hi, &y_width](double y) { return y_hi + y_width * 0.5 * (y - 1); });

    /*Integrate over the y component per value of x, collate into a single
     * vector */
    Eigen::Vector<double, 32> int_yvals =
        xvals
            /*This unaryExpr evaluates the integral over y keeping x constant*/
            .unaryExpr([this, &yvals, &y_width, &function](double x) {
                return 0.5 * y_width *
                       this->GL32_weights
                           .cwiseProduct(yvals.unaryExpr(
                               std::bind(function, x, std::placeholders::_1)))
                           .sum();
            });

    return 0.5 * x_width * this->GL32_weights.cwiseProduct(int_yvals).sum();
}

double TwoDimQuadrature::CircleIntegratorPolar(
    const std::function<double(const double &x, const double &y)> function,
    const double &x, const double &y, const double &R) {
    double theta_hi = 2 * 3.141592653589;
    double theta_width = theta_hi;

    /* Scale -1,1 to 0, 2 pi */
    Eigen::Vector<double, 32> theta_vals =
        this->GL32_abscissae.unaryExpr([&theta_hi, &theta_width](double x) {
            return theta_hi + theta_width * 0.5 * (x - 1);
        });

    /* Represents the integration fixing theta, and varying R */
    Eigen::Vector<double, 32> intr_vals =
        theta_vals.unaryExpr([this, &function, &R](double theta) {
            std::function<double(const double &)> fnc{
                [&theta, &function](double r) {
                    return function(r * cos(theta), r * sin(theta));
                }};

            return this->GaussLegendre32Point(fnc, 0, R);
        });

    return 0.5 * theta_width * this->GL32_weights.cwiseProduct(intr_vals).sum();
}

double TwoDimQuadrature::CircleIntegrator(
    const std::function<double(const double &x, const double &y)> function,
    const double &x, const double &y, const double &R) {
    const double x_width = 2 * R;
    const double x_hi = x + R;

    /*Scaled sweep on the x-range*/
    Eigen::Vector<double, 32> xvals =
        this->GL32_abscissae.unaryExpr([&x_hi, &x_width](const double x) {
            return x_hi + x_width * 0.5 * (x - 1);
        });

    Eigen::Vector<double, 32> int_yvals =
        xvals
            /*This unaryExpr evaluates the integral over y keeping x constant */
            .unaryExpr([this, &function, &x, &y, &R](double xx) {
                double y_halfwidth =
                    sqrt(R * R - (xx - x) * (xx - x));  // half width
                return this->GaussLegendre32Point(
                    std::bind(function, xx, std::placeholders::_1),
                    y - y_halfwidth, y + y_halfwidth);
            });

    return 0.5 * x_width * this->GL32_weights.cwiseProduct(int_yvals).sum();
}

