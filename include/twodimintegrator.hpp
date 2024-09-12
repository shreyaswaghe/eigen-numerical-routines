#ifndef TWODIMINTEGRATOR_HPP_C3VVYDQU
#define TWODIMINTEGRATOR_HPP_C3VVYDQU

#include "onedimintegrator.hpp"

class TwoDimQuadrature : private OneDimQuadrature {
   public:
    /**
     * @brief Integrates the function on the rectangle domain
     * x -> x + L_x,
     * y -> y + L_y
     *
     * Through a one-dimensional 32 point Gauss-Legendre subroutine.
     *
     * @param function: f(x,y) -> double
     * @param x - for the bottom left corner of the domain
     * @param y - for the bottom left corner of the domain
     * @param L_x - the size in the x-domain
     * @param L_y - the size in the y-domain
     *
     * @return
     */
    double RectangleIntegrator(
        const std::function<double(const double &, const double &)> function,
        const double &x, const double &y, const double &L_x, const double &L_y);

    /**
     * @brief Integrates the function on the circular domain
     *  of radius R centered at (x,y)
     *
     * @param function - f(x,y) -> double
     * @param x - for the center of the circular domain
     * @param y - for the center of the circular domain
     * @param R - the radius of the circular domain
     *
     * @return
     */
    double CircleIntegrator(
        const std::function<double(const double &, const double &)> function,
        const double &x, const double &y, const double &R);

    double CircleIntegratorPolar(
        const std::function<double(const double &x, const double &y)> function,
        const double &x, const double &y, const double &R);
};

#endif /* end of include guard: TWODIMINTEGRATOR_HPP_C3VVYDQU */
