#include <cmath>
#include <iostream>
#include <iomanip>

#include "include/onedimintegrator.hpp"
#include "include/twodimintegrator.hpp"

double x(double y) { return exp(y); }

double id(double x, double y) { return 1.0; }

int main(int argc, char *argv[]) {
    const double pi = 3.14159265389;

    std::cout << std::fixed << std::showpoint;
    std::cout << std::setprecision(10);

    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << std::endl;

    TwoDimQuadrature b;
    /* The Rectangular grid 2D integration */
    std::cout << "Area of circle of radius 1.0: "
              << b.CircleIntegrator(id, 0.0, 0.0, 1.0) << ". Should be " << pi
              << std::endl;
    std::cout << "Area of circle of radius 5.0: "
              << b.CircleIntegrator(id, 0.0, 0.0, 5.0) << ". Should be "
              << pi * 25 << std::endl;
    std::cout << "Area of circle of radius 10.0: "
              << b.CircleIntegrator(id, 0.0, 0.0, 10.0) << ". Should be "
              << pi * 100 << std::endl;
    std::cout << "Area of circle of radius 20.0: "
              << b.CircleIntegrator(id, 0.0, 0.0, 20.0) << ". Should be "
              << pi * 400 << std::endl;

    std::cout << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << std::endl;
        
    /* The polar grid 2D integration */
    std::cout << "(Polar) Area of circle of radius 1.0: "
              << b.CircleIntegratorPolar(id, 0.0, 0.0, 1.0) << ". Should be " << pi
              << std::endl;
    std::cout << "(Polar) Area of circle of radius 5.0: "
              << b.CircleIntegratorPolar(id, 0.0, 0.0, 5.0) << ". Should be "
              << pi * 25 << std::endl;
    std::cout << "(Polar) Area of circle of radius 10.0: "
              << b.CircleIntegratorPolar(id, 0.0, 0.0, 10.0) << ". Should be "
              << pi * 100 << std::endl;
    std::cout << "(Polar) Area of circle of radius 20.0: "
              << b.CircleIntegratorPolar(id, 0.0, 0.0, 20.0) << ". Should be "
              << pi * 400 << std::endl;

    std::cout << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << std::endl;
     
    OneDimQuadrature a{};
    /* 1D Integration */
    std::cout << "int_0^1 e^x dx = " << a.GaussLegendre32Point(x, 0.0, 1.0)
              << ". Should be " << exp(1.0) - 1.0 << std::endl;
    std::cout << "int_0^5 e^x dx = " << a.GaussLegendre32Point(x, 0.0, 5.0)
              << ". Should be " << exp(5.0) - 1.0 << std::endl;
    std::cout << "int_0^10 e^x dx = " << a.GaussLegendre32Point(x, 0.0, 10.0)
              << ". Should be " << exp(10.0) - 1.0 << std::endl;
    std::cout << "int_0^20 e^x dx = " << a.GaussLegendre32Point(x, 0.0, 20.0)
              << ". Should be " << exp(20.0) - 1.0 << std::endl;

    return 0;
}
