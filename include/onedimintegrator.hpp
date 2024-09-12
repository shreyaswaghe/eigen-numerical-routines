#ifndef QUADRATURE_HPP_L7VWULHA
#define QUADRATURE_HPP_L7VWULHA

#include <Eigen/Core>

/**
 * @brief This class houses two numerical integration methods: adaptive
 * trapezoidal rule and
 */
class OneDimQuadrature {
   public:
    /**
     * @brief This function performs an adaptive trapezoidal integration for the
     * 1D function function on [low, high]
     *
     * @param function: f(double) -> double
     * @param low: the lower range of integration
     * @param high: the higher range of integration
     *
     * @return int_a^b f(x) d x
     */
    double AdaptiveTrapezoidal(const std::function<double(double)> &function,
                               const double &low, const double &high);

    /**
     * @brief This function performs a 32 Point Gauss Legendre quadrature
     * routine
     *
     * @param function: f(double) -> double
     * @param low: the lower range of integration
     * @param high: the higher range of integration
     *
     * @return int_a^b f(x) dx
     */
    double GaussLegendre32Point(const std::function<double(double)> &function,
                                const double &low, const double &high);

    OneDimQuadrature() {
        this->GL32_weights << 0.0965400885147278, 0.0965400885147278,
            0.0956387200792749, 0.0956387200792749, 0.0938443990808046,
            0.0938443990808046, 0.0911738786957639, 0.0911738786957639,
            0.0876520930044038, 0.0876520930044038, 0.0833119242269467,
            0.0833119242269467, 0.0781938957870703, 0.0781938957870703,
            0.0723457941088485, 0.0723457941088485, 0.0658222227763618,
            0.0658222227763618, 0.0586840934785355, 0.0586840934785355,
            0.0509980592623762, 0.0509980592623762, 0.0428358980222267,
            0.0428358980222267, 0.0342738629130214, 0.0342738629130214,
            0.0253920653092621, 0.0253920653092621, 0.0162743947309057,
            0.0162743947309057, 0.0070186100094701, 0.0070186100094701;

        this->GL32_abscissae << -0.0483076656877383, 0.0483076656877383,
            -0.1444719615827965, 0.1444719615827965, -0.2392873622521371,
            0.2392873622521371, -0.3318686022821277, 0.3318686022821277,
            -0.4213512761306353, 0.4213512761306353, -0.5068999089322294,
            0.5068999089322294, -0.5877157572407623, 0.5877157572407623,
            -0.6630442669302152, 0.6630442669302152, -0.7321821187402897,
            0.7321821187402897, -0.7944837959679424, 0.7944837959679424,
            -0.8493676137325700, 0.8493676137325700, -0.8963211557660521,
            0.8963211557660521, -0.9349060759377397, 0.9349060759377397,
            -0.9647622555875064, 0.9647622555875064, -0.9856115115452684,
            0.9856115115452684, -0.9972638618494816, 0.9972638618494816;
    }

   protected:
    Eigen::Vector<double, 32> GL32_weights;
    Eigen::Vector<double, 32> GL32_abscissae;

   private:
    /*For Adaptive Trapezoidal rule*/
    /**
     * @brief This stores the trapezoidal step.
     */
    double s;

    /**
     * @brief Computes the trapezoidal approximation to the the integral int_a^b
     * f(x) dx
     *
     * @param function: f(double) -> double
     * @param low: the lower range of integration
     * @param high: the higher range of integration
     * @param N: The number of points to include in the trapezoidal integration.
     *
     * @return int_a^b f(x) dx
     */
    double trapzd(const std::function<double(const double &)> &function,
                  const double &low, const double &high, const int &N);
};

#endif /* end of include guard: QUADRATURE_HPP_L7VWULHA */
