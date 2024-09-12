#ifndef RANDOMGEN_HPP_OKGSH9FE
#define RANDOMGEN_HPP_OKGSH9FE

#include <RngStream.h>
#include <Eigen/Core>

namespace Random {
    inline void N01(RngStream& rng1, RngStream& rng2, double& out1, double& out2)
    {
        static double u1, u2;
        u1 = rng1.RandU01();
        u2 = rng2.RandU01();

        u1 = sqrt(-2 * log(u1));
        u2 = u2 * 2 * 3.141592653589;
        
        out1 = u1 * cos(u2);
        out2 = u2 * sin(u2);
    }
}

#endif /* end of include guard: RANDOMGEN_HPP_OKGSH9FE */
