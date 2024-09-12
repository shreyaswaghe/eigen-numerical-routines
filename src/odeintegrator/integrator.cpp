
/**
State gradient(const State &x, double mass);

State gradient(const State& x, double mass){
    int dimensions = constants::dimensions;
    int num_cells = constants::num_cells;
    double alpha = constants::alpha;

    State grad(dimensions, num_cells);
    Eigen::Vector2d dxi, xi, xj, diff;

    double norm_diff;
    grad.setZero();
    for(int i = 0; i < num_cells; ++i){
        dxi.setZero();
        for(int j = 0; j < num_cells; j++){
            if (i != j){
                xi = x.col(i);
                xj = x.col(j);
                diff = xi - xj;
                norm_diff =  diff.norm();
                dxi += (*kernel)(norm_diff)*(diff / norm_diff);
            }
        }
        dxi *= -mass;
        grad.col(i) = dxi;
    }

    grad.colwise() += alpha * 0 * Eigen::Vector2d::UnitX();

    return grad;
}
*/
/**
State RK4(State& x, double t0, double timestep){
    t0 = t0 + timestep;
    double mass = constants::cell_mass;

    State k1 = timestep * gradient(x, mass);
    State k2 = timestep * gradient(x + k1 * 0.5, mass);
    State k3 = timestep * gradient(x + k2 * 0.5, mass);
    State k4 = timestep * gradient(x + k3, mass);

    x = x + (k1/6.0) + (k2/3.0) + (k3/3.0) + (k4/6.0);
    return x;
}
*/
