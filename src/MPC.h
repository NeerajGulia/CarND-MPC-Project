#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class Constants
{
  public:
    // TODO: Set the timestep length and duration
    static const size_t N = 8;
    static constexpr  double DT = 0.15; //seconds
    static const int POLY_ORDER = 3;

    // This value assumes the model presented in the classroom is used.
    //
    // It was obtained by measuring the radius formed by running the vehicle in the
    // simulator around in a circle with a constant steering angle and velocity on a
    // flat terrain.
    //
    // Lf was tuned until the the radius formed by the simulating the model
    // presented in the classroom matched the previous radius.
    //
    // This is the length from front to CoG that has a similar radius.
    static constexpr  double LF = 2.67;

    static constexpr double FACTOR_MILES_TO_METER_PER_SEC = 0.44704;
    // Both the reference cross track and orientation errors are 0.
    // The reference velocity is set to 40 mph.
    static constexpr  double REF_V = 70 * FACTOR_MILES_TO_METER_PER_SEC; //convert from m/h to m/s
    
    // The solver takes all the state variables and actuator
    // variables in a singular vector. Thus, we should to establish
    // when one variable starts and another ends to make our lifes easier.
    static const size_t X_START = 0;
    static const size_t Y_START = X_START + N;
    static const size_t PSI_START = Y_START + N;
    static const size_t V_START = PSI_START + N;
    static const size_t CTE_START = V_START + N;
    static const size_t EPSI_START = CTE_START + N;
    static const size_t DELTA_START = EPSI_START + N;
    static const size_t A_START = DELTA_START + N - 1;
};

class Solution
{
  public:

    Solution();

    double steer_angle;
    double throttle;
    vector<double> mpc_x_vals;
    vector<double> mpc_y_vals;
};

class MPC {
 public:

  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
 Solution Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
