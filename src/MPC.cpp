#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;
    size_t t;

    // The part of the cost based on the reference state.
    for (t = 0; t < Constants::N; t++) {
      fg[0] += 2 * CppAD::pow(vars[Constants::CTE_START+ t], 2);
      fg[0] += 20 * CppAD::pow(vars[Constants::EPSI_START + t], 2);
      fg[0] += 20 * CppAD::pow(vars[Constants::V_START + t] - Constants::REF_V, 2);
    }

    // Minimize the use of actuators.
    for (t = 0; t < Constants::N - 1; t++) {
      fg[0] += 450 * CppAD::pow(vars[Constants::DELTA_START + t], 2);
      fg[0] += 20 * CppAD::pow(vars[Constants::A_START + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (t = 0; t < Constants::N - 2; t++) {
      fg[0] += 450 * CppAD::pow(vars[Constants::DELTA_START + t + 1] - vars[Constants::DELTA_START + t], 2);
      fg[0] += 20 * CppAD::pow(vars[Constants::A_START + t + 1] - vars[Constants::A_START + t], 2);
    }

    //
    // Setup Constraints
    //
    // NOTE: In this section you'll setup the model constraints.

    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + Constants::X_START] = vars[Constants::X_START];
    fg[1 + Constants::Y_START] = vars[Constants::Y_START];
    fg[1 + Constants::PSI_START] = vars[Constants::PSI_START];
    fg[1 + Constants::V_START] = vars[Constants::V_START];
    fg[1 + Constants::CTE_START] = vars[Constants::CTE_START];
    fg[1 + Constants::EPSI_START] = vars[Constants::EPSI_START];

    // The rest of the constraints
    for (t = 1; t < Constants::N; t++) {
      // The state at time t+1 .
      AD<double> x1 = vars[Constants::X_START + t];
      AD<double> y1 = vars[Constants::Y_START + t];
      AD<double> psi1 = vars[Constants::PSI_START + t];
      AD<double> v1 = vars[Constants::V_START + t];
      AD<double> cte1 = vars[Constants::CTE_START + t];
      AD<double> epsi1 = vars[Constants::EPSI_START + t];

      // The state at time t.
      AD<double> x0 = vars[Constants::X_START + t - 1];
      AD<double> y0 = vars[Constants::Y_START + t - 1];
      AD<double> psi0 = vars[Constants::PSI_START + t - 1];
      AD<double> v0 = vars[Constants::V_START + t - 1];
      AD<double> cte0 = vars[Constants::CTE_START + t - 1];
      AD<double> epsi0 = vars[Constants::EPSI_START + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[Constants::DELTA_START + t - 1];
      AD<double> a0 = vars[Constants::A_START + t - 1];

      AD<double> f0;
      AD<double> sum;
      AD<double> power = 1.0;

      //degree polynomial
      for(int i = 0; i <= Constants::POLY_ORDER; i++)
      {
        f0 += coeffs[i] * power;
        if(i + 1 <= Constants::POLY_ORDER)
          sum += ( (i+1) * coeffs(i + 1) * power);
        power *= x0;
      }
      // AD<double> f0 = coeffs[0] + coeffs[1]*x0 + coeffs[2]*x02 + coeffs[3]*x02*x0;
      // AD<double> psides0 = CppAD::atan(coeffs[1] + 2*coeffs[2]*x0  + 3*coeffs[3]*x02);
      AD<double> psides0 = CppAD::atan(sum);
      // Here's `x` to get you started.
      // The idea here is to constraint this value to be 0.
      //
      // Recall the equations for the model:
      // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
      // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
      // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
      // v_[t+1] = v[t] + a[t] * dt
      // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
      // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
      AD<double> cal = delta0 / Constants::LF * Constants::DT;
      fg[1 + Constants::X_START + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * Constants::DT);
      fg[1 + Constants::Y_START + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * Constants::DT);
      fg[1 + Constants::PSI_START + t] = psi1 - (psi0 + v0 * cal);
      fg[1 + Constants::V_START + t] = v1 - (v0 + a0 * Constants::DT);
      fg[1 + Constants::CTE_START + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * Constants::DT));
      fg[1 + Constants::EPSI_START + t] = epsi1 - ((psi0 - psides0) + v0 * cal);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

Solution::Solution(){}

Solution MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = Constants::N * 6 + (Constants::N - 1) * 2;
  // TODO: Set the number of constraints
  size_t n_constraints = Constants::N * 6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (i = 0; i < n_vars; i++) {
    vars[i] = 0.0;
  }

  // Set the initial variable values
  vars[Constants::X_START] = x;
  vars[Constants::Y_START] = y;
  vars[Constants::PSI_START] = psi;
  vars[Constants::V_START] = v;
  vars[Constants::CTE_START] = cte;
  vars[Constants::EPSI_START] = epsi;

  // Lower and upper limits for x
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // TODO: Set lower and upper limits for variables.

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (i = 0; i < Constants::DELTA_START; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (i = Constants::DELTA_START; i < Constants::A_START; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (i = Constants::A_START; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[Constants::X_START] = x;
  constraints_lowerbound[Constants::Y_START] = y;
  constraints_lowerbound[Constants::PSI_START] = psi;
  constraints_lowerbound[Constants::V_START] = v;
  constraints_lowerbound[Constants::CTE_START] = cte;
  constraints_lowerbound[Constants::EPSI_START] = epsi;

  constraints_upperbound[Constants::X_START] = x;
  constraints_upperbound[Constants::Y_START] = y;
  constraints_upperbound[Constants::PSI_START] = psi;
  constraints_upperbound[Constants::V_START] = v;
  constraints_upperbound[Constants::CTE_START] = cte;
  constraints_upperbound[Constants::EPSI_START] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  Solution sol;
  sol.steer_angle = solution.x[Constants::DELTA_START];
  sol.throttle = solution.x[Constants::A_START];

  //rest all contains the x and y coordinates
  for(i=0; i<Constants::N; i++)
  {
    sol.mpc_x_vals.push_back(solution.x[Constants::X_START + i]);
    sol.mpc_y_vals.push_back(solution.x[Constants::Y_START + i]);
  }
  return sol;
}
