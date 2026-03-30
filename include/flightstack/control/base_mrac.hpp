/***********************************************************************************************************************
 * Copyright (c) 2024 Mattia Gramuglia, Giri M. Kumar, Andrea L'Afflitto. All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
 * following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following
 *    disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the
 *    following disclaimer in the documentation and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote
 *    products derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **********************************************************************************************************************/

/***********************************************************************************************************************
 * File:        base_mrac.hpp
 * Author:      Mattia Gramuglia
 * Date:        June 17, 2025
 * For info:    Andrea L'Afflitto 
 *              a.lafflitto@vt.edu
 * 
 * Description: Common/base functions for MRAC controllers.
 * 
 * GitHub:    https://github.com/andrealaffly/ACSL-flightstack.git
 **********************************************************************************************************************/

#pragma once

#include <cmath>
#include <utility>
#include <type_traits>

#include <Eigen/Dense>

namespace base_mrac {

namespace funnel {

// Compute the H function
template <typename VectorType, typename MatrixType>
inline double computeHfunctionFunnel(double eta_max, double eta_funnel, const VectorType& e, const MatrixType& M)
{
  static_assert(static_cast<int>(VectorType::ColsAtCompileTime) == 1, "e must be a column vector");
  static_assert(
    static_cast<int>(MatrixType::RowsAtCompileTime) == static_cast<int>(MatrixType::ColsAtCompileTime),
    "M must be a square matrix"
  );
  static_assert(
    static_cast<int>(VectorType::RowsAtCompileTime) == static_cast<int>(MatrixType::RowsAtCompileTime),
    "Dimension mismatch between e and M"
  );

  // eᵀ * M * e
  double eT_M_e = e.transpose() * M * e;

  // funnel_diameter = ηₘₐₓ - η²
  double funnel_diameter = eta_max - std::pow(eta_funnel, 2);

  // H = funnel_diameter - eᵀMe
  double H_function = funnel_diameter - eT_M_e;

  return H_function;
}

// Compute the V_e function
template <typename VectorType, typename MatrixType>
inline double computeVeFunnel(const VectorType& e, const MatrixType& P, double H_function)
{
  static_assert(static_cast<int>(VectorType::ColsAtCompileTime) == 1, "e must be a column vector");
  static_assert(
    static_cast<int>(MatrixType::RowsAtCompileTime) == static_cast<int>(MatrixType::ColsAtCompileTime),
    "P must be a square matrix"
  );
  static_assert(
    static_cast<int>(VectorType::RowsAtCompileTime) == static_cast<int>(MatrixType::RowsAtCompileTime),
    "Dimension mismatch between e and P"
  );

  // Compute eᵀ * P * e
  double eTPe = (e.transpose() * P * e).value();

  // Ve_funnel = eᵀ * P * e / H
  double Ve_funnel = eTPe / H_function;
  
  return Ve_funnel;
}

/*
  Compute the lambda saturation value using the less conservative method:

  \\lambda_{\\rm sat}(t, e) &= \\max \\left( 0, \\frac{H(t,e) [\\lambda_{\\rm min}(Q + V_e(t,e)Q_{\\rm M}
  - \\frac{2\\overline{\\xi}_{\\rm d}}{\\|e\\|} (P + V_e(t,e)M)) - \\nu]}{2 \\lambda_{\\rm max}(P) \\eta(t)} \\right)
*/
template <typename MatrixType>
inline double computeLambdaSatFunnelFromMatrixEigenvalue(
  const MatrixType& Q,
  double Ve_function,
  const MatrixType& Q_M,
  double xi_bar,
  double e_norm,
  const MatrixType& P,
  const MatrixType& M,
  double H_function,
  double nu,
  double lambda_max_P,
  double eta_funnel
) {
  static_assert(
    static_cast<int>(MatrixType::RowsAtCompileTime) == static_cast<int>(MatrixType::ColsAtCompileTime),
    "All matrices must be square and of the same size"
  );

  // matrix_saturation = Q + Ve * Q_M - (2 * xi_bar / ||e||) * (P + Ve * M)
  MatrixType matrix_saturation = Q
                               + Ve_function * Q_M
                               - ((2.0 * xi_bar) / e_norm) * (P + Ve_function * M);
                          
  if (!matrix_saturation.isApprox(matrix_saturation.transpose(), 1e-8)) {
    std::cerr << "[WARNING] matrix_saturation is not symmetric!" << std::endl;
  }

  // Compute the smallest eigenvalue of matrix_saturation
  Eigen::SelfAdjointEigenSolver<MatrixType> eigen_solver(matrix_saturation);
  if (eigen_solver.info() != Eigen::Success) {
    // throw std::runtime_error("Eigen decomposition failed in computeLambdaSatFunnelFromMatrixEigenvalue");
    std::cerr << "[WARNING] Eigen decomposition failed in computeLambdaSatFunnelFromMatrixEigenvalue.\n";
    std::cerr << " Returning lambda_sat = 0.\n";
    return 0.0;
  }

  const auto& eigenvalues = eigen_solver.eigenvalues();
  double lambda_min_matrix_saturation = eigenvalues.minCoeff();

  // Compute lambda_sat
  double lambda_sat = std::max(
    0.0,
    H_function * (lambda_min_matrix_saturation - nu) / (2.0 * lambda_max_P * eta_funnel)
  );

  return lambda_sat;
}

/*
  Compute the xi value for funnel
*/
template <typename InputType>
inline double computeXiFunnel(
  const InputType& control_input,
  double control_input_max,
  double control_input_min,
  double Delta_control_input_min)
{
  double control_input_magnitude;

  if constexpr (std::is_same<InputType, double>::value) {
    control_input_magnitude = control_input;
  } else {
    static_assert(static_cast<int>(InputType::ColsAtCompileTime) == 1, "control_input must be a column vector");
    control_input_magnitude = control_input.norm();
  }

  double numerator = control_input_max - control_input_magnitude;
  double denominator = std::max(control_input_magnitude - control_input_min, Delta_control_input_min);

  double xi = numerator / denominator;

  return xi;
}

/*
  Compute sigma_ideal and sigma_nom
*/
inline std::pair<double, double> computeSigmasFunnel(double xi, double eta, double lambda_sat)
{
  double sigma_ideal = eta * xi;
  double sigma_nom = std::min(sigma_ideal, lambda_sat);
  return std::make_pair(sigma_ideal, sigma_nom);
}

/*
  Computes the time derivative of the funnel variable η (eta).
*/
template <typename VectorType, typename MatrixType>
inline std::pair<Eigen::Matrix<double, 1, 1>, int> computeEtaDotFunnel(
  const VectorType& e,
  const VectorType& e_dot,
  double e_norm,
  double eta,
  double H_function,
  double sigma_nom,
  const MatrixType& M,
  double e_min,
  double eta_max,
  double delta_1,
  double delta_2,
  double delta_3,
  bool print_flag = false
) {
  static_assert(static_cast<int>(VectorType::ColsAtCompileTime) == 1, "e and e_dot must be a column vector");
  static_assert(
    static_cast<int>(MatrixType::RowsAtCompileTime) == static_cast<int>(MatrixType::ColsAtCompileTime),
    "M must be a square matrix"
  );
  static_assert(
    static_cast<int>(VectorType::RowsAtCompileTime) == static_cast<int>(MatrixType::RowsAtCompileTime),
    "Dimension mismatch between e/e_dot and M"
  );

  const double eta_lower_bound = std::sqrt(delta_1);
  const double eta_upper_bound = std::sqrt(eta_max - delta_3);
  double eta_dot = 0.0;
  int case_eta_dot = 0;

  if (e_norm > e_min && H_function <= delta_2 && eta > eta_lower_bound && eta < eta_upper_bound)
  {
    double arg = - (e.transpose() * M * e_dot).value() / eta;
    eta_dot = std::min(arg, sigma_nom);
    case_eta_dot = 1;
    if (print_flag) {std::cout << "CASE 1" << std::endl;}
  }
  else if (e_norm > e_min && H_function > delta_2 && eta > eta_lower_bound && eta < eta_upper_bound)
  {
    eta_dot = sigma_nom;
    case_eta_dot = 2;
    if (print_flag) {std::cout << "CASE 2" << std::endl;}
  }
  else if (e_norm > e_min && eta <= eta_lower_bound)
  {
    eta_dot = std::max(0.0, sigma_nom);
    case_eta_dot = 3;
    if (print_flag) {std::cout << "CASE 3" << std::endl;}
  }
  else if (e_norm > e_min && eta >= eta_upper_bound)
  {   
    eta_dot = std::min(0.0, sigma_nom);
    case_eta_dot = 4;
    if (print_flag) {std::cout << "CASE 4" << std::endl;}
  }
  else
  { // e_norm <= e_min  
    eta_dot = 0.0;
    case_eta_dot = 5;
    if (print_flag) {std::cout << "CASE 5" << std::endl;}
  }
  return std::make_pair(Eigen::Matrix<double, 1, 1>(eta_dot), case_eta_dot);
}

} // namespace funnel


namespace hybrid {

/*
  Store the previous iteration trajectory tracking error.
*/
template <typename VectorType>
inline VectorType computePreviousTrajectoryTrackingError(const VectorType& e)
{
  static_assert(
    static_cast<int>(VectorType::ColsAtCompileTime) == 1,
    "e must be a column vector"
  );
  VectorType e_previous = e;

  return e_previous;
}

/*
  Evaluate the s-th element of the series used for reference resetting events.
*/
inline double evaluateSeriesElement(double s, double alpha)
{
  return 1.0 / std::pow(s, alpha);
}

/*
  Compute the next index s in the convergent series.
  weighted_e_squared = eᵀ(t_j) P e(t_j)
*/
inline int findSseries(int s_previous, double weighted_e_squared, double alpha)
{
  // Compute ceil(weighted_e_squared^(-1/alpha))
  double exponent = -1.0 / alpha;
  int candidate = static_cast<int>(
    std::ceil(std::pow(weighted_e_squared, exponent))
  );
  int s = std::max(candidate, s_previous + 1);
  return s;
}

/*
  Compute the quadratic form: vectorᵀ * matrix * vector
*/
template <typename VectorType, typename MatrixType>
inline double computeQuadraticForm(const VectorType& vector, const MatrixType& matrix)
{
  static_assert(static_cast<int>(VectorType::ColsAtCompileTime) == 1, "vector must be a column vector");

  static_assert(static_cast<int>(MatrixType::RowsAtCompileTime) == static_cast<int>(MatrixType::ColsAtCompileTime),
    "matrix must be square"
  );

  static_assert(static_cast<int>(VectorType::RowsAtCompileTime) == static_cast<int>(MatrixType::RowsAtCompileTime),
    "Dimension mismatch between vector and matrix"
  );

  double result = (vector.transpose() * matrix * vector).value();

  return result;
}

/*
  Update the hybrid summation term.
*/
inline double updateSummationHybridP(double ePe, double ePe_previous, double summation)
{
  if (ePe > ePe_previous)
  {
    summation += (ePe - ePe_previous);
  }
  return summation;
}

/*
  Reset the hybrid series index when either:
  - it grows beyond a prescribed threshold, or
  - sufficient time has elapsed since the last trajectory reset.
*/
inline int resetSeriesIfNeeded(
  int s_hybrid,
  double time_now,
  double time_of_last_trajectory_reset,
  double tolerance_time_reset_series
)
{
  constexpr int S_MAX_THRESHOLD = 1000000000;
  const double time_since_last_trajectory_reset = time_now - time_of_last_trajectory_reset;

  if (s_hybrid > S_MAX_THRESHOLD || (s_hybrid > 0 && time_since_last_trajectory_reset > tolerance_time_reset_series))
  {
    s_hybrid = 0;
  }
  return s_hybrid;
}

/*
  Check the hybrid reset condition and, if satisfied, update:
  - reset time,
  - series index s,
  - reference trajectory,
  - integral term storage.
*/
template <typename VectorType, typename MatrixType, typename DerivedRef>
inline std::pair<int, double> checkAndTriggerTrajectoryReset(
  double summation_hybrid_P,
  double time_now,
  double time_of_last_trajectory_reset,
  double ePe,
  double alpha_hybrid_series,
  int s_hybrid,
  const VectorType& e,
  const VectorType& e_previous,
  const MatrixType& Q,
  Eigen::MatrixBase<DerivedRef>& reference_model_state_map,
  Eigen::Ref<Eigen::Matrix<double,1,1>> integral_eQe_map
)
{
  static_assert(VectorType::ColsAtCompileTime == 1, "Vectors must be column vectors");
  static_assert(MatrixType::RowsAtCompileTime == MatrixType::ColsAtCompileTime, "Q must be square");

  double integral_eQe = integral_eQe_map(0);

  // std::cout << "integral_eQe: " << integral_eQe << std::endl;
  // std::cout << "summation_hybrid_P: " << summation_hybrid_P << std::endl;

  if (integral_eQe >= summation_hybrid_P)
  {
    // Update reset time
    time_of_last_trajectory_reset = time_now;

    // Update series index
    s_hybrid = findSseries(s_hybrid, ePe, alpha_hybrid_series);

    // Update reference trajectory (jump)
    const double series_element = evaluateSeriesElement(static_cast<double>(s_hybrid), alpha_hybrid_series);
    const double delta_series = ePe - series_element;
    const double jump_factor = 1.0 - std::sqrt(delta_series / ePe);
    const VectorType jump_reference_trajectory = jump_factor * e;
    reference_model_state_map +=  jump_reference_trajectory;

    // Update integral storage term related to (e^T * Q * e)
    const double eQe = computeQuadraticForm(e, Q);
    const double eQe_previous = computeQuadraticForm(e_previous, Q);
    const double delta_integral_eQe = eQe - eQe_previous;
    integral_eQe_map(0) += delta_integral_eQe;

    std::cout << "checkAndTriggerTrajectoryReset() cond. verified [integral_eQe >= summation_hybrid_P]" << std::endl;
  }

  return {s_hybrid, time_of_last_trajectory_reset};
}

/*
  Execute one hybrid MRAC step:
  - Update quadratic forms and summation term
  - Handle series reset logic
  - Trigger trajectory reset if condition is met
*/
template <typename VectorType, typename MatrixType, typename DerivedRef>
inline std::tuple<Eigen::Matrix<double,1,1>, double, int, double> runHybridStep(
  const VectorType& e,
  const VectorType& e_previous,
  const MatrixType& P,
  const MatrixType& Q,
  double summation_hybrid_P,
  int s_hybrid,
  double time_of_last_trajectory_reset,
  double alpha_hybrid_series,
  double time_now,
  double tolerance_time_reset_series_hybrid,
  Eigen::MatrixBase<DerivedRef>& reference_model_state_map,
  Eigen::Ref<Eigen::Matrix<double,1,1>> integral_eQe_map,
  bool& first_controller_loop
)
{
  static_assert(VectorType::ColsAtCompileTime == 1, "Vectors must be column vectors");
  static_assert(MatrixType::RowsAtCompileTime == MatrixType::ColsAtCompileTime, "Matrices must be square");

  // Compute quadratic terms
  Eigen::Matrix<double,1,1> eQe;
  eQe(0,0) = computeQuadraticForm(e, Q);
  const double ePe = computeQuadraticForm(e, P);
  const double ePe_previous = computeQuadraticForm(e_previous, P);
  // Update hybrid summation term
  summation_hybrid_P = updateSummationHybridP(ePe, ePe_previous, summation_hybrid_P);

  if (!first_controller_loop)
  {
    // Reset series index if required
    s_hybrid = resetSeriesIfNeeded(
      s_hybrid,
      time_now,
      time_of_last_trajectory_reset,
      tolerance_time_reset_series_hybrid
    );

    // Check reset condition
    std::tie(s_hybrid, time_of_last_trajectory_reset) = checkAndTriggerTrajectoryReset(
      summation_hybrid_P,
      time_now,
      time_of_last_trajectory_reset,
      ePe,
      alpha_hybrid_series,
      s_hybrid,
      e,
      e_previous,
      Q,
      reference_model_state_map,
      integral_eQe_map               
    );
  }

  first_controller_loop = false;

  return {eQe, summation_hybrid_P, s_hybrid, time_of_last_trajectory_reset};
}


} // namespace hybrid


namespace ebci {

/*
  Compute the error bounding control input (non-adaptive EBCI).
*/
template <typename MatrixBType, typename MatrixPType, typename VectorType>
inline Eigen::Matrix<typename MatrixBType::Scalar, MatrixBType::ColsAtCompileTime, 1> computeErrorBoundingControlInput(
  double xi_bar_d,
  double lambda_bar,
  double delta_ebci,
  const MatrixBType& B,
  const MatrixPType& P,
  const VectorType& e,
  bool use_ebci
)
{
  using Scalar = typename MatrixBType::Scalar;

  constexpr int n = MatrixBType::RowsAtCompileTime;
  constexpr int m = MatrixBType::ColsAtCompileTime;

  static_assert(VectorType::ColsAtCompileTime == 1, "e must be a column vector");
  static_assert(MatrixPType::RowsAtCompileTime == MatrixPType::ColsAtCompileTime, "P must be square");
  static_assert(MatrixPType::RowsAtCompileTime == n, "Dimension mismatch: P and B must have same row size");
  static_assert(VectorType::RowsAtCompileTime == n, "Dimension mismatch: e and B must have same row size");

  using InputVector = Eigen::Matrix<Scalar, m, 1>;

  if (!use_ebci)
  {
    return InputVector::Zero();
  }

  InputVector BPe = B.transpose() * P * e;
  Scalar BPe_norm = BPe.norm();

  // Dead-zone
  if (BPe_norm < delta_ebci)
  {
    return InputVector::Zero();
  }

  Eigen::Matrix<Scalar, n, 1> Pe = P * e;
  Scalar sum_Pe = Pe.cwiseAbs().sum();
  InputVector control_input = - (xi_bar_d / lambda_bar) * (BPe / BPe_norm) * sum_Pe;

  return control_input;
}

} // namespace ebci


} // namespace base_mrac