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
 * File:        hybrid_two_layer_mrac_gains.hpp
 * Author:      Mattia Gramuglia
 * Date:        February 25, 2026
 * For info:    Andrea L'Afflitto 
 *              a.lafflitto@vt.edu
 * 
 * Description: Tuning gains of the HybridTwoLayerMRAC controller.
 * 
 * GitHub:    https://github.com/andrealaffly/ACSL-flightstack.git
 **********************************************************************************************************************/

#pragma once

#include <cmath>
#include <algorithm>

#include <Eigen/Dense>

#include "json_parser.hpp"


// Struct containing the HybridTwoLayerMRAC tuning gains coming from the .json file 
struct GainsHybridTwoLayerMRAC
{
  /*****************************************************************************************
   *  BASELINE parameters
  ******************************************************************************************/
  // TRANSLATIONAL gains to let the reference model follow the user-defined trajectory
  Eigen::Matrix3d KP_refmod_translational;
  Eigen::Matrix3d KD_refmod_translational;
  Eigen::Matrix3d KI_refmod_translational;

  // TRANSLATIONAL gains for the BASELINE controller to follow the reference model
  Eigen::Matrix3d KP_translational;
  Eigen::Matrix3d KD_translational;
  Eigen::Matrix3d KI_translational;

  // ROTATIONAL gains to generate a command angular velocity to let the current attitude
  // track the desired attitude
  Eigen::Matrix3d KP_omega_cmd_rotational;
  Eigen::Matrix3d KI_omega_cmd_rotational;

  // ROTATIONAL gains to let the reference model track the command angular velocity
  Eigen::Matrix3d KP_omega_ref_rotational;
  Eigen::Matrix3d KI_omega_ref_rotational;

  // ROTATIONAL gains for the BASELINE controller to follow the reference model
  Eigen::Matrix3d KP_rotational;
  Eigen::Matrix3d KD_rotational;
  Eigen::Matrix3d KI_rotational;

  /*****************************************************************************************
   *  ADAPTIVE parameters
  ******************************************************************************************/
  // TRANSLATIONAL gains ADAPTIVE controller
  Eigen::Matrix<double, 6, 6> Gamma_x_translational;
  Eigen::Matrix<double, 3, 3> Gamma_r_translational;
  Eigen::Matrix<double, 6, 6> Gamma_Theta_translational;
  Eigen::Matrix<double, 6, 6> Gamma_g_translational;

  // TRANSLATIONAL parameters LYAPUNOV equation
  Eigen::Matrix<double, 6, 6> Q_translational;
  Eigen::Matrix<double, 6, 6> P_translational;

  // TRANSLATIONAL parameters ROBUST MRAC (e-modification)
  double sigma_x_translational;
  double sigma_r_translational; 
  double sigma_Theta_translational;
  double sigma_g_translational;

  // TRANSLATIONAL parameters ROBUST MRAC (dead-zone modification)
  double dead_zone_delta_translational;
  double dead_zone_e0_translational;

  // TRANSLATIONAL parameters ROBUST MRAC (ellipsoid projection operator)
  Eigen::Matrix<double, 18, 1> x_e_x_translational;
  Eigen::Matrix<double, 18, 1> S_diagonal_x_translational;
  Eigen::Matrix<double, 18, 18> S_x_translational;
  double alpha_x_translational;
  double epsilon_x_translational;
  Eigen::Matrix<double, 9, 1> x_e_r_translational;
  Eigen::Matrix<double, 9, 1> S_diagonal_r_translational;
  Eigen::Matrix<double, 9, 9> S_r_translational;
  double alpha_r_translational;
  double epsilon_r_translational;
  Eigen::Matrix<double, 18, 1> x_e_Theta_translational;
  Eigen::Matrix<double, 18, 1> S_diagonal_Theta_translational;
  Eigen::Matrix<double, 18, 18> S_Theta_translational;
  double alpha_Theta_translational;
  double epsilon_Theta_translational;
  Eigen::Matrix<double, 18, 1> x_e_g_translational;
  Eigen::Matrix<double, 18, 1> S_diagonal_g_translational;
  Eigen::Matrix<double, 18, 18> S_g_translational;
  double alpha_g_translational;
  double epsilon_g_translational;

  // ROTATIONAL gains ADAPTIVE controller
  Eigen::Matrix<double, 3, 3> Gamma_x_rotational;
  Eigen::Matrix<double, 3, 3> Gamma_r_rotational;
  Eigen::Matrix<double, 6, 6> Gamma_Theta_rotational;
  Eigen::Matrix<double, 3, 3> Gamma_g_rotational;

  // ROTATIONAL parameters LYAPUNOV equation
  Eigen::Matrix<double, 3, 3> Q_rotational;
  Eigen::Matrix<double, 3, 3> P_rotational;

  // ROTATIONAL parameters ROBUST MRAC (e-modification)
  double sigma_x_rotational;
  double sigma_r_rotational; 
  double sigma_Theta_rotational;
  double sigma_g_rotational;

  // ROTATIONAL parameters ROBUST MRAC (dead-zone modification)
  double dead_zone_delta_rotational;
  double dead_zone_e0_rotational;

  // ROTATIONAL parameters ROBUST MRAC (ellipsoid projection operator)
  Eigen::Matrix<double, 9, 1> x_e_x_rotational;
  Eigen::Matrix<double, 9, 1> S_diagonal_x_rotational;
  Eigen::Matrix<double, 9, 9> S_x_rotational;
  double alpha_x_rotational;
  double epsilon_x_rotational;
  Eigen::Matrix<double, 9, 1> x_e_r_rotational;
  Eigen::Matrix<double, 9, 1> S_diagonal_r_rotational;
  Eigen::Matrix<double, 9, 9> S_r_rotational;
  double alpha_r_rotational;
  double epsilon_r_rotational;
  Eigen::Matrix<double, 18, 1> x_e_Theta_rotational;
  Eigen::Matrix<double, 18, 1> S_diagonal_Theta_rotational;
  Eigen::Matrix<double, 18, 18> S_Theta_rotational;
  double alpha_Theta_rotational;
  double epsilon_Theta_rotational;
  Eigen::Matrix<double, 9, 1> x_e_g_rotational;
  Eigen::Matrix<double, 9, 1> S_diagonal_g_rotational;
  Eigen::Matrix<double, 9, 9> S_g_rotational;
  double alpha_g_rotational;
  double epsilon_g_rotational;

  /*****************************************************************************************
   *  PLANT & REFERENCE MODEL parameters
  ******************************************************************************************/
  // TRANSLATIONAL dynamics plant parameters
  Eigen::Matrix<double, 6, 6> A_translational;
  Eigen::Matrix<double, 6, 3> B_translational;

  // TRANSLATIONAL reference model parameters
  Eigen::Matrix<double, 6, 6> A_ref_translational;
  Eigen::Matrix<double, 6, 3> B_ref_translational;

  // TRANSLATIONAL two-layer auxillary model parameters
  Eigen::Matrix<double, 6, 6> A_transient_translational;

  // ROTATIONAL dynamics plant parameters
  Eigen::Matrix<double, 3, 3> A_rotational;
  Eigen::Matrix<double, 3, 3> B_rotational;

  // ROTATIONAL reference model parameters
  Eigen::Matrix<double, 3, 3> A_ref_rotational;
  Eigen::Matrix<double, 3, 3> B_ref_rotational;

  // ROTATIONAL two-layer auxillary model parameters
  Eigen::Matrix<double, 3, 3> A_transient_rotational;

  /*****************************************************************************************
   *  HYBRID parameters
  ******************************************************************************************/
  // TRANSLATIONAL gains HYBRID controller
  bool use_hybrid;
  double alpha_hybrid_series_translational;
  double tolerance_time_reset_series_hybrid_translational;

   /*****************************************************************************************
   *  ERROR BOUNDING CONTROL INPUT (EBCI) parameters
  ******************************************************************************************/
  // TRANSLATIONAL gains EBCI controller
  bool use_error_bounding_control_input;
  double xi_bar_d_ebci_translational;
  double lambda_bar_ebci_translational;
  double delta_ebci_translational;

};

namespace nlohmann
{

inline void to_json(nlohmann::json& j, const GainsHybridTwoLayerMRAC& g)
{
  j["A_translational"] = g.A_translational;
  j["B_translational"] = g.B_translational;
  j["A_ref_translational"] = g.A_ref_translational;
  j["B_ref_translational"] = g.B_ref_translational;

  j["A_rotational"] = g.A_rotational;
  j["B_rotational"] = g.B_rotational;
  j["A_ref_rotational"] = g.A_ref_rotational;
  j["B_ref_rotational"] = g.B_ref_rotational;

  j["P_translational"] = g.P_translational;
  j["P_rotational"] = g.P_rotational;

  j["S_x_translational"] = g.S_x_translational;
  j["S_r_translational"] = g.S_r_translational;
  j["S_Theta_translational"] = g.S_Theta_translational;
  j["S_g_translational"] = g.S_g_translational;

  j["S_x_rotational"] = g.S_x_rotational;
  j["S_r_rotational"] = g.S_r_rotational;
  j["S_Theta_rotational"] = g.S_Theta_rotational;
  j["S_g_rotational"] = g.S_g_rotational;

  j["epsilon_x_translational"] = g.epsilon_x_translational;
  j["epsilon_r_translational"] = g.epsilon_r_translational;
  j["epsilon_Theta_translational"] = g.epsilon_Theta_translational;
  j["epsilon_g_translational"] = g.epsilon_g_translational;

  j["epsilon_x_rotational"] = g.epsilon_x_rotational;
  j["epsilon_r_rotational"] = g.epsilon_r_rotational;
  j["epsilon_Theta_rotational"] = g.epsilon_Theta_rotational;
  j["epsilon_g_rotational"] = g.epsilon_g_rotational;
}

} // namspace nlohmann