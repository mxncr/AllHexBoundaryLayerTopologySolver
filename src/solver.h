// Use of this source code is governed by an MIT-style license that can be found
// in the LICENSE file or at https://opensource.org/licenses/MIT. 
//
// Copyright (c) 2022 Maxence Reberol

#pragma once

#include <vector>
#include <array>
#include <cstdint>

namespace hbl {

  /* Hardcoded maximal range for the integer values, 
   * can be changed at compile time */
  constexpr int VAL_MIN = 1;
  constexpr int VAL_MAX = 4;

  /**
   * @brief Search for an optimal integer solution to the all-hexahedral
   *        boundary layer problem. In practice, the problem should have hundred
   *        unkwowns at most for quick results. This function does not perform problem
   *        decomposition, which should be done beforehand before calling this one.
   *        See the paper "All-hexahedral boundary layer meshing" for more context.
   *        Internally, this function uses the Gecode library (https://www.gecode.org).
   *
   * @param[in] N number of integer variables
   * @param[in] polygons the boundaries of disk triangulations (vertex indexes)
   * @param[in] x_i ideal values for each variable
   * @param[in] n_min_max allowed range for each variable
   * @param[in] n_priority variable priority in the branch and bound (highest value selected first)
   * @param[out] n the best solution found
   * @param[out] stopped true if search stopped due timeout limit
   * @param[out] iterTime one pair (iter,time) for each solution encountered during the search
   * @param[in] timeMaxInit maximum time budget to found the initial solution (in miliseconds)
   * @param[in] timeMaxImprove maximum time budget to improve after initial solution found (in miliseconds)
   *
   * @return true if integer solution n found
   */
  bool solveAllHexLayerTopology(
      size_t N,
      const std::vector<std::vector<int> >& polygons,
      const std::vector<double>& x_i,
      const std::vector<std::array<int,2> >& n_min_max,
      const std::vector<int>& n_priority,
      std::vector<int>& n,
      bool& stopped,
      std::vector<double>& iterTime,
      double timeMaxInit = 100e3,
      double timeMaxImprove = 10e3);

  /**
   * @brief Search all solutions in the given ranges.  
   *        Useful to validate the solver, not for practical meshing.
   *
   * @param[in] N number of integer variables
   * @param[in] polygons the boundaries of disk triangulations (vertex indexes)
   * @param[in] x_i ideal values for each variable
   * @param[in] n_min_max allowed range for each variable
   * @param[in] n_priority variable priority in the branch and bound (highest value selected first)
   * @param[out] solutions list of solutions, each solution is a vector of integer values
   *
   * @return true if at least one solution found
   */
  bool solveAllHexLayerTopologyAllSolutions(
      size_t N,
      const std::vector<std::vector<int> >& polygons,
      const std::vector<double>& x_i,
      const std::vector<std::array<int,2> >& n_min_max,
      const std::vector<int>& n_priority,
      std::vector<std::vector<int> >& solutions);
}
