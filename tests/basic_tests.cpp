// Use of this source code is governed by an MIT-style license that can be found
// in the LICENSE file or at https://opensource.org/licenses/MIT. 
//
// Copyright (c) 2022 Maxence Reberol

#define CATCH_CONFIG_MAIN 
#include "catch.hpp"

#include "solver.h"

using namespace hbl;

bool boundary_disk_triangulation_exists(const std::vector<int>& n) {
  size_t N = n.size();
  std::vector<std::vector<int> > polygons(1);
  std::vector<double> x_i(n.size(),0.);
  std::vector<std::array<int,2> > n_min_max(n.size());
  std::vector<int> n_priority(n.size(),1);
  for (int i = 0; (size_t) i < N; ++i) {
    polygons[0].push_back(i);
    x_i[i] = double(n[i]);
    n_min_max[i] = {VAL_MIN,VAL_MAX};
  };

  bool stopped;
  std::vector<double> iterTime;
  double timeMaxInit = 100e3;
  double timeMaxImprove = 10e3;

  std::vector<int> slt(N,0);
  bool found = solveAllHexLayerTopology(N, polygons, x_i, n_min_max, n_priority, slt,
      stopped, iterTime, timeMaxInit, timeMaxImprove);
  if (found && slt == n) {
    return true;
  }
  return false;
}

TEST_CASE("Check small boundary disk triangulations exist", "[basic]") {
  REQUIRE(boundary_disk_triangulation_exists({1,1,1}));
  REQUIRE(boundary_disk_triangulation_exists({2,2,2}));
  REQUIRE(boundary_disk_triangulation_exists({3,3,3}));
  REQUIRE(boundary_disk_triangulation_exists({1,2,1,2}));
  REQUIRE(boundary_disk_triangulation_exists({2,3,3,2}));
}

TEST_CASE("Check small boundary triangulations do not exist", "[basic]") {
  REQUIRE(!boundary_disk_triangulation_exists({1,2,2}));
  REQUIRE(!boundary_disk_triangulation_exists({1,1,1,1}));
  REQUIRE(!boundary_disk_triangulation_exists({2,2,2,1}));
  REQUIRE(!boundary_disk_triangulation_exists({2,2,2,3}));
  REQUIRE(!boundary_disk_triangulation_exists({2,2,2,1,1,2}));
}


