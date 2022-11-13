// Use of this source code is governed by an MIT-style license that can be found
// in the LICENSE file or at https://opensource.org/licenses/MIT. 
//
// Copyright (c) 2022 Maxence Reberol

#include "solver.h"
#include "disk_triangulations.h"

#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>

using namespace hbl;

std::string directory_name(const std::string& fname) {
  size_t pos = fname.find_last_of("\\/");
  return (std::string::npos == pos) ? "" : fname.substr(0, pos);
}

template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& values) { 
  os << "[";
  for (size_t i = 0; i < values.size(); ++i) {
    const  T & x = values[i];
    os << x;
    if (i != values.size() - 1) {
      os << ", ";
    }
  }
  os << "]";
  return os;
}

template<class T> 
void sort_unique(std::vector<T>& vec) {
    std::sort( vec.begin(), vec.end() );
    vec.erase( std::unique( vec.begin(), vec.end() ), vec.end() );
}

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


bool verify_all_disk_triangulation_boundaries_are_found_by_solver(
    hbl::DiskTriangulations& dts,
    size_t nVerticesOnBoundaryMax) {

  size_t nFound = 0;
  size_t nTry = 0;
  for (size_t b = 3; b <= nVerticesOnBoundaryMax; ++b) {
    for (auto& kv: dts.nbv_sign_to_trgl[b]) {
      bdrsign_t bdrValences =  kv.first;

      bool okVmax = true;
      std::vector<int> n;
      for (size_t k = 0; k < b; ++k) {
        n.push_back((int)bdrValences[k]);
        if ((int)bdrValences[k] > VAL_MAX) {
          okVmax = false;
        }
      }
      if (!okVmax) continue;

      nTry += 1;
      bool found = boundary_disk_triangulation_exists(n);
      if (found) {
        nFound += 1;
      }
    }
  }
  if (nTry == nFound) {
    printf("\n");
    printf("Verify all disk triangulation boundaries are found by integer solver ...\n");
    printf("up to %li vertices on polygon boundary\n",nVerticesOnBoundaryMax);
    printf("- found %li/%li, OK\n", nFound, nTry);
    printf("\n");
  } else {
    printf("\n");
    printf("Verify all disk triangulation boundaries are found by integer solver ...\n");
    printf("up to %li vertices on polygon boundary\n",nVerticesOnBoundaryMax);
    printf("- found %li/%li, NOT OK\n", nFound, nTry);
    printf("\n");
  }

  return nTry == nFound;
}

bool compare_disk_triangulation_boundaries_and_integer_sets_generated_by_solver(
    hbl::DiskTriangulations& dts,
    size_t nVerticesOnBoundaryMax) {

  bool okProof = true;
  for (size_t b = 3; b <= nVerticesOnBoundaryMax; ++b) {

    /* Collect disk triangulation boundaries from dataset */
    std::vector<std::vector<int> > diskTriangulationBoundaries;
    for (auto& kv: dts.nbv_sign_to_trgl[b]) {
      bdrsign_t bdrValences =  kv.first;
      bool okVmax = true;
      std::vector<int> n;
      for (size_t k = 0; k < b; ++k) {
        n.push_back((int)bdrValences[k]);
        if ((int)bdrValences[k] > VAL_MAX) {
          okVmax = false;
        }
      }
      if (!okVmax) continue;
      diskTriangulationBoundaries.push_back(n);
    }
    sort_unique(diskTriangulationBoundaries);

    /* Collect integer sets of boundary valences generated
     * by the integer solver */
    std::vector<std::vector<int> > polygons(1);
    std::vector<double> x_i(b,0.);
    std::vector<std::array<int,2> > n_min_max(b);
    std::vector<int> n_priority(b,1);
    for (int i = 0; (size_t) i < b; ++i) {
      polygons[0].push_back(i);
      x_i[i] = double(0.);
      n_min_max[i] = {VAL_MIN,VAL_MAX};
    };
    std::vector<std::vector<int> > solutions;
    bool found = solveAllHexLayerTopologyAllSolutions(b, polygons, x_i, n_min_max, n_priority, solutions);

    /* - get the unique solutions by removing the rotations/symmetries */
    std::vector<std::vector<int> > uniqueSolutions;
    for (size_t k = 0; k < solutions.size(); ++k) {
      std::vector<int> canonical = get_smallest_rotation(solutions[k]);
      uniqueSolutions.push_back(canonical);
    }
    sort_unique(uniqueSolutions);

    /* Comparison */
    printf("\n");
    printf("Compare disk triangulation boundary and integer sets generted by solver ...\n");
    printf("- polygon size: %li\n",b);
    printf("- from dataset: %li boundaries\n",diskTriangulationBoundaries.size());
    printf("- from solver: %li unique solutions\n",uniqueSolutions.size());
    printf("- verify each solution match ... \n");
    for (size_t k = 0; k < uniqueSolutions.size(); ++k) {
      std::vector<int> slt = uniqueSolutions[k];
      auto it = std::find(diskTriangulationBoundaries.begin(),diskTriangulationBoundaries.end(),slt);
      if (it == diskTriangulationBoundaries.end()) {
        std::cout << "-- slt from solver not found in dataset: " << slt << std::endl;
      }
    }
    printf("\n");

  }

  return okProof;
}

int main (int argc, char *argv[]) { 
  printf("Experimental proof ... \n");

  std::string dt_path = directory_name(std::string(argv[0])) + "/assets/disk_trgls_3.data";
  hbl::DiskTriangulations dts;
  bool okl = hbl::load_disk_triangulations(dt_path, dts);
  if (!okl) {
    printf("error: failed to load disk triangulations from path: %s\n",dt_path.c_str());
    return 1;
  }

  size_t nVerticesOnBoundaryMax = 6;

  verify_all_disk_triangulation_boundaries_are_found_by_solver(dts, nVerticesOnBoundaryMax);

  compare_disk_triangulation_boundaries_and_integer_sets_generated_by_solver(dts, nVerticesOnBoundaryMax);

  printf("end of program\n");
  return 0;
}

