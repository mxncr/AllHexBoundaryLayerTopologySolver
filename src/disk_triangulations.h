// Use of this source code is governed by an MIT-style license that can be found
// in the LICENSE file or at https://opensource.org/licenses/MIT. 
//
// Copyright (c) 2022 Maxence Reberol

#pragma once

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include <iostream>
#include <fstream>

namespace hbl {

  using bdrsign_t = std::array<unsigned char,12>;

  struct bdrsign_t_hash {
    /* from https://en.wikipedia.org/wiki/Jenkins_hash_function */
    size_t operator()(bdrsign_t p) const noexcept {
      uint32_t hash = 0;
      for (size_t i = 0; i < p.size(); ++i) {
        hash += p[i];
        hash += hash << 10;
        hash ^= hash >> 6;
      }
      hash += hash << 3;
      hash ^= hash >> 11;
      hash += hash << 15;
      return hash;
    }
  };

  using int2 = std::array<int,2>;
  using int3 = std::array<int,3>;

  struct DTriangulation {
    int nb = 0;
    int nv = 0;
    std::vector<int3> triangles; /* vertex ids */
    std::vector<bool> vOnBdr;
    std::vector<int> vValence;
    std::vector<int2> edges;
    std::vector<bool> eOnBdr;

    bool fill_struct_from_triangles();
    bool load_triangulation(const std::vector<int>& tri_vertices);
    bool transform_to_canonical_boundary_loop();
  };

  bool compute_triangulation_boundary_valence_loop(
      const std::vector<int2>& bdrEdges,
      const std::vector<int>& vValence,
      std::vector<int>& bdr_val_loop);

  std::vector<int> get_smallest_rotation(const std::vector<int>& vec);
  bool get_smallest_rotation(const std::vector<int>& bdrVal, const std::vector<int>& bdrLoop, 
      std::vector<int>& smallestBdrVal, std::vector<int>& smallestBdrLoop);

  bdrsign_t signature_from_boundary(const std::vector<int>& bdr_val_loop);

  struct DiskTriangulations {
    std::vector<std::unique_ptr<DTriangulation> > trgls;
    std::vector<std::unordered_map<bdrsign_t,std::vector<DTriangulation*>,bdrsign_t_hash> > nbv_sign_to_trgl;
    bool write(const std::string& path);
    bool read(const std::string& path);
  };

  bool load_disk_triangulations(const std::string& path, DiskTriangulations& trgls);

}
