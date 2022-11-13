// Use of this source code is governed by an MIT-style license that can be found
// in the LICENSE file or at https://opensource.org/licenses/MIT. 
//
// Copyright (c) 2022 Maxence Reberol

#include "disk_triangulations.h"

#include <map>
#include <unordered_map>
#include <algorithm>

#include <sstream>
#include <iostream>
#include <fstream>

using std::vector;
using std::map;
using std::unordered_map;
using std::array;
using std::pair;

namespace hbl {
  constexpr int NO_INT = -1;

  bool trgl_is_file(const char *fileName) {
    std::ifstream in(fileName);
    return in.good();
  }

  std::vector<std::string> split_string(const std::string& str, char delim) {
    std::vector<std::string> strings;
    size_t start;
    size_t end = 0;
    while ((start = str.find_first_not_of(delim, end)) != std::string::npos) {
      end = str.find(delim, start);
      strings.push_back(str.substr(start, end - start));
    }
    return strings;
  }

  inline int2 sorted(int v1, int v2) { if (v1 < v2) { return {v1,v2}; } else { return {v2,v1}; } }
  inline int2 sorted(int2 e) { if (e[0] < e[1]) { return {e[0],e[1]}; } else { return {e[1],e[0]}; } }
  struct int2Hash {
    size_t operator()(int2 p) const noexcept {
      return size_t(p[0]) << 32 | p[1];
    }
  };

  vector<int> get_smallest_rotation(const vector<int>& vec) {
    vector<int> smallest = vec;
    vector<int> rot_l = vec;
    vector<int> rot_r = vec;
    std::reverse(rot_r.begin(),rot_r.end());
    for (size_t i = 0; i < vec.size(); ++i) {
      std::rotate(rot_l.begin(),rot_l.begin()+1,rot_l.end());
      if (rot_l < smallest) smallest = rot_l;
      std::rotate(rot_r.begin(),rot_r.begin()+1,rot_r.end());
      if (rot_r < smallest) smallest = rot_r;
    }
    return smallest;
  }

  bool get_smallest_rotation(const vector<int>& bdrVal, const vector<int>& bdrLoop, 
      vector<int>& smallestBdrVal, vector<int>& smallestBdrLoop) {
    smallestBdrVal = bdrVal;
    smallestBdrLoop = bdrLoop;
    vector<int> rot = bdrVal;
    vector<int> ids = bdrLoop;
    for (size_t i = 0; i < bdrVal.size(); ++i) {
      std::rotate(rot.begin(),rot.begin()+1,rot.end());
      std::rotate(ids.begin(),ids.begin()+1,ids.end());
      if (rot < smallestBdrVal) {
        smallestBdrVal = rot;
        smallestBdrLoop = ids;
      }
    }
    std::reverse(rot.begin(),rot.end());
    std::reverse(ids.begin(),ids.end());
    for (size_t i = 0; i < bdrVal.size(); ++i) {
      std::rotate(rot.begin(),rot.begin()+1,rot.end());
      std::rotate(ids.begin(),ids.begin()+1,ids.end());
      if (rot < smallestBdrVal) {
        smallestBdrVal = rot;
        smallestBdrLoop = ids;
      }
    }
    return true;
  }

  bool getOrderedVerticesFromEdges(int vStart, const vector<int2>& edges, vector<int>& orderedVertices) {
    orderedVertices.clear();
    int eStart = NO_INT;
    int lvStart = NO_INT;
    for (int e = 0; e < edges.size(); ++e) for (int lv = 0; lv < 2; ++lv) {
      if (edges[e][lv] == vStart) {
        eStart = e;
        lvStart = lv;
        break;
      }
    }
    if (eStart == NO_INT) return false;

    size_t iter = 0;
    int e = eStart;
    int v = vStart;
    while (true) {
      iter += 1;
      if (iter > 100000) {
        printf("error: infinite loop ? iter = %li\n", iter);
        return false;
      }
      orderedVertices.push_back(v);
      int next_v = (edges[e][0]  != v) ? edges[e][0] : edges[e][1];
      if (next_v == vStart) { /* closed chain */
        break;
      }
      int next_e = NO_INT;
      for (int ee = 0; ee < edges.size(); ++ee) if (ee != e) for (int lv = 0; lv < 2; ++lv) {
        if (edges[ee][lv] == next_v) {
          next_e = ee;
          break;
        }
      }
      if (next_e == NO_INT) { /* open chain */
        orderedVertices.push_back(next_v);
        break;
      }
      v = next_v;
      e = next_e;
    }
    return true;
  }

  bool DTriangulation::fill_struct_from_triangles() {
    vOnBdr.clear();
    vOnBdr.resize(nv,false);

    vValence.clear();
    vValence.resize(nv,0);
    std::unordered_map<int2,int,int2Hash> sedge_val;
    for (size_t f = 0; f < triangles.size(); ++f) {
      for (size_t le = 0; le < 3; ++le) {
        int v1 = triangles[f][le];
        int v2 = triangles[f][(le+1)%3];
        int2 sedge = sorted(v1,v2);
        sedge_val[sedge] += 1;
        vValence[v1] += 1;
      }
    }

    edges.clear();
    edges.reserve(sedge_val.size());

    eOnBdr.clear();
    eOnBdr.reserve(sedge_val.size());

    for (const auto& kv : sedge_val) {
      edges.push_back(kv.first);
      if (kv.second == 1) {
        eOnBdr.push_back(true);
        vOnBdr[kv.first[0]] = true;
        vOnBdr[kv.first[1]] = true;
      } else if (kv.second == 2) {
        eOnBdr.push_back(false);
      } else {
        printf("error: load_triangulation, non manifold edge\n");
        return false;
      }
    }
    nb = 0;
    for (size_t v = 0; v < nv; ++v) if (vOnBdr[v]) {
      nb += 1;
    }
    return true;
  }

  bool DTriangulation::load_triangulation(const std::vector<int>& tri_vertices) {
    if(tri_vertices.size() == 0 || tri_vertices.size() % 3 != 0) {
      printf("error: load_triangulation | wrong input\n");
      return false;
    } 
    triangles.resize(tri_vertices.size()/3);
    nv = 0;
    for (size_t i = 0; i < triangles.size(); ++i) {
      for (size_t j = 0; j < 3; ++j) {
        int v = tri_vertices[3*i+j];
        triangles[i][j] = v;
        nv = std::max(nv,v+1);
      }
    }
    bool okf = fill_struct_from_triangles();
    if (!okf) {
      printf("error: failed to build the DTriangulation\n");
      return false;
    }
    return true;
  }

  bool DTriangulation::transform_to_canonical_boundary_loop() {
    vector<int> old2new(nv,NO_INT);

    /* Get vertex boundary loop */
    vector<int2> bdrEdges;
    for (size_t le = 0; le < edges.size(); ++le) if (eOnBdr[le]) {
      bdrEdges.push_back(edges[le]);
    }
    std::vector<int> overt;
    bool ok = getOrderedVerticesFromEdges(bdrEdges[0][0], bdrEdges, overt);
    if (!ok) {
      printf("error: failed to get ordered vertices from edges\n");
      return false;
    }
    vector<int> bdr_val_loop(overt.size());
    for (size_t i = 0; i < bdr_val_loop.size(); ++i) {
      bdr_val_loop[i] = vValence[overt[i]];
    }

    /* Smallest 'boundary val' loop */
    vector<int> s_bdrVal(overt.size());
    vector<int> s_overt(overt.size());
    get_smallest_rotation(bdr_val_loop, overt, s_bdrVal, s_overt);

    /* Re-numbering map */
    int cv = 0;
    for (size_t i = 0; i < s_overt.size(); ++i) {
      old2new[s_overt[i]] = cv;
      cv += 1;
    }
    for (size_t i = 0; i < nv; ++i) if (old2new[i] == NO_INT) {
      old2new[i] = cv;
      cv += 1;
    }

    /* Apply */
    for (size_t i = 0; i < triangles.size(); ++i) {
      for (size_t j = 0; j < 3; ++j) {
        triangles[i][j] = old2new[triangles[i][j]];
      }
    }

    bool okf = fill_struct_from_triangles();
    if (!okf) {
      printf("error: failed to build the DTriangulation after transform_to_canonical_boundary_loop\n");
    }

    return true;
  }

  bool load_disk_triangulations(const std::string& path, DiskTriangulations& data) {
    if(!trgl_is_file(path.c_str())) {
      printf("error: file not found: %s\n", path.c_str());
      return false;
    }

    std::ifstream input(path);
    std::string nb_trgl_str;
    std::getline(input, nb_trgl_str, '\n');
    size_t nbTriangulations = std::stoi(nb_trgl_str);
    printf("loading %li disk triangulations from disk ...\n", nbTriangulations);
    int nbmax = 0;
    for (size_t i = 0; i < nbTriangulations; ++i) {
      std::string str_line;
      std::getline(input, str_line, '\n');
      std::vector<std::string> str_vertices = split_string(str_line, ' ');

      if(str_vertices.size() % 3 != 0) {
        printf("error: load_disk_triangulations | trgl %li: vertices are not a multiple of 3\n", i);
        return false;
      } 
      std::vector<int> tri_vertices(str_vertices.size());
      for (size_t j = 0; j < tri_vertices.size(); ++j) {
        tri_vertices[j] = std::stoi(str_vertices[j]);
      }
      std::unique_ptr<DTriangulation> uptr(new DTriangulation);
      bool okl = uptr->load_triangulation(tri_vertices);
      if (!okl) {
        printf("failed to load triangulation with vertices\n");
        return false;
      }
      bool okc = uptr->transform_to_canonical_boundary_loop();
      if (!okc) {
        printf("failed to transform to canonical representation\n");
        return false;
      }
      data.trgls.push_back(std::move(uptr));
      nbmax = std::max(nbmax,data.trgls.back()->nb);
    }

    /* Create the hash mapping */
    printf("computing signatures and redirections ...\n");
    data.nbv_sign_to_trgl.clear();
    data.nbv_sign_to_trgl.resize(nbmax+1);

    for (size_t i = 0; i < data.trgls.size();++i) {
      const DTriangulation& dt = *(data.trgls[i]);
      if (((size_t) dt.nb >= data.nbv_sign_to_trgl.size())) {
        printf("error nbv\n");
        return false;
      }
      vector<int2> bdrEdges;
      for (size_t le = 0; le < dt.edges.size(); ++le) if (dt.eOnBdr[le]) {
        bdrEdges.push_back(dt.edges[le]);
      }
      vector<int> bdr_val_loop;
      for (size_t j = 0; j < dt.nv; ++j) if (dt.vOnBdr[j]) {
        bdr_val_loop.push_back(dt.vValence[j]);
      }
      bdrsign_t bsign = signature_from_boundary(bdr_val_loop);
      data.nbv_sign_to_trgl[dt.nb][bsign].push_back(data.trgls[i].get());
    }

    return true;
  }

  bool compute_triangulation_boundary_valence_loop(
      const std::vector<int2>& bdrEdges,
      const std::vector<int>& vValence,
      std::vector<int>& bdr_val_loop) {
    std::vector<int> overt;
    bool ok = getOrderedVerticesFromEdges(bdrEdges[0][0], bdrEdges, overt);
    if (!ok) {
      printf("failed to get ordered vertices from edges\n");
      return false;
    }
    bdr_val_loop.clear();
    bdr_val_loop.resize(overt.size());
    for (size_t i = 0; i < bdr_val_loop.size(); ++i) {
      bdr_val_loop[i] = vValence[overt[i]];
    }
    return true;
  }

  bdrsign_t signature_from_boundary(const std::vector<int>& bdr_val_loop) {
    bdrsign_t sign;
    if (bdr_val_loop.size() >= sign.size()) {
      printf("error: too much vertices in boundary: %li", bdr_val_loop.size());
    }
    std::fill(sign.begin(),sign.end(),0);
    for (size_t i = 0; i < bdr_val_loop.size(); ++i) {
      sign[i] = (unsigned char) bdr_val_loop[i];
    }
    return sign;
  }
}
