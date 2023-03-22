// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PARABOLOID_RECONSTRUCTION_PARAMETRIZED_SURFACE_TPP_
#define IRL_PARABOLOID_RECONSTRUCTION_PARAMETRIZED_SURFACE_TPP_

#include <fstream>
#include <iomanip>

#define IRL_USE_EARCUT
// #define IRL_USE_TRIANGLE
// #define IRL_USE_CGAL

#include "external/NumericalIntegration/NumericalIntegration.h"

namespace IRL {

template <class List>
bool noDuplicates(List& list) {
  if (list.size() == 0) {
    return true;
  }
  for (UnsignedIndex_t i = 0; i < list.size() - 1; i++) {
    for (UnsignedIndex_t j = i + 1; j < list.size(); j++) {
      if (list[i] == list[j]) {
        return false;
      }
    }
  }
  return true;
}

template <class VertexList, class TriList, class EdgeList,
          class TriConnectivity, class TriEdges, class VertConnectivity,
          class VertexValence>
void splitEdge(const UnsignedIndex_t edge_id, VertexList& vertices,
               TriList& triangles, EdgeList& edges, TriConnectivity& tri_neigh,
               TriEdges& tri_edges, VertConnectivity& vert_tri,
               VertexValence& vert_neigh) {
  // Start and end point
  const int v0 = edges[edge_id][0];
  const int v1 = edges[edge_id][1];
  assert(v0 != v1);
  assert(v0 >= 0 && v0 < vertices.size());
  assert(v1 >= 0 && v1 < vertices.size());

  // Find opposite triagnles and vertices
  const int tp = edges[edge_id][2];
  assert(tp >= 0);
  int vp = -1;
  int dp = -1;
  int ep = -1;
  int tpneigh = -1;
  const int tm = edges[edge_id][3];
  int vm = -1;
  int dm = -1;
  int em = -1;
  int tmneigh = -1;
  assert(tp != tm);

  if (tp >= 0) {
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tp + d] != v0 && triangles[3 * tp + d] != v1) {
        vp = triangles[3 * tp + d];
        break;
      }
    }
    for (int d = 0; d < 3; d++) {
      if ((triangles[3 * tp + d] == v1 &&
           triangles[3 * tp + (d + 1) % 3] == vp) ||
          (triangles[3 * tp + d] == vp &&
           triangles[3 * tp + (d + 1) % 3] == v1)) {
        dp = d;
        tpneigh = tri_neigh[tp][d];
        ep = tri_edges[tp][d];
        assert(ep >= 0);
        assert((edges[ep][0] == v1 && edges[ep][1] == vp) ||
               (edges[ep][1] == v1 && edges[ep][0] == vp));
        assert((edges[ep][2] == tp && edges[ep][3] == tpneigh) ||
               (edges[ep][3] == tp && edges[ep][2] == tpneigh));
        break;
      }
    }
    assert(vp >= 0);
  }
  if (tm >= 0) {
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tm + d] != v0 && triangles[3 * tm + d] != v1) {
        vm = triangles[3 * tm + d];
        break;
      }
    }
    for (int d = 0; d < 3; d++) {
      if ((triangles[3 * tm + d] == v1 &&
           triangles[3 * tm + (d + 1) % 3] == vm) ||
          (triangles[3 * tm + d] == vm &&
           triangles[3 * tm + (d + 1) % 3] == v1)) {
        dm = d;
        tmneigh = tri_neigh[tm][d];
        em = tri_edges[tm][d];
        assert(em >= 0);
        assert((edges[em][0] == v1 && edges[em][1] == vm) ||
               (edges[em][1] == v1 && edges[em][0] == vm));
        assert((edges[em][2] == tm && edges[em][3] == tmneigh) ||
               (edges[em][3] == tm && edges[em][2] == tmneigh));
        break;
      }
    }

    assert(vm >= 0);
  }
  assert(vm != vp);

  // Create new mid-point
  const int vnew = vertices.size();
  vertices.push_back(Pt(0.5 * (vertices[v0] + vertices[v1])));
  vert_tri.resize(vert_tri.size() + 1);
  vert_neigh.resize(vert_neigh.size() + 1);
  vert_tri.back().resize(0);
  vert_neigh.back().resize(0);

  // Create new triangle(s)
  int tpnew = -1;
  if (tp >= 0) {
    tpnew = triangles.size() / 3;
    tri_neigh.push_back({-1, -1, -1});
    tri_edges.push_back({-1, -1, -1});
    triangles.push_back(vnew);
    triangles.push_back(v1);
    triangles.push_back(vp);
    // std::cout << "adding triangle " << tpnew << " = " << vnew << ", " << v1
    //           << ", " << vp << std::endl;
  }
  int tmnew = -1;
  if (tm >= 0) {
    tmnew = triangles.size() / 3;
    tri_neigh.push_back({-1, -1, -1});
    tri_edges.push_back({-1, -1, -1});
    triangles.push_back(vnew);
    triangles.push_back(v1);
    triangles.push_back(vm);
    // std::cout << "adding triangle " << tmnew << " = " << vnew << ", " << v1
    //           << ", " << vm << std::endl;
  }
  assert(tpnew != tmnew);

  // Create new edges
  edges[edge_id][1] = vnew;
  const int enew = edges.size();
  if (tp >= 0) {
    edges.push_back({vnew, v1, tpnew, tmnew});
  } else {
    edges.push_back({vnew, v1, tmnew, tpnew});
  }
  int epnew = -1;
  if (tp >= 0) {
    epnew = edges.size();
    edges.push_back({vp, vnew, tp, tpnew});
  }
  int emnew = -1;
  if (tm >= 0) {
    emnew = edges.size();
    edges.push_back({vm, vnew, tm, tmnew});
  }

  // Update old and new triangles
  if (tp >= 0) {
    bool found_vertex = false;
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tp + d] == v1) {
        triangles[3 * tp + d] = vnew;
        found_vertex = true;
        break;
      }
    }
    assert(found_vertex);
    tri_neigh[tp][dp] = tpnew;
    tri_edges[tp][dp] = epnew;
    tri_neigh[tpnew][0] = tmnew;
    tri_neigh[tpnew][1] = tpneigh;
    tri_neigh[tpnew][2] = tp;
    tri_edges[tpnew][0] = enew;
    tri_edges[tpnew][1] = ep;
    tri_edges[tpnew][2] = epnew;
    if (edges[ep][2] == tp) {
      edges[ep][2] = tpnew;
    } else if (edges[ep][3] == tp) {
      edges[ep][3] = tpnew;
    } else {
      assert(false);
    }
    if (tpneigh >= 0) {
      bool found_neigh = false;
      for (int d = 0; d < 3; d++) {
        if (tri_neigh[tpneigh][d] == tp) {
          tri_neigh[tpneigh][d] = tpnew;
          found_neigh = true;
          break;
        }
      }
      assert(found_neigh);
    }
  }
  if (tm >= 0) {
    bool found_vertex = false;
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tm + d] == v1) {
        triangles[3 * tm + d] = vnew;
        found_vertex = true;
        break;
      }
    }
    assert(found_vertex);
    tri_neigh[tm][dm] = tmnew;
    tri_edges[tm][dm] = emnew;
    tri_neigh[tmnew][0] = tpnew;
    tri_neigh[tmnew][1] = tmneigh;
    tri_neigh[tmnew][2] = tm;
    tri_edges[tmnew][0] = enew;
    tri_edges[tmnew][1] = em;
    tri_edges[tmnew][2] = emnew;
    if (edges[em][2] == tm) {
      edges[em][2] = tmnew;
    } else if (edges[em][3] == tm) {
      edges[em][3] = tmnew;
    } else {
      assert(false);
    }
    if (tmneigh >= 0) {
      bool found_neigh = false;
      for (int d = 0; d < 3; d++) {
        if (tri_neigh[tmneigh][d] == tm) {
          tri_neigh[tmneigh][d] = tmnew;
          found_neigh = true;
          break;
        }
      }
      assert(found_neigh);
    }
  }

  // Update vertice neighbours
  if (vp >= 0) {
    vert_neigh[vp].push_back(vnew);
  }
  if (vm >= 0) {
    vert_neigh[vm].push_back(vnew);
  }
  bool found_neigh = false;
  for (int d = 0; d < vert_neigh[v0].size(); d++) {
    if (vert_neigh[v0][d] == v1) {
      vert_neigh[v0][d] = vnew;
      found_neigh = true;
      break;
    }
  }
  assert(found_neigh);
  found_neigh = false;
  for (int d = 0; d < vert_neigh[v1].size(); d++) {
    if (vert_neigh[v1][d] == v0) {
      vert_neigh[v1][d] = vnew;
      found_neigh = true;
      break;
    }
  }
  assert(found_neigh);
  vert_neigh[vnew].push_back(v0);
  vert_neigh[vnew].push_back(v1);
  if (vp >= 0) {
    vert_neigh[vnew].push_back(vp);
  }
  if (vm >= 0) {
    vert_neigh[vnew].push_back(vm);
  }

  // Update vertex->tri connectivity
  if (tp >= 0) {
    bool found_tri = false;
    for (int d = 0; d < vert_tri[v1].size(); d++) {
      if (vert_tri[v1][d] == tp) {
        vert_tri[v1][d] = tpnew;
        found_tri = true;
        break;
      }
    }
    assert(found_tri);
    vert_tri[vp].push_back(tpnew);
    vert_tri[vnew].push_back(tp);
    vert_tri[vnew].push_back(tpnew);
  }
  if (tm >= 0) {
    bool found_tri = false;
    for (int d = 0; d < vert_tri[v1].size(); d++) {
      if (vert_tri[v1][d] == tm) {
        vert_tri[v1][d] = tmnew;
        found_tri = true;
        break;
      }
    }
    assert(found_tri);
    vert_tri[vm].push_back(tmnew);
    vert_tri[vnew].push_back(tm);
    vert_tri[vnew].push_back(tmnew);
  }

  if (vp >= 0) {
    assert(noDuplicates(vert_tri[vp]));
    assert(noDuplicates(vert_neigh[vp]));
  }
  if (vm >= 0) {
    assert(noDuplicates(vert_tri[vm]));
    assert(noDuplicates(vert_neigh[vm]));
  }
  assert(noDuplicates(vert_tri[v0]));
  assert(noDuplicates(vert_tri[v1]));
  assert(noDuplicates(vert_tri[vnew]));
  assert(noDuplicates(vert_neigh[v0]));
  assert(noDuplicates(vert_neigh[v1]));
  assert(noDuplicates(vert_neigh[vnew]));
}

template <class VertexList, class TriList, class EdgeList,
          class TriConnectivity, class TriEdges, class VertConnectivity,
          class VertexValence>
void collapseEdge(const UnsignedIndex_t edge_id, VertexList& vertices,
                  TriList& triangles, EdgeList& edges,
                  TriConnectivity& tri_neigh, TriEdges& tri_edges,
                  VertConnectivity& vert_tri, VertexValence& vert_neigh) {
  // Start and end point
  const int v0 = edges[edge_id][0];
  const int v1 = edges[edge_id][1];
  assert(v0 != v1);
  assert(v0 >= 0 && v0 < vertices.size());
  assert(v1 >= 0 && v1 < vertices.size());

  // Find opposite triagnles and vertices
  const int tp = edges[edge_id][2];
  assert(tp >= 0);
  int vp = -1;
  int dp0 = -1;
  int dp1 = -1;
  int ep0 = -1;
  int ep1 = -1;
  int tpneigh0 = -1;
  int tpneigh1 = -1;
  const int tm = edges[edge_id][3];
  int vm = -1;
  int dm0 = -1;
  int dm1 = -1;
  int em0 = -1;
  int em1 = -1;
  int tmneigh0 = -1;
  int tmneigh1 = -1;
  assert(tp != tm);

  if (tp >= 0) {
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tp + d] != v0 && triangles[3 * tp + d] != v1) {
        vp = triangles[3 * tp + d];
        break;
      }
    }
    for (int d = 0; d < 3; d++) {
      if ((triangles[3 * tp + d] == v1 &&
           triangles[3 * tp + (d + 1) % 3] == vp) ||
          (triangles[3 * tp + d] == vp &&
           triangles[3 * tp + (d + 1) % 3] == v1)) {
        dp1 = d;
        tpneigh1 = tri_neigh[tp][d];
        ep1 = tri_edges[tp][d];
        assert(ep1 >= 0);
        assert((edges[ep1][2] == tp && edges[ep1][3] == tpneigh1) ||
               (edges[ep1][3] == tp && edges[ep1][2] == tpneigh1));
      } else if ((triangles[3 * tp + d] == v0 &&
                  triangles[3 * tp + (d + 1) % 3] == vp) ||
                 (triangles[3 * tp + d] == vp &&
                  triangles[3 * tp + (d + 1) % 3] == v0)) {
        dp0 = d;
        tpneigh0 = tri_neigh[tp][d];
        ep0 = tri_edges[tp][d];
        assert(ep0 >= 0);
        assert((edges[ep0][2] == tp && edges[ep0][3] == tpneigh0) ||
               (edges[ep0][3] == tp && edges[ep0][2] == tpneigh0));
      }
    }
    assert(vp >= 0);
  }
  if (tm >= 0) {
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tm + d] != v0 && triangles[3 * tm + d] != v1) {
        vm = triangles[3 * tm + d];
        break;
      }
    }
    for (int d = 0; d < 3; d++) {
      if ((triangles[3 * tm + d] == v1 &&
           triangles[3 * tm + (d + 1) % 3] == vm) ||
          (triangles[3 * tm + d] == vm &&
           triangles[3 * tm + (d + 1) % 3] == v1)) {
        dm1 = d;
        tmneigh1 = tri_neigh[tm][d];
        em1 = tri_edges[tm][d];
        assert(em1 >= 0);
        assert((edges[em1][2] == tm && edges[em1][3] == tmneigh1) ||
               (edges[em1][3] == tm && edges[em1][2] == tmneigh1));
      } else if ((triangles[3 * tm + d] == v0 &&
                  triangles[3 * tm + (d + 1) % 3] == vm) ||
                 (triangles[3 * tm + d] == vm &&
                  triangles[3 * tm + (d + 1) % 3] == v0)) {
        dm0 = d;
        tmneigh0 = tri_neigh[tm][d];
        em0 = tri_edges[tm][d];
        assert(em0 >= 0);
        assert((edges[em0][2] == tm && edges[em0][3] == tmneigh0) ||
               (edges[em0][3] == tm && edges[em0][2] == tmneigh0));
      }
    }
    assert(vm >= 0);
  }
  assert(vm != vp);

  // Update start triangles
  for (int i = 0; i < vert_tri[v0].size(); i++) {
    const int tri = vert_tri[v0][i];
    bool found_vertex = false;
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tri + d] == v0) {
        triangles[3 * tri + d] = v1;
        found_vertex = true;
      }
      const int edge = tri_edges[tri][d];
      if (edges[edge][0] == v0) {
        edges[edge][0] = v1;
      } else if (edges[edge][1] == v0) {
        edges[edge][1] = v1;
      }
    }
    assert(found_vertex);
    bool tri_not_here = true;
    for (int j = 0; j < vert_tri[v1].size(); j++) {
      if (vert_tri[v1][j] == tri) {
        tri_not_here = false;
        break;
      }
    }
    if (tri_not_here) {
      // if (tri != tp && tri != tm) {
      vert_tri[v1].push_back(tri);
    }
  }
  vert_tri[v0].resize(0);

  for (int i = 0; i < vert_neigh[v0].size(); i++) {
    const int neigh = vert_neigh[v0][i];
    if (neigh != v1) {
      bool neigh_has_v1 = false;
      int index_v0_neigh = -1;
      int count_v0 = 0;
      for (int j = 0; j < vert_neigh[neigh].size(); j++) {
        if (vert_neigh[neigh][j] == v1) {
          neigh_has_v1 = true;
        } else if (vert_neigh[neigh][j] == v0) {
          index_v0_neigh = j;
          count_v0++;
        }
      }
      assert(count_v0 == 1);
      assert(index_v0_neigh >= 0);
      if (!neigh_has_v1) {
        vert_neigh[neigh][index_v0_neigh] = v1;
        bool v1_had_neigh = false;
        for (int k = 0; k < vert_neigh[v1].size(); k++) {
          if (vert_neigh[v1][k] == neigh) {
            v1_had_neigh = true;
            break;
          }
        }
        assert(!v1_had_neigh);
        vert_neigh[v1].push_back(neigh);
      } else {
        vert_neigh[neigh].erase(vert_neigh[neigh].begin() + index_v0_neigh);
      }
    } else {
      bool found_vertex = false;
      for (int j = 0; j < vert_neigh[v1].size(); j++) {
        if (vert_neigh[v1][j] == v0) {
          vert_neigh[v1].erase(vert_neigh[v1].begin() + j);
          found_vertex = true;
          break;
        }
      }
      assert(found_vertex);
    }
  }
  vert_neigh[v0].resize(0);

  if (tpneigh0 >= 0) {
    bool found_edge = false;
    bool found_neigh = false;
    for (int d = 0; d < 3; d++) {
      if (tri_edges[tpneigh0][d] == ep0) {
        tri_edges[tpneigh0][d] = ep1;
        found_edge = true;
      }
      if (tri_neigh[tpneigh0][d] == tp) {
        tri_neigh[tpneigh0][d] = tpneigh1;
        found_neigh = true;
      }
    }
    assert(found_edge);
    assert(found_neigh);
  }
  if (tpneigh1 >= 0) {
    bool found_neigh = false;
    for (int d = 0; d < 3; d++) {
      if (tri_neigh[tpneigh1][d] == tp) {
        tri_neigh[tpneigh1][d] = tpneigh0;
        found_neigh = true;
        break;
      }
    }
    assert(found_neigh);
  }
  if (tmneigh0 >= 0) {
    bool found_edge = false;
    bool found_neigh = false;
    for (int d = 0; d < 3; d++) {
      if (tri_edges[tmneigh0][d] == em0) {
        found_edge = true;
        tri_edges[tmneigh0][d] = em1;
      }
      if (tri_neigh[tmneigh0][d] == tm) {
        tri_neigh[tmneigh0][d] = tmneigh1;
        found_neigh = true;
      }
    }
    assert(found_edge);
    assert(found_neigh);
  }
  if (tmneigh1 >= 0) {
    bool found_neigh = false;
    for (int d = 0; d < 3; d++) {
      if (tri_neigh[tmneigh1][d] == tm) {
        tri_neigh[tmneigh1][d] = tmneigh0;
        found_neigh = true;
        break;
      }
    }
    assert(found_neigh);
  }
  if (ep1 >= 0) {
    if (edges[ep1][2] == tp) {
      assert(edges[ep1][3] == tpneigh1);
      edges[ep1][2] = tpneigh0;
    } else if (edges[ep1][3] == tp) {
      assert(edges[ep1][2] == tpneigh1);
      edges[ep1][3] = tpneigh0;
    } else {
      assert(false);
    }
  }
  if (em1 >= 0) {
    if (edges[em1][2] == tm) {
      assert(edges[em1][3] == tmneigh1);
      edges[em1][2] = tmneigh0;
    } else if (edges[em1][3] == tm) {
      assert(edges[em1][2] == tmneigh1);
      edges[em1][3] = tmneigh0;
    } else {
      assert(false);
    }
  }

  // Remove triangles adjacent to edge
  if (tp >= 0) {
    for (int d = 0; d < 3; d++) {
      triangles[3 * tp + d] = -1;
    }
  }
  if (tm >= 0) {
    for (int d = 0; d < 3; d++) {
      triangles[3 * tm + d] = -1;
    }
  }
  if (tp >= 0) {
    bool found = false;
    for (int i = 0; i < vert_tri[v1].size(); i++) {
      if (vert_tri[v1][i] == tp) {
        vert_tri[v1].erase(vert_tri[v1].begin() + i);
        found = true;
        break;
      }
    }
    assert(found);
    found = false;
    for (int i = 0; i < vert_tri[vp].size(); i++) {
      if (vert_tri[vp][i] == tp) {
        vert_tri[vp].erase(vert_tri[vp].begin() + i);
        found = true;
        break;
      }
    }
    assert(found);
  }
  if (tm >= 0) {
    bool found = false;
    for (int i = 0; i < vert_tri[v1].size(); i++) {
      if (vert_tri[v1][i] == tm) {
        vert_tri[v1].erase(vert_tri[v1].begin() + i);
        found = true;
        break;
      }
    }
    assert(found);
    found = false;
    for (int i = 0; i < vert_tri[vm].size(); i++) {
      if (vert_tri[vm][i] == tm) {
        vert_tri[vm].erase(vert_tri[vm].begin() + i);
        found = true;
        break;
      }
    }
    assert(found);
  }

  // Remove edge
  edges[edge_id] = {-1, -1, -1, -1};
  if (ep0 >= 0) {
    edges[ep0] = {-1, -1, -1, -1};
  }
  if (em0 >= 0) {
    edges[em0] = {-1, -1, -1, -1};
  }
}

template <class VertexList, class TriList, class EdgeList,
          class TriConnectivity, class TriEdges, class VertConnectivity,
          class VertexValence>
void flipEdge(const UnsignedIndex_t edge_id, VertexList& vertices,
              TriList& triangles, EdgeList& edges, TriConnectivity& tri_neigh,
              TriEdges& tri_edges, VertConnectivity& vert_tri,
              VertexValence& vert_neigh) {
  // Start and end point
  const int v0 = edges[edge_id][0];
  const int v1 = edges[edge_id][1];
  assert(v0 != v1);
  assert(v0 >= 0 && v0 < vertices.size());
  assert(v1 >= 0 && v1 < vertices.size());

  // Find opposite triagnles and vertices
  const int tp = edges[edge_id][2];
  assert(tp >= 0);
  int vp = -1;
  int dp0 = -1;
  int dp1 = -1;
  int ep0 = -1;
  int ep1 = -1;
  int tpneigh0 = -1;
  int tpneigh1 = -1;
  const int tm = edges[edge_id][3];
  assert(tm >= 0);
  int vm = -1;
  int dm0 = -1;
  int dm1 = -1;
  int em0 = -1;
  int em1 = -1;
  int tmneigh0 = -1;
  int tmneigh1 = -1;
  assert(tp != tm);

  // std::cout << "Flipping egde " << edge_id << " = " << v0 << " " << v1 << " "
  //           << tp << " " << tm << std::endl;
  for (int d = 0; d < 3; d++) {
    if (triangles[3 * tp + d] != v0 && triangles[3 * tp + d] != v1) {
      vp = triangles[3 * tp + d];
      break;
    }
  }
  for (int d = 0; d < 3; d++) {
    if ((triangles[3 * tp + d] == v1 &&
         triangles[3 * tp + (d + 1) % 3] == vp) ||
        (triangles[3 * tp + d] == vp &&
         triangles[3 * tp + (d + 1) % 3] == v1)) {
      dp1 = d;
      tpneigh1 = tri_neigh[tp][d];
      ep1 = tri_edges[tp][d];
      assert(ep1 >= 0);
      assert((edges[ep1][0] == v1 && edges[ep1][1] == vp) ||
             (edges[ep1][1] == v1 && edges[ep1][0] == vp));
      assert((edges[ep1][2] == tp && edges[ep1][3] == tpneigh1) ||
             (edges[ep1][3] == tp && edges[ep1][2] == tpneigh1));
    } else if ((triangles[3 * tp + d] == v0 &&
                triangles[3 * tp + (d + 1) % 3] == vp) ||
               (triangles[3 * tp + d] == vp &&
                triangles[3 * tp + (d + 1) % 3] == v0)) {
      dp0 = d;
      tpneigh0 = tri_neigh[tp][d];
      ep0 = tri_edges[tp][d];
      assert(ep0 >= 0);
      assert((edges[ep0][0] == v0 && edges[ep0][1] == vp) ||
             (edges[ep0][1] == v0 && edges[ep0][0] == vp));
      assert((edges[ep0][2] == tp && edges[ep0][3] == tpneigh0) ||
             (edges[ep0][3] == tp && edges[ep0][2] == tpneigh0));
    }
  }
  assert(vp >= 0);
  for (int d = 0; d < 3; d++) {
    if (triangles[3 * tm + d] != v0 && triangles[3 * tm + d] != v1) {
      vm = triangles[3 * tm + d];
      break;
    }
  }
  for (int d = 0; d < 3; d++) {
    if ((triangles[3 * tm + d] == v1 &&
         triangles[3 * tm + (d + 1) % 3] == vm) ||
        (triangles[3 * tm + d] == vm &&
         triangles[3 * tm + (d + 1) % 3] == v1)) {
      dm1 = d;
      tmneigh1 = tri_neigh[tm][d];
      em1 = tri_edges[tm][d];
      assert(em1 >= 0);
      assert((edges[em1][0] == v1 && edges[em1][1] == vm) ||
             (edges[em1][1] == v1 && edges[em1][0] == vm));
      assert((edges[em1][2] == tm && edges[em1][3] == tmneigh1) ||
             (edges[em1][3] == tm && edges[em1][2] == tmneigh1));
    } else if ((triangles[3 * tm + d] == v0 &&
                triangles[3 * tm + (d + 1) % 3] == vm) ||
               (triangles[3 * tm + d] == vm &&
                triangles[3 * tm + (d + 1) % 3] == v0)) {
      dm0 = d;
      tmneigh0 = tri_neigh[tm][d];
      em0 = tri_edges[tm][d];
      assert(em0 >= 0);
      assert((edges[em0][0] == v0 && edges[em0][1] == vm) ||
             (edges[em0][1] == v0 && edges[em0][0] == vm));
      assert((edges[em0][2] == tm && edges[em0][3] == tmneigh0) ||
             (edges[em0][3] == tm && edges[em0][2] == tmneigh0));
    }
  }
  assert(vm >= 0);
  assert(vm != vp);

  bool do_flip = true;
  if (tpneigh0 >= 0 && tmneigh0 >= 0) {
    if (tpneigh0 == tmneigh0) {
      do_flip = false;
    }
  }
  if (tpneigh1 >= 0 && tmneigh1 >= 0) {
    if (tpneigh1 == tmneigh1) {
      do_flip = false;
    }
  }

  if (do_flip) {
    // Flip edge
    edges[edge_id][0] = vp;
    edges[edge_id][1] = vm;

    bool found_vertex = false;
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tp + d] == v0) {
        triangles[3 * tp + (d + 1) % 3] = vm;
        triangles[3 * tp + (d + 2) % 3] = vp;
        tri_neigh[tp][d] = tmneigh0;
        tri_neigh[tp][(d + 1) % 3] = tm;
        tri_neigh[tp][(d + 2) % 3] = tpneigh0;
        tri_edges[tp][d] = em0;
        tri_edges[tp][(d + 1) % 3] = edge_id;
        tri_edges[tp][(d + 2) % 3] = ep0;
        found_vertex = true;
        break;
      }
    }
    found_vertex = false;
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tm + d] == v1) {
        triangles[3 * tm + (d + 1) % 3] = vm;
        triangles[3 * tm + (d + 2) % 3] = vp;
        tri_neigh[tm][d] = tmneigh1;
        tri_neigh[tm][(d + 1) % 3] = tp;
        tri_neigh[tm][(d + 2) % 3] = tpneigh1;
        tri_edges[tm][d] = em1;
        tri_edges[tm][(d + 1) % 3] = edge_id;
        tri_edges[tm][(d + 2) % 3] = ep1;
        found_vertex = true;
        break;
      }
    }
    assert(found_vertex);
    if (tpneigh1 >= 0) {
      bool found_tri = false;
      for (int d = 0; d < 3; d++) {
        if (tri_neigh[tpneigh1][d] == tp) {
          tri_neigh[tpneigh1][d] = tm;
          found_tri = true;
          break;
        }
      }
      assert(found_tri);
    }
    assert(found_vertex);
    if (tmneigh0 >= 0) {
      bool found_tri = false;
      for (int d = 0; d < 3; d++) {
        if (tri_neigh[tmneigh0][d] == tm) {
          tri_neigh[tmneigh0][d] = tp;
          found_tri = true;
          break;
        }
      }
      assert(found_tri);
    }
    if (edges[ep1][2] == tp) {
      edges[ep1][2] = tm;
    } else if (edges[ep1][3] == tp) {
      edges[ep1][3] = tm;
    } else {
      assert(false);
    }
    if (edges[em0][2] == tm) {
      edges[em0][2] = tp;
    } else if (edges[em0][3] == tm) {
      edges[em0][3] = tp;
    } else {
      assert(false);
    }

    vert_tri[vp].push_back(tm);
    vert_tri[vm].push_back(tp);
    for (int i = 0; i < vert_tri[v0].size(); i++) {
      if (vert_tri[v0][i] == tm) {
        vert_tri[v0].erase(vert_tri[v0].begin() + i);
        break;
      }
    }
    for (int i = 0; i < vert_tri[v1].size(); i++) {
      if (vert_tri[v1][i] == tp) {
        vert_tri[v1].erase(vert_tri[v1].begin() + i);
        break;
      }
    }

    bool has_vm = false;
    for (int i = 0; i < vert_neigh[vp].size(); i++) {
      if (vert_neigh[vp][i] == vm) {
        has_vm = true;
        break;
      }
    }
    if (!has_vm) {
      vert_neigh[vp].push_back(vm);
    }
    bool has_vp = false;
    for (int i = 0; i < vert_neigh[vm].size(); i++) {
      if (vert_neigh[vm][i] == vp) {
        has_vp = true;
        break;
      }
    }
    if (!has_vp) {
      vert_neigh[vm].push_back(vp);
    }
    for (int i = 0; i < vert_neigh[v0].size(); i++) {
      if (vert_neigh[v0][i] == v1) {
        vert_neigh[v0].erase(vert_neigh[v0].begin() + i);
        break;
      }
    }
    for (int i = 0; i < vert_neigh[v1].size(); i++) {
      if (vert_neigh[v1][i] == v0) {
        vert_neigh[v1].erase(vert_neigh[v1].begin() + i);
        break;
      }
    }

    assert(noDuplicates(vert_tri[vp]));
    assert(noDuplicates(vert_tri[vm]));
    assert(noDuplicates(vert_tri[v0]));
    assert(noDuplicates(vert_tri[v1]));
    assert(noDuplicates(vert_neigh[vp]));
    assert(noDuplicates(vert_neigh[vm]));
    assert(noDuplicates(vert_neigh[v0]));
    assert(noDuplicates(vert_neigh[v1]));
    assert((edges[em0][2] == tmneigh0 && edges[em0][3] == tp) ||
           (edges[em0][3] == tmneigh0 && edges[em0][2] == tp));
    assert((edges[ep0][2] == tpneigh0 && edges[ep0][3] == tp) ||
           (edges[ep0][3] == tpneigh0 && edges[ep0][2] == tp));
    assert((edges[em1][2] == tmneigh1 && edges[em1][3] == tm) ||
           (edges[em1][3] == tmneigh1 && edges[em1][2] == tm));
    assert((edges[ep1][2] == tpneigh1 && edges[ep1][3] == tm) ||
           (edges[ep1][3] == tpneigh1 && edges[ep1][2] == tm));
  }
}

template <class VertexList, class EdgeList, class VertexValence>
void relaxVertices(VertexList& vertices, EdgeList& edges,
                   VertexValence& vert_neigh,
                   const UnsignedIndex_t fixed_vertices,
                   const UnsignedIndex_t iterations,
                   const double relax_factor) {
  if (vertices.size() > fixed_vertices) {
    std::vector<bool> fixed;
    fixed.resize(vertices.size(), false);
    for (UnsignedIndex_t i = 0; i < edges.size(); i++) {
      if ((edges[i][0] >= 0 || edges[i][1] >= 0) &&
          (edges[i][2] < 0 || edges[i][3] < 0)) {
        fixed[edges[i][0]] = true;
        fixed[edges[i][1]] = true;
      }
    }
    assert(vert_neigh.size() == vertices.size());
    std::vector<Pt> barycenters;
    barycenters.resize(vertices.size() - fixed_vertices);
    for (UnsignedIndex_t it = 0; it < iterations; it++) {
      for (UnsignedIndex_t i = fixed_vertices; i < vertices.size(); i++) {
        barycenters[i - fixed_vertices] = Pt(0, 0, 0);
        if (!fixed[i] && vert_neigh[i].size() > 0) {
          assert(noDuplicates(vert_neigh[i]));
          for (UnsignedIndex_t j = 0; j < vert_neigh[i].size(); j++) {
            assert(vert_neigh[i][j] >= 0);
            assert(vert_neigh[i][j] < vertices.size());
            barycenters[i - fixed_vertices] += vertices[vert_neigh[i][j]];
          }
          barycenters[i - fixed_vertices] /=
              static_cast<double>(vert_neigh[i].size());
        }
      }
      for (UnsignedIndex_t i = fixed_vertices; i < vertices.size(); i++) {
        if (!fixed[i] && vert_neigh[i].size() > 0) {
          vertices[i] = (1.0 - relax_factor) * vertices[i] +
                        relax_factor * barycenters[i - fixed_vertices];
        }
      }
    }
  }
}

template <class VertexList>
void projectOnSurface(VertexList& vertices, const AlignedParaboloid paraboloid,
                      const UnsignedIndex_t fixed_vertices) {
  if (vertices.size() > fixed_vertices) {
    for (UnsignedIndex_t i = fixed_vertices; i < vertices.size(); i++) {
      const double x = vertices[i][0];
      const double y = vertices[i][1];
      vertices[i][2] = -paraboloid.a() * x * x - paraboloid.b() * y * y;
    }
  }
}

template <class VertexList, class TriList, class EdgeList,
          class TriConnectivity, class TriEdges, class VertConnectivity,
          class VertexValence>
void splitLongEdges(VertexList& vertices, TriList& triangles, EdgeList& edges,
                    TriConnectivity& tri_neigh, TriEdges& tri_edges,
                    VertConnectivity& vert_tri, VertexValence& vert_neigh,
                    const double high) {
  const double high_sq = high * high;
  bool remains_long_edges = true;
  UnsignedIndex_t it = 0;
  while (it < 1 && remains_long_edges) {
    it++;
    const int n_edges = edges.size();
    remains_long_edges = false;
    for (int i = 0; i < n_edges; ++i) {
      const int v0 = edges[i][0];
      const int v1 = edges[i][1];
      if (v0 < 0 || v1 < 0) {
        continue;
      }
      const double edge_length_sq =
          squaredMagnitude(vertices[v0] - vertices[v1]);
      if (edge_length_sq > high_sq) {
        remains_long_edges = true;
        splitEdge(i, vertices, triangles, edges, tri_neigh, tri_edges, vert_tri,
                  vert_neigh);
      }
    }
  }
}

template <class VertexList, class TriList, class EdgeList,
          class TriConnectivity, class TriEdges, class VertConnectivity,
          class VertexValence>
void collapseShortEdges(VertexList& vertices, TriList& triangles,
                        EdgeList& edges, TriConnectivity& tri_neigh,
                        TriEdges& tri_edges, VertConnectivity& vert_tri,
                        VertexValence& vert_neigh, const double low,
                        const double high,
                        const UnsignedIndex_t fixed_vertices) {
  const double low_sq = low * low;
  const double high_sq = high * high;
  bool remains_short_edges = true;
  UnsignedIndex_t it = 0;
  while (it < 1 && remains_short_edges) {
    it++;
    const int n_edges = edges.size();
    remains_short_edges = false;
    for (int i = 0; i < n_edges; ++i) {
      const int v0 = edges[i][0];
      const int v1 = edges[i][1];
      if (v0 < fixed_vertices || v1 < 0 || vert_tri[v0].size() < 5 ||
          vert_tri[v1].size() < 5) {
        continue;
      }
      const double edge_length_sq =
          squaredMagnitude(vertices[v0] - vertices[v1]);
      if (edge_length_sq < low_sq) {
        bool collapse = true;
        for (int j = 0; j < vert_neigh[v0].size(); j++) {
          const int neigh = vert_neigh[v0][j];
          const double test_edge_length_sq =
              squaredMagnitude(vertices[neigh] - vertices[v1]);
          if (test_edge_length_sq > high_sq) {
            collapse = false;
            break;
          }
        }
        if (collapse) {
          std::cout << "Collapsing edge " << i << " " << v0 << " -- " << v1
                    << " -- " << edges[i][2] << " -- " << edges[i][3]
                    << std::endl;
          remains_short_edges = true;
          collapseEdge(i, vertices, triangles, edges, tri_neigh, tri_edges,
                       vert_tri, vert_neigh);
        }
        // else {
        //   break;
        // }
      }
    }
  }
}

template <class VertexList, class TriList, class EdgeList,
          class TriConnectivity, class TriEdges, class VertConnectivity,
          class VertexValence>
void equalizeValence(VertexList& vertices, TriList& triangles, EdgeList& edges,
                     TriConnectivity& tri_neigh, TriEdges& tri_edges,
                     VertConnectivity& vert_tri, VertexValence& vert_neigh,
                     const UnsignedIndex_t fixed_vertices) {
  const int n_edges = edges.size();
  for (int i = 0; i < n_edges; ++i) {
    const int v0 = edges[i][0];
    const int v1 = edges[i][1];
    const int tp = edges[i][2];
    const int tm = edges[i][3];
    if (v0 < 0 || v1 < 0 || vert_neigh[v0].size() == 0 ||
        vert_neigh[v1].size() == 0 || tp < 0 || tm < 0) {
      continue;
    }
    int vp = -1;
    int vm = -1;
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tp + d] != v0 && triangles[3 * tp + d] != v1) {
        vp = triangles[3 * tp + d];
        break;
      }
    }
    assert(vp >= 0);
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tm + d] != v0 && triangles[3 * tm + d] != v1) {
        vm = triangles[3 * tm + d];
        break;
      }
    }
    assert(vm >= 0);
    int val0 = vert_neigh[v0].size();
    int val1 = vert_neigh[v1].size();
    int valp = vert_neigh[vp].size();
    int valm = vert_neigh[vm].size();
    int target0 = v0 < fixed_vertices ? 4 : 6;
    int target1 = v1 < fixed_vertices ? 4 : 6;
    int targetp = vp < fixed_vertices ? 4 : 6;
    int targetm = vm < fixed_vertices ? 4 : 6;
    int deviation_pre = std::abs(val0 - target0) + std::abs(val1 - target1) +
                        std::abs(valp - targetp) + std::abs(valm - targetm);
    flipEdge(i, vertices, triangles, edges, tri_neigh, tri_edges, vert_tri,
             vert_neigh);
    val0 = vert_neigh[v0].size();
    val1 = vert_neigh[v1].size();
    valp = vert_neigh[vp].size();
    valm = vert_neigh[vm].size();
    int deviation_post = std::abs(val0 - target0) + std::abs(val1 - target1) +
                         std::abs(valp - targetp) + std::abs(valm - targetm);
    if (deviation_pre <= deviation_post) {
      flipEdge(i, vertices, triangles, edges, tri_neigh, tri_edges, vert_tri,
               vert_neigh);
    }
  }
}

template <class VertexList, class EdgeList, class TriList>
void reMeshPolygon(VertexList& vertices, EdgeList& edges, TriList& triangles,
                   const double length_scale,
                   const AlignedParaboloid paraboloid) {
  const double low = 4.0 / 5.0 * length_scale;
  const double high = 4.0 / 3.0 * length_scale;
  const int fixed_vertices = vertices.size();
  const int original_tris = triangles.size() / 3;

  // Construct edges and connectivity
  std::vector<std::array<int, 3>> tri_neigh;
  tri_neigh.resize(original_tris, std::array<int, 3>({-1, -1, -1}));
  std::vector<std::array<int, 3>> tri_edges;
  tri_edges.resize(original_tris, std::array<int, 3>({-1, -1, -1}));
  std::vector<std::vector<int>> vert_tri;
  vert_tri.resize(fixed_vertices);
  std::vector<std::vector<int>> vert_neigh;
  vert_neigh.resize(fixed_vertices);

  for (int i = 0; i < fixed_vertices; i++) {
    vert_tri[i].resize(0);
    vert_neigh[i].resize(0);
  }

  // List of triangles linked to a vertex
  for (int i = 0; i < original_tris; ++i) {
    for (int d = 0; d < 3; ++d) {
      vert_tri[triangles[3 * i + d]].push_back(i);
      tri_neigh[i][d] = -1;
      tri_edges[i][d] = -1;
    }
  }

  // Triangle neighbours
  for (int i = 0; i < original_tris; ++i) {
    // Loop over the 3 edges
    for (int d = 0; d < 3; ++d) {
      // tri_neigh[i][d] = -1;
      const int v0 = triangles[3 * i + d];
      const int v1 = triangles[3 * i + (d + 1) % 3];
      // Loop over triangles attached to start point
      bool found = false;
      for (int j = 0; j < vert_tri[v0].size(); ++j) {
        const int neigh = vert_tri[v0][j];
        if (neigh != i) {
          for (int k = 0; k < 3; ++k) {
            if (triangles[3 * neigh + k] == v1) {
              tri_neigh[i][d] = neigh;
              found = true;
              assert(triangles[3 * neigh + (k + 1) % 3] == v0 ||
                     triangles[3 * neigh + (k + 2) % 3] == v0);
              break;
            }
          }
        }
        if (found) {
          break;
        }
      }
      // if (!found) {
      //   std::cout << "Triangle " << i << " has boundary" << std::endl;
      // }
      assert(tri_neigh[i][d] != i);
    }
  }

  // Build edges
  for (int i = 0; i < original_tris; ++i) {
    // Loop over the 3 edges
    for (int d = 0; d < 3; ++d) {
      const int neigh = tri_neigh[i][d];
      if (i > neigh) {
        tri_edges[i][d] = edges.size();
        if (neigh >= 0) {
          bool found_edge = false;
          for (int k = 0; k < 3; k++) {
            if (tri_neigh[neigh][k] == i) {
              tri_edges[neigh][k] = edges.size();
              found_edge = true;
              break;
            }
          }
          // assert(found_edge);
        }
        edges.push_back({static_cast<int>(triangles[3 * i + d]),
                         static_cast<int>(triangles[3 * i + (d + 1) % 3]),
                         static_cast<int>(i), tri_neigh[i][d]});
      }
    }
  }

  bool edges_are_valid = true;
  for (int i = 0; i < original_tris; ++i) {
    for (int d = 0; d < 3; ++d) {
      if (tri_edges[i][d] < 0) {
        edges_are_valid = false;
      }
      // assert(tri_edges[i][d] >= 0);
    }
  }

  if (edges_are_valid) {
    // Build valence
    const int original_edges = edges.size();
    for (int i = 0; i < original_edges; ++i) {
      vert_neigh[edges[i][0]].push_back(edges[i][1]);
      vert_neigh[edges[i][1]].push_back(edges[i][0]);
    }

    bool no_duplicates = true;
    for (int i = 0; i < fixed_vertices; i++) {
      if (!noDuplicates(vert_neigh[i])) {
        no_duplicates = false;
        // std::cout << "Duplicate neighbours" << std::endl;
        break;
      }
      if (!noDuplicates(vert_tri[i])) {
        no_duplicates = false;
        // std::cout << "Duplicate triangles" << std::endl;
        break;
      }
    }

    if (no_duplicates) {
      for (int i = 0; i < 5; ++i) {
        splitLongEdges(vertices, triangles, edges, tri_neigh, tri_edges,
                       vert_tri, vert_neigh, high);
        // TODO: fix bug in edge collapse
        // collapseShortEdges(vertices, triangles, edges, tri_neigh, tri_edges,
        //                    vert_tri, vert_neigh, low, high, fixed_vertices);
        equalizeValence(vertices, triangles, edges, tri_neigh, tri_edges,
                        vert_tri, vert_neigh, fixed_vertices);
        relaxVertices(vertices, edges, vert_neigh, fixed_vertices, 10, 0.5);
        projectOnSurface(vertices, paraboloid, fixed_vertices);
      }
    }
  }
}

template <class MomentType, class SurfaceType>
MomentType& AddSurfaceOutput<MomentType, SurfaceType>::getMoments(void) {
  return volume_moments_m;
}

template <class MomentType, class SurfaceType>
const MomentType& AddSurfaceOutput<MomentType, SurfaceType>::getMoments(
    void) const {
  return volume_moments_m;
}

template <class MomentType, class SurfaceType>
SurfaceType& AddSurfaceOutput<MomentType, SurfaceType>::getSurface(void) {
  return surface_m;
}

template <class MomentType, class SurfaceType>
const SurfaceType& AddSurfaceOutput<MomentType, SurfaceType>::getSurface(
    void) const {
  return surface_m;
}

inline ParametrizedSurfaceOutput::ParametrizedSurfaceOutput()
    : knows_surface_area_m{false},
      knows_avg_normal_m{false},
      knows_int_mean_curv_m{false},
      knows_int_gaussian_curv_m{false},
      length_scale_m{-1.0} {}

inline ParametrizedSurfaceOutput::ParametrizedSurfaceOutput(
    const Paraboloid& a_paraboloid)
    : paraboloid_m{a_paraboloid},
      knows_surface_area_m{false},
      knows_avg_normal_m{false},
      knows_int_mean_curv_m{false},
      knows_int_gaussian_curv_m{false},
      length_scale_m{-1.0} {}

inline ParametrizedSurfaceOutput::ParametrizedSurfaceOutput(
    ParametrizedSurfaceOutput&& a_rhs)
    : paraboloid_m(a_rhs.paraboloid_m),
      arc_list_m(std::move(a_rhs.arc_list_m)),
      knows_surface_area_m(a_rhs.knows_surface_area_m),
      surface_area_m(a_rhs.surface_area_m),
      knows_avg_normal_m(a_rhs.knows_avg_normal_m),
      avg_normal_m(a_rhs.avg_normal_m),
      knows_int_mean_curv_m(a_rhs.knows_int_mean_curv_m),
      int_mean_curv_m(a_rhs.int_mean_curv_m),
      knows_int_gaussian_curv_m(a_rhs.knows_int_gaussian_curv_m),
      int_gaussian_curv_m(a_rhs.int_gaussian_curv_m),
      length_scale_m{a_rhs.length_scale_m},
      pt_from_bezier_split_m(std::move(a_rhs.pt_from_bezier_split_m)) {}

inline ParametrizedSurfaceOutput::ParametrizedSurfaceOutput(
    const ParametrizedSurfaceOutput& a_rhs)
    : paraboloid_m(a_rhs.paraboloid_m),
      arc_list_m(a_rhs.arc_list_m),
      knows_surface_area_m(a_rhs.knows_surface_area_m),
      surface_area_m(a_rhs.surface_area_m),
      knows_avg_normal_m(a_rhs.knows_avg_normal_m),
      avg_normal_m(a_rhs.avg_normal_m),
      knows_int_mean_curv_m(a_rhs.knows_int_mean_curv_m),
      int_mean_curv_m(a_rhs.int_mean_curv_m),
      knows_int_gaussian_curv_m(a_rhs.knows_int_gaussian_curv_m),
      int_gaussian_curv_m(a_rhs.int_gaussian_curv_m),
      length_scale_m{a_rhs.length_scale_m} {}

inline ParametrizedSurfaceOutput& ParametrizedSurfaceOutput::operator=(
    ParametrizedSurfaceOutput&& a_rhs) {
  if (this != &a_rhs) {
    paraboloid_m = a_rhs.paraboloid_m;
    arc_list_m = std::move(a_rhs.arc_list_m);
    knows_surface_area_m = a_rhs.knows_surface_area_m;
    surface_area_m = a_rhs.surface_area_m;
    knows_avg_normal_m = a_rhs.knows_avg_normal_m;
    avg_normal_m = a_rhs.avg_normal_m;
    knows_int_mean_curv_m = a_rhs.knows_int_mean_curv_m;
    int_mean_curv_m = a_rhs.int_mean_curv_m;
    knows_int_gaussian_curv_m = a_rhs.knows_int_gaussian_curv_m;
    int_gaussian_curv_m = a_rhs.int_gaussian_curv_m;
    length_scale_m = a_rhs.length_scale_m;
    pt_from_bezier_split_m = std::move(a_rhs.pt_from_bezier_split_m);
  }
  return *this;
}

inline ParametrizedSurfaceOutput& ParametrizedSurfaceOutput::operator=(
    const ParametrizedSurfaceOutput& a_rhs) {
  if (this != &a_rhs) {
    paraboloid_m = a_rhs.paraboloid_m;
    arc_list_m = a_rhs.arc_list_m;
    knows_surface_area_m = a_rhs.knows_surface_area_m;
    surface_area_m = a_rhs.surface_area_m;
    knows_avg_normal_m = a_rhs.knows_avg_normal_m;
    avg_normal_m = a_rhs.avg_normal_m;
    knows_int_mean_curv_m = a_rhs.knows_int_mean_curv_m;
    int_mean_curv_m = a_rhs.int_mean_curv_m;
    knows_int_gaussian_curv_m = a_rhs.knows_int_gaussian_curv_m;
    int_gaussian_curv_m = a_rhs.int_gaussian_curv_m;
    length_scale_m = a_rhs.length_scale_m;
  }
  return *this;
}

inline void ParametrizedSurfaceOutput::setLengthScale(
    const double a_length_scale) {
  length_scale_m = a_length_scale;
}

inline void ParametrizedSurfaceOutput::setParaboloid(
    const Paraboloid& a_paraboloid) {
  paraboloid_m = a_paraboloid;
}

inline RationalBezierArc& ParametrizedSurfaceOutput::operator[](
    const UnsignedIndex_t a_index) {
  return arc_list_m[a_index];
}

inline const RationalBezierArc& ParametrizedSurfaceOutput::operator[](
    const UnsignedIndex_t a_index) const {
  return arc_list_m[a_index];
}

inline const Paraboloid& ParametrizedSurfaceOutput::getParaboloid(void) const {
  return paraboloid_m;
}

inline std::vector<RationalBezierArc>& ParametrizedSurfaceOutput::getArcs(
    void) {
  return arc_list_m;
}

inline std::vector<Pt*>& ParametrizedSurfaceOutput::getPts(void) {
  return pt_from_bezier_split_m;
}

inline void ParametrizedSurfaceOutput::addArc(
    const RationalBezierArc& a_rational_bezier_arc) {
  arc_list_m.push_back(a_rational_bezier_arc);
}

inline void ParametrizedSurfaceOutput::addPt(Pt* a_pt) {
  pt_from_bezier_split_m.push_back(a_pt);
}

inline const std::vector<RationalBezierArc>::size_type
ParametrizedSurfaceOutput::size(void) const {
  return arc_list_m.size();
}

inline void ParametrizedSurfaceOutput::clearArcs(void) { arc_list_m.clear(); }

inline void ParametrizedSurfaceOutput::clearPts(void) {
  for (auto& elem : pt_from_bezier_split_m) {
    delete elem;
  }
  pt_from_bezier_split_m.clear();
}

inline void ParametrizedSurfaceOutput::clear(void) {
  this->clearArcs();
  this->clearPts();
}

inline ParametrizedSurfaceOutput::~ParametrizedSurfaceOutput(void) {
  for (auto elem : pt_from_bezier_split_m) {
    delete elem;
  }
}

class ArcContributionToSurfaceArea_Functor {
 public:
  ArcContributionToSurfaceArea_Functor(const RationalBezierArc& a_arc,
                                       const AlignedParaboloid& a_paraboloid)
      : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      if (std::fabs(a) < 10.0 * DBL_EPSILON &&
          std::fabs(b) < 10.0 * DBL_EPSILON) {
        return 0.0;
      } else if (std::fabs(a) > std::fabs(b)) {
        const double primitive0 =
            (2. * pt0[0] *
                 std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0]) +
                           4. * (b * b) * (pt0[1] * pt0[1])) -
             (1. + 4. * (b * b) * (pt0[1] * pt0[1])) *
                 std::log(-2. * std::fabs(a) * pt0[0] +
                          std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0]) +
                                    4. * (b * b) * (pt0[1] * pt0[1]))) /
                 std::fabs(a)) /
            4.;
        const double primitive1 =
            (2. * pt1[0] *
                 std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0]) +
                           4. * (b * b) * (pt1[1] * pt1[1])) -
             (1. + 4. * (b * b) * (pt1[1] * pt1[1])) *
                 std::log(-2. * std::fabs(a) * pt1[0] +
                          std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0]) +
                                    4. * (b * b) * (pt1[1] * pt1[1]))) /
                 std::fabs(a)) /
            4.;
        return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);
      } else {
        const double primitive0 =
            -(2. * pt0[1] *
                  std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0]) +
                            4. * (b * b) * (pt0[1] * pt0[1])) -
              (1. + 4. * (a * a) * (pt0[0] * pt0[0])) *
                  std::log(-2. * std::fabs(b) * pt0[1] +
                           std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0]) +
                                     4. * (b * b) * (pt0[1] * pt0[1]))) /
                  std::fabs(b)) /
            (4.);
        const double primitive1 =
            -(2. * pt1[1] *
                  std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0]) +
                            4. * (b * b) * (pt1[1] * pt1[1])) -
              (1. + 4. * (a * a) * (pt1[0] * pt1[0])) *
                  std::log(-2. * std::fabs(b) * pt1[1] +
                           std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0]) +
                                     4. * (b * b) * (pt1[1] * pt1[1]))) /
                  std::fabs(b)) /
            (4.);
        return 0.5 * (primitive0 * der0[0] + primitive1 * der1[0]);
      }
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      if (std::fabs(a) < 10.0 * DBL_EPSILON &&
          std::fabs(b) < 10.0 * DBL_EPSILON) {
        return pt[0] * der[1];
      } else if (std::fabs(a) > std::fabs(b)) {
        const double primitive =
            (2. * pt[0] *
                 std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0]) +
                           4. * (b * b) * (pt[1] * pt[1])) -
             (1. + 4. * (b * b) * (pt[1] * pt[1])) *
                 std::log(-2. * std::fabs(a) * pt[0] +
                          std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0]) +
                                    4. * (b * b) * (pt[1] * pt[1]))) /
                 std::fabs(a)) /
            4.;
        if (std::isnan(primitive)) {
          std::cout << "Pr = " << pt << std::endl;
          std::cout << "Der = " << der << std::endl;
          std::cout << "a = " << a << std::endl;
          std::cout << "b = " << b << std::endl;
          std::cout << "Arc: weight = " << arc_m.weight() << std::endl;
          std::cout << "Arc: start = " << arc_m.start_point() << std::endl;
          std::cout << "Arc: ctrl  = " << arc_m.control_point() << std::endl;
          std::cout << "Arc: end   = " << arc_m.end_point() << std::endl;
          std::cout << "Primitive is NaN" << std::endl;
          exit(1);
        }
        if (std::isnan(der[1])) {
          std::cout << "der[1] is NaN" << std::endl;
          exit(1);
        }

        return primitive * der[1];
      } else {
        const double primitive =
            -(2. * pt[1] *
                  std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0]) +
                            4. * (b * b) * (pt[1] * pt[1])) -
              (1. + 4. * (a * a) * (pt[0] * pt[0])) *
                  std::log(-2. * std::fabs(b) * pt[1] +
                           std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0]) +
                                     4. * (b * b) * (pt[1] * pt[1]))) /
                  std::fabs(b)) /
            (4.);
        if (std::isnan(primitive)) {
          std::cout << "Pr = " << pt << std::endl;
          std::cout << "Der = " << der << std::endl;
          std::cout << "a = " << a << std::endl;
          std::cout << "b = " << b << std::endl;
          std::cout << "Arc: weight = " << arc_m.weight() << std::endl;
          std::cout << "Arc: start = " << arc_m.start_point() << std::endl;
          std::cout << "Arc: ctrl  = " << arc_m.control_point() << std::endl;
          std::cout << "Arc: end   = " << arc_m.end_point() << std::endl;
          std::cout << "Primitive is NaN" << std::endl;
          exit(1);
        }
        if (std::isnan(der[0])) {
          std::cout << "der[0] is NaN" << std::endl;
          exit(1);
        }

        return primitive * der[0];
      }
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedParaboloid& paraboloid_m;
};

class ArcContributionToNormalX_Functor {
 public:
  ArcContributionToNormalX_Functor(const RationalBezierArc& a_arc,
                                   const AlignedParaboloid& a_paraboloid)
      : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      const double primitive0 = a * pt0[0] * pt0[0];
      const double primitive1 = a * pt1[0] * pt1[0];
      return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      const double primitive = a * pt[0] * pt[0];
      return primitive * der[1];
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedParaboloid& paraboloid_m;
};

class ArcContributionToNormalY_Functor {
 public:
  ArcContributionToNormalY_Functor(const RationalBezierArc& a_arc,
                                   const AlignedParaboloid& a_paraboloid)
      : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      const double primitive0 = -b * pt0[1] * pt0[1];
      const double primitive1 = -b * pt1[1] * pt1[1];
      return 0.5 * (primitive0 * der0[0] + primitive1 * der1[0]);
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      const double primitive = -b * pt[1] * pt[1];
      return primitive * der[0];
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedParaboloid& paraboloid_m;
};

class ArcContributionToNormalZ_Functor {
 public:
  ArcContributionToNormalZ_Functor(const RationalBezierArc& a_arc,
                                   const AlignedParaboloid& a_paraboloid)
      : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double primitive0 = pt0[0];
      const double primitive1 = pt1[0];
      return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double primitive = pt[0];
      return primitive * der[1];
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedParaboloid& paraboloid_m;
};

class ArcContributionToMeanCurvature_Functor {
 public:
  ArcContributionToMeanCurvature_Functor(const RationalBezierArc& a_arc,
                                         const AlignedParaboloid& a_paraboloid)
      : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      if (std::fabs(a) < 10.0 * DBL_EPSILON &&
          std::fabs(b) < 10.0 * DBL_EPSILON) {
        return 0.0;
      } else if (std::fabs(a) > std::fabs(b)) {
        const double primitive0 =
            2. * b * pt0[0] +
            ((a + 4. * a * (b * b) * (pt0[1] * pt0[1]) -
              4. * (b * b * b) * (pt0[1] * pt0[1])) *
             std::atan((2. * a * pt0[0]) /
                       std::sqrt(1. + 4. * (b * b) * (pt0[1] * pt0[1])))) /
                (a * std::sqrt(1. + 4. * (b * b) * (pt0[1] * pt0[1])));
        const double primitive1 =
            2. * b * pt1[0] +
            ((a + 4. * a * (b * b) * (pt1[1] * pt1[1]) -
              4. * (b * b * b) * (pt1[1] * pt1[1])) *
             std::atan((2. * a * pt1[0]) /
                       std::sqrt(1. + 4. * (b * b) * (pt1[1] * pt1[1])))) /
                (a * std::sqrt(1. + 4. * (b * b) * (pt1[1] * pt1[1])));
        return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);
      } else {
        const double primitive0 =
            -2. * a * pt0[1] -
            ((b - 4. * (a * a * a) * (pt0[0] * pt0[0]) +
              4. * (a * a) * b * (pt0[0] * pt0[0])) *
             std::atan((2. * b * pt0[1]) /
                       std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0])))) /
                (b * std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0])));
        const double primitive1 =
            -2. * a * pt1[1] -
            ((b - 4. * (a * a * a) * (pt1[0] * pt1[0]) +
              4. * (a * a) * b * (pt1[0] * pt1[0])) *
             std::atan((2. * b * pt1[1]) /
                       std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0])))) /
                (b * std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0])));
        return 0.5 * (primitive0 * der0[0] + primitive1 * der1[0]);
      }
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      if (std::fabs(a) < 10.0 * DBL_EPSILON &&
          std::fabs(b) < 10.0 * DBL_EPSILON) {
        return 0.0;
      } else if (std::fabs(a) > std::fabs(b)) {
        const double primitive =
            2. * b * pt[0] +
            ((a + 4. * a * (b * b) * (pt[1] * pt[1]) -
              4. * (b * b * b) * (pt[1] * pt[1])) *
             std::atan((2. * a * pt[0]) /
                       std::sqrt(1. + 4. * (b * b) * (pt[1] * pt[1])))) /
                (a * std::sqrt(1. + 4. * (b * b) * (pt[1] * pt[1])));
        return primitive * der[1];
      } else {
        const double primitive =
            -2. * a * pt[1] -
            ((b - 4. * (a * a * a) * (pt[0] * pt[0]) +
              4. * (a * a) * b * (pt[0] * pt[0])) *
             std::atan((2. * b * pt[1]) /
                       std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0])))) /
                (b * std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0])));
        return primitive * der[0];
      }
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedParaboloid& paraboloid_m;
};  // namespace IRL

class ArcContributionToGaussianCurvature_Functor {
 public:
  ArcContributionToGaussianCurvature_Functor(
      const RationalBezierArc& a_arc, const AlignedParaboloid& a_paraboloid)
      : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      if (std::fabs(a) < 10.0 * DBL_EPSILON &&
          std::fabs(b) < 10.0 * DBL_EPSILON) {
        return 0.0;
      } else {
        const double primitive0 =
            4.0 * a * b * pt0[0] /
            ((1.0 + 4.0 * b * b * pt0[1] * pt0[1]) *
             std::sqrt(1.0 + 4.0 * a * a * pt0[0] * pt0[0] +
                       4.0 * b * b * pt0[1] * pt0[1]));
        const double primitive1 =
            4.0 * a * b * pt1[0] /
            ((1.0 + 4.0 * b * b * pt1[1] * pt1[1]) *
             std::sqrt(1.0 + 4.0 * a * a * pt1[0] * pt1[0] +
                       4.0 * b * b * pt1[1] * pt1[1]));
        return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);
      }
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      if (std::fabs(a) < 10.0 * DBL_EPSILON &&
          std::fabs(b) < 10.0 * DBL_EPSILON) {
        return 0.0;
      } else {
        const double primitive = 4.0 * a * b * pt[0] /
                                 ((1.0 + 4.0 * b * b * pt[1] * pt[1]) *
                                  std::sqrt(1.0 + 4.0 * a * a * pt[0] * pt[0] +
                                            4.0 * b * b * pt[1] * pt[1]));
        return primitive * der[1];
      }
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedParaboloid& paraboloid_m;
};  // namespace IRL

inline double ParametrizedSurfaceOutput::getSurfaceArea(void) {
  if (!knows_surface_area_m) {
    const UnsignedIndex_t nArcs = this->size();
    surface_area_m = 0.0;
    size_t limit = 128;

    const double epsabs = 10.0 * DBL_EPSILON;
    const double epsrel = 0.0;
    auto& aligned_paraboloid = paraboloid_m.getAlignedParaboloid();
    for (std::size_t t = 0; t < nArcs; ++t) {
      // Define the functor
      ArcContributionToSurfaceArea_Functor functor(arc_list_m[t],
                                                   aligned_paraboloid);

      // Define the integrator.
      Eigen::Integrator<double> integrator(limit);

      // Define a quadrature rule.
      Eigen::Integrator<double>::QuadratureRule quadrature_rule =
          Eigen::Integrator<double>::GaussKronrod61;

      // Integrate.
      surface_area_m += integrator.quadratureAdaptive(functor, 0.0, 1.0, epsabs,
                                                      epsrel, quadrature_rule);
    }
    knows_surface_area_m = true;
  }
  return surface_area_m;
}

inline Normal ParametrizedSurfaceOutput::getAverageNormal(void) {
  if (!knows_avg_normal_m) {
    const UnsignedIndex_t nArcs = this->size();
    avg_normal_m = Normal();
    size_t limit = 128;

    const double epsabs = 10.0 * DBL_EPSILON;
    const double epsrel = 0.0;
    auto& aligned_paraboloid = paraboloid_m.getAlignedParaboloid();
    for (std::size_t t = 0; t < nArcs; ++t) {
      // Define the functor
      ArcContributionToNormalX_Functor functorx(arc_list_m[t],
                                                aligned_paraboloid);
      ArcContributionToNormalY_Functor functory(arc_list_m[t],
                                                aligned_paraboloid);
      ArcContributionToNormalZ_Functor functorz(arc_list_m[t],
                                                aligned_paraboloid);

      // Define the integrator.
      Eigen::Integrator<double> integrator(limit);

      // Define a quadrature rule.
      Eigen::Integrator<double>::QuadratureRule quadrature_rule =
          Eigen::Integrator<double>::GaussKronrod61;

      // Integrate.
      avg_normal_m[0] += integrator.quadratureAdaptive(
          functorx, 0.0, 1.0, epsabs, epsrel, quadrature_rule);
      avg_normal_m[1] += integrator.quadratureAdaptive(
          functory, 0.0, 1.0, epsabs, epsrel, quadrature_rule);
      avg_normal_m[2] += integrator.quadratureAdaptive(
          functorz, 0.0, 1.0, epsabs, epsrel, quadrature_rule);
    }
    avg_normal_m.normalize();
    knows_avg_normal_m = true;
  }
  return avg_normal_m;
}

inline Normal ParametrizedSurfaceOutput::getAverageNormalNonAligned(void) {
  auto aligned_normal = this->getAverageNormal();
  const auto& ref_frame = this->getParaboloid().getReferenceFrame();
  auto normal = Normal();
  for (std::size_t d = 0; d < 3; ++d) {
    for (std::size_t n = 0; n < 3; ++n) {
      normal[n] += ref_frame[d][n] * aligned_normal[d];
    }
  }
  return normal;
}

inline double ParametrizedSurfaceOutput::getMeanCurvatureIntegral(void) {
  if (!knows_int_mean_curv_m) {
    const UnsignedIndex_t nArcs = this->size();
    int_mean_curv_m = 0.0;
    size_t limit = 128;

    const double epsabs = 10.0 * DBL_EPSILON;
    const double epsrel = 0.0;
    auto& aligned_paraboloid = paraboloid_m.getAlignedParaboloid();
    for (std::size_t t = 0; t < nArcs; ++t) {
      // Define the functor
      ArcContributionToMeanCurvature_Functor functor(arc_list_m[t],
                                                     aligned_paraboloid);

      // Define the integrator.
      Eigen::Integrator<double> integrator(limit);

      // Define a quadrature rule.
      Eigen::Integrator<double>::QuadratureRule quadrature_rule =
          Eigen::Integrator<double>::GaussKronrod61;

      // Integrate.
      int_mean_curv_m += integrator.quadratureAdaptive(
          functor, 0.0, 1.0, epsabs, epsrel, quadrature_rule);
    }
    knows_int_mean_curv_m = true;
  }
  return int_mean_curv_m;
}

inline double ParametrizedSurfaceOutput::getAverageMeanCurvature(void) {
  return this->getMeanCurvatureIntegral() / safelyTiny(this->getSurfaceArea());
}

inline double ParametrizedSurfaceOutput::getGaussianCurvatureIntegral(void) {
  if (!knows_int_gaussian_curv_m) {
    const UnsignedIndex_t nArcs = this->size();
    int_gaussian_curv_m = 0.0;
    size_t limit = 128;

    const double epsabs = 10.0 * DBL_EPSILON;
    const double epsrel = 0.0;
    auto& aligned_paraboloid = paraboloid_m.getAlignedParaboloid();
    for (std::size_t t = 0; t < nArcs; ++t) {
      // Define the functor
      ArcContributionToGaussianCurvature_Functor functor(arc_list_m[t],
                                                         aligned_paraboloid);

      // Define the integrator.
      Eigen::Integrator<double> integrator(limit);

      // Define a quadrature rule.
      Eigen::Integrator<double>::QuadratureRule quadrature_rule =
          Eigen::Integrator<double>::GaussKronrod61;

      // Integrate.
      int_gaussian_curv_m += integrator.quadratureAdaptive(
          functor, 0.0, 1.0, epsabs, epsrel, quadrature_rule);
    }
    knows_int_gaussian_curv_m = true;
  }
  return int_gaussian_curv_m;
}

inline double ParametrizedSurfaceOutput::getAverageGaussianCurvature(void) {
  return this->getGaussianCurvatureIntegral() /
         safelyTiny(this->getSurfaceArea());
}

inline Normal ParametrizedSurfaceOutput::getNormalAligned(const Pt a_pt) {
  auto& aligned_paraboloid = this->getParaboloid().getAlignedParaboloid();
  auto aligned_normal = getParaboloidSurfaceNormal(aligned_paraboloid, a_pt);
  aligned_normal.normalize();
  return aligned_normal;
}

inline Normal ParametrizedSurfaceOutput::getNormalNonAligned(const Pt a_pt) {
  const auto& datum = this->getParaboloid().getDatum();
  const auto& ref_frame = this->getParaboloid().getReferenceFrame();
  // assert(ref_frame.isOrthonormalBasis());
  const Pt original_pt = a_pt - datum;
  auto aligned_pt = a_pt;
  for (std::size_t n = 0; n < 3; ++n) {
    aligned_pt[n] = ref_frame[n] * original_pt;
  }
  auto aligned_normal = this->getNormalAligned(aligned_pt);
  auto normal = Normal();
  for (std::size_t d = 0; d < 3; ++d) {
    for (std::size_t n = 0; n < 3; ++n) {
      normal[n] += ref_frame[d][n] * aligned_normal[d];
    }
  }
  return normal;
}

inline double ParametrizedSurfaceOutput::getMeanCurvatureAligned(
    const Pt a_pt) {
  auto& aligned_paraboloid = this->getParaboloid().getAlignedParaboloid();
  return (2. * (aligned_paraboloid.a() + aligned_paraboloid.b() +
                4. * (aligned_paraboloid.a() * aligned_paraboloid.a()) *
                    aligned_paraboloid.b() * (a_pt[0] * a_pt[0]) +
                4. * aligned_paraboloid.a() *
                    (aligned_paraboloid.b() * aligned_paraboloid.b()) *
                    (a_pt[1] * a_pt[1]))) /
         std::pow(1. +
                      4. * (aligned_paraboloid.a() * aligned_paraboloid.a()) *
                          (a_pt[0] * a_pt[0]) +
                      4. * (aligned_paraboloid.b() * aligned_paraboloid.b()) *
                          (a_pt[1] * a_pt[1]),
                  1.5);
}

inline double ParametrizedSurfaceOutput::getMeanCurvatureNonAligned(
    const Pt a_pt) {
  const auto& datum = this->getParaboloid().getDatum();
  const auto& ref_frame = this->getParaboloid().getReferenceFrame();
  // assert(ref_frame.isOrthonormalBasis());
  const Pt original_pt = a_pt - datum;
  auto aligned_pt = a_pt;
  for (std::size_t n = 0; n < 3; ++n) {
    aligned_pt[n] = ref_frame[n] * original_pt;
  }
  return this->getMeanCurvatureAligned(aligned_pt);
}

inline double ParametrizedSurfaceOutput::getGaussianCurvatureAligned(
    const Pt a_pt) {
  auto& aligned_paraboloid = this->getParaboloid().getAlignedParaboloid();
  return 4. * aligned_paraboloid.a() * aligned_paraboloid.b() /
         ((1. +
           4. * (aligned_paraboloid.a() * aligned_paraboloid.a()) *
               (a_pt[0] * a_pt[0]) +
           4. * (aligned_paraboloid.b() * aligned_paraboloid.b()) *
               (a_pt[1] * a_pt[1])) *
          (1. +
           4. * (aligned_paraboloid.a() * aligned_paraboloid.a()) *
               (a_pt[0] * a_pt[0]) +
           4. * (aligned_paraboloid.b() * aligned_paraboloid.b()) *
               (a_pt[1] * a_pt[1])));
}

inline double ParametrizedSurfaceOutput::getGaussianCurvatureNonAligned(
    const Pt a_pt) {
  const auto& datum = this->getParaboloid().getDatum();
  const auto& ref_frame = this->getParaboloid().getReferenceFrame();
  // assert(ref_frame.isOrthonormalBasis());
  const Pt original_pt = a_pt - datum;
  auto aligned_pt = a_pt;
  for (std::size_t n = 0; n < 3; ++n) {
    aligned_pt[n] = ref_frame[n] * original_pt;
  }
  return this->getGaussianCurvatureAligned(aligned_pt);
}

inline TriangulatedSurfaceOutput ParametrizedSurfaceOutput::triangulate(
    const double a_length_scale, const UnsignedIndex_t a_nsplit) const {
  TriangulatedSurfaceOutput returned_surface;
  this->triangulate_fromPtr(a_length_scale, a_nsplit, &returned_surface);
  return returned_surface;
}

inline void ParametrizedSurfaceOutput::triangulate_fromPtr(
    const double a_length_scale, const UnsignedIndex_t a_nsplit,
    TriangulatedSurfaceOutput* returned_surface) const {
  const UnsignedIndex_t nArcs = this->size();
  double length_scale, length_scale_ref = length_scale_m;
  if (a_length_scale > 0.0) {
    length_scale_ref = a_length_scale;
  } else if (length_scale_ref <= 0.0) {
    auto surf = (*this);
    const double avg_length = std::sqrt(std::abs(surf.getSurfaceArea())) / 3.0;
    const double curv = std::fabs(surf.getAverageMeanCurvature());
    length_scale_ref = std::min(0.1 / curv, avg_length);
  }
  const auto& aligned_paraboloid = paraboloid_m.getAlignedParaboloid();

  std::vector<std::vector<RationalBezierArc>> list_of_closed_curves;
  std::vector<bool> visited(nArcs, false);

  // First, we need to order the arcs so as to form closed
  // curves
  double min_arc_length = DBL_MAX;
  bool valid_curves = true;
  for (std::size_t t = 0; t < nArcs; ++t) {
    if (visited[t]) {
      continue;
    }
    visited[t] = true;
    // Start with next available arc
    if (arc_list_m[t].weight() > 1.0e15) {
      const Pt p0 = arc_list_m[t].start_point();
      const Pt p1 = arc_list_m[t].control_point();
      const Pt p2 = arc_list_m[t].end_point();
      list_of_closed_curves.push_back(std::vector<RationalBezierArc>(
          {RationalBezierArc(p0, 0.5 * (p0 + p1), p1, 0.0),
           RationalBezierArc(p1, 0.5 * (p1 + p2), p2, 0.0)}));
    } else {
      list_of_closed_curves.push_back(
          std::vector<RationalBezierArc>({arc_list_m[t]}));
    }
    const std::uintptr_t start_id = arc_list_m[t].start_point_id();
    std::uintptr_t end_id = arc_list_m[t].end_point_id();
    int counter = 0;
    while (end_id != start_id) {
      for (std::size_t e = t + 1; e < nArcs; ++e) {
        if (arc_list_m[e].start_point_id() == end_id) {
          visited[e] = true;
          if (arc_list_m[e].weight() > 1.0e15) {
            const Pt p0 = arc_list_m[e].start_point();
            const Pt p1 = arc_list_m[e].control_point();
            const Pt p2 = arc_list_m[e].end_point();
            list_of_closed_curves.back().push_back(
                RationalBezierArc(p0, 0.5 * (p0 + p1), p1, 0.0));
            list_of_closed_curves.back().push_back(
                RationalBezierArc(p1, 0.5 * (p1 + p2), p2, 0.0));
          } else {
            list_of_closed_curves.back().push_back(arc_list_m[e]);
          }
          end_id = arc_list_m[e].end_point_id();
          break;
        }
      }
      if (++counter > nArcs) {
        valid_curves = false;
        break;
      }
    }
  }

  returned_surface->clearAll();

  if (valid_curves) {
#ifdef IRL_USE_EARCUT
    // The number type to use for tessellation
    using Coord = double;
    // The index type. Defaults to uint32_t, but you can also
    // pass uint16_t if you know that your data won't have
    // more than 65536 vertices.

    length_scale = DBL_MAX;
    const UnsignedIndex_t nCurves =
        static_cast<UnsignedIndex_t>(list_of_closed_curves.size());

    // Create array
    using Point = std::array<Coord, 2>;
    std::vector<std::vector<Point>> polygon;
    polygon.resize(nCurves);

    for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
      const UnsignedIndex_t nLocalArcs = list_of_closed_curves[i].size();
      // Loop over arcs of curve
      UnsignedIndex_t added_points = 0;
      double signed_area = 0.0;
      for (UnsignedIndex_t j = 0; j < nLocalArcs; ++j) {
        // Compute approximate arc length
        const RationalBezierArc& arc = list_of_closed_curves[i][j];
        const double arc_length = arc.arc_length();
        // Split arc
        UnsignedIndex_t nSplit = a_nsplit <= 0 ? 1 : a_nsplit;
        if (length_scale_ref > 0.0) {
          nSplit = static_cast<UnsignedIndex_t>(arc_length / length_scale_ref);
          nSplit = nSplit < a_nsplit ? a_nsplit : nSplit;
        }
        const double step = 1.0 / static_cast<double>(nSplit);
        length_scale = std::min(length_scale, step * arc_length);
        if (length_scale_ref > 0.0) length_scale = length_scale_ref;
        Pt previous_pt = arc.point(0.0);
        for (UnsignedIndex_t k = 1; k <= nSplit; ++k) {
          const double t = static_cast<double>(k) * step;
          const auto pt = arc.point(t);
          polygon[i].push_back({pt[0], pt[1]});
          signed_area +=
              0.5 * (previous_pt[0] * pt[1] - pt[0] * previous_pt[1]);
          previous_pt = pt;
        }
      }
    }

    if (a_length_scale > 0.0) {
      length_scale = a_length_scale;
    }
    // Run tessellation
    // Returns array of indices that refer to the vertices of
    // the input polygon. e.g: the index 6 would refer to {25,
    // 75} in this example. Three subsequent indices form a
    // triangle. Output triangles are clockwise.
    std::vector<int> indices = mapbox::earcut<int>(polygon);

    auto& vlist = returned_surface->getVertexList();
    auto& tlist = returned_surface->getTriangleList();
    std::vector<std::array<int, 4>> elist;

    UnsignedIndex_t count = 0;
    for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
      count += polygon[i].size();
    }
    vlist.resize(count);

    count = 0;
    for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
      for (UnsignedIndex_t j = 0; j < polygon[i].size(); ++j) {
        double x = polygon[i][j][0];
        double y = polygon[i][j][1];
        double z =
            -aligned_paraboloid.a() * x * x - aligned_paraboloid.b() * y * y;
        vlist[count++] = Pt(x, y, z);
      }
    }

    reMeshPolygon(vlist, elist, indices, length_scale, aligned_paraboloid);

    count = 0;
    for (UnsignedIndex_t i = 0; i < indices.size() / 3; ++i) {
      if (indices[3 * i + 0] >= 0 && indices[3 * i + 1] >= 0 &&
          indices[3 * i + 2] >= 0) {
        count++;
      }
    }
    // assert(count == indices.size() / 3);
    tlist.resize(count, TriangulatedSurfaceOutput::TriangleStorage::value_type::
                            fromNoExistencePlane(vlist, {0, 0, 0}));
    count = 0;
    for (UnsignedIndex_t i = 0; i < indices.size() / 3; ++i) {
      if (indices[3 * i + 0] >= 0 && indices[3 * i + 1] >= 0 &&
          indices[3 * i + 2] >= 0) {
        assert(indices[3 * i + 0] != indices[3 * i + 1]);
        assert(indices[3 * i + 1] != indices[3 * i + 2]);
        assert(indices[3 * i + 2] != indices[3 * i + 0]);
        tlist[count++] = TriangulatedSurfaceOutput::TriangleStorage::
            value_type::fromNoExistencePlane(
                vlist, {static_cast<UnsignedIndex_t>(indices[3 * i]),
                        static_cast<UnsignedIndex_t>(indices[3 * i + 1]),
                        static_cast<UnsignedIndex_t>(indices[3 * i + 2])});
      }
    }

    for (UnsignedIndex_t i = 0; i < elist.size(); ++i) {
      if ((elist[i][0] >= 0 || elist[i][1] >= 0) &&
          (elist[i][2] < 0 || elist[i][3] < 0)) {
        returned_surface->addBoundaryEdge(elist[i][0], elist[i][1]);
      }
    }

    // Translate and rotate triangulated surface vertices
    const auto& datum = paraboloid_m.getDatum();
    const auto& ref_frame = paraboloid_m.getReferenceFrame();
    for (auto& vertex : vlist) {
      const Pt base_pt = vertex;
      vertex = Pt(0.0, 0.0, 0.0);
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        for (UnsignedIndex_t n = 0; n < 3; ++n) {
          vertex[n] += ref_frame[d][n] * base_pt[d];
        }
      }
      vertex += datum;
    }

#elif defined IRL_USE_CGAL
    typedef CGAL::Exact_predicates_exact_constructions_kernel K;
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kexact;
    typedef CGAL::Delaunay_mesh_vertex_base_2<K> Vb;
    typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
    typedef CGAL::Exact_predicates_tag Itag;
    typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
    typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
    typedef CDT::Vertex_handle Vertex_handle;
    typedef CDT::Point Point;
    typedef CGAL::Arr_segment_traits_2<Kexact> Traits_2;
    typedef Traits_2::Curve_2 SegmentExact;
    typedef Kexact::Point_2 PointExact;

    CDT cdt;

    std::ofstream myfile;
    myfile.open("triangulation_log.txt");
    myfile << "Starting triangulating surface.\n";
    myfile << std::setprecision(16) << std::scientific
           << "Paraboloid: " << aligned_paraboloid << "\n";

    // Create boundaries
    std::vector<Point> points;
    std::list<Point> list_of_seeds;
    const UnsignedIndex_t nCurves =
        static_cast<UnsignedIndex_t>(list_of_closed_curves.size());
    UnsignedIndex_t start_points = 0;
    double total_signed_area = 0.0;
    double xmin = DBL_MAX, xmax = -DBL_MAX;
    double ymin = DBL_MAX, ymax = -DBL_MAX;
    UnsignedIndex_t vertex_count = 0;
    bool previous_valid = false;
    for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
      points.resize(0);
      const UnsignedIndex_t nLocalArcs = list_of_closed_curves[i].size();
      // Loop over arcs of curve
      UnsignedIndex_t added_points = 0;
      double signed_area = 0.0;
      for (UnsignedIndex_t j = 0; j < nLocalArcs; ++j) {
        // Compute approximate arc length
        const RationalBezierArc& arc = list_of_closed_curves[i][j];
        const double arc_length = arc.arc_length();
        myfile << std::setprecision(16) << std::scientific << "Curve " << i
               << " has arc: " << arc << "\n";
        // Split arc
        UnsignedIndex_t nSplit = a_nsplit <= 0 ? 1 : a_nsplit;
        if (length_scale_ref > 0.0) {
          nSplit = static_cast<UnsignedIndex_t>(arc_length / length_scale_ref);
          nSplit = nSplit < a_nsplit ? a_nsplit : nSplit;
        }
        const double step = 1.0 / static_cast<double>(nSplit);
        length_scale = std::min(length_scale, step * arc_length);
        if (length_scale_ref > 0.0) length_scale = length_scale_ref;
        Pt previous_pt = arc.point(0.0);
        for (UnsignedIndex_t k = 1; k <= nSplit; ++k) {
          const double t = static_cast<double>(k) * step;
          const auto pt = arc.point(t);
          myfile << std::setprecision(16) << std::scientific << "Adding vertex "
                 << vertex_count++ << " at " << pt[0] << ", " << pt[1] << ".\n";
          points.push_back(Point(pt[0], pt[1]));
          signed_area +=
              0.5 * (previous_pt[0] * pt[1] - pt[0] * previous_pt[1]);
          previous_pt = pt;
        }
      }

      /* Remove duplicates */
      UnsignedIndex_t id0 = 0;
      do {
        UnsignedIndex_t id1 = (id0 + 1) % points.size();
        if ((points[id1].x() - points[id0].x()) *
                    (points[id1].x() - points[id0].x()) +
                (points[id1].y() - points[id0].y()) *
                    (points[id1].y() - points[id0].y()) <
            1.0e12 * DBL_EPSILON * DBL_EPSILON) {
          myfile << std::setprecision(16) << std::scientific
                 << "Removing duplicate " << id1 << " at " << points[id1].x()
                 << ", " << points[id1].y() << " too close to " << id0 << " at "
                 << points[id0].x() << ", " << points[id0].y() << ".\n";
          points.erase(points.begin() + id1);
          continue;
        } else {
          id0 = id1;
        }
      } while (id0 != 0);

      /* Create constraints */
      if (points.size() >= 3 &&
          std::fabs(signed_area) >
              std::max(1.0e-4 * length_scale * length_scale, 1.0e-14)) {
        // Construct the input segments.
        std::vector<SegmentExact> segments;
        segments.resize(points.size());
        segments[0] = SegmentExact(PointExact(points[points.size() - 1].x(),
                                              points[points.size() - 1].y()),
                                   PointExact(points[0].x(), points[0].y()));
        for (UnsignedIndex_t j = 0; j < points.size() - 1; ++j) {
          segments[j + 1] =
              SegmentExact(PointExact(points[j].x(), points[j].y()),
                           PointExact(points[j + 1].x(), points[j + 1].y()));
        }

        if (!CGAL::do_curves_intersect(segments.begin(), segments.end())) {
          if (nCurves > 1 && signed_area < 0.0) {
            // Add hole
            const auto p1x = CGAL::to_double(points[0].x());
            const auto p1y = CGAL::to_double(points[0].y());
            const auto p2x = CGAL::to_double(points[1].x());
            const auto p2y = CGAL::to_double(points[1].y());
            std::array<double, 2> hole_location{
                {0.5 * (p1x + p2x), 0.5 * (p1y + p2y)}};
            Normal shift_dir = Normal(p2y - p1y, p1x - p2x, 0.0);
            shift_dir.normalize();
            myfile << std::setprecision(16) << std::scientific << "Adding hole "
                   << hole_location[0] + (1.0e3 * DBL_EPSILON) * shift_dir[0]
                   << ", "
                   << hole_location[1] + (1.0e3 * DBL_EPSILON) * shift_dir[1]
                   << ".\n";
            list_of_seeds.push_back(
                Point(hole_location[0] + (1.0e3 * DBL_EPSILON) * shift_dir[0],
                      hole_location[1] + (1.0e3 * DBL_EPSILON) * shift_dir[1]));
          }

          // Create segments
          myfile << "Adding constraint " << points.size() - 1 << " -- " << 0
                 << ".\n";
          cdt.insert_constraint(points[points.size() - 1], points[0]);

          for (UnsignedIndex_t j = 0; j < points.size() - 1; ++j) {
            myfile << "Adding constraint " << j << " -- " << j + 1 << ".\n";
            cdt.insert_constraint(points[j], points[j + 1]);
          }
        }
        start_points += added_points;
        total_signed_area += 0.5 * signed_area;
      }
    }

    myfile << "Mesh has " << cdt.number_of_vertices() << " vertices.\n";
    myfile << "Refining with length-scale " << length_scale << ".\n";
    // sleep(1.0e-4);
    CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(),
                                 list_of_seeds.end(),
                                 Criteria(0.15, length_scale), false);
    myfile << "Mesh has " << cdt.number_of_vertices() << " vertices.\n";
    myfile << "Mesh has " << cdt.number_of_faces() << " faces.\n";
    // CGAL::lloyd_optimize_mesh_2(cdt,
    //                             CGAL::parameters::max_iteration_number
    //                             = 20);
    auto& vlist = returned_surface->getVertexList();
    auto& tlist = returned_surface->getTriangleList();
    UnsignedIndex_t count = 0;
    CDT::Finite_faces_iterator face;
    myfile << "Counting faces.\n";
    for (face = cdt.finite_faces_begin(); face != cdt.finite_faces_end();
         face++) {
      if (face->is_in_domain()) {
        count++;
      }
    }
    vlist.resize(3 * count);
    tlist.resize(count, TriangulatedSurfaceOutput::TriangleStorage::value_type::
                            fromNoExistencePlane(vlist, {0, 0, 0}));
    count = 0;
    myfile << "Adding faces and vertices.\n";
    for (face = cdt.finite_faces_begin(); face != cdt.finite_faces_end();
         face++) {
      if (face->is_in_domain()) {
        myfile << "Adding face " << count << ".\n";
        tlist[count] = TriangulatedSurfaceOutput::TriangleStorage::value_type::
            fromNoExistencePlane(vlist,
                                 {3 * count, 3 * count + 1, 3 * count + 2});
        for (UnsignedIndex_t d = 0; d < 3; d++) {
          const double x = CGAL::to_double(face->vertex(d)->point().x());
          const double y = CGAL::to_double(face->vertex(d)->point().y());
          const double z =
              -aligned_paraboloid.a() * x * x - aligned_paraboloid.b() * y * y;
          vlist[3 * count + d] = Pt(x, y, z);
          auto neigh = face->neighbor(d);
          if (!neigh->is_in_domain()) {
            returned_surface->addBoundaryEdge(3 * count + (d + 1) % 3,
                                              3 * count + (d + 2) % 3);
          }
        }
        count++;
      }
    }

    // Translate and rotate triangulated surface vertices
    myfile << "Moving vertices in canonical frame.\n";
    const auto& datum = paraboloid_m.getDatum();
    const auto& ref_frame = paraboloid_m.getReferenceFrame();
    for (auto& vertex : vlist) {
      const Pt base_pt = vertex;
      vertex = Pt(0.0, 0.0, 0.0);
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        for (UnsignedIndex_t n = 0; n < 3; ++n) {
          vertex[n] += ref_frame[d][n] * base_pt[d];
        }
      }
      vertex += datum;
    }

    myfile << "Finished triangulating surface.\n";
    myfile.close();
#elif defined IRL_USE_TRIANGLE
    // Second, we approximate the arc length of the arc, so as
    // to know how many times it needs to be split
    std::vector<REAL> input_points;
    std::vector<REAL> input_holes;
    std::vector<int> input_segments;
    const UnsignedIndex_t nCurves =
        static_cast<UnsignedIndex_t>(list_of_closed_curves.size());
    // Loop over curves
    UnsignedIndex_t start_points = 0;
    double total_signed_area = 0.0;
    for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
      const UnsignedIndex_t nLocalArcs = list_of_closed_curves[i].size();
      // Loop over arcs of curve
      UnsignedIndex_t added_points = 0;
      double signed_area = 0.0;
      for (UnsignedIndex_t j = 0; j < nLocalArcs; ++j) {
        // Compute approximate arc length
        const RationalBezierArc& arc = list_of_closed_curves[i][j];
        // const auto& sp = arc.start_point();
        // const auto& ep = arc.start_point();
        // signed_area += (sp[0] * ep[1] - ep[0] * sp[1]);
        const double arc_length = arc.arc_length();

        // Split arc
        UnsignedIndex_t nSplit = a_nsplit <= 0 ? 1 : a_nsplit;
        if (length_scale_ref > 0.0) {
          nSplit = static_cast<UnsignedIndex_t>(arc_length / length_scale_ref);
          nSplit = nSplit < a_nsplit ? a_nsplit : nSplit;
        }
        const double step = 1.0 / static_cast<double>(nSplit);
        length_scale = std::min(length_scale, step * arc_length);
        if (length_scale_ref > 0.0) length_scale = length_scale_ref;
        // added_points += nSplit;
        // const auto start_ind = input_points.size();
        // input_points.resize(start_ind + 2 * nSplit);
        // auto loc = input_points.begin() + start_ind;
        Pt previous_pt = arc.point(0.0);
        for (UnsignedIndex_t k = 1; k <= nSplit; ++k) {
          const double t = static_cast<double>(k) * step;
          const auto pt = arc.point(t);
          if (squaredMagnitude(pt - previous_pt) >
              1.0e8 * DBL_EPSILON * DBL_EPSILON) {
            input_points.push_back(pt[0]);
            input_points.push_back(pt[1]);
            previous_pt = pt;
            added_points++;
          }
        }
      }

      if (added_points >= 3) {
        signed_area += (input_points[start_points + 2 * added_points - 2] *
                            input_points[start_points + 1] -
                        input_points[start_points + 0] *
                            input_points[start_points + 2 * added_points - 1]);
        for (UnsignedIndex_t j = 0; j < added_points - 1; ++j) {
          signed_area += (input_points[start_points + 2 * j + 0] *
                              input_points[start_points + 2 * j + 3] -
                          input_points[start_points + 2 * j + 2] *
                              input_points[start_points + 2 * j + 1]);
        }

        if (nCurves > 1 && signed_area < 0.0) {
          // Add hole
          const auto p1x = input_points[start_points];
          const auto p1y = input_points[start_points + 1];
          const auto p2x = input_points[start_points + 2];
          const auto p2y = input_points[start_points + 3];
          std::array<double, 2> hole_location{
              {0.5 * (p1x + p2x), 0.5 * (p1y + p2y)}};
          Normal shift_dir = Normal(p2y - p1y, p1x - p2x, 0.0);
          shift_dir.normalize();
          const auto start_ind = input_holes.size();
          input_holes.resize(start_ind + 2);
          input_holes[start_ind] = 0.0;
          // hole_location[0] - (500.0 * DBL_EPSILON) *
          // shift_dir[0];
          input_holes[start_ind + 1] = 0.0;
          // hole_location[1] - (500.0 * DBL_EPSILON) *
          // shift_dir[1];
        }

        // Create segments
        const int seg_size = input_segments.size();
        input_segments.resize(seg_size + 2 * (added_points));
        auto seg_loc = input_segments.begin() + seg_size;
        *(seg_loc++) = start_points + added_points - 1;
        *(seg_loc++) = start_points;
        for (UnsignedIndex_t j = start_points;
             j < start_points + added_points - 1; ++j) {
          *(seg_loc++) = j;
          *(seg_loc++) = j + 1;
        }
        start_points += added_points;
        total_signed_area += 0.5 * signed_area;
      }
    }

    // Below section is for Triangle library
    if (input_points.size() > 0) {
      // std::cout << " Total area = " << total_signed_area <<
      // " compared to "
      //           << 2.0 * length_scale * length_scale <<
      //           std::endl;
      if (std::fabs(total_signed_area) > length_scale * length_scale) {
        // Calling triangulation library
        struct triangulateio in = {0}, out = {0};
        in.numberofpoints = input_points.size() / 2;
        in.pointlist = input_points.data();

        std::vector<int> pointmarkerlist(in.numberofpoints, 1);
        in.pointmarkerlist = pointmarkerlist.data();

        in.numberofsegments = input_segments.size() / 2;
        in.segmentlist = input_segments.data();
        std::vector<int> segmentmarkerlist(in.numberofsegments, 1);
        in.segmentmarkerlist = segmentmarkerlist.data();

        in.numberofholes = input_holes.size() / 2;
        if (in.numberofholes > 0) {
          in.holelist = input_holes.data();
        }

        char flags[50];
        sprintf(flags, "pza%.15feiQ", 0.5 * length_scale * length_scale);

        // std::cout << "Calling triangle with flags " <<
        // flags << " and with "
        //           << in.numberofpoints << " points and " <<
        //           in.numberofsegments
        //           << " segments and " << in.numberofholes
        //           << " holes and max area = "
        //           << 0.5 * length_scale * length_scale <<
        //           std::endl;

        // for (UnsignedIndex_t i = 0; i < in.numberofpoints;
        // ++i) {
        //   const double x = in.pointlist[2 * i + 0];
        //   const double y = in.pointlist[2 * i + 1];
        //   std::cout << "Point " << i << " = (" << x << ", "
        //   << y << ")"
        //             << std::endl;
        // }
        // for (UnsignedIndex_t i = 0; i <
        // in.numberofsegments; ++i) {
        //   const int j = in.segmentlist[2 * i + 0];
        //   const int k = in.segmentlist[2 * i + 1];
        //   std::cout << "Segment " << i << " = (" << j << ",
        //   " << k << ")"
        //             << std::endl;
        // }

        try {
          triangulate_from_lib(flags, &in, &out, (struct triangulateio*)NULL);
          // std::cout << "Triangle finished" << std::endl;

        } catch (std::runtime_error& e) {
          // std::cerr << e.what() << std::endl;
          // free(in.pointlist);
          free(in.pointattributelist);
          // free(in.pointmarkerlist);
          free(in.trianglelist);
          free(in.triangleattributelist);
          free(in.trianglearealist);
          free(in.neighborlist);
          // free(in.segmentlist);
          // free(in.segmentmarkerlist);
          // free(in.holelist);
          free(in.regionlist);
          free(in.edgelist);
          free(in.edgemarkerlist);
          free(in.normlist);
          free(out.pointlist);
          free(out.pointattributelist);
          free(out.pointmarkerlist);
          free(out.trianglelist);
          free(out.triangleattributelist);
          free(out.trianglearealist);
          free(out.neighborlist);
          free(out.segmentlist);
          free(out.segmentmarkerlist);
          free(out.regionlist);
          free(out.edgelist);
          free(out.edgemarkerlist);
          free(out.normlist);
        }

        auto& vlist = returned_surface->getVertexList();
        vlist.resize(out.numberofpoints);
        for (UnsignedIndex_t i = 0; i < out.numberofpoints; ++i) {
          const double x = out.pointlist[2 * i + 0];
          const double y = out.pointlist[2 * i + 1];
          const double z =
              -aligned_paraboloid.a() * x * x - aligned_paraboloid.b() * y * y;
          vlist[i] = Pt(x, y, z);
        }

        // Translate and rotate triangulated surface vertices
        const auto& datum = paraboloid_m.getDatum();
        const auto& ref_frame = paraboloid_m.getReferenceFrame();
        for (auto& vertex : vlist) {
          const Pt base_pt = vertex;
          vertex = Pt(0.0, 0.0, 0.0);
          for (UnsignedIndex_t d = 0; d < 3; ++d) {
            for (UnsignedIndex_t n = 0; n < 3; ++n) {
              vertex[n] += ref_frame[d][n] * base_pt[d];
            }
          }
          vertex += datum;
        }

        for (UnsignedIndex_t i = 0; i < out.numberofedges; ++i) {
          if (out.edgemarkerlist[i] == 1) {
            returned_surface->addBoundaryEdge(out.edgelist[2 * i],
                                              out.edgelist[2 * i + 1]);
          }
        }

        auto& tlist = returned_surface->getTriangleList();
        tlist.resize(out.numberoftriangles,
                     TriangulatedSurfaceOutput::TriangleStorage::value_type::
                         fromNoExistencePlane(vlist, {0, 0, 0}));
        for (UnsignedIndex_t i = 0; i < out.numberoftriangles; ++i) {
          tlist[i] = TriangulatedSurfaceOutput::TriangleStorage::value_type::
              fromNoExistencePlane(
                  vlist,
                  {static_cast<UnsignedIndex_t>(out.trianglelist[3 * i + 0]),
                   static_cast<UnsignedIndex_t>(out.trianglelist[3 * i + 1]),
                   static_cast<UnsignedIndex_t>(out.trianglelist[3 * i + 2])});
        }

        /* free all allocated arrays, including those
         * allocated by Triangle.
         */
        free(out.pointlist);
        free(out.pointattributelist);
        free(out.pointmarkerlist);
        free(out.trianglelist);
        free(out.triangleattributelist);
        free(out.trianglearealist);
        free(out.neighborlist);
        free(out.segmentlist);
        free(out.segmentmarkerlist);
        free(out.regionlist);
        free(out.edgelist);
        free(out.edgemarkerlist);
        free(out.normlist);
        // free(in.pointlist);
        free(in.pointattributelist);
        // free(in.pointmarkerlist);
        free(in.trianglelist);
        free(in.triangleattributelist);
        free(in.trianglearealist);
        free(in.neighborlist);
        // free(in.segmentlist);
        // free(in.segmentmarkerlist);
        // free(in.holelist);
        free(in.regionlist);
        free(in.edgelist);
        free(in.edgemarkerlist);
        free(in.normlist);
      } else {  // Triangulate by hand
        auto& vlist = returned_surface->getVertexList();
        vlist.resize(input_points.size() / 2);
        for (UnsignedIndex_t i = 0; i < input_points.size() / 2; ++i) {
          const double x = input_points[2 * i + 0];
          const double y = input_points[2 * i + 1];
          const double z =
              -aligned_paraboloid.a() * x * x - aligned_paraboloid.b() * y * y;
          vlist[i] = Pt(x, y, z);
        }

        // Translate and rotate triangulated surface vertices
        const auto& datum = paraboloid_m.getDatum();
        const auto& ref_frame = paraboloid_m.getReferenceFrame();
        for (auto& vertex : vlist) {
          const Pt base_pt = vertex;
          vertex = Pt(0.0, 0.0, 0.0);
          for (UnsignedIndex_t d = 0; d < 3; ++d) {
            for (UnsignedIndex_t n = 0; n < 3; ++n) {
              vertex[n] += ref_frame[d][n] * base_pt[d];
            }
          }
          vertex += datum;
        }

        returned_surface->addBoundaryEdge(input_points.size() / 2 - 1, 0);
        for (UnsignedIndex_t i = 0; i < input_points.size() / 2 - 1; ++i) {
          returned_surface->addBoundaryEdge(i, i + 1);
        }

        auto& tlist = returned_surface->getTriangleList();
        tlist.resize(input_points.size() / 2 - 2,
                     TriangulatedSurfaceOutput::TriangleStorage::value_type::
                         fromNoExistencePlane(vlist, {0, 0, 0}));
        for (UnsignedIndex_t i = 0; i < input_points.size() / 2 - 2; ++i) {
          tlist[i] = TriangulatedSurfaceOutput::TriangleStorage::value_type::
              fromNoExistencePlane(vlist, {0, i + 1, i + 2});
        }
      }
    }
#endif
  }
}

inline std::ostream& operator<<(
    std::ostream& out,
    const ParametrizedSurfaceOutput& a_parametrized_surface) {
  const auto& aligned_paraboloid =
      a_parametrized_surface.getParaboloid().getAlignedParaboloid();
  out.precision(16);
  out << std::scientific << aligned_paraboloid.a() << " "
      << aligned_paraboloid.b() << std::endl;
  for (UnsignedIndex_t i = 0; i < a_parametrized_surface.size(); ++i) {
    out << a_parametrized_surface[i];
    if (i < a_parametrized_surface.size() - 1) out << std::endl;
  }
  return out;
}

}  // namespace IRL

#endif  // IRL_PARABOLOID_RECONSTRUCTION_PARAMETRIZED_SURFACE_TPP_
