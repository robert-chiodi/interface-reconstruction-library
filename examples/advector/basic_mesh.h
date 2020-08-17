// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_ADVECTOR_BASIC_MESH_H_
#define EXAMPLES_ADVECTOR_BASIC_MESH_H_

#include <cstring>

#include "src/geometry/general/pt.h"

/// \brief This is a basic mesh class that stores
/// the number of cells in the mesh, the number of
/// ghost cells, and provides access to the underlying
/// cell face locations. These cell face locations are
/// known for -ngc_m to n+1+ngc_m.
class BasicMesh {
 public:
  /// \brief Default constructor
  BasicMesh(void)
      : x_m(nullptr),
        y_m(nullptr),
        z_m(nullptr),
        dx_m(0.0),
        dy_m(0.0),
        dz_m(0.0),
        nx_m(0),
        ny_m(0),
        nz_m(0),
        ngc_m(0),
        nxo_m(0),
        nyo_m(0),
        nzo_m(0),
        imin_m(0),
        jmin_m(0),
        kmin_m(0),
        imax_m(0),
        jmax_m(0),
        kmax_m(0),
        imino_m(0),
        jmino_m(0),
        kmino_m(0),
        imaxo_m(0),
        jmaxo_m(0),
        kmaxo_m(0) {}

  /// \brief Constructor with allocation for mesh
  BasicMesh(const int a_nx, const int a_ny, const int a_nz, const int a_ngc)
      : dx_m(0.0),
        dy_m(0.0),
        dz_m(0.0),
        nx_m(a_nx),
        ny_m(a_ny),
        nz_m(a_nz),
        ngc_m(a_ngc),
        nxo_m(a_nx + 2 * a_ngc),
        nyo_m(a_ny + 2 * a_ngc),
        nzo_m(a_nz + 2 * a_ngc),
        imin_m(0),
        jmin_m(0),
        kmin_m(0),
        imax_m(nx_m - 1),
        jmax_m(ny_m - 1),
        kmax_m(nz_m - 1),
        imino_m(-ngc_m),
        jmino_m(-ngc_m),
        kmino_m(-ngc_m),
        imaxo_m(nx_m - 1 + ngc_m),
        jmaxo_m(ny_m - 1 + ngc_m),
        kmaxo_m(nz_m - 1 + ngc_m) {
    x_m = new double[a_nx + 1 + 2 * ngc_m];
    y_m = new double[a_ny + 1 + 2 * ngc_m];
    z_m = new double[a_nz + 1 + 2 * ngc_m];
  }

  /// \brief Copy construction
  BasicMesh(const BasicMesh& other) {
    delete[] x_m;
    delete[] y_m;
    delete[] z_m;

    dx_m = other.dx_m;
    dy_m = other.dy_m;
    dz_m = other.dz_m;
    nx_m = other.nx_m;
    ny_m = other.ny_m;
    nz_m = other.nz_m;
    ngc_m = other.ngc_m;
    nxo_m = other.nxo_m;
    nyo_m = other.nyo_m;
    nzo_m = other.nzo_m;
    imin_m = other.imin_m;
    jmin_m = other.jmin_m;
    kmin_m = other.kmin_m;
    imax_m = other.imax_m;
    jmax_m = other.jmax_m;
    kmax_m = other.kmax_m;
    imino_m = other.imin_m;
    jmino_m = other.jmino_m;
    kmino_m = other.kmino_m;
    imaxo_m = other.imaxo_m;
    jmaxo_m = other.jmaxo_m;
    kmaxo_m = other.kmaxo_m;
    x_m = new double[nxo_m + 1];
    y_m = new double[nyo_m + 1];
    z_m = new double[nzo_m + 1];
    std::memcpy(x_m, other.x_m,
                static_cast<std::size_t>(sizeof(double)) *
                        static_cast<std::size_t>(nxo_m) +
                    1);
    std::memcpy(y_m, other.y_m,
                static_cast<std::size_t>(sizeof(double)) *
                        static_cast<std::size_t>(nyo_m) +
                    1);
    std::memcpy(z_m, other.z_m,
                static_cast<std::size_t>(sizeof(double)) *
                        static_cast<std::size_t>(nzo_m) +
                    1);
  }

  /// \brief Move construction
  BasicMesh(BasicMesh&& other) {
    delete[] x_m;
    delete[] y_m;
    delete[] z_m;

    dx_m = other.dx_m;
    dy_m = other.dy_m;
    dz_m = other.dz_m;
    nx_m = other.nx_m;
    ny_m = other.ny_m;
    nz_m = other.nz_m;
    ngc_m = other.ngc_m;
    nxo_m = other.nxo_m;
    nyo_m = other.nyo_m;
    nzo_m = other.nzo_m;
    imin_m = other.imin_m;
    jmin_m = other.jmin_m;
    kmin_m = other.kmin_m;
    imax_m = other.imax_m;
    jmax_m = other.jmax_m;
    kmax_m = other.kmax_m;
    imino_m = other.imin_m;
    jmino_m = other.jmino_m;
    kmino_m = other.kmino_m;
    imaxo_m = other.imaxo_m;
    jmaxo_m = other.jmaxo_m;
    kmaxo_m = other.kmaxo_m;
    x_m = other.x_m;
    y_m = other.y_m;
    z_m = other.z_m;
    other.x_m = nullptr;
    other.y_m = nullptr;
    other.z_m = nullptr;
  }

  /// \brief Copy assignment
  BasicMesh& operator=(const BasicMesh& other) {
    if (this != &other) {
      delete[] x_m;
      delete[] y_m;
      delete[] z_m;

      dx_m = other.dx_m;
      dy_m = other.dy_m;
      dz_m = other.dz_m;
      nx_m = other.nx_m;
      ny_m = other.ny_m;
      nz_m = other.nz_m;
      ngc_m = other.ngc_m;
      nxo_m = other.nxo_m;
      nyo_m = other.nyo_m;
      nzo_m = other.nzo_m;
      imin_m = other.imin_m;
      jmin_m = other.jmin_m;
      kmin_m = other.kmin_m;
      imax_m = other.imax_m;
      jmax_m = other.jmax_m;
      kmax_m = other.kmax_m;
      imino_m = other.imin_m;
      jmino_m = other.jmino_m;
      kmino_m = other.kmino_m;
      imaxo_m = other.imaxo_m;
      jmaxo_m = other.jmaxo_m;
      kmaxo_m = other.kmaxo_m;
      x_m = new double[nxo_m + 1];
      y_m = new double[nyo_m + 1];
      z_m = new double[nzo_m + 1];
      std::memcpy(x_m, other.x_m,
                  static_cast<std::size_t>(sizeof(double)) *
                          static_cast<std::size_t>(nxo_m) +
                      1);
      std::memcpy(y_m, other.y_m,
                  static_cast<std::size_t>(sizeof(double)) *
                          static_cast<std::size_t>(nyo_m) +
                      1);
      std::memcpy(z_m, other.z_m,
                  static_cast<std::size_t>(sizeof(double)) *
                          static_cast<std::size_t>(nzo_m) +
                      1);
    }
    return (*this);
  }

  /// \brief Move assignment
  BasicMesh& operator=(BasicMesh&& other) {
    if (this != &other) {
      delete[] x_m;
      delete[] y_m;
      delete[] z_m;

      dx_m = other.dx_m;
      dy_m = other.dy_m;
      dz_m = other.dz_m;
      nx_m = other.nx_m;
      ny_m = other.ny_m;
      nz_m = other.nz_m;
      ngc_m = other.ngc_m;
      nxo_m = other.nxo_m;
      nyo_m = other.nyo_m;
      nzo_m = other.nzo_m;
      imin_m = other.imin_m;
      jmin_m = other.jmin_m;
      kmin_m = other.kmin_m;
      imax_m = other.imax_m;
      jmax_m = other.jmax_m;
      kmax_m = other.kmax_m;
      imino_m = other.imin_m;
      jmino_m = other.jmino_m;
      kmino_m = other.kmino_m;
      imaxo_m = other.imaxo_m;
      jmaxo_m = other.jmaxo_m;
      kmaxo_m = other.kmaxo_m;
      x_m = other.x_m;
      y_m = other.y_m;
      z_m = other.z_m;
      other.x_m = nullptr;
      other.y_m = nullptr;
      other.z_m = nullptr;
    }
    return (*this);
  }

  /// \brief Set x/y/z cell.
  void setCellBoundaries(const IRL::Pt& a_bottom_bounding_box,
                         const IRL::Pt& a_top_bounding_box);

  int getNx(void) { return nx_m; }
  int getNy(void) { return ny_m; }
  int getNz(void) { return nz_m; }
  int getNgc(void) { return ngc_m; }
  int getNxo(void) { return nxo_m; }
  int getNyo(void) { return nyo_m; }
  int getNzo(void) { return nzo_m; }
  int imin(void) { return imin_m; }
  int jmin(void) { return jmin_m; }
  int kmin(void) { return kmin_m; }
  int imax(void) { return imax_m; }
  int jmax(void) { return jmax_m; }
  int kmax(void) { return kmax_m; }
  int imino(void) { return imino_m; }
  int jmino(void) { return jmino_m; }
  int kmino(void) { return kmino_m; }
  int imaxo(void) { return imaxo_m; }
  int jmaxo(void) { return jmaxo_m; }
  int kmaxo(void) { return kmaxo_m; }

  int getNx(void) const { return nx_m; }
  int getNy(void) const { return ny_m; }
  int getNz(void) const { return nz_m; }
  int getNgc(void) const { return ngc_m; }
  int getNxo(void) const { return nxo_m; }
  int getNyo(void) const { return nyo_m; }
  int getNzo(void) const { return nzo_m; }
  int imin(void) const { return imin_m; }
  int jmin(void) const { return jmin_m; }
  int kmin(void) const { return kmin_m; }
  int imax(void) const { return imax_m; }
  int jmax(void) const { return jmax_m; }
  int kmax(void) const { return kmax_m; }
  int imino(void) const { return imino_m; }
  int jmino(void) const { return jmino_m; }
  int kmino(void) const { return kmino_m; }
  int imaxo(void) const { return imaxo_m; }
  int jmaxo(void) const { return jmaxo_m; }
  int kmaxo(void) const { return kmaxo_m; }

  /// \brief Access to x cell locations
  double& x(const int i) {
    assert(i >= this->imino() && i <= this->imaxo() + 1);
    return x_m[i + getNgc()];
  }

  /// \brief Access to y cell locations
  double& y(const int j) {
    assert(j >= this->jmino() && j <= this->jmaxo() + 1);
    return y_m[j + getNgc()];
  }

  /// \brief Access to z cell locations.
  double& z(const int k) {
    assert(k >= this->kmino() && k <= this->kmaxo() + 1);
    return z_m[k + getNgc()];
  }

  /// \brief Access to xm cell center locations
  double xm(const int i) { return 0.5 * (this->x(i) + this->x(i + 1)); }

  /// \brief Access to ym cell center locations
  double ym(const int j) { return 0.5 * (this->y(j) + this->y(j + 1)); }

  /// \brief Access to zm cell center locations
  double zm(const int k) { return 0.5 * (this->z(k) + this->z(k + 1)); }

  /// \brief Access to x cell locations
  double& x(const int i) const {
    assert(i >= this->imino() && i <= this->imaxo() + 1);
    return x_m[i + getNgc()];
  }

  /// \brief Access to y cell locations
  double& y(const int j) const {
    assert(j >= this->jmino() && j <= this->jmaxo() + 1);
    return y_m[j + getNgc()];
  }

  /// \brief Access to z cell locations.
  double& z(const int k) const {
    assert(k >= this->kmino() && k <= this->kmaxo() + 1);
    return z_m[k + getNgc()];
  }

  /// \brief Access to xm cell center locations
  double xm(const int i) const { return 0.5 * (this->x(i) + this->x(i + 1)); }

  /// \brief Access to ym cell center locations
  double ym(const int j) const { return 0.5 * (this->y(j) + this->y(j + 1)); }

  /// \brief Access to zm cell center locations
  double zm(const int k) const { return 0.5 * (this->z(k) + this->z(k + 1)); }

  double dx(void) { return dx_m; }
  double dy(void) { return dy_m; }
  double dz(void) { return dz_m; }
  double dx(void) const { return dx_m; }
  double dy(void) const { return dy_m; }
  double dz(void) const { return dz_m; }

  double lx(void) const {
    return this->dx() * static_cast<double>(this->getNx());
  }
  double ly(void) const {
    return this->dy() * static_cast<double>(this->getNy());
  }
  double lz(void) const {
    return this->dz() * static_cast<double>(this->getNz());
  }

  /// \brief Const version of getIndices.
  void getIndices(const IRL::Pt& a_location, int* a_indices) const;

  /// \brief Return size of mesh, the amount that needs to
  /// be allocated for data arrays.
  std::size_t size(void) {
    return static_cast<std::size_t>(this->getNxo() * this->getNyo() *
                                    this->getNzo());
  }

  /// \brief Const version of size.
  std::size_t size(void) const {
    return static_cast<std::size_t>(this->getNxo() * this->getNyo() *
                                    this->getNzo());
  }

  /// \brief Default destructor
  ~BasicMesh(void) {
    delete[] x_m;
    delete[] y_m;
    delete[] z_m;
  }

 private:
  double* x_m;
  double* y_m;
  double* z_m;
  double dx_m, dy_m, dz_m;
  int nx_m, ny_m, nz_m;
  int ngc_m;
  int nxo_m, nyo_m, nzo_m;
  int imin_m, jmin_m, kmin_m;
  int imax_m, jmax_m, kmax_m;
  int imino_m, jmino_m, kmino_m;
  int imaxo_m, jmaxo_m, kmaxo_m;
};

#endif  // EXAMPLES_ADVECTOR_BASIC_MESH_H_
