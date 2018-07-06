// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_ADVECTOR_DATA_H_
#define EXAMPLES_ADVECTOR_DATA_H_

#include <algorithm>
#include <cstring>

#include "src/geometry/general/pt.h"

#include "examples/advector/basic_mesh.h"

/// \brief A basic multi-dimensional data container.
template <class ContainedType>
class Data {
 public:
  using value_t = ContainedType;

  /// \brief Default constructor
  Data(void) = default;

  /// \brief Construct and set initial size.
  explicit Data(const BasicMesh* a_mesh_ptr) : mesh_m(a_mesh_ptr) {
    data_m = new ContainedType[this->getMesh().size()];
  }

  // Copy constructor
  Data(const Data& other) {
    delete[] this->data_m;
    const BasicMesh& mesh = other.getMesh();
    mesh_m = other.mesh_m;
    data_m = new ContainedType[mesh.size()];
    std::memcpy(data_m, other.data_m, sizeof(ContainedType) * mesh.size());
  }

  // Move constructor
  Data(Data&& other) {
    delete[] data_m;
    mesh_m = other.mesh_m;
    data_m = other.data_m;
    other.mesh_m = nullptr;
    other.data_m = nullptr;
  }

  // Copy assignment
  Data& operator=(const Data& other) {
    if (this != &other) {
      delete[] data_m;
      const BasicMesh& mesh = other.getMesh();
      mesh_m = other.mesh_m;
      data_m = new ContainedType[mesh.size()];
      std::memcpy(data_m, other.data_m, sizeof(ContainedType) * mesh.size());
    }
    return (*this);
  }

  // Move assignment
  Data& operator=(Data&& other) {
    if (this != &other) {
      delete[] data_m;
      mesh_m = other.mesh_m;
      data_m = other.data_m;
      other.mesh_m = nullptr;
      other.data_m = nullptr;
    }
    return (*this);
  }

  /// \brief Return the pointer to the mesh.
  const BasicMesh& getMesh(void) { return *mesh_m; }
  const BasicMesh& getMesh(void) const { return *mesh_m; }

  /// \brief Provide access to data.
  ContainedType& operator()(const int i, const int j, const int k) {
    return data_m[this->calculateIndex(i, j, k)];
  }

  /// \brief Provide const access to data.
  const ContainedType& operator()(const int i, const int j, const int k) const {
    return data_m[this->calculateIndex(i, j, k)];
  }

  /// \brief Update the border for periodicity.
  void updateBorder(void);
  void updateLowerX(void);
  void updateUpperX(void);
  void updateLowerY(void);
  void updateUpperY(void);
  void updateLowerZ(void);
  void updateUpperZ(void);

  /// \brief Tri-linearly interpolate data to a point.
  ContainedType interpolate(const IRL::Pt& a_location);

  /// \brief Const version of interpolate.
  ContainedType interpolate(const IRL::Pt& a_location) const;

  /// \brief Destructor to delete memory allocated during construction.
  ~Data(void) { delete[] data_m; }

 private:
  /// \brief Calculate the index, where the first real (non-ghost) cell is at 0
  /// and the fastest changing indices are k, j, i.
  std::size_t calculateIndex(const int i, const int j, const int k) {
    assert(i >= this->getMesh().imino() && i <= this->getMesh().imaxo());
    assert(j >= this->getMesh().jmino() && j <= this->getMesh().jmaxo());
    assert(k >= this->getMesh().kmino() && k <= this->getMesh().kmaxo());
    assert(static_cast<std::size_t>(
               (i + this->getMesh().getNgc()) * this->getMesh().getNyo() *
                   this->getMesh().getNzo() +
               (j + this->getMesh().getNgc()) * this->getMesh().getNzo() +
               (k + this->getMesh().getNgc())) < this->getMesh().size());
    return static_cast<std::size_t>(
        (i + this->getMesh().getNgc()) * this->getMesh().getNyo() *
            this->getMesh().getNzo() +
        (j + this->getMesh().getNgc()) * this->getMesh().getNzo() +
        (k + this->getMesh().getNgc()));
  }

  /// \brief Const getInd
  std::size_t calculateIndex(const int i, const int j, const int k) const {
    assert(i >= this->getMesh().imino() && i <= this->getMesh().imaxo());
    assert(j >= this->getMesh().jmino() && j <= this->getMesh().jmaxo());
    assert(k >= this->getMesh().kmino() && k <= this->getMesh().kmaxo());
    assert(static_cast<std::size_t>(
               (i + this->getMesh().getNgc()) * this->getMesh().getNyo() *
                   this->getMesh().getNzo() +
               (j + this->getMesh().getNgc()) * this->getMesh().getNzo() +
               (k + this->getMesh().getNgc())) < this->getMesh().size());
    return static_cast<std::size_t>(
        (i + this->getMesh().getNgc()) * this->getMesh().getNyo() *
            this->getMesh().getNzo() +
        (j + this->getMesh().getNgc()) * this->getMesh().getNzo() +
        (k + this->getMesh().getNgc()));
  }

  ContainedType* data_m = nullptr;
  const BasicMesh* mesh_m = nullptr;
};

template <class ContainedType>
void Data<ContainedType>::updateBorder(void) {
  this->updateLowerX();
  this->updateUpperX();
  this->updateLowerY();
  this->updateUpperY();
  this->updateLowerZ();
  this->updateUpperZ();
}

template <class ContainedType>
void Data<ContainedType>::updateLowerX(void) {
  const BasicMesh& mesh = this->getMesh();
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*this)(i, j, k) = (*this)(mesh.imax() + 1 + i, j, k);
      }
    }
  }
}
template <class ContainedType>
void Data<ContainedType>::updateUpperX(void) {
  const BasicMesh& mesh = this->getMesh();
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*this)(i, j, k) = (*this)(i - mesh.getNx(), j, k);
      }
    }
  }
}

template <class ContainedType>
void Data<ContainedType>::updateLowerY(void) {
  const BasicMesh& mesh = this->getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*this)(i, j, k) = (*this)(i, mesh.jmax() + 1 + j, k);
      }
    }
  }
}
template <class ContainedType>
void Data<ContainedType>::updateUpperY(void) {
  const BasicMesh& mesh = this->getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*this)(i, j, k) = (*this)(i, j - mesh.getNy(), k);
      }
    }
  }
}
template <class ContainedType>
void Data<ContainedType>::updateLowerZ(void) {
  const BasicMesh& mesh = this->getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k < mesh.kmin(); ++k) {
        (*this)(i, j, k) = (*this)(i, j, mesh.kmax() + 1 + k);
      }
    }
  }
}
template <class ContainedType>
void Data<ContainedType>::updateUpperZ(void) {
  const BasicMesh& mesh = this->getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmax() + 1; k <= mesh.kmaxo(); ++k) {
        (*this)(i, j, k) = (*this)(i, j, k - mesh.getNz());
      }
    }
  }
}

template <class ContainedType>
ContainedType Data<ContainedType>::interpolate(const IRL::Pt& a_location) {
  const BasicMesh& mesh = this->getMesh();
  static int indices[3];
  mesh.getIndices(a_location, indices);
  int i = std::max(this->getMesh().imino(),
                   std::min(indices[0], this->getMesh().imaxo() - 1));
  int j = std::max(this->getMesh().jmino(),
                   std::min(indices[1], this->getMesh().jmaxo() - 1));
  int k = std::max(this->getMesh().kmino(),
                   std::min(indices[2], this->getMesh().kmaxo() - 1));
  double wx1 = (a_location[0] - mesh.xm(i)) / (mesh.xm(i + 1) - mesh.xm(i));
  double wy1 = (a_location[1] - mesh.ym(j)) / (mesh.ym(j + 1) - mesh.ym(j));
  double wz1 = (a_location[2] - mesh.zm(k)) / (mesh.zm(k + 1) - mesh.zm(k));
  double wx2 = 1.0 - wx1;
  double wy2 = 1.0 - wy1;
  double wz2 = 1.0 - wz1;
  return wz1 * (wy1 * (wx1 * (*this)(i + 1, j + 1, k + 1) +
                       wx2 * (*this)(i, j + 1, k + 1)) +
                wy2 * (wx1 * (*this)(i + 1, j, k + 1) +
                       wx2 * (*this)(i, j, k + 1))) +
         wz2 * (wy1 * (wx1 * (*this)(i + 1, j + 1, k) +
                       wx2 * (*this)(i, j + 1, k)) +
                wy2 * (wx1 * (*this)(i + 1, j, k) + wx2 * (*this)(i, j, k)));
}

template <class ContainedType>
ContainedType Data<ContainedType>::interpolate(
    const IRL::Pt& a_location) const {
  const BasicMesh& mesh = this->getMesh();
  static int indices[3];
  mesh.getIndices(a_location, indices);
  const int i = std::max(this->getMesh().imino(),
                         std::min(indices[0], this->getMesh().imaxo() - 1));
  const int j = std::max(this->getMesh().jmino(),
                         std::min(indices[1], this->getMesh().jmaxo() - 1));
  const int k = std::max(this->getMesh().kmino(),
                         std::min(indices[2], this->getMesh().kmaxo() - 1));
  const double wx1 =
      (a_location[0] - mesh.xm(i)) / (mesh.xm(i + 1) - mesh.xm(i));
  const double wy1 =
      (a_location[1] - mesh.ym(j)) / (mesh.ym(j + 1) - mesh.ym(j));
  const double wz1 =
      (a_location[2] - mesh.zm(k)) / (mesh.zm(k + 1) - mesh.zm(k));
  const double wx2 = 1.0 - wx1;
  const double wy2 = 1.0 - wy1;
  const double wz2 = 1.0 - wz1;
  return wz1 * (wy1 * (wx1 * (*this)(i + 1, j + 1, k + 1) +
                       wx2 * (*this)(i, j + 1, k + 1)) +
                wy2 * (wx1 * (*this)(i + 1, j, k + 1) +
                       wx2 * (*this)(i, j, k + 1))) +
         wz2 * (wy1 * (wx1 * (*this)(i + 1, j + 1, k) +
                       wx2 * (*this)(i, j + 1, k)) +
                wy2 * (wx1 * (*this)(i + 1, j, k) + wx2 * (*this)(i, j, k)));
}

#endif  // EXAMPLES_ADVECTOR_DATA_H_
