// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <mpi.h>
#include <iostream>

#include "examples/paraboloid_advector/reconstruction_types.h"
#include "examples/paraboloid_advector/vof_advection.h"
#include "examples/paraboloid_advector/vtk.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/paraboloid_intersection/paraboloid_intersection_amr.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"

void resetCentroids(const Data<IRL::LocalizedParaboloidLink<double>>&
                        a_link_localized_paraboloid,
                    Data<IRL::Pt>* a_liquid_centroid,
                    Data<IRL::Pt>* a_gas_centroid) {
  const BasicMesh& mesh = a_link_localized_paraboloid.getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        auto moments = IRL::getNormalizedVolumeMoments<
            IRL::SeparatedMoments<IRL::VolumeMoments>>(
            cell, a_link_localized_paraboloid(i, j, k).getNextReconstruction());
        (*a_liquid_centroid)(i, j, k) = moments[0].centroid();
        (*a_gas_centroid)(i, j, k) = moments[1].centroid();
      }
    }
  }
}

void connectMesh(
    const BasicMesh& a_mesh,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid) {
  IRL::LocalizedParaboloidLink<double>* neighbor_ptr;
  // Provide mesh connectivity information.
  IRL::UnsignedIndex_t unique_id = 0;

  for (int i = a_mesh.imino(); i <= a_mesh.imaxo(); ++i) {
    for (int j = a_mesh.jmino(); j <= a_mesh.jmaxo(); ++j) {
      for (int k = a_mesh.kmino(); k <= a_mesh.kmaxo(); ++k) {
        (*a_link_localized_paraboloid)(i, j, k).setId(unique_id);
        neighbor_ptr = i - 1 < a_mesh.imino()
                           ? nullptr
                           : &(*a_link_localized_paraboloid)(i - 1, j, k);
        (*a_link_localized_paraboloid)(i, j, k).setEdgeConnectivity(
            0, neighbor_ptr);
        neighbor_ptr = i + 1 > a_mesh.imaxo()
                           ? nullptr
                           : &(*a_link_localized_paraboloid)(i + 1, j, k);
        (*a_link_localized_paraboloid)(i, j, k).setEdgeConnectivity(
            1, neighbor_ptr);
        neighbor_ptr = j - 1 < a_mesh.jmino()
                           ? nullptr
                           : &(*a_link_localized_paraboloid)(i, j - 1, k);
        (*a_link_localized_paraboloid)(i, j, k).setEdgeConnectivity(
            2, neighbor_ptr);
        neighbor_ptr = j + 1 > a_mesh.jmaxo()
                           ? nullptr
                           : &(*a_link_localized_paraboloid)(i, j + 1, k);
        (*a_link_localized_paraboloid)(i, j, k).setEdgeConnectivity(
            3, neighbor_ptr);
        neighbor_ptr = k - 1 < a_mesh.kmino()
                           ? nullptr
                           : &(*a_link_localized_paraboloid)(i, j, k - 1);
        (*a_link_localized_paraboloid)(i, j, k).setEdgeConnectivity(
            4, neighbor_ptr);
        neighbor_ptr = k + 1 > a_mesh.kmaxo()
                           ? nullptr
                           : &(*a_link_localized_paraboloid)(i, j, k + 1);
        (*a_link_localized_paraboloid)(i, j, k).setEdgeConnectivity(
            5, neighbor_ptr);
        ++unique_id;
      }
    }
  }
}

std::array<int, 3> getIndexFromTag(const BasicMesh& a_mesh,
                                   const IRL::UnsignedIndex_t a_tag) {
  std::array<int, 3> indices;
  auto int_tag = static_cast<int64_t>(a_tag);
  indices[0] = static_cast<int>(int_tag / (a_mesh.getNzo() * a_mesh.getNyo()));
  indices[1] = static_cast<int>(
      (int_tag - indices[0] * (a_mesh.getNzo() * a_mesh.getNyo())) /
      a_mesh.getNzo());
  indices[2] =
      static_cast<int>(int_tag - a_mesh.getNzo() * indices[1] -
                       (indices[0] * (a_mesh.getNzo() * a_mesh.getNyo())));

  return {indices[0] + a_mesh.imino(), indices[1] + a_mesh.jmino(),
          indices[2] + a_mesh.kmino()};
}

void advectVOF(
    const std::string& a_advection_method,
    const std::string& a_reconstruction_method, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid,
    Data<double>* a_liquid_volume_fraction, Data<IRL::Pt>* a_liquid_centroid,
    Data<IRL::Pt>* a_gas_centroid, Data<IRL::Paraboloid>* a_interface) {
  if (a_advection_method == "FullLagrangian") {
    FullLagrangian::advectVOF(a_reconstruction_method, a_dt, a_U, a_V, a_W,
                              a_link_localized_paraboloid,
                              a_liquid_volume_fraction, a_liquid_centroid,
                              a_gas_centroid);
  } else if (a_advection_method == "SemiLagrangian") {
    SemiLagrangian::advectVOF(a_reconstruction_method, a_dt, a_U, a_V, a_W,
                              a_link_localized_paraboloid,
                              a_liquid_volume_fraction, a_liquid_centroid,
                              a_gas_centroid);
  } else if (a_advection_method == "SemiLagrangianCorrected") {
    SemiLagrangianCorrected::advectVOF(a_reconstruction_method, a_dt, a_U, a_V,
                                       a_W, a_link_localized_paraboloid,
                                       a_liquid_volume_fraction,
                                       a_liquid_centroid, a_gas_centroid);
  } else if (a_advection_method == "Split") {
    Split::advectVOF(a_reconstruction_method, a_dt, a_U, a_V, a_W,
                     a_link_localized_paraboloid, a_liquid_volume_fraction,
                     a_liquid_centroid, a_gas_centroid, a_interface);
  } else {
    std::cout << "Unknown advection method of : " << a_advection_method << '\n';
    std::cout << "Value entries are: FullLagrangian, SemiLagrangian, "
                 "SemiLagrangianCorrected, Split. \n";
    std::exit(-1);
  }
}

void Split::advectVOF(
    const std::string& a_reconstruction_method, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid,
    Data<double>* a_liquid_volume_fraction, Data<IRL::Pt>* a_liquid_centroid,
    Data<IRL::Pt>* a_gas_centroid, Data<IRL::Paraboloid>* a_interface) {
  const BasicMesh& mesh = a_liquid_volume_fraction->getMesh();
  Data<IRL::VolumeMoments> face_flux[3] = {Data<IRL::VolumeMoments>(&mesh),
                                           Data<IRL::VolumeMoments>(&mesh),
                                           Data<IRL::VolumeMoments>(&mesh)};
  Data<double> U_face[3] = {Data<double>(&mesh), Data<double>(&mesh),
                            Data<double>(&mesh)};
  auto cc = Data<double>(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        U_face[0](i, j, k) = 0.5 * (a_U(i, j, k) + a_U(i - 1, j, k));
        U_face[1](i, j, k) = 0.5 * (a_V(i, j, k) + a_V(i, j - 1, k));
        U_face[2](i, j, k) = 0.5 * (a_W(i, j, k) + a_W(i, j, k - 1));
      }
    }
  }
  U_face[0].updateBorder();
  U_face[1].updateBorder();
  U_face[2].updateBorder();
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        cc(i, j, k) = (*a_liquid_volume_fraction)(i, j, k) > 0.5 ? 1.0 : 0.0;
      }
    }
  }

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int split_proc = static_cast<int>(std::cbrt(static_cast<double>(size)));
  int imin, nx, jmin, ny, kmin, nz;
  const int NX = mesh.getNx(), NY = mesh.getNy(), NZ = mesh.getNz();

  if (size > 1) {
    if (size == split_proc * split_proc * split_proc) {
      for (int i = 0; i < split_proc; i++) {
        for (int j = 0; j < split_proc; j++) {
          for (int k = 0; k < split_proc; k++) {
            if (i + split_proc * j + split_proc * split_proc * k == rank) {
              imin = i * (NX / split_proc);
              nx = std::min((i + 1) * (NX / split_proc), NX) - imin;
              jmin = j * (NY / split_proc);
              ny = std::min((j + 1) * (NY / split_proc), NY) - jmin;
              kmin = k * (NZ / split_proc);
              nz = std::min((k + 1) * (NZ / split_proc), NZ) - kmin;
            }
          }
        }
      }
    } else {
      imin = rank * (NX / size);
      nx = std::min((rank + 1) * (NX / size), NX) - imin;
      jmin = 0;
      ny = NZ;
      kmin = 0;
      nz = NZ;
    }
  } else {
    imin = 0, nx = NX, jmin = 0, ny = NY, kmin = 0, nz = NZ;
  }

  // For now, naively advect everywhere in domain
  const IRL::Pt vec_dt[3] = {IRL::Pt(a_dt, 0.0, 0.0), IRL::Pt(0.0, a_dt, 0.0),
                             IRL::Pt(0.0, 0.0, a_dt)};
  // Direction sweep
  for (int dim = 0; dim < 3; ++dim) {
    getReconstruction(a_reconstruction_method, *a_liquid_volume_fraction,
                      *a_liquid_centroid, *a_gas_centroid,
                      *a_link_localized_paraboloid, a_dt, a_U, a_V, a_W,
                      a_interface);

    // for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    //   for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
    //     for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int i = imin; i < imin + nx; ++i) {
      for (int j = jmin; j < jmin + ny; ++j) {
        for (int k = kmin; k < kmin + nz; ++k) {
          int shift[3] = {0, 0, 0};
          double sign = -1.0;
          auto x0 = IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k));
          auto x1 = IRL::Pt(mesh.x(i + (dim == 0 ? 0 : 1)),
                            mesh.y(j + (dim == 1 ? 0 : 1)),
                            mesh.z(k + (dim == 2 ? 0 : 1)));
          if ((U_face[dim])(i, j, k) > 0.0) {
            x0 -= (U_face[dim])(i, j, k) * vec_dt[dim];
            shift[dim] = -1;
            sign = 1.0;
          } else {
            x1 += -(U_face[dim])(i, j, k) * vec_dt[dim];
          }
          auto face_cell = IRL::RectangularCuboid::fromBoundingPts(x0, x1);
          IRL::Paraboloid local_paraboloid = (*a_interface)(i, j, k);
          (face_flux[dim])(i, j, k) =
              sign * IRL::getVolumeMoments<IRL::VolumeMoments>(
                         face_cell, (*a_interface)(i + shift[0], j + shift[1],
                                                   k + shift[2]));
        }
      }
    }

    int vector_size = mesh.size();
    std::vector<double> face_flux_dim_local, face_flux_dim;
    face_flux_dim_local.resize(vector_size);
    face_flux_dim.resize(vector_size);
    std::fill(face_flux_dim_local.begin(), face_flux_dim_local.end(), 0.0);
    std::fill(face_flux_dim.begin(), face_flux_dim.end(), 0.0);

    for (int i = imin; i < imin + nx; ++i) {
      for (int j = jmin; j < jmin + ny; ++j) {
        for (int k = kmin; k < kmin + nz; ++k) {
          int address = (i + mesh.getNgc()) * mesh.getNyo() * mesh.getNzo() +
                        (j + mesh.getNgc()) * mesh.getNzo() +
                        (k + mesh.getNgc());
          face_flux_dim_local[address] = (face_flux[dim](i, j, k)).volume();
        }
      }
    }

    MPI_Allreduce(face_flux_dim_local.data(), face_flux_dim.data(), vector_size,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
      for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
        for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
          int address = (i + mesh.getNgc()) * mesh.getNyo() * mesh.getNzo() +
                        (j + mesh.getNgc()) * mesh.getNzo() +
                        (k + mesh.getNgc());
          (face_flux[dim](i, j, k)).volume() = face_flux_dim[address];
        }
      }
    }

    face_flux[dim].updateBorder();
    for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
      for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
        for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
          auto cell = IRL::RectangularCuboid::fromBoundingPts(
              IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
              IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
          const double cell_volume = cell.calculateVolume();
          const double previous_liquid_volume_fraction =
              (*a_liquid_volume_fraction)(i, j, k);

          // Update VOF
          (*a_liquid_volume_fraction)(i, j, k) =
              (previous_liquid_volume_fraction * cell_volume +
               +(face_flux[dim])(i, j, k).volume() -
               (face_flux[dim])(i + (dim == 0 ? 1 : 0), j + (dim == 1 ? 1 : 0),
                                k + (dim == 2 ? 1 : 0))
                   .volume()) /
                  (cell_volume) +
              a_dt * cc(i, j, k) *
                  ((U_face[dim])(i + (dim == 0 ? 1 : 0), j + (dim == 1 ? 1 : 0),
                                 k + (dim == 2 ? 1 : 0)) -
                   (U_face[dim])(i, j, k)) /
                  mesh.dx();
        }
      }
    }
    a_liquid_volume_fraction->updateBorder();
  }

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        if ((*a_liquid_volume_fraction)(i, j, k) < 0.0) {
          (*a_liquid_volume_fraction)(i, j, k) = 0.0;
        } else if ((*a_liquid_volume_fraction)(i, j, k) > 1.0) {
          (*a_liquid_volume_fraction)(i, j, k) = 1.0;
        }
      }
    }
  }

  a_liquid_volume_fraction->updateBorder();
}

void FullLagrangian::advectVOF(
    const std::string& a_reconstruction_method, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid,
    Data<double>* a_liquid_volume_fraction, Data<IRL::Pt>* a_liquid_centroid,
    Data<IRL::Pt>* a_gas_centroid) {
  const BasicMesh& mesh = a_liquid_volume_fraction->getMesh();
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int split_proc = static_cast<int>(std::cbrt(static_cast<double>(size)));
  int imin, nx, jmin, ny, kmin, nz;
  const int NX = mesh.getNx(), NY = mesh.getNy(), NZ = mesh.getNz();

  if (size > 1) {
    if (size == split_proc * split_proc * split_proc) {
      for (int i = 0; i < split_proc; i++) {
        for (int j = 0; j < split_proc; j++) {
          for (int k = 0; k < split_proc; k++) {
            if (i + split_proc * j + split_proc * split_proc * k == rank) {
              imin = i * (NX / split_proc);
              nx = std::min((i + 1) * (NX / split_proc), NX) - imin;
              jmin = j * (NY / split_proc);
              ny = std::min((j + 1) * (NY / split_proc), NY) - jmin;
              kmin = k * (NZ / split_proc);
              nz = std::min((k + 1) * (NZ / split_proc), NZ) - kmin;
            }
          }
        }
      }
    } else {
      imin = rank * (NX / size);
      nx = std::min((rank + 1) * (NX / size), NX) - imin;
      jmin = 0;
      ny = NZ;
      kmin = 0;
      nz = NZ;
    }
  } else {
    imin = 0, nx = NX, jmin = 0, ny = NY, kmin = 0, nz = NZ;
  }

  // For now, naively advect everywhere in domain
  // for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
  //   for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
  //     for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
  for (int i = imin; i < imin + nx; ++i) {
    for (int j = jmin; j < jmin + ny; ++j) {
      for (int k = kmin; k < kmin + nz; ++k) {
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        // Get the back project CappedDodecahedron.
        IRL::Dodecahedron transported_cell;
        for (IRL::UnsignedIndex_t n = 0; n < 8; ++n) {
          transported_cell[n] =
              back_project_vertex(cell[n], -a_dt, a_U, a_V, a_W);
        }

        // Now perform the actual cutting.
        IRL::SeparatedMoments<IRL::VolumeMoments> cell_volume_moments =
            IRL::getNormalizedVolumeMoments<
                IRL::SeparatedMoments<IRL::VolumeMoments>>(
                transported_cell, (*a_link_localized_paraboloid)(i, j, k));
        const double cell_volume = cell.calculateVolume();
        (*a_liquid_volume_fraction)(i, j, k) =
            cell_volume_moments[0].volume() / (cell_volume);
        (*a_liquid_centroid)(i, j, k) = cell_volume_moments[0].centroid();
        (*a_gas_centroid)(i, j, k) = cell_volume_moments[1].centroid();

        if ((*a_liquid_volume_fraction)(i, j, k) < 0.0) {
          (*a_liquid_volume_fraction)(i, j, k) = 0.0;
        } else if ((*a_liquid_volume_fraction)(i, j, k) > 1.0) {
          (*a_liquid_volume_fraction)(i, j, k) = 1.0;
        }

        // if ((*a_liquid_volume_fraction)(i, j, k) <
        //     IRL::global_constants::VF_LOW) {
        //   // (*a_liquid_volume_fraction)(i, j, k) = 0.0;
        //   (*a_liquid_centroid)(i, j, k) = cell.calculateCentroid();
        //   (*a_gas_centroid)(i, j, k) = (*a_liquid_centroid)(i, j, k);
        // } else if ((*a_liquid_volume_fraction)(i, j, k) >
        //            IRL::global_constants::VF_HIGH) {
        //   // (*a_liquid_volume_fraction)(i, j, k) = 1.0;
        //   (*a_liquid_centroid)(i, j, k) = cell.calculateCentroid();
        //   (*a_gas_centroid)(i, j, k) = (*a_liquid_centroid)(i, j, k);
        // } else {
        //   (*a_liquid_centroid)(i, j, k) = back_project_vertex(
        //       (*a_liquid_centroid)(i, j, k), a_dt, a_U, a_V, a_W);
        //   (*a_gas_centroid)(i, j, k) = back_project_vertex(
        //       (*a_gas_centroid)(i, j, k), a_dt, a_U, a_V, a_W);
        // }
      }
    }
  }

  int vector_size = mesh.size();
  std::vector<double> vfrac_local, vfrac;
  vfrac_local.resize(vector_size);
  vfrac.resize(vector_size);
  std::fill(vfrac_local.begin(), vfrac_local.end(), 0.0);
  std::fill(vfrac.begin(), vfrac.end(), 0.0);

  for (int i = imin; i < imin + nx; ++i) {
    for (int j = jmin; j < jmin + ny; ++j) {
      for (int k = kmin; k < kmin + nz; ++k) {
        int address = (i + mesh.getNgc()) * mesh.getNyo() * mesh.getNzo() +
                      (j + mesh.getNgc()) * mesh.getNzo() + (k + mesh.getNgc());
        vfrac_local[address] = (*a_liquid_volume_fraction)(i, j, k);
      }
    }
  }

  MPI_Allreduce(vfrac_local.data(), vfrac.data(), vector_size, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        int address = (i + mesh.getNgc()) * mesh.getNyo() * mesh.getNzo() +
                      (j + mesh.getNgc()) * mesh.getNzo() + (k + mesh.getNgc());
        (*a_liquid_volume_fraction)(i, j, k) = vfrac[address];
      }
    }
  }

  a_liquid_volume_fraction->updateBorder();

  // Technically wrong below, need to move to new reference frame for periodic
  a_liquid_centroid->updateBorder();
  a_gas_centroid->updateBorder();
  correctCentroidLocation(a_liquid_centroid, a_gas_centroid);
}

void SemiLagrangian::advectVOF(
    const std::string& a_reconstruction_method, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid,
    Data<double>* a_liquid_volume_fraction, Data<IRL::Pt>* a_liquid_centroid,
    Data<IRL::Pt>* a_gas_centroid) {
  const BasicMesh& mesh = a_liquid_volume_fraction->getMesh();
  Data<IRL::SeparatedMoments<IRL::VolumeMoments>> face_flux[3] = {
      Data<IRL::SeparatedMoments<IRL::VolumeMoments>>(&mesh),
      Data<IRL::SeparatedMoments<IRL::VolumeMoments>>(&mesh),
      Data<IRL::SeparatedMoments<IRL::VolumeMoments>>(&mesh)};

  resetCentroids(*a_link_localized_paraboloid, a_liquid_centroid,
                 a_gas_centroid);

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int split_proc = static_cast<int>(std::cbrt(static_cast<double>(size)));
  int imin, nx, jmin, ny, kmin, nz;
  const int NX = mesh.getNx(), NY = mesh.getNy(), NZ = mesh.getNz();

  if (size > 1) {
    if (size == split_proc * split_proc * split_proc) {
      for (int i = 0; i < split_proc; i++) {
        for (int j = 0; j < split_proc; j++) {
          for (int k = 0; k < split_proc; k++) {
            if (i + split_proc * j + split_proc * split_proc * k == rank) {
              imin = i * (NX / split_proc);
              nx = std::min((i + 1) * (NX / split_proc), NX) - imin;
              jmin = j * (NY / split_proc);
              ny = std::min((j + 1) * (NY / split_proc), NY) - jmin;
              kmin = k * (NZ / split_proc);
              nz = std::min((k + 1) * (NZ / split_proc), NZ) - kmin;
            }
          }
        }
      }
    } else {
      imin = rank * (NX / size);
      nx = std::min((rank + 1) * (NX / size), NX) - imin;
      jmin = 0;
      ny = NZ;
      kmin = 0;
      nz = NZ;
    }
  } else {
    imin = 0, nx = NX, jmin = 0, ny = NY, kmin = 0, nz = NZ;
  }

  // For now, naively advect everywhere in domain
  // for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
  //   for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
  //     for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
  for (int i = imin; i < imin + nx; ++i) {
    for (int j = jmin; j < jmin + ny; ++j) {
      for (int k = kmin; k < kmin + nz; ++k) {
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        // Get the back project CappedDodecahedron.
        IRL::Dodecahedron transported_cell;
        for (IRL::UnsignedIndex_t n = 0; n < 8; ++n) {
          transported_cell[n] =
              back_project_vertex(cell[n], -a_dt, a_U, a_V, a_W);
        }
        IRL::Dodecahedron face_cell[3];
        face_cell[0] = IRL::Dodecahedron(
            {cell[7], cell[4], cell[5], cell[6], transported_cell[7],
             transported_cell[4], transported_cell[5], transported_cell[6]});
        face_cell[1] = IRL::Dodecahedron(
            {cell[0], cell[4], cell[7], cell[3], transported_cell[0],
             transported_cell[4], transported_cell[7], transported_cell[3]});
        face_cell[2] = IRL::Dodecahedron(
            {cell[5], cell[4], cell[0], cell[1], transported_cell[5],
             transported_cell[4], transported_cell[0], transported_cell[1]});

        for (int dim = 0; dim < 3; ++dim) {
          // Store face flux
          (face_flux[dim])(i, j, k) =
              IRL::getVolumeMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(
                  face_cell[dim], (*a_link_localized_paraboloid)(i, j, k));
        }
      }
    }
  }

  int vector_size = mesh.size();
  std::vector<double> face_flux_x_local, face_flux_y_local, face_flux_z_local;
  std::vector<double> face_flux_x, face_flux_y, face_flux_z;
  face_flux_x_local.resize(vector_size);
  face_flux_y_local.resize(vector_size);
  face_flux_z_local.resize(vector_size);
  face_flux_x.resize(vector_size);
  face_flux_y.resize(vector_size);
  face_flux_z.resize(vector_size);
  std::fill(face_flux_x_local.begin(), face_flux_x_local.end(), 0.0);
  std::fill(face_flux_y_local.begin(), face_flux_y_local.end(), 0.0);
  std::fill(face_flux_z_local.begin(), face_flux_z_local.end(), 0.0);
  std::fill(face_flux_x.begin(), face_flux_x.end(), 0.0);
  std::fill(face_flux_y.begin(), face_flux_y.end(), 0.0);
  std::fill(face_flux_z.begin(), face_flux_z.end(), 0.0);

  for (int i = imin; i < imin + nx; ++i) {
    for (int j = jmin; j < jmin + ny; ++j) {
      for (int k = kmin; k < kmin + nz; ++k) {
        int address = (i + mesh.getNgc()) * mesh.getNyo() * mesh.getNzo() +
                      (j + mesh.getNgc()) * mesh.getNzo() + (k + mesh.getNgc());
        face_flux_x_local[address] = (face_flux[0](i, j, k))[0].volume();
        face_flux_y_local[address] = (face_flux[1](i, j, k))[0].volume();
        face_flux_z_local[address] = (face_flux[2](i, j, k))[0].volume();
      }
    }
  }

  MPI_Allreduce(face_flux_x_local.data(), face_flux_x.data(), vector_size,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(face_flux_y_local.data(), face_flux_y.data(), vector_size,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(face_flux_z_local.data(), face_flux_z.data(), vector_size,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        int address = (i + mesh.getNgc()) * mesh.getNyo() * mesh.getNzo() +
                      (j + mesh.getNgc()) * mesh.getNzo() + (k + mesh.getNgc());
        (face_flux[0](i, j, k))[0].volume() = face_flux_x[address];
        (face_flux[1](i, j, k))[0].volume() = face_flux_y[address];
        (face_flux[2](i, j, k))[0].volume() = face_flux_z[address];
      }
    }
  }

  std::fill(face_flux_x_local.begin(), face_flux_x_local.end(), 0.0);
  std::fill(face_flux_y_local.begin(), face_flux_y_local.end(), 0.0);
  std::fill(face_flux_z_local.begin(), face_flux_z_local.end(), 0.0);
  std::fill(face_flux_x.begin(), face_flux_x.end(), 0.0);
  std::fill(face_flux_y.begin(), face_flux_y.end(), 0.0);
  std::fill(face_flux_z.begin(), face_flux_z.end(), 0.0);

  for (int i = imin; i < imin + nx; ++i) {
    for (int j = jmin; j < jmin + ny; ++j) {
      for (int k = kmin; k < kmin + nz; ++k) {
        int address = (i + mesh.getNgc()) * mesh.getNyo() * mesh.getNzo() +
                      (j + mesh.getNgc()) * mesh.getNzo() + (k + mesh.getNgc());
        face_flux_x_local[address] = (face_flux[0](i, j, k))[1].volume();
        face_flux_y_local[address] = (face_flux[1](i, j, k))[1].volume();
        face_flux_z_local[address] = (face_flux[2](i, j, k))[1].volume();
      }
    }
  }

  MPI_Allreduce(face_flux_x_local.data(), face_flux_x.data(), vector_size,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(face_flux_y_local.data(), face_flux_y.data(), vector_size,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(face_flux_z_local.data(), face_flux_z.data(), vector_size,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        int address = (i + mesh.getNgc()) * mesh.getNyo() * mesh.getNzo() +
                      (j + mesh.getNgc()) * mesh.getNzo() + (k + mesh.getNgc());
        (face_flux[0](i, j, k))[1].volume() = face_flux_x[address];
        (face_flux[1](i, j, k))[1].volume() = face_flux_y[address];
        (face_flux[2](i, j, k))[1].volume() = face_flux_z[address];
      }
    }
  }

  face_flux[0].updateBorder();
  face_flux[1].updateBorder();
  face_flux[2].updateBorder();

  // Now calculate VOF from the face fluxes.
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        const double cell_volume = cell.calculateVolume();
        const double previous_liquid_volume_fraction =
            (*a_liquid_volume_fraction)(i, j, k);

        // Update VOF
        (*a_liquid_volume_fraction)(i, j, k) =
            (previous_liquid_volume_fraction * cell_volume +
             (face_flux[0])(i, j, k)[0].volume() -
             (face_flux[0])(i + 1, j, k)[0].volume() +
             (face_flux[1])(i, j, k)[0].volume() -
             (face_flux[1])(i, j + 1, k)[0].volume() +
             (face_flux[2])(i, j, k)[0].volume() -
             (face_flux[2])(i, j, k + 1)[0].volume()) /
            (cell_volume + (face_flux[0])(i, j, k)[0].volume() -
             (face_flux[0])(i + 1, j, k)[0].volume() +
             (face_flux[1])(i, j, k)[0].volume() -
             (face_flux[1])(i, j + 1, k)[0].volume() +
             (face_flux[2])(i, j, k)[0].volume() -
             (face_flux[2])(i, j, k + 1)[0].volume() +
             (face_flux[0])(i, j, k)[1].volume() -
             (face_flux[0])(i + 1, j, k)[1].volume() +
             (face_flux[1])(i, j, k)[1].volume() -
             (face_flux[1])(i, j + 1, k)[1].volume() +
             (face_flux[2])(i, j, k)[1].volume() -
             (face_flux[2])(i, j, k + 1)[1].volume());

        if ((*a_liquid_volume_fraction)(i, j, k) < 0.0) {
          (*a_liquid_volume_fraction)(i, j, k) = 0.0;
        } else if ((*a_liquid_volume_fraction)(i, j, k) > 1.0) {
          (*a_liquid_volume_fraction)(i, j, k) = 1.0;
        }

        if ((*a_liquid_volume_fraction)(i, j, k) <
            IRL::global_constants::VF_LOW) {
          // (*a_liquid_volume_fraction)(i, j, k) = 0.0;
          (*a_liquid_centroid)(i, j, k) = cell.calculateCentroid();
          (*a_gas_centroid)(i, j, k) = cell.calculateCentroid();
        } else if ((*a_liquid_volume_fraction)(i, j, k) >
                   IRL::global_constants::VF_HIGH) {
          // (*a_liquid_volume_fraction)(i, j, k) = 1.0;
          (*a_liquid_centroid)(i, j, k) = cell.calculateCentroid();
          (*a_gas_centroid)(i, j, k) = cell.calculateCentroid();
        } else {
          // Update liquid centroid, .centroid() is un-normalized
          (*a_liquid_centroid)(i, j, k) =
              IRL::Pt(previous_liquid_volume_fraction * cell_volume *
                          (*a_liquid_centroid)(i, j, k) +
                      (face_flux[0])(i, j, k)[0].centroid() -
                      (face_flux[0])(i + 1, j, k)[0].centroid() +
                      (face_flux[1])(i, j, k)[0].centroid() -
                      (face_flux[1])(i, j + 1, k)[0].centroid() +
                      (face_flux[2])(i, j, k)[0].centroid() -
                      (face_flux[2])(i, j, k + 1)[0].centroid()) /
              (previous_liquid_volume_fraction * cell_volume +
               (face_flux[0])(i, j, k)[0].volume() -
               (face_flux[0])(i + 1, j, k)[0].volume() +
               (face_flux[1])(i, j, k)[0].volume() -
               (face_flux[1])(i, j + 1, k)[0].volume() +
               (face_flux[2])(i, j, k)[0].volume() -
               (face_flux[2])(i, j, k + 1)[0].volume());

          // Update gas centroid, .centroid() is un-normalized
          (*a_gas_centroid)(i, j, k) =
              IRL::Pt((1.0 - previous_liquid_volume_fraction) * cell_volume *
                          (*a_gas_centroid)(i, j, k) +
                      (face_flux[0])(i, j, k)[1].centroid() -
                      (face_flux[0])(i + 1, j, k)[1].centroid() +
                      (face_flux[1])(i, j, k)[1].centroid() -
                      (face_flux[1])(i, j + 1, k)[1].centroid() +
                      (face_flux[2])(i, j, k)[1].centroid() -
                      (face_flux[2])(i, j, k + 1)[1].centroid()) /
              ((1.0 - previous_liquid_volume_fraction) * cell_volume +
               (face_flux[0])(i, j, k)[1].volume() -
               (face_flux[0])(i + 1, j, k)[1].volume() +
               (face_flux[1])(i, j, k)[1].volume() -
               (face_flux[1])(i, j + 1, k)[1].volume() +
               (face_flux[2])(i, j, k)[1].volume() -
               (face_flux[2])(i, j, k + 1)[1].volume());

          (*a_liquid_centroid)(i, j, k) = back_project_vertex(
              (*a_liquid_centroid)(i, j, k), a_dt, a_U, a_V, a_W);
          (*a_gas_centroid)(i, j, k) = back_project_vertex(
              (*a_gas_centroid)(i, j, k), a_dt, a_U, a_V, a_W);
        }
      }
    }
  }
  a_liquid_volume_fraction->updateBorder();
  // Technically wrong below, need to move to new reference frame for periodic
  a_liquid_centroid->updateBorder();
  a_gas_centroid->updateBorder();
  correctCentroidLocation(a_liquid_centroid, a_gas_centroid);
}

void SemiLagrangianCorrected::advectVOF(
    const std::string& a_reconstruction_method, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid,
    Data<double>* a_liquid_volume_fraction, Data<IRL::Pt>* a_liquid_centroid,
    Data<IRL::Pt>* a_gas_centroid) {
  const BasicMesh& mesh = a_liquid_volume_fraction->getMesh();

  resetCentroids(*a_link_localized_paraboloid, a_liquid_centroid,
                 a_gas_centroid);

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int split_proc = static_cast<int>(std::cbrt(static_cast<double>(size)));
  int imin, nx, jmin, ny, kmin, nz;
  const int NX = mesh.getNx(), NY = mesh.getNy(), NZ = mesh.getNz();

  if (size > 1) {
    if (size == split_proc * split_proc * split_proc) {
      for (int i = 0; i < split_proc; i++) {
        for (int j = 0; j < split_proc; j++) {
          for (int k = 0; k < split_proc; k++) {
            if (i + split_proc * j + split_proc * split_proc * k == rank) {
              imin = i * (NX / split_proc);
              nx = std::min((i + 1) * (NX / split_proc), NX) - imin;
              jmin = j * (NY / split_proc);
              ny = std::min((j + 1) * (NY / split_proc), NY) - jmin;
              kmin = k * (NZ / split_proc);
              nz = std::min((k + 1) * (NZ / split_proc), NZ) - kmin;
            }
          }
        }
      }
    } else {
      imin = rank * (NX / size);
      nx = std::min((rank + 1) * (NX / size), NX) - imin;
      jmin = 0;
      ny = NZ;
      kmin = 0;
      nz = NZ;
    }
  } else {
    imin = 0, nx = NX, jmin = 0, ny = NY, kmin = 0, nz = NZ;
  }

  Data<double> U_face(&mesh);
  Data<double> V_face(&mesh);
  Data<double> W_face(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax() + 1; ++k) {
        // for (int i = imin; i < imin + nx + 1; ++i) {
        //   for (int j = jmin; j < jmin + ny + 1; ++j) {
        //     for (int k = kmin; k < kmin + nz + 1; ++k) {
        U_face(i, j, k) = 0.5 * (a_U(i, j, k) + a_U(i - 1, j, k));
        V_face(i, j, k) = 0.5 * (a_V(i, j, k) + a_V(i, j - 1, k));
        W_face(i, j, k) = 0.5 * (a_W(i, j, k) + a_W(i, j, k - 1));
      }
    }
  }

  // Allocate storage for face fluxes
  Data<IRL::SeparatedMoments<IRL::VolumeMoments>> face_flux[3] = {
      Data<IRL::SeparatedMoments<IRL::VolumeMoments>>(&mesh),
      Data<IRL::SeparatedMoments<IRL::VolumeMoments>>(&mesh),
      Data<IRL::SeparatedMoments<IRL::VolumeMoments>>(&mesh)};

  // For now, naively advect everywhere in domain
  // for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
  //   for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
  //     for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
  for (int i = imin; i < imin + nx; ++i) {
    for (int j = jmin; j < jmin + ny; ++j) {
      for (int k = kmin; k < kmin + nz; ++k) {
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        // Get the back projected Dodecahedron.
        IRL::Dodecahedron transported_cell;
        for (IRL::UnsignedIndex_t n = 0; n < 8; ++n) {
          transported_cell[n] =
              back_project_vertex(cell[n], -a_dt, a_U, a_V, a_W);
        }
        // Initialize face flux hexahedra
        IRL::CappedDodecahedron face_cell[3];
        IRL::Pt face_center_pt = 0.25 * (cell[7] + cell[4] + cell[5] + cell[6]);
        IRL::Pt correction_pt =
            back_project_vertex(face_center_pt, -a_dt, a_U, a_V, a_W);
        face_cell[0] = IRL::CappedDodecahedron(
            {cell[7], cell[4], cell[5], cell[6], transported_cell[7],
             transported_cell[4], transported_cell[5], transported_cell[6],
             correction_pt});
        face_center_pt = 0.25 * (cell[0] + cell[4] + cell[7] + cell[3]);
        correction_pt =
            back_project_vertex(face_center_pt, -a_dt, a_U, a_V, a_W);
        face_cell[1] = IRL::CappedDodecahedron(
            {cell[0], cell[4], cell[7], cell[3], transported_cell[0],
             transported_cell[4], transported_cell[7], transported_cell[3],
             correction_pt});
        face_center_pt = 0.25 * (cell[5] + cell[4] + cell[0] + cell[1]);
        correction_pt =
            back_project_vertex(face_center_pt, -a_dt, a_U, a_V, a_W);
        face_cell[2] = IRL::CappedDodecahedron(
            {cell[5], cell[4], cell[0], cell[1], transported_cell[5],
             transported_cell[4], transported_cell[0], transported_cell[1],
             correction_pt});
        // Tack on the corrective 9th vertex for each face hexahedra
        face_cell[0].adjustCapToMatchVolume(a_dt * U_face(i, j, k) * mesh.dy() *
                                            mesh.dz());
        face_cell[1].adjustCapToMatchVolume(a_dt * V_face(i, j, k) * mesh.dx() *
                                            mesh.dz());
        face_cell[2].adjustCapToMatchVolume(a_dt * W_face(i, j, k) * mesh.dx() *
                                            mesh.dy());

        for (int dim = 0; dim < 3; ++dim) {
          // Store face flux
          (face_flux[dim])(i, j, k) =
              IRL::getVolumeMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(
                  face_cell[dim], (*a_link_localized_paraboloid)(i, j, k));
        }
      }
    }
  }

  int vector_size = mesh.size();
  std::vector<double> face_flux_x_local, face_flux_y_local, face_flux_z_local;
  std::vector<double> face_flux_x, face_flux_y, face_flux_z;
  face_flux_x_local.resize(vector_size);
  face_flux_y_local.resize(vector_size);
  face_flux_z_local.resize(vector_size);
  face_flux_x.resize(vector_size);
  face_flux_y.resize(vector_size);
  face_flux_z.resize(vector_size);
  std::fill(face_flux_x_local.begin(), face_flux_x_local.end(), 0.0);
  std::fill(face_flux_y_local.begin(), face_flux_y_local.end(), 0.0);
  std::fill(face_flux_z_local.begin(), face_flux_z_local.end(), 0.0);
  std::fill(face_flux_x.begin(), face_flux_x.end(), 0.0);
  std::fill(face_flux_y.begin(), face_flux_y.end(), 0.0);
  std::fill(face_flux_z.begin(), face_flux_z.end(), 0.0);

  for (int i = imin; i < imin + nx; ++i) {
    for (int j = jmin; j < jmin + ny; ++j) {
      for (int k = kmin; k < kmin + nz; ++k) {
        int address = (i + mesh.getNgc()) * mesh.getNyo() * mesh.getNzo() +
                      (j + mesh.getNgc()) * mesh.getNzo() + (k + mesh.getNgc());
        face_flux_x_local[address] = (face_flux[0](i, j, k))[0].volume();
        face_flux_y_local[address] = (face_flux[1](i, j, k))[0].volume();
        face_flux_z_local[address] = (face_flux[2](i, j, k))[0].volume();
      }
    }
  }

  MPI_Allreduce(face_flux_x_local.data(), face_flux_x.data(), vector_size,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(face_flux_y_local.data(), face_flux_y.data(), vector_size,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(face_flux_z_local.data(), face_flux_z.data(), vector_size,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        int address = (i + mesh.getNgc()) * mesh.getNyo() * mesh.getNzo() +
                      (j + mesh.getNgc()) * mesh.getNzo() + (k + mesh.getNgc());
        (face_flux[0](i, j, k))[0].volume() = face_flux_x[address];
        (face_flux[1](i, j, k))[0].volume() = face_flux_y[address];
        (face_flux[2](i, j, k))[0].volume() = face_flux_z[address];
      }
    }
  }

  std::fill(face_flux_x_local.begin(), face_flux_x_local.end(), 0.0);
  std::fill(face_flux_y_local.begin(), face_flux_y_local.end(), 0.0);
  std::fill(face_flux_z_local.begin(), face_flux_z_local.end(), 0.0);
  std::fill(face_flux_x.begin(), face_flux_x.end(), 0.0);
  std::fill(face_flux_y.begin(), face_flux_y.end(), 0.0);
  std::fill(face_flux_z.begin(), face_flux_z.end(), 0.0);

  for (int i = imin; i < imin + nx; ++i) {
    for (int j = jmin; j < jmin + ny; ++j) {
      for (int k = kmin; k < kmin + nz; ++k) {
        int address = (i + mesh.getNgc()) * mesh.getNyo() * mesh.getNzo() +
                      (j + mesh.getNgc()) * mesh.getNzo() + (k + mesh.getNgc());
        face_flux_x_local[address] = (face_flux[0](i, j, k))[1].volume();
        face_flux_y_local[address] = (face_flux[1](i, j, k))[1].volume();
        face_flux_z_local[address] = (face_flux[2](i, j, k))[1].volume();
      }
    }
  }

  MPI_Allreduce(face_flux_x_local.data(), face_flux_x.data(), vector_size,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(face_flux_y_local.data(), face_flux_y.data(), vector_size,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(face_flux_z_local.data(), face_flux_z.data(), vector_size,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        int address = (i + mesh.getNgc()) * mesh.getNyo() * mesh.getNzo() +
                      (j + mesh.getNgc()) * mesh.getNzo() + (k + mesh.getNgc());
        (face_flux[0](i, j, k))[1].volume() = face_flux_x[address];
        (face_flux[1](i, j, k))[1].volume() = face_flux_y[address];
        (face_flux[2](i, j, k))[1].volume() = face_flux_z[address];
      }
    }
  }

  face_flux[0].updateBorder();
  face_flux[1].updateBorder();
  face_flux[2].updateBorder();

  // Now calculate VOF from the face fluxes.
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        const double cell_volume = cell.calculateVolume();
        const double previous_liquid_volume_fraction =
            (*a_liquid_volume_fraction)(i, j, k);

        // Update VOF
        (*a_liquid_volume_fraction)(i, j, k) =
            (previous_liquid_volume_fraction * cell_volume +
             (face_flux[0])(i, j, k)[0].volume() -
             (face_flux[0])(i + 1, j, k)[0].volume() +
             (face_flux[1])(i, j, k)[0].volume() -
             (face_flux[1])(i, j + 1, k)[0].volume() +
             (face_flux[2])(i, j, k)[0].volume() -
             (face_flux[2])(i, j, k + 1)[0].volume()) /
            (cell_volume + (face_flux[0])(i, j, k)[0].volume() -
             (face_flux[0])(i + 1, j, k)[0].volume() +
             (face_flux[1])(i, j, k)[0].volume() -
             (face_flux[1])(i, j + 1, k)[0].volume() +
             (face_flux[2])(i, j, k)[0].volume() -
             (face_flux[2])(i, j, k + 1)[0].volume() +
             (face_flux[0])(i, j, k)[1].volume() -
             (face_flux[0])(i + 1, j, k)[1].volume() +
             (face_flux[1])(i, j, k)[1].volume() -
             (face_flux[1])(i, j + 1, k)[1].volume() +
             (face_flux[2])(i, j, k)[1].volume() -
             (face_flux[2])(i, j, k + 1)[1].volume());

        if ((*a_liquid_volume_fraction)(i, j, k) < 0.0) {
          (*a_liquid_volume_fraction)(i, j, k) = 0.0;
        } else if ((*a_liquid_volume_fraction)(i, j, k) > 1.0) {
          (*a_liquid_volume_fraction)(i, j, k) = 1.0;
        }

        if ((*a_liquid_volume_fraction)(i, j, k) <
            IRL::global_constants::VF_LOW) {
          // (*a_liquid_volume_fraction)(i, j, k) = 0.0;
          (*a_liquid_centroid)(i, j, k) = cell.calculateCentroid();
          (*a_gas_centroid)(i, j, k) = cell.calculateCentroid();
        } else if ((*a_liquid_volume_fraction)(i, j, k) >
                   IRL::global_constants::VF_HIGH) {
          // (*a_liquid_volume_fraction)(i, j, k) = 1.0;
          (*a_liquid_centroid)(i, j, k) = cell.calculateCentroid();
          (*a_gas_centroid)(i, j, k) = cell.calculateCentroid();
        } else {
          // Update liquid centroid, .centroid() is un-normalized
          (*a_liquid_centroid)(i, j, k) =
              IRL::Pt(previous_liquid_volume_fraction * cell_volume *
                          (*a_liquid_centroid)(i, j, k) +
                      (face_flux[0])(i, j, k)[0].centroid() -
                      (face_flux[0])(i + 1, j, k)[0].centroid() +
                      (face_flux[1])(i, j, k)[0].centroid() -
                      (face_flux[1])(i, j + 1, k)[0].centroid() +
                      (face_flux[2])(i, j, k)[0].centroid() -
                      (face_flux[2])(i, j, k + 1)[0].centroid()) /
              (previous_liquid_volume_fraction * cell_volume +
               (face_flux[0])(i, j, k)[0].volume() -
               (face_flux[0])(i + 1, j, k)[0].volume() +
               (face_flux[1])(i, j, k)[0].volume() -
               (face_flux[1])(i, j + 1, k)[0].volume() +
               (face_flux[2])(i, j, k)[0].volume() -
               (face_flux[2])(i, j, k + 1)[0].volume());

          // Update gas centroid, .centroid() is un-normalized
          (*a_gas_centroid)(i, j, k) =
              IRL::Pt((1.0 - previous_liquid_volume_fraction) * cell_volume *
                          (*a_gas_centroid)(i, j, k) +
                      (face_flux[0])(i, j, k)[1].centroid() -
                      (face_flux[0])(i + 1, j, k)[1].centroid() +
                      (face_flux[1])(i, j, k)[1].centroid() -
                      (face_flux[1])(i, j + 1, k)[1].centroid() +
                      (face_flux[2])(i, j, k)[1].centroid() -
                      (face_flux[2])(i, j, k + 1)[1].centroid()) /
              ((1.0 - previous_liquid_volume_fraction) * cell_volume +
               (face_flux[0])(i, j, k)[1].volume() -
               (face_flux[0])(i + 1, j, k)[1].volume() +
               (face_flux[1])(i, j, k)[1].volume() -
               (face_flux[1])(i, j + 1, k)[1].volume() +
               (face_flux[2])(i, j, k)[1].volume() -
               (face_flux[2])(i, j, k + 1)[1].volume());

          (*a_liquid_centroid)(i, j, k) = back_project_vertex(
              (*a_liquid_centroid)(i, j, k), a_dt, a_U, a_V, a_W);
          (*a_gas_centroid)(i, j, k) = back_project_vertex(
              (*a_gas_centroid)(i, j, k), a_dt, a_U, a_V, a_W);
        }
      }
    }
  }
  a_liquid_volume_fraction->updateBorder();
  a_liquid_centroid->updateBorder();
  a_gas_centroid->updateBorder();
  correctCentroidLocation(a_liquid_centroid, a_gas_centroid);
}

void correctCentroidLocation(Data<IRL::Pt>* a_liquid_centroid,
                             Data<IRL::Pt>* a_gas_centroid) {
  const BasicMesh& mesh = (*a_liquid_centroid).getMesh();
  // Fix distance to recreate volume fraction

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_liquid_centroid)(i, j, k) =
            (*a_liquid_centroid)(i, j, k)[0] - mesh.lx();
        (*a_gas_centroid)(i, j, k) = (*a_gas_centroid)(i, j, k)[0] - mesh.lx();
      }
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_liquid_centroid)(i, j, k) =
            (*a_liquid_centroid)(i, j, k)[0] + mesh.lx();
        (*a_gas_centroid)(i, j, k) = (*a_gas_centroid)(i, j, k)[0] + mesh.lx();
      }
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_liquid_centroid)(i, j, k) =
            (*a_liquid_centroid)(i, j, k)[1] - mesh.ly();
        (*a_gas_centroid)(i, j, k) = (*a_gas_centroid)(i, j, k)[1] - mesh.ly();
      }
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_liquid_centroid)(i, j, k) =
            (*a_liquid_centroid)(i, j, k)[1] + mesh.ly();
        (*a_gas_centroid)(i, j, k) = (*a_gas_centroid)(i, j, k)[1] + mesh.ly();
      }
    }
  }

  // z- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k < mesh.kmin(); ++k) {
        (*a_liquid_centroid)(i, j, k) =
            (*a_liquid_centroid)(i, j, k)[2] - mesh.lz();
        (*a_gas_centroid)(i, j, k) = (*a_gas_centroid)(i, j, k)[2] - mesh.lz();
      }
    }
  }

  // z+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmax() + 1; k <= mesh.kmaxo(); ++k) {
        (*a_liquid_centroid)(i, j, k) =
            (*a_liquid_centroid)(i, j, k)[2] + mesh.lz();
        (*a_gas_centroid)(i, j, k) = (*a_gas_centroid)(i, j, k)[2] + mesh.lz();
      }
    }
  }
}
