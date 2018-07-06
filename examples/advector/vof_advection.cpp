// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <iostream>

#include "examples/advector/vof_advection.h"

#include "src/generic_cutting/generic_cutting.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"

void resetCentroids(
    const Data<IRL::LocalizedSeparatorLink>& a_link_localized_separator,
    Data<IRL::Pt>* a_liquid_centroid, Data<IRL::Pt>* a_gas_centroid) {
  const BasicMesh& mesh = a_link_localized_separator.getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        auto moments = IRL::getNormalizedVolumeMoments<
            IRL::SeparatedMoments<IRL::VolumeMoments>>(
            cell, a_link_localized_separator(i, j, k).getNextReconstruction());
        (*a_liquid_centroid)(i, j, k) = moments[0].centroid();
        (*a_gas_centroid)(i, j, k) = moments[1].centroid();
      }
    }
  }
}

void connectMesh(
    const BasicMesh& a_mesh,
    Data<IRL::LocalizedSeparatorLink>* a_link_localized_separator) {
  IRL::LocalizedSeparatorLink* neighbor_ptr;
  // Provide mesh connectivity information.
  IRL::UnsignedIndex_t unique_id = 0;

  for (int i = a_mesh.imino(); i <= a_mesh.imaxo(); ++i) {
    for (int j = a_mesh.jmino(); j <= a_mesh.jmaxo(); ++j) {
      for (int k = a_mesh.kmino(); k <= a_mesh.kmaxo(); ++k) {
        (*a_link_localized_separator)(i, j, k).setId(unique_id);
        neighbor_ptr = i - 1 < a_mesh.imino()
                           ? nullptr
                           : &(*a_link_localized_separator)(i - 1, j, k);
        (*a_link_localized_separator)(i, j, k).setEdgeConnectivity(
            0, neighbor_ptr);
        neighbor_ptr = i + 1 > a_mesh.imaxo()
                           ? nullptr
                           : &(*a_link_localized_separator)(i + 1, j, k);
        (*a_link_localized_separator)(i, j, k).setEdgeConnectivity(
            1, neighbor_ptr);
        neighbor_ptr = j - 1 < a_mesh.jmino()
                           ? nullptr
                           : &(*a_link_localized_separator)(i, j - 1, k);
        (*a_link_localized_separator)(i, j, k).setEdgeConnectivity(
            2, neighbor_ptr);
        neighbor_ptr = j + 1 > a_mesh.jmaxo()
                           ? nullptr
                           : &(*a_link_localized_separator)(i, j + 1, k);
        (*a_link_localized_separator)(i, j, k).setEdgeConnectivity(
            3, neighbor_ptr);
        neighbor_ptr = k - 1 < a_mesh.kmino()
                           ? nullptr
                           : &(*a_link_localized_separator)(i, j, k - 1);
        (*a_link_localized_separator)(i, j, k).setEdgeConnectivity(
            4, neighbor_ptr);
        neighbor_ptr = k + 1 > a_mesh.kmaxo()
                           ? nullptr
                           : &(*a_link_localized_separator)(i, j, k + 1);
        (*a_link_localized_separator)(i, j, k).setEdgeConnectivity(
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

void advectVOF(const std::string& a_advection_method, const double a_dt,
               const Data<double>& a_U, const Data<double>& a_V,
               const Data<double>& a_W,
               Data<IRL::LocalizedSeparatorLink>* a_link_localized_separator,
               Data<double>* a_liquid_volume_fraction,
               Data<IRL::Pt>* a_liquid_centroid,
               Data<IRL::Pt>* a_gas_centroid) {
  if (a_advection_method == "FullLagrangian") {
    FullLagrangian::advectVOF(a_dt, a_U, a_V, a_W, a_link_localized_separator,
                              a_liquid_volume_fraction, a_liquid_centroid,
                              a_gas_centroid);
  } else if (a_advection_method == "SemiLagrangian") {
    SemiLagrangian::advectVOF(a_dt, a_U, a_V, a_W, a_link_localized_separator,
                              a_liquid_volume_fraction, a_liquid_centroid,
                              a_gas_centroid);
  } else if (a_advection_method == "SemiLagrangianCorrected") {
    SemiLagrangianCorrected::advectVOF(
        a_dt, a_U, a_V, a_W, a_link_localized_separator,
        a_liquid_volume_fraction, a_liquid_centroid, a_gas_centroid);
  } else {
    std::cout << "Unknown advection method of : " << a_advection_method << '\n';
    std::cout << "Value entries are: FullLagrangian, SemiLagrangian, "
                 "SemiLagrangianCorrected. \n";
    std::exit(-1);
  }
}

void FullLagrangian::advectVOF(
    const double a_dt, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W,
    Data<IRL::LocalizedSeparatorLink>* a_link_localized_separator,
    Data<double>* a_liquid_volume_fraction, Data<IRL::Pt>* a_liquid_centroid,
    Data<IRL::Pt>* a_gas_centroid) {
  const BasicMesh& mesh = a_liquid_volume_fraction->getMesh();
  // For now, naively advect everywhere in domain
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        // Get the back project CappedDodecahedron.
        IRL::Dodecahedron transported_cell;
        for (IRL::UnsignedIndex_t n = 0; n < 8; ++n) {
          transported_cell[n] =
              back_project_vertex(cell[n], -a_dt, a_U, a_V, a_W);
        }
        // Determine the bounding volume that the transported Dodecahedron
        // covers.

        // Now perform the actual cutting.
        IRL::SeparatedMoments<IRL::VolumeMoments> cell_volume_moments =
            IRL::getNormalizedVolumeMoments<
                IRL::SeparatedMoments<IRL::VolumeMoments>>(
                transported_cell, (*a_link_localized_separator)(i, j, k));
        (*a_liquid_volume_fraction)(i, j, k) =
            cell_volume_moments[0].volume() /
            (cell_volume_moments[0].volume() + cell_volume_moments[1].volume());
        (*a_liquid_centroid)(i, j, k) = cell_volume_moments[0].centroid();
        (*a_gas_centroid)(i, j, k) = cell_volume_moments[1].centroid();

        if ((*a_liquid_volume_fraction)(i, j, k) <
            IRL::global_constants::VF_LOW) {
          (*a_liquid_volume_fraction)(i, j, k) = 0.0;
          (*a_liquid_centroid)(i, j, k) = cell.calculateCentroid();
          (*a_gas_centroid)(i, j, k) = (*a_liquid_centroid)(i, j, k);
        } else if ((*a_liquid_volume_fraction)(i, j, k) >
                   IRL::global_constants::VF_HIGH) {
          (*a_liquid_volume_fraction)(i, j, k) = 1.0;
          (*a_liquid_centroid)(i, j, k) = cell.calculateCentroid();
          (*a_gas_centroid)(i, j, k) = (*a_liquid_centroid)(i, j, k);
        } else {
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

void SemiLagrangian::advectVOF(
    const double a_dt, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W,
    Data<IRL::LocalizedSeparatorLink>* a_link_localized_separator,
    Data<double>* a_liquid_volume_fraction, Data<IRL::Pt>* a_liquid_centroid,
    Data<IRL::Pt>* a_gas_centroid) {
  const BasicMesh& mesh = a_liquid_volume_fraction->getMesh();
  Data<IRL::SeparatedMoments<IRL::VolumeMoments>> face_flux[3] = {
      Data<IRL::SeparatedMoments<IRL::VolumeMoments>>(&mesh),
      Data<IRL::SeparatedMoments<IRL::VolumeMoments>>(&mesh),
      Data<IRL::SeparatedMoments<IRL::VolumeMoments>>(&mesh)};

  resetCentroids(*a_link_localized_separator, a_liquid_centroid,
                 a_gas_centroid);

  // For now, naively advect everywhere in domain
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
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
                  face_cell[dim], (*a_link_localized_separator)(i, j, k));
        }
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
        if ((*a_liquid_volume_fraction)(i, j, k) <
            IRL::global_constants::VF_LOW) {
          (*a_liquid_volume_fraction)(i, j, k) = 0.0;
          (*a_liquid_centroid)(i, j, k) = cell.calculateCentroid();
          (*a_gas_centroid)(i, j, k) = cell.calculateCentroid();
        } else if ((*a_liquid_volume_fraction)(i, j, k) >
                   IRL::global_constants::VF_HIGH) {
          (*a_liquid_volume_fraction)(i, j, k) = 1.0;
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
    const double a_dt, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W,
    Data<IRL::LocalizedSeparatorLink>* a_link_localized_separator,
    Data<double>* a_liquid_volume_fraction, Data<IRL::Pt>* a_liquid_centroid,
    Data<IRL::Pt>* a_gas_centroid) {
  const BasicMesh& mesh = a_liquid_volume_fraction->getMesh();

  resetCentroids(*a_link_localized_separator, a_liquid_centroid,
                 a_gas_centroid);

  Data<double> U_face(&mesh);
  Data<double> V_face(&mesh);
  Data<double> W_face(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax() + 1; ++k) {
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
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
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
                  face_cell[dim], (*a_link_localized_separator)(i, j, k));
        }
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
        if ((*a_liquid_volume_fraction)(i, j, k) <
            IRL::global_constants::VF_LOW) {
          (*a_liquid_volume_fraction)(i, j, k) = 0.0;
          (*a_liquid_centroid)(i, j, k) = cell.calculateCentroid();
          (*a_gas_centroid)(i, j, k) = cell.calculateCentroid();
        } else if ((*a_liquid_volume_fraction)(i, j, k) >
                   IRL::global_constants::VF_HIGH) {
          (*a_liquid_volume_fraction)(i, j, k) = 1.0;
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
