!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_cutpolygon.f90
!!
!! This file deals with intersecting polygons
!! and generating polygons corresponding
!! to planar reconstructions.

!> This module contains mappings to the
!! IRL C interface that deal with intersecting planes
!! to generate polygons and creating polygons
!! that are representative of planar reconstructions
!! in given cells.
module f_CutPolygon
  use f_DefinedTypes
  use f_RectCub_class
  use f_Hex_class
  use f_Tet_class
  use f_Poly_class
  use f_DivPoly_class
  use f_PlanarSep_class
  implicit none

  interface getPoly
    module procedure getPoly_RectCub_Poly
    module procedure getPoly_Tet_Poly
    module procedure getPoly_Hex_Poly         
    module procedure getPoly_RectCub_DivPoly
  end interface getPoly

  interface getSA
    module procedure getSA_RectCub
  end interface getSA

  interface
    subroutine F_getPoly_RectCub_Poly(a_rectangular_cuboid, a_planar_separator, a_plane_index, a_polygon) &
    bind(C, name="c_getPoly_RectCub_Poly")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_RectCub) :: a_rectangular_cuboid ! Pointer to RectCub object
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
      integer(C_INT), intent(in) :: a_plane_index ! Plane to get polygon for
      type(c_Poly) :: a_polygon ! Pointer to Polyg object
    end subroutine F_getPoly_RectCub_Poly
  end interface

  interface
    subroutine F_getPoly_Tet_Poly(a_tet, a_planar_separator, a_plane_index, a_polygon) &
    bind(C, name="c_getPoly_Tet_Poly")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Tet) :: a_tet ! Pointer to Tet object
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
      integer(C_INT), intent(in) :: a_plane_index ! Plane to get polygon for
      type(c_Poly) :: a_polygon ! Pointer to Polyg object
    end subroutine F_getPoly_Tet_Poly
  end interface

  interface
    subroutine F_getPoly_Hex_Poly(a_hexahedron, a_planar_separator, a_plane_index, a_polygon) &
    bind(C, name="c_getPoly_Hex_Poly")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Hex) :: a_hexahedron ! Pointer to Hex object
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
      integer(C_INT), intent(in) :: a_plane_index ! Plane to get polygon for
      type(c_Poly) :: a_polygon ! Pointer to Polyg object
    end subroutine F_getPoly_Hex_Poly
  end interface  

  interface
    subroutine F_getPoly_RectCub_DivPoly(a_rectangular_cuboid, a_planar_separator, a_plane_index, a_divided_polygon) &
    bind(C, name="c_getPoly_RectCub_DivPoly")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_RectCub) :: a_rectangular_cuboid ! Pointer to RectCub object
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
      integer(C_INT), intent(in) :: a_plane_index ! Plane to get polygon for
      type(c_DivPoly) :: a_divided_polygon ! Pointer to DivPoly object
    end subroutine F_getPoly_RectCub_DivPoly
  end interface

  interface
    function F_getSA_RectCub(a_rectangular_cuboid, a_planar_separator) result(a_surface_area) &
    bind(C, name="c_getSA_RectCub")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_RectCub) :: a_rectangular_cuboid ! Pointer to RectCub object
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
      real(C_DOUBLE) :: a_surface_area
    end function F_getSA_RectCub
  end interface




contains

  subroutine getPoly_RectCub_Poly(a_rectangular_cuboid, a_planar_separator, a_plane_index, a_polygon)
    use, intrinsic :: iso_c_binding
    implicit none
      type(RectCub_type) :: a_rectangular_cuboid
      type(PlanarSep_type) :: a_planar_separator
      integer(IRL_UnsignedIndex_t), intent(in) :: a_plane_index
      type(Poly_type) :: a_polygon
      call F_getPoly_RectCub_Poly &
          (a_rectangular_cuboid%c_object, a_planar_separator%c_object, a_plane_index, a_polygon%c_object)
  end subroutine getPoly_RectCub_Poly

  subroutine getPoly_Tet_Poly(a_tet, a_planar_separator, a_plane_index, a_polygon)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Tet_type) :: a_tet
      type(PlanarSep_type) :: a_planar_separator
      integer(IRL_UnsignedIndex_t), intent(in) :: a_plane_index
      type(Poly_type) :: a_polygon
      call F_getPoly_Tet_Poly &
          (a_tet%c_object, a_planar_separator%c_object, a_plane_index, a_polygon%c_object)
  end subroutine getPoly_Tet_Poly

  subroutine getPoly_Hex_Poly(a_hexahedron, a_planar_separator, a_plane_index, a_polygon)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Hex_type) :: a_hexahedron
      type(PlanarSep_type) :: a_planar_separator
      integer(IRL_UnsignedIndex_t), intent(in) :: a_plane_index
      type(Poly_type) :: a_polygon
      call F_getPoly_Hex_Poly &
          (a_hexahedron%c_object, a_planar_separator%c_object, a_plane_index, a_polygon%c_object)
  end subroutine getPoly_Hex_Poly
  
  subroutine getPoly_RectCub_DivPoly(a_rectangular_cuboid, a_planar_separator, a_plane_index, a_divided_polygon)
    use, intrinsic :: iso_c_binding
    implicit none
      type(RectCub_type) :: a_rectangular_cuboid
      type(PlanarSep_type) :: a_planar_separator
      integer(IRL_UnsignedIndex_t), intent(in) :: a_plane_index
      type(DivPoly_type) :: a_divided_polygon
      call F_getPoly_RectCub_DivPoly &
          (a_rectangular_cuboid%c_object, a_planar_separator%c_object, a_plane_index, a_divided_polygon%c_object)
  end subroutine getPoly_RectCub_DivPoly

  function getSA_RectCub(a_rectangular_cuboid, a_planar_separator) result(a_surface_area)
    use, intrinsic :: iso_c_binding
    implicit none
      type(RectCub_type) :: a_rectangular_cuboid
      type(PlanarSep_type) :: a_planar_separator
      real(IRL_double) :: a_surface_area
      a_surface_area = F_getSA_RectCub(a_rectangular_cuboid%c_object, a_planar_separator%c_object)
  end function getSA_RectCub

end module f_CutPolygon
