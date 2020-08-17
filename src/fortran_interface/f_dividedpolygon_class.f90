!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_DividedPolygon_class.f90
!!
!! This file contains the Fortran interface for the
!! DividedPolygon class.

!> \brief A fortran type class that allows the creation of
!! IRL's DividedPolygon class along with enabling
!! some of its methods.
module f_DivPoly_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  use f_Poly_class
  use f_Tri_class
  use f_PlanarLoc_class
  implicit none

  type, public, bind(C) :: c_DivPoly
    private
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_DivPoly

  type, public :: DivPoly_type
    type(c_DivPoly) :: c_object
  contains
    final :: DivPoly_class_delete
  end type DivPoly_type

  interface new
    module procedure DivPoly_class_new
  end interface

  interface construct
    module procedure DivPoly_class_construct
  end interface

  interface constructFromPolygon
    module procedure DivPoly_class_constructFromPolygon
  end interface

  interface resetCentroid
    module procedure DivPoly_class_resetCentroid
  end interface

  interface getNumberOfSimplicesInDecomposition
    module procedure DivPoly_class_getNumberOfSimplicesInDecomposition
  end interface

  interface getSimplexFromDecomposition
    module procedure DivPoly_class_getSimplexFromDecomposition
  end interface

  interface getLocalizer
    module procedure DivPoly_class_getLocalizer
  end interface

  interface calculateNormal
    module procedure DivPoly_class_calculateNormal
  end interface

  interface reversePtOrdering
    module procedure DivPoly_class_reversePtOrdering
  end interface

  interface getBoundingPts
    module procedure DivPoly_class_getBoundingPts
  end interface

  interface getNumberOfVertices
    module procedure DivPoly_class_getNumberOfPts
  end interface

  interface getPt
    module procedure DivPoly_class_getPt
  end interface

  interface zeroPolygon
    module procedure DivPoly_class_zeroPolygon
  end interface

  interface calculateSurfaceArea
    module procedure DivPoly_class_calculateSurfaceArea
  end interface

  interface calculateSign
    module procedure DivPoly_class_calculateSign
  end interface

  interface setPlaneOfExistence
    module procedure DivPoly_class_setPlaneOfExistence
  end interface

  interface calculateAndSetPlaneOfExistence
    module procedure DivPoly_class_calculateAndSetPlaneOfExistence
  end interface

  interface getPlaneOfExistence
    module procedure DivPoly_class_getPlaneOfExistence
  end interface

  interface printToScreen
    module procedure DivPoly_class_printToScreen
  end interface

  interface

    subroutine F_DivPoly_new(this) &
      bind(C, name="c_DivPoly_new")
      import
      implicit none
      type(c_DivPoly) :: this
    end subroutine F_DivPoly_new

    subroutine F_DivPoly_delete(this) &
      bind(C, name="c_DivPoly_delete")
      import
      implicit none
      type(c_DivPoly) :: this
    end subroutine F_DivPoly_delete

    subroutine F_DivPoly_construct(this, a_npts, a_pts) &
      bind(C, name="c_DivPoly_construct")
      import
      implicit none
      type(c_DivPoly) :: this
      integer(C_INT), intent(in) :: a_npts
      real(C_DOUBLE), dimension(*), intent(in) :: a_pts ! dimension(1:3,1:npts)
    end subroutine F_DivPoly_construct

    subroutine F_DivPoly_constructFromPolygon(this, a_polygon) &
      bind(C, name="c_DivPoly_constructFromPolygon")
      import
      implicit none
      type(c_DivPoly) :: this
      type(c_Poly) :: a_polygon
    end subroutine F_DivPoly_constructFromPolygon

    subroutine F_DivPoly_resetCentroid(this) &
      bind(C, name="c_DivPoly_resetCentroid")
      import
      implicit none
      type(c_DivPoly) :: this
    end subroutine F_DivPoly_resetCentroid

    function F_DivPoly_getNumberOfSimplicesInDecomposition(this) result(a_number_of_triangles) &
      bind(C, name="c_DivPoly_getNumberOfSimplicesInDecomposition")
      import
      implicit none
      type(c_DivPoly) :: this
      integer(C_INT) :: a_number_of_triangles
    end function F_DivPoly_getNumberOfSimplicesInDecomposition

    subroutine F_DivPoly_getSimplexFromDecomposition(this, a_tri_number_to_get, a_triangle_in_decomposition) &
      bind(C, name="c_DivPoly_getSimplexFromDecomposition")
      import
      implicit none
      type(c_DivPoly) :: this
      integer(C_INT) :: a_tri_number_to_get
      type(c_Tri) :: a_triangle_in_decomposition
    end subroutine F_DivPoly_getSimplexFromDecomposition

    subroutine F_DivPoly_getLocalizer(this, a_planar_localizer) &
      bind(C, name="c_DivPoly_getLocalizer")
      import
      implicit none
      type(c_DivPoly) :: this
      type(c_PlanarLoc) :: a_planar_localizer ! Ptr to PlanarLoc object
    end subroutine F_DivPoly_getLocalizer

    subroutine F_DivPoly_calculateNormal(this, a_normal) &
      bind(C, name="c_DivPoly_calculateNormal")
      import
      implicit none
      type(c_DivPoly) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_normal ! dimension(1:3)
    end subroutine F_DivPoly_calculateNormal

    subroutine F_DivPoly_reversePtOrdering(this) &
      bind(C, name="c_DivPoly_reversePtOrdering")
      import
      implicit none
      type(c_DivPoly) :: this
    end subroutine F_DivPoly_reversePtOrdering

    subroutine F_DivPoly_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_DivPoly_getBoundingPts")
      import
      implicit none
      type(c_DivPoly) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_DivPoly_getBoundingPts

    function F_DivPoly_getNumberOfPts(this) result(a_number_of_pts) &
      bind(C, name="c_DivPoly_getNumberOfPts")
      import
      implicit none
      type(c_DivPoly) :: this
      integer(C_INT) :: a_number_of_pts
    end function F_DivPoly_getNumberOfPts

    subroutine F_DivPoly_getPt(this, a_index, a_pt) &
      bind(C, name="c_DivPoly_getPt")
      import
      implicit none
      type(c_DivPoly) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(out) :: a_pt ! dimension(1:3)
    end subroutine F_DivPoly_getPt

    subroutine F_DivPoly_zeroPolygon(this) &
      bind(C, name="c_DivPoly_zeroPolygon")
      import
      implicit none
      type(c_DivPoly) :: this
    end subroutine F_DivPoly_zeroPolygon

    function F_DivPoly_calculateSurfaceArea(this) result(a_surface_area) &
      bind(C, name="c_DivPoly_calculateSurfaceArea")
      import
      implicit none
      type(c_DivPoly) :: this
      real(C_DOUBLE) :: a_surface_area
    end function F_DivPoly_calculateSurfaceArea

    function F_DivPoly_calculateSign(this) result(a_sign) &
      bind(C, name="c_DivPoly_calculateSign")
      import
      implicit none
      type(c_DivPoly) :: this
      real(C_DOUBLE) :: a_sign
    end function F_DivPoly_calculateSign

    subroutine F_DivPoly_setPlaneOfExistence(this, a_plane) &
      bind(C, name="c_DivPoly_setPlaneOfExistence")
      import
      implicit none
      type(c_DivPoly) :: this
      real(C_DOUBLE), dimension(*) :: a_plane ! dimension(1:4)
    end subroutine F_DivPoly_setPlaneOfExistence

    subroutine F_DivPoly_calculateAndSetPlaneOfExistence(this) &
      bind(C, name="c_DivPoly_calculateAndSetPlaneOfExistence")
      import
      implicit none
      type(c_DivPoly) :: this
    end subroutine F_DivPoly_calculateAndSetPlaneOfExistence

    subroutine F_DivPoly_getPlaneOfExistence(this, a_plane) &
      bind(C, name="c_DivPoly_getPlaneOfExistence")
      import
      implicit none
      type(c_DivPoly) :: this
      real(C_DOUBLE), dimension(*) :: a_plane ! dimension(1:4)
    end subroutine F_DivPoly_getPlaneOfExistence

    subroutine F_DivPoly_printToScreen(this) &
      bind(C, name="c_DivPoly_printToScreen")
      import
      implicit none
      type(c_DivPoly) :: this
    end subroutine F_DivPoly_printToScreen

  end interface


  contains

    subroutine DivPoly_class_new(this)
      implicit none
      type(DivPoly_type), intent(inout) :: this
      call F_DivPoly_new(this%c_object)
    end subroutine DivPoly_class_new

    impure elemental subroutine DivPoly_class_delete(this)
      implicit none
      type(DivPoly_type), intent(in) :: this
      call F_DivPoly_delete(this%c_object)
    end subroutine DivPoly_class_delete

    subroutine DivPoly_class_construct(this, a_npts, a_pts)
      implicit none
      type(DivPoly_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_npts
      real(IRL_double), dimension(1:3,1:a_npts), intent(in) :: a_pts
      call F_DivPoly_construct(this%c_object, a_npts, a_pts)
    end subroutine DivPoly_class_construct

    subroutine DivPoly_class_constructFromPolygon(this, a_polygon)
      implicit none
      type(DivPoly_type), intent(inout) :: this
      type(Poly_type), intent(in) :: a_polygon
      call F_DivPoly_constructFromPolygon(this%c_object, a_polygon%c_object)
    end subroutine DivPoly_class_constructFromPolygon

    subroutine DivPoly_class_resetCentroid(this)
      implicit none
      type(DivPoly_type), intent(inout) :: this
      call F_DivPoly_resetCentroid(this%c_object)
    end subroutine DivPoly_class_resetCentroid

    function DivPoly_class_getNumberOfSimplicesInDecomposition(this) result(a_number_of_triangles)
      implicit none
      type(DivPoly_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_number_of_triangles
      a_number_of_triangles = F_DivPoly_getNumberOfSimplicesInDecomposition(this%c_object)
    end function DivPoly_class_getNumberOfSimplicesInDecomposition

    subroutine DivPoly_class_getSimplexFromDecomposition(this, a_tri_number_to_get, a_tri_in_decomposition)
      implicit none
      type(DivPoly_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_tri_number_to_get
      type(Tri_type), intent(inout) :: a_tri_in_decomposition
      call F_DivPoly_getSimplexFromDecomposition(this%c_object, a_tri_number_to_get, a_tri_in_decomposition%c_object)
    end subroutine DivPoly_class_getSimplexFromDecomposition

    subroutine DivPoly_class_getLocalizer(this, a_planar_localizer)
      implicit none
      type(DivPoly_type), intent(in) :: this
      type(PlanarLoc_type) :: a_planar_localizer
      call F_DivPoly_getLocalizer(this%c_object, a_planar_localizer%c_object)
      return
    end subroutine DivPoly_class_getLocalizer

    function DivPoly_class_calculateNormal(this) result(a_normal)
      implicit none
      type(DivPoly_type), intent(in) :: this
      real(IRL_double), dimension(1:3) :: a_normal
      call F_DivPoly_calculateNormal(this%c_object, a_normal)
      return
    end function DivPoly_class_calculateNormal

    subroutine DivPoly_class_reversePtOrdering(this)
      implicit none
      type(DivPoly_type), intent(in) :: this
      call F_DivPoly_reversePtOrdering(this%c_object)
      return
    end subroutine DivPoly_class_reversePtOrdering

    subroutine DivPoly_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(DivPoly_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_DivPoly_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine DivPoly_class_getBoundingPts

    function DivPoly_class_getNumberOfPts(this) result(a_number_of_pts)
      implicit none
      type(DivPoly_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t) :: a_number_of_pts
      a_number_of_pts = F_DivPoly_getNumberOfPts(this%c_object)
      return
    end function DivPoly_class_getNumberOfPts

    function DivPoly_class_getPt(this, a_index) result(a_pt)
      implicit none
      type(DivPoly_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(3) :: a_pt
      call F_DivPoly_getPt(this%c_object, a_index, a_pt)
      return
    end function DivPoly_class_getPt

    subroutine DivPoly_class_zeroPolygon(this)
      implicit none
      type(DivPoly_type), intent(in) :: this
      call F_DivPoly_zeroPolygon(this%c_object)
    end subroutine DivPoly_class_zeroPolygon

    function DivPoly_class_calculateSurfaceArea(this) result(a_surface_area)
      implicit none
      type(DivPoly_type), intent(in) :: this
      real(IRL_double) :: a_surface_area
      a_surface_area = F_DivPoly_calculateSurfaceArea(this%c_object)
      return
    end function DivPoly_class_calculateSurfaceArea

    function DivPoly_class_calculateSign(this) result(a_sign)
      implicit none
      type(DivPoly_type), intent(in) :: this
      real(IRL_double) :: a_sign
      a_sign = F_DivPoly_calculateSign(this%c_object)
      return
    end function DivPoly_class_calculateSign

    subroutine DivPoly_class_setPlaneOfExistence(this, a_plane)
      implicit none
      type(DivPoly_type), intent(in) :: this
      real(IRL_double), dimension(4) :: a_plane
      call F_DivPoly_setPlaneOfExistence(this%c_object, a_plane)
      return
    end subroutine DivPoly_class_setPlaneOfExistence

    subroutine DivPoly_class_calculateAndSetPlaneOfExistence(this)
      implicit none
      type(DivPoly_type), intent(in) :: this
      call F_DivPoly_calculateAndSetPlaneOfExistence(this%c_object)
      return
    end subroutine DivPoly_class_calculateAndSetPlaneOfExistence

    function DivPoly_class_getPlaneOfExistence(this) result(a_plane)
      implicit none
      type(DivPoly_type), intent(in) :: this
      real(IRL_double), dimension(4) :: a_plane
      call F_DivPoly_getPlaneOfExistence(this%c_object, a_plane)
      return
    end function DivPoly_class_getPlaneOfExistence

    subroutine DivPoly_class_printToScreen(this)
      implicit none
      type(DivPoly_type), intent(in) :: this
      call F_DivPoly_printToScreen(this%c_object)
    end subroutine DivPoly_class_printToScreen

end module f_DivPoly_class
