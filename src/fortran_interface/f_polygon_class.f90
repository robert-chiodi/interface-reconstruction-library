!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_Polygon_class.f90
!!
!! This file contains the Fortran interface for the
!! Polygon class.

!> \brief A fortran type class that allows the creation of
!! IRL's Polygon class along with enabling
!! some of its methods.
module f_Poly_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  use f_PlanarLoc_class
  use f_Tri_class
  implicit none

  type, public, bind(C) :: c_Poly
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_Poly

  type, public :: Poly_type
    type(c_Poly) :: c_object
  contains
    final :: Poly_class_delete
  end type Poly_type

  interface new
    module procedure Poly_class_new
  end interface
  interface construct
    module procedure Poly_class_construct
  end interface
  interface calculateNormal
    module procedure Poly_class_calculateNormal
  end interface
  interface getLocalizer
    module procedure Poly_class_getLocalizer
  end interface
  interface reversePtOrdering
    module procedure Poly_class_reversePtOrdering
  end interface
  interface getBoundingPts
    module procedure Poly_class_getBoundingPts
  end interface
  interface getNumberOfVertices
    module procedure Poly_class_getNumberOfPts
  end interface
  interface getPt
    module procedure Poly_class_getPt
  end interface
  interface getNumberOfSimplicesInDecomposition
    module procedure Poly_class_getNumberOfSimplicesInDecomposition
  end interface
  interface getSimplexFromDecomposition
    module procedure Poly_class_getSimplexFromDecomposition
  end interface
  interface zeroPolygon
    module procedure Poly_class_zeroPolygon
  end interface
  interface calculateNearestPtOnSurface
    module procedure Poly_class_calculateNearestPtOnSurface
  end interface
  interface calculateVolume
    module procedure Poly_class_calculateVolume
  end interface
  interface calculateSign
    module procedure Poly_class_calculateSign
  end interface
  interface setPlaneOfExistence
    module procedure Poly_class_setPlaneOfExistence
  end interface
  interface calculateAndSetPlaneOfExistence
    module procedure Poly_class_calculateAndSetPlaneOfExistence
  end interface
  interface calculateCentroid
    module procedure Poly_class_calculateCentroid
  end interface
  interface getPlaneOfExistence
    module procedure Poly_class_getPlaneOfExistence
  end interface
  interface printToScreen
    module procedure Poly_class_printToScreen
  end interface

  interface

    subroutine F_Poly_new(this) &
      bind(C, name="c_Poly_new")
      import
      implicit none
      type(c_Poly) :: this
    end subroutine F_Poly_new

    subroutine F_Poly_delete(this) &
      bind(C, name="c_Poly_delete")
      import
      implicit none
      type(c_Poly) :: this
    end subroutine F_Poly_delete

    subroutine F_Poly_construct(this, a_npts, a_pts) &
      bind(C, name="c_Poly_construct")
      import
      implicit none
      type(c_Poly) :: this
      integer(C_INT), intent(in) :: a_npts
      real(C_DOUBLE), dimension(*), intent(in) :: a_pts ! dimension(1:3,1:npts)
    end subroutine F_Poly_construct

    subroutine F_Poly_calculateNormal(this, a_normal) &
      bind(C, name="c_Poly_calculateNormal")
      import
      implicit none
      type(c_Poly) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_normal ! dimension(1:3)
    end subroutine F_Poly_calculateNormal

    subroutine F_Poly_getLocalizer(this, a_planar_localizer) &
      bind(C, name="c_Poly_getLocalizer")
      import
      implicit none
      type(c_Poly) :: this
      type(c_PlanarLoc) :: a_planar_localizer ! Ptr to PlanarLocalizer object
    end subroutine F_Poly_getLocalizer

    subroutine F_Poly_reversePtOrdering(this) &
      bind(C, name="c_Poly_reversePtOrdering")
      import
      implicit none
      type(c_Poly) :: this
    end subroutine F_Poly_reversePtOrdering

    subroutine F_Poly_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_Poly_getBoundingPts")
      import
      implicit none
      type(c_Poly) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_Poly_getBoundingPts

    function F_Poly_getNumberOfPts(this) result(a_number_of_pts) &
      bind(C, name="c_Poly_getNumberOfPts")
      import
      implicit none
      type(c_Poly) :: this
      integer(C_INT) :: a_number_of_pts
    end function F_Poly_getNumberOfPts

    subroutine F_Poly_getPt(this, a_index, a_pt) &
      bind(C, name="c_Poly_getPt")
      import
      implicit none
      type(c_Poly) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(out) :: a_pt ! dimension(1:3)
    end subroutine F_Poly_getPt

    function F_Poly_getNumberOfSimplicesInDecomposition(this) result(a_number_of_triangles) &
      bind(C, name="c_Poly_getNumberOfSimplicesInDecomposition")
      import
      implicit none
      type(c_Poly) :: this
      integer(C_INT) :: a_number_of_triangles
    end function F_Poly_getNumberOfSimplicesInDecomposition

    subroutine F_Poly_getSimplexFromDecomposition(this, a_tri_number_to_get, a_triangle_in_decomposition) &
      bind(C, name="c_Poly_getSimplexFromDecomposition")
      import
      implicit none
      type(c_Poly) :: this
      integer(C_INT) :: a_tri_number_to_get
      type(c_Tri) :: a_triangle_in_decomposition
    end subroutine F_Poly_getSimplexFromDecomposition

    subroutine F_Poly_zeroPolygon(this) &
      bind(C, name="c_Poly_zeroPolygon")
      import
      implicit none
      type(c_Poly) :: this
    end subroutine F_Poly_zeroPolygon

    subroutine F_Poly_calculateNearestPtOnSurface(this, a_pt, a_pt_on_polygon) &
      bind(C, name="c_Poly_calculateNearestPtOnSurface")
      import
      implicit none
      type(c_Poly) :: this
      real(C_DOUBLE), dimension(*) :: a_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*) :: a_pt_on_polygon ! dimension(1:3)
    end subroutine F_Poly_calculateNearestPtOnSurface

    function F_Poly_calculateVolume(this) result(a_surface_area) &
      bind(C, name="c_Poly_calculateVolume")
      import
      implicit none
      type(c_Poly) :: this
      real(C_DOUBLE) :: a_surface_area
    end function F_Poly_calculateVolume

    function F_Poly_calculateSign(this) result(a_sign) &
      bind(C, name="c_Poly_calculateSign")
      import
      implicit none
      type(c_Poly) :: this
      real(C_DOUBLE) :: a_sign
    end function F_Poly_calculateSign

    subroutine F_Poly_setPlaneOfExistence(this, a_plane) &
      bind(C, name="c_Poly_setPlaneOfExistence")
      import
      implicit none
      type(c_Poly) :: this
      real(C_DOUBLE), dimension(*) :: a_plane ! dimension(1:4)
    end subroutine F_Poly_setPlaneOfExistence

    subroutine F_Poly_calculateAndSetPlaneOfExistence(this) &
      bind(C, name="c_Poly_calculateAndSetPlaneOfExistence")
      import
      implicit none
      type(c_Poly) :: this
    end subroutine F_Poly_calculateAndSetPlaneOfExistence

    subroutine F_Poly_getPlaneOfExistence(this, a_plane) &
      bind(C, name="c_Poly_getPlaneOfExistence")
      import
      implicit none
      type(c_Poly) :: this
      real(C_DOUBLE), dimension(*) :: a_plane ! dimension(1:4)
    end subroutine F_Poly_getPlaneOfExistence

    subroutine F_Poly_calculateCentroid(this, a_centroid) &
      bind(C, name="c_Poly_calculateCentroid")
      import
      implicit none
      type(c_Poly) :: this
      real(C_DOUBLE), dimension(*) :: a_centroid ! dimension(1:3)
    end subroutine F_Poly_calculateCentroid

    subroutine F_Poly_printToScreen(this) &
      bind(C, name="c_Poly_printToScreen")
      import
      implicit none
      type(c_Poly) :: this
    end subroutine F_Poly_printToScreen


  end interface


  contains

    subroutine Poly_class_new(this)
      implicit none
      type(Poly_type), intent(inout) :: this
      call F_Poly_new(this%c_object)
    end subroutine Poly_class_new

    impure elemental subroutine Poly_class_delete(this)
      implicit none
      type(Poly_type), intent(in) :: this
      call F_Poly_delete(this%c_object)
    end subroutine Poly_class_delete

    subroutine Poly_class_construct(this, a_npts, a_pts)
      implicit none
      type(Poly_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_npts
      real(IRL_double), dimension(1:3,1:a_npts), intent(in) :: a_pts
      call F_Poly_construct(this%c_object, a_npts, a_pts)
    end subroutine Poly_class_construct

    function Poly_class_calculateNormal(this) result(a_normal)
      implicit none
      type(Poly_type), intent(in) :: this
      real(IRL_double), dimension(1:3) :: a_normal
      call F_Poly_calculateNormal(this%c_object, a_normal)
      return
    end function Poly_class_calculateNormal

    subroutine Poly_class_getLocalizer(this, a_planar_localizer)
      implicit none
      type(Poly_type), intent(in) :: this
      type(PlanarLoc_type) :: a_planar_localizer
      call F_Poly_getLocalizer(this%c_object, a_planar_localizer%c_object)
      return
    end subroutine Poly_class_getLocalizer

    subroutine Poly_class_reversePtOrdering(this)
      implicit none
      class(Poly_type), intent(in) :: this
      call F_Poly_reversePtOrdering(this%c_object)
      return
    end subroutine Poly_class_reversePtOrdering

    subroutine Poly_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(Poly_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_Poly_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine Poly_class_getBoundingPts

    function Poly_class_getNumberOfPts(this) result(a_number_of_pts)
      implicit none
      type(Poly_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t) :: a_number_of_pts
      a_number_of_pts = F_Poly_getNumberOfPts(this%c_object)
      return
    end function Poly_class_getNumberOfPts

    function Poly_class_getPt(this, a_index) result(a_pt)
      implicit none
      type(Poly_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(3) :: a_pt
      call F_Poly_getPt(this%c_object, a_index, a_pt)
      return
    end function Poly_class_getPt

    function Poly_class_getNumberOfSimplicesInDecomposition(this) result(a_number_of_triangles)
      implicit none
      type(Poly_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_number_of_triangles
      a_number_of_triangles = F_Poly_getNumberOfSimplicesInDecomposition(this%c_object)
    end function Poly_class_getNumberOfSimplicesInDecomposition

    subroutine Poly_class_getSimplexFromDecomposition(this, a_tri_number_to_get, a_tri_in_decomposition)
      implicit none
      type(Poly_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_tri_number_to_get
      type(Tri_type), intent(inout) :: a_tri_in_decomposition
      call F_Poly_getSimplexFromDecomposition(this%c_object, a_tri_number_to_get, a_tri_in_decomposition%c_object)
    end subroutine Poly_class_getSimplexFromDecomposition

    subroutine Poly_class_zeroPolygon(this)
      implicit none
      class(Poly_type), intent(in) :: this
      call F_Poly_zeroPolygon(this%c_object)
    end subroutine Poly_class_zeroPolygon

    function Poly_class_calculateNearestPtOnSurface(this, a_pt) result(a_pt_on_polygon)
      implicit none
      type(Poly_type), intent(in) :: this
      real(IRL_double), dimension(3) :: a_pt
      real(IRL_double), dimension(3) :: a_pt_on_polygon
      call F_Poly_calculateNearestPtOnSurface(this%c_object, a_pt, a_pt_on_polygon)
      return
    end function Poly_class_calculateNearestPtOnSurface

    function Poly_class_calculateVolume(this) result(a_surface_area)
      implicit none
      type(Poly_type), intent(in) :: this
      real(IRL_double) :: a_surface_area
      a_surface_area = F_Poly_calculateVolume(this%c_object)
      return
    end function Poly_class_calculateVolume

    function Poly_class_calculateSign(this) result(a_sign)
      implicit none
      type(Poly_type), intent(in) :: this
      real(IRL_double) :: a_sign
      a_sign = F_Poly_calculateSign(this%c_object)
      return
    end function Poly_class_calculateSign

    subroutine Poly_class_setPlaneOfExistence(this, a_plane)
      implicit none
      type(Poly_type), intent(in) :: this
      real(IRL_double), dimension(4) :: a_plane
      call F_Poly_setPlaneOfExistence(this%c_object, a_plane)
      return
    end subroutine Poly_class_setPlaneOfExistence

    subroutine Poly_class_calculateAndSetPlaneOfExistence(this)
      implicit none
      type(Poly_type), intent(in) :: this
      call F_Poly_calculateAndSetPlaneOfExistence(this%c_object)
      return
    end subroutine Poly_class_calculateAndSetPlaneOfExistence

    function Poly_class_getPlaneOfExistence(this) result(a_plane)
      implicit none
      type(Poly_type), intent(in) :: this
      real(IRL_double), dimension(4) :: a_plane
      call F_Poly_getPlaneOfExistence(this%c_object, a_plane)
      return
    end function Poly_class_getPlaneOfExistence

    function Poly_class_calculateCentroid(this) result(a_centroid)
      implicit none
      type(Poly_type), intent(in) :: this
      real(IRL_double), dimension(3) :: a_centroid
      call F_Poly_calculateCentroid(this%c_object, a_centroid)
      return
    end function Poly_class_calculateCentroid

    subroutine Poly_class_printToScreen(this)
      implicit none
      type(Poly_type), intent(in) :: this
      call F_Poly_printToScreen(this%c_object)
    end subroutine Poly_class_printToScreen

end module f_Poly_class
