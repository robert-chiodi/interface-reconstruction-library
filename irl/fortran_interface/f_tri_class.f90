!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_Tri_class.f90
!!
!! This file contains the Fortran interface for the
!! Tri class.

!> \brief A fortran type class that allows the creation of
!! IRL's Tri class along with enabling
!! some of its methods.
module f_Tri_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  use f_PlanarLoc_class
  implicit none

  type, public, bind(C) :: c_Tri
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_Tri
  type, public :: Tri_type
    type(c_Tri) :: c_object
  contains
    final :: Tri_class_delete
  end type Tri_type

  interface new
    module procedure Tri_class_new
  end interface
  interface construct
    module procedure Tri_class_construct
  end interface
  interface getVertices
    module procedure Tri_class_getVertices
  end interface
  interface calculateVolume
    module procedure Tri_class_calculateVolume
  end interface
  interface calculateCentroid
    module procedure Tri_class_calculateCentroid
  end interface
  interface calculateNormal
    module procedure Tri_class_calculateNormal
  end interface
  interface getLocalizer
    module procedure Tri_class_getLocalizer
  end interface
  interface reversePtOrdering
    module procedure Tri_class_reversePtOrdering
  end interface
  interface getBoundingPts
    module procedure Tri_class_getBoundingPts
  end interface
  interface calculateSign
    module procedure Tri_class_calculateSign
  end interface
  interface setPlaneOfExistence
    module procedure Tri_class_setPlaneOfExistence
  end interface
  interface calculateAndSetPlaneOfExistence
    module procedure Tri_class_calculateAndSetPlaneOfExistence
  end interface
  interface getPlaneOfExistence
    module procedure Tri_class_getPlaneOfExistence
  end interface


  interface

    subroutine F_Tri_new(this) &
      bind(C, name="c_Tri_new")
      import
      implicit none
      type(c_Tri) :: this
    end subroutine F_Tri_new

    subroutine F_Tri_delete(this) &
      bind(C, name="c_Tri_delete")
      import
      implicit none
      type(c_Tri) :: this
    end subroutine F_Tri_delete

    subroutine F_Tri_construct(this, a_pts) &
      bind(C, name="c_Tri_construct")
      import
      implicit none
      type(c_Tri) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_pts ! dimension(1:3,1:3)
    end subroutine F_Tri_construct

    subroutine F_Tri_getVertices(this, a_pts) &
      bind(C, name="c_Tri_getVertices")
      import
      implicit none
      type(c_Tri) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_pts ! dimension(1:3,1:3)
    end subroutine F_Tri_getVertices

    function F_Tri_calculateVolume(this) result(a_volume) &
      bind(C, name="c_Tri_calculateVolume")
      import
      implicit none
      type(c_Tri) :: this
      real(C_DOUBLE) :: a_volume
    end function F_Tri_calculateVolume

    subroutine F_Tri_calculateCentroid(this, a_centroid) &
      bind(C, name="c_Tri_calculateCentroid")
      import
      implicit none
      type(c_Tri) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_centroid ! dimension(1:3)
    end subroutine F_Tri_calculateCentroid

    subroutine F_Tri_calculateNormal(this, a_normal) &
      bind(C, name="c_Tri_calculateNormal")
      import
      implicit none
      type(c_Tri) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_normal ! dimension(1:3)
    end subroutine F_Tri_calculateNormal

    subroutine F_Tri_getLocalizer(this, a_planar_localizer) &
      bind(C, name="c_Tri_getLocalizer")
      import
      implicit none
      type(c_Tri) :: this
      type(c_PlanarLoc) :: a_planar_localizer ! Ptr to PlanarLoc object
    end subroutine F_Tri_getLocalizer

    subroutine F_Tri_reversePtOrdering(this) &
      bind(C, name="c_Tri_reversePtOrdering")
      import
      implicit none
      type(c_Tri) :: this
    end subroutine F_Tri_reversePtOrdering

    subroutine F_Tri_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_Tri_getBoundingPts")
      import
      implicit none
      type(c_Tri) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_Tri_getBoundingPts

    function F_Tri_calculateSign(this) result(a_sign) &
      bind(C, name="c_Tri_calculateSign")
      import
      implicit none
      type(c_Tri) :: this
      real(C_DOUBLE) :: a_sign
    end function F_Tri_calculateSign

    subroutine F_Tri_setPlaneOfExistence(this, a_plane) &
      bind(C, name="c_Tri_setPlaneOfExistence")
      import
      implicit none
      type(c_Tri) :: this
      real(C_DOUBLE), dimension(*) :: a_plane ! dimension(1:4)
    end subroutine F_Tri_setPlaneOfExistence

    subroutine F_Tri_calculateAndSetPlaneOfExistence(this) &
      bind(C, name="c_Tri_calculateAndSetPlaneOfExistence")
      import
      implicit none
      type(c_Tri) :: this
    end subroutine F_Tri_calculateAndSetPlaneOfExistence

    subroutine F_Tri_getPlaneOfExistence(this, a_plane) &
      bind(C, name="c_Tri_getPlaneOfExistence")
      import
      implicit none
      type(c_Tri) :: this
      real(C_DOUBLE), dimension(*) :: a_plane ! dimension(1:4)
    end subroutine F_Tri_getPlaneOfExistence

  end interface


  contains

    subroutine Tri_class_new(this)
      implicit none
      type(Tri_type), intent(inout) :: this
      call F_Tri_new(this%c_object)
    end subroutine Tri_class_new

    impure elemental subroutine Tri_class_delete(this)
      implicit none
      type(Tri_type), intent(in) :: this
      call F_Tri_delete(this%c_object)
    end subroutine Tri_class_delete

    subroutine Tri_class_construct(this, a_pts)
      implicit none
      type(Tri_type), intent(inout) :: this
      real(IRL_double), dimension(1:3,1:3), intent(in) :: a_pts
      call F_Tri_construct(this%c_object, a_pts)
    end subroutine Tri_class_construct

    function Tri_class_getVertices(this) result(a_pts)
      implicit none
      type(Tri_type), intent(inout) :: this
      real(IRL_double), dimension(1:3,1:3) :: a_pts
      call F_Tri_getVertices(this%c_object, a_pts)
      return
    end function Tri_class_getVertices

    function Tri_class_calculateVolume(this) result(a_volume)
      implicit none
      type(Tri_type), intent(in) :: this
      real(IRL_double) :: a_volume
      a_volume = F_Tri_calculateVolume(this%c_object)
      return
    end function Tri_class_calculateVolume

    function Tri_class_calculateCentroid(this) result(a_centroid)
      implicit none
      type(Tri_type), intent(in) :: this
      real(IRL_double), dimension(1:3) :: a_centroid
      call F_Tri_calculateCentroid(this%c_object, a_centroid)
      return
    end function Tri_class_calculateCentroid

    function Tri_class_calculateNormal(this) result(a_normal)
      implicit none
      type(Tri_type), intent(in) :: this
      real(IRL_double), dimension(1:3) :: a_normal
      call F_Tri_calculateNormal(this%c_object, a_normal)
      return
    end function Tri_class_calculateNormal

    subroutine Tri_class_getLocalizer(this, a_planar_localizer)
      implicit none
      type(Tri_type), intent(in) :: this
      type(PlanarLoc_type) :: a_planar_localizer
      call F_Tri_getLocalizer(this%c_object, a_planar_localizer%c_object)
      return
    end subroutine Tri_class_getLocalizer

    subroutine Tri_class_reversePtOrdering(this)
      implicit none
      type(Tri_type), intent(in) :: this
      call F_Tri_reversePtOrdering(this%c_object)
      return
    end subroutine Tri_class_reversePtOrdering

    subroutine Tri_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(Tri_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_Tri_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine Tri_class_getBoundingPts

    function Tri_class_calculateSign(this) result(a_sign)
      implicit none
      type(Tri_type), intent(in) :: this
      real(IRL_double) :: a_sign
      a_sign = F_Tri_calculateSign(this%c_object)
      return
    end function Tri_class_calculateSign

    subroutine Tri_class_setPlaneOfExistence(this, a_plane)
      implicit none
      type(Tri_type), intent(in) :: this
      real(IRL_double), dimension(4) :: a_plane
      call F_Tri_setPlaneOfExistence(this%c_object, a_plane)
      return
    end subroutine Tri_class_setPlaneOfExistence

    subroutine Tri_class_calculateAndSetPlaneOfExistence(this)
      implicit none
      type(Tri_type), intent(in) :: this
      call F_Tri_calculateAndSetPlaneOfExistence(this%c_object)
      return
    end subroutine Tri_class_calculateAndSetPlaneOfExistence

    function Tri_class_getPlaneOfExistence(this) result(a_plane)
      implicit none
      type(Tri_type), intent(in) :: this
      real(IRL_double), dimension(4) :: a_plane
      call F_Tri_getPlaneOfExistence(this%c_object, a_plane)
      return
    end function Tri_class_getPlaneOfExistence


end module f_Tri_class
