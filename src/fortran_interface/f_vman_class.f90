!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_VMAN_class.f90
!!
!! This file contains the Fortran interface for volume moments classes.

!> \brief A fortran type class that allows the creation of
!! IRL's AccumulatedListedVolumeMomentsM<VMAN>
!! class along with enabling some of its methods.
module f_VMAN_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_VMAN
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_VMAN

  type, public :: VMAN_type
    type(c_VMAN) :: c_object
  contains
    final :: VMAN_class_delete
  end type VMAN_type

  interface new
    module procedure VMAN_class_new
  end interface
  interface getVolume
    module procedure VMAN_class_getVolume
  end interface
  interface getCentroid
    module procedure VMAN_class_getCentroid
  end interface
  interface getNormal
    module procedure VMAN_class_getNormal
  end interface
  interface normalizeByVolume
    module procedure VMAN_class_normalizeByVolume
  end interface
  interface multiplyByVolume
    module procedure VMAN_class_multiplyByVolume
  end interface

  interface

    subroutine F_VMAN_new(this) &
      bind(C, name="c_VMAN_new")
      import
      implicit none
      type(c_VMAN) :: this
    end subroutine F_VMAN_new

    subroutine F_VMAN_delete(this) &
      bind(C, name="c_VMAN_delete")
      import
      implicit none
      type(c_VMAN) :: this
    end subroutine F_VMAN_delete

    subroutine F_VMAN_getVolume(this, a_volume) &
      bind(C, name="c_VMAN_getVolume")
      import
      implicit none
      type(c_VMAN) :: this
      real(C_DOUBLE) :: a_volume
    end subroutine F_VMAN_getVolume

    subroutine F_VMAN_getCentroid(this, a_centroid) &
      bind(C, name="c_VMAN_getCentroid")
      import
      implicit none
      type(c_VMAN) :: this
      real(C_DOUBLE), dimension(*) :: a_centroid
    end subroutine F_VMAN_getCentroid

    subroutine F_VMAN_getNormal(this, a_normal) &
      bind(C, name="c_VMAN_getNormal")
      import
      implicit none
      type(c_VMAN) :: this
      real(C_DOUBLE), dimension(*) :: a_normal
    end subroutine F_VMAN_getNormal

    subroutine F_VMAN_normalizeByVolume(this) &
      bind(C, name="c_VMAN_normalizeByVolume")
      import
      implicit none
      type(c_VMAN) :: this
    end subroutine F_VMAN_normalizeByVolume

    subroutine F_VMAN_multiplyByVolume(this) &
      bind(C, name="c_VMAN_multiplyByVolume")
      import
      implicit none
      type(c_VMAN) :: this
    end subroutine F_VMAN_multiplyByVolume

  end interface


  contains

    subroutine VMAN_class_new(this)
      implicit none
      type(VMAN_type), intent(inout) :: this
      call F_VMAN_new(this%c_object)
    end subroutine VMAN_class_new

    impure elemental subroutine VMAN_class_delete(this)
      implicit none
      type(VMAN_type), intent(in) :: this
      call F_VMAN_delete(this%c_object)
    end subroutine VMAN_class_delete

    function VMAN_class_getVolume(this) result(a_volume)
      implicit none
      type(VMAN_type), intent(in) :: this
      real(IRL_double) :: a_volume
      call F_VMAN_getVolume(this%c_object, a_volume)
      return
    end function VMAN_class_getVolume

    function VMAN_class_getCentroid(this) result(a_centroid)
      implicit none
      type(VMAN_type), intent(in) :: this
      real(IRL_double), dimension(3) :: a_centroid
      call F_VMAN_getCentroid(this%c_object, a_centroid)
      return
    end function VMAN_class_getCentroid

    function VMAN_class_getNormal(this) result(a_normal)
      implicit none
      type(VMAN_type), intent(in) :: this
      real(IRL_double), dimension(3) :: a_normal
      call F_VMAN_getNormal(this%c_object, a_normal)
      return
    end function VMAN_class_getNormal

    subroutine VMAN_class_normalizeByVolume(this)
      implicit none
      type(VMAN_type), intent(inout) :: this
      call F_VMAN_normalizeByVolume(this%c_object)
    end subroutine VMAN_class_normalizeByVolume

    subroutine VMAN_class_multiplyByVolume(this)
      implicit none
      type(VMAN_type), intent(inout) :: this
      call F_VMAN_multiplyByVolume(this%c_object)
    end subroutine VMAN_class_multiplyByVolume


end module f_VMAN_class
