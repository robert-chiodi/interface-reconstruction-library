!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_RectangularCuboid_class.f90
!!
!! This file contains the Fortran interface for the
!! RectangularCuboid class.

!> \brief A fortran type class that allows the creation of
!! IRL's RectangularCuboid class along with enabling
!! some of its methods.
module f_RectCub_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_RectCub
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_RectCub

  type, public :: RectCub_type
    type(c_RectCub) :: c_object
  contains
    final :: RectCub_class_delete
  end type RectCub_type

  interface new
    module procedure RectCub_class_new
  end interface
  interface construct
    module procedure RectCub_class_construct
  end interface
  interface construct_2pt
    module procedure RectCub_class_construct_2pt
  end interface
  interface calculateVolume
    module procedure RectCub_class_calculateVolume
  end interface
  interface getBoundingPts
    module procedure RectCub_class_getBoundingPts
  end interface

  interface

    subroutine F_RectCub_new(this) &
      bind(C, name="c_RectCub_new")
      import
      implicit none
      type(c_RectCub) :: this
    end subroutine F_RectCub_new

    subroutine F_RectCub_delete(this) &
      bind(C, name="c_RectCub_delete")
      import
      implicit none
      type(c_RectCub) :: this
    end subroutine F_RectCub_delete

    subroutine F_RectCub_construct(this, a_transported_cell) &
      bind(C, name="c_RectCub_construct")
      import
      implicit none
      type(c_RectCub) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_transported_cell ! dimension(1:3,1:8)
    end subroutine F_RectCub_construct

    subroutine F_RectCub_construct_2pt(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_RectCub_construct_2pt")
      import
      implicit none
      type(c_RectCub) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(in) :: a_upper_pt ! dimension(1:3)
    end subroutine F_RectCub_construct_2pt

    function F_RectCub_calculateVolume(this) result(a_cuboid_volume) &
      bind(C, name="c_RectCub_calculateVolume")
      import
      implicit none
      type(c_RectCub) :: this
      real(C_DOUBLE) :: a_cuboid_volume
    end function F_RectCub_calculateVolume

    subroutine F_RectCub_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_RectCub_getBoundingPts")
      import
      implicit none
      type(c_RectCub) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_RectCub_getBoundingPts

  end interface


  contains

    subroutine RectCub_class_new(this)
      implicit none
      type(RectCub_type), intent(inout) :: this
      call F_RectCub_new(this%c_object)
    end subroutine RectCub_class_new

    impure elemental subroutine RectCub_class_delete(this)
      use, intrinsic :: iso_c_binding
      implicit none
      type(RectCub_type), intent(in) :: this
      call F_RectCub_delete(this%c_object)
    end subroutine RectCub_class_delete

    subroutine RectCub_class_construct(this, a_transported_cell)
      implicit none
      type(RectCub_type), intent(inout) :: this
      real(IRL_double), dimension(1:3,1:8), intent(in) :: a_transported_cell
      call F_RectCub_construct(this%c_object, a_transported_cell)
    end subroutine RectCub_class_construct

    subroutine RectCub_class_construct_2pt(this, a_lower_pt, a_upper_pt)
      implicit none
      type(RectCub_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(in) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(in) :: a_upper_pt
      call F_RectCub_construct_2pt(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine RectCub_class_construct_2pt

    function RectCub_class_calculateVolume(this) result(a_cuboid_volume)
      implicit none
      type(RectCub_type), intent(inout) :: this
      real(IRL_double) :: a_cuboid_volume
      a_cuboid_volume = F_RectCub_calculateVolume(this%c_object)
    end function RectCub_class_calculateVolume

    subroutine RectCub_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(RectCub_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_RectCub_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine RectCub_class_getBoundingPts


end module f_RectCub_class
