!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_pyramid_class.f90
!!
!! This file contains the Fortran interface for the
!! Pyramid class.

!> \brief A fortran type class that allows the creation of
!! IRL's Pyramid class along with enabling
!! some of its methods.
module f_Pyrmd_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_Pyrmd
    private
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_Pyrmd

  type, public :: Pyrmd_type
    type(c_Pyrmd) :: c_object
  contains
    final :: Pyrmd_class_delete
  end type Pyrmd_type

  interface new
    module procedure Pyrmd_class_new
  end interface
  interface construct
    module procedure Pyrmd_class_construct
  end interface
  interface calculateVolume
    module procedure Pyrmd_class_calculateVolume
  end interface
  interface printToScreen
    module procedure Pyrmd_class_printToScreen
  end interface
  interface getBoundingPts
    module procedure Pyrmd_class_getBoundingPts
  end interface

  interface

    subroutine F_Pyrmd_new(this) &
      bind(C, name="c_Pyrmd_new")
      import
      implicit none
      type(c_Pyrmd) :: this
    end subroutine F_Pyrmd_new

    subroutine F_Pyrmd_delete(this) &
      bind(C, name="c_Pyrmd_delete")
      import
      implicit none
      type(c_Pyrmd) :: this
    end subroutine F_Pyrmd_delete

    subroutine F_Pyrmd_construct(this, a_transported_cell) &
      bind(C, name="c_Pyrmd_construct")
      import
      implicit none
      type(c_Pyrmd) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_transported_cell ! dimension(1:3,1:5)
    end subroutine F_Pyrmd_construct

    function F_Pyrmd_calculateVolume(this) result(a_octahedron_volume) &
      bind(C, name="c_Pyrmd_calculateVolume")
      import
      implicit none
      type(c_Pyrmd) :: this
      real(C_DOUBLE) :: a_octahedron_volume
    end function F_Pyrmd_calculateVolume

    subroutine F_Pyrmd_printToScreen(this) &
      bind(C, name="c_Pyrmd_printToScreen")
      import
      implicit none
      type(c_Pyrmd) :: this
    end subroutine F_Pyrmd_printToScreen

    subroutine F_Pyrmd_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_Pyrmd_getBoundingPts")
      import
      implicit none
      type(c_Pyrmd) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_Pyrmd_getBoundingPts

  end interface

  contains

    subroutine Pyrmd_class_new(this)
      implicit none
      type(Pyrmd_type), intent(inout) :: this
      call F_Pyrmd_new(this%c_object)
    end subroutine Pyrmd_class_new

    impure elemental subroutine Pyrmd_class_delete(this)
      implicit none
      type(Pyrmd_type), intent(in) :: this
      call F_Pyrmd_delete(this%c_object)
    end subroutine Pyrmd_class_delete

    subroutine Pyrmd_class_construct(this, a_transported_cell)
      implicit none
      type(Pyrmd_type), intent(inout) :: this
      real(IRL_double), dimension(1:3,1:5), intent(in) :: a_transported_cell
      call F_Pyrmd_construct(this%c_object, a_transported_cell)
    end subroutine Pyrmd_class_construct

    function Pyrmd_class_calculateVolume(this) result(a_octahedron_volume)
      implicit none
      type(Pyrmd_type), intent(inout) :: this
      real(IRL_double) :: a_octahedron_volume
      a_octahedron_volume = F_Pyrmd_calculateVolume(this%c_object)
    end function Pyrmd_class_calculateVolume

    subroutine Pyrmd_class_printToScreen(this)
      implicit none
      type(Pyrmd_type), intent(in) :: this
      call F_Pyrmd_printToScreen(this%c_object)
    end subroutine Pyrmd_class_printToScreen

    subroutine Pyrmd_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(Pyrmd_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_Pyrmd_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine Pyrmd_class_getBoundingPts


end module f_Pyrmd_class
