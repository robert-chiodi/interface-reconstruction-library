!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_triangularprism_class.f90
!!
!! This file contains the Fortran interface for the
!! TriangularPrism class.

!> \brief A fortran type class that allows the creation of
!! IRL's TriangularPrism class along with enabling
!! some of its methods.
module f_TriPrism_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_TriPrism
    private
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_TriPrism

  type, public :: TriPrism_type
    type(c_TriPrism) :: c_object
  contains
    final :: TriPrism_class_delete
  end type TriPrism_type

  interface new
    module procedure TriPrism_class_new
  end interface
  interface construct
    module procedure TriPrism_class_construct
  end interface
  interface calculateVolume
    module procedure TriPrism_class_calculateVolume
  end interface
  interface printToScreen
    module procedure TriPrism_class_printToScreen
  end interface
  interface getBoundingPts
    module procedure TriPrism_class_getBoundingPts
  end interface

  interface

    subroutine F_TriPrism_new(this) &
      bind(C, name="c_TriPrism_new")
      import
      implicit none
      type(c_TriPrism) :: this
    end subroutine F_TriPrism_new

    subroutine F_TriPrism_delete(this) &
      bind(C, name="c_TriPrism_delete")
      import
      implicit none
      type(c_TriPrism) :: this
    end subroutine F_TriPrism_delete

    subroutine F_TriPrism_construct(this, a_transported_cell) &
      bind(C, name="c_TriPrism_construct")
      import
      implicit none
      type(c_TriPrism) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_transported_cell ! dimension(1:3,1:6)
    end subroutine F_TriPrism_construct

    function F_TriPrism_calculateVolume(this) result(a_triprism_volume) &
      bind(C, name="c_TriPrism_calculateVolume")
      import
      implicit none
      type(c_TriPrism) :: this
      real(C_DOUBLE) :: a_triprism_volume
    end function F_TriPrism_calculateVolume

    subroutine F_TriPrism_printToScreen(this) &
      bind(C, name="c_TriPrism_printToScreen")
      import
      implicit none
      type(c_TriPrism) :: this
    end subroutine F_TriPrism_printToScreen

    subroutine F_TriPrism_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_TriPrism_getBoundingPts")
      import
      implicit none
      type(c_TriPrism) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_TriPrism_getBoundingPts

  end interface

  contains

    subroutine TriPrism_class_new(this)
      implicit none
      type(TriPrism_type), intent(inout) :: this
      call F_TriPrism_new(this%c_object)
    end subroutine TriPrism_class_new

    impure elemental subroutine TriPrism_class_delete(this)
      implicit none
      type(TriPrism_type), intent(in) :: this
      call F_TriPrism_delete(this%c_object)
    end subroutine TriPrism_class_delete

    subroutine TriPrism_class_construct(this, a_transported_cell)
      implicit none
      type(TriPrism_type), intent(inout) :: this
      real(IRL_double), dimension(1:3,1:6), intent(in) :: a_transported_cell
      call F_TriPrism_construct(this%c_object, a_transported_cell)
    end subroutine TriPrism_class_construct

    function TriPrism_class_calculateVolume(this) result(a_triprism_volume)
      implicit none
      type(TriPrism_type), intent(inout) :: this
      real(IRL_double) :: a_triprism_volume
      a_triprism_volume = F_TriPrism_calculateVolume(this%c_object)
    end function TriPrism_class_calculateVolume

    subroutine TriPrism_class_printToScreen(this)
      implicit none
      type(TriPrism_type), intent(in) :: this
      call F_TriPrism_printToScreen(this%c_object)
    end subroutine TriPrism_class_printToScreen

    subroutine TriPrism_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(TriPrism_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_TriPrism_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine TriPrism_class_getBoundingPts


end module f_TriPrism_class
