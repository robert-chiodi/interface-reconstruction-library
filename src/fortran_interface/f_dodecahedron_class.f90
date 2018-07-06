!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_Dodecahedron_class.f90
!!
!! This file contains the Fortran interface for the
!! Dodecahedron class.

!> \brief A fortran type class that allows the creation of
!! IRL's Dodecahedron class along with enabling
!! some of its methods.
module f_Dod_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_Dod
    private
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_Dod

  type, public :: Dod_type
    type(c_Dod) :: c_object
  contains
    final :: Dod_class_delete
  end type Dod_type

  interface new
    module procedure Dod_class_new
  end interface
  interface construct
    module procedure Dod_class_construct
  end interface
  interface calculateVolume
    module procedure Dod_class_calculateVolume
  end interface
  interface printToScreen
    module procedure Dod_class_printToScreen
  end interface
  interface getBoundingPts
    module procedure Dod_class_getBoundingPts
  end interface

  interface

    subroutine F_Dod_new(this) &
      bind(C, name="c_Dod_new")
      import
      implicit none
      type(c_Dod) :: this
    end subroutine F_Dod_new

    subroutine F_Dod_delete(this) &
      bind(C, name="c_Dod_delete")
      import
      implicit none
      type(c_Dod) :: this
    end subroutine F_Dod_delete

    subroutine F_Dod_construct(this, a_transported_cell) &
      bind(C, name="c_Dod_construct")
      import
      implicit none
      type(c_Dod) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_transported_cell ! dimension(1:3,1:8)
    end subroutine F_Dod_construct

    function F_Dod_calculateVolume(this) result(a_dodecahedron_volume) &
      bind(C, name="c_Dod_calculateVolume")
      import
      implicit none
      type(c_Dod) :: this
      real(C_DOUBLE) :: a_dodecahedron_volume
    end function F_Dod_calculateVolume

    subroutine F_Dod_printToScreen(this) &
      bind(C, name="c_Dod_printToScreen")
      import
      implicit none
      type(c_Dod) :: this
    end subroutine F_Dod_printToScreen

    subroutine F_Dod_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_Dod_getBoundingPts")
      import
      implicit none
      type(c_Dod) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_Dod_getBoundingPts

  end interface

  contains

    subroutine Dod_class_new(this)
      implicit none
      type(Dod_type), intent(inout) :: this
      call F_Dod_new(this%c_object)
    end subroutine Dod_class_new

    impure elemental subroutine Dod_class_delete(this)
      implicit none
      type(Dod_type), intent(in) :: this
      call F_Dod_delete(this%c_object)
    end subroutine Dod_class_delete

    subroutine Dod_class_construct(this, a_transported_cell)
      implicit none
      type(Dod_type), intent(inout) :: this
      real(IRL_double), dimension(1:3,1:8), intent(in) :: a_transported_cell
      call F_Dod_construct(this%c_object, a_transported_cell)
    end subroutine Dod_class_construct

    function Dod_class_calculateVolume(this) result(a_dodecahedron_volume)
      implicit none
      type(Dod_type), intent(inout) :: this
      real(IRL_double) :: a_dodecahedron_volume
      a_dodecahedron_volume = F_Dod_calculateVolume(this%c_object)
    end function Dod_class_calculateVolume

    subroutine Dod_class_printToScreen(this)
      implicit none
      type(Dod_type), intent(in) :: this
      call F_Dod_printToScreen(this%c_object)
    end subroutine Dod_class_printToScreen

    subroutine Dod_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(Dod_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_Dod_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine Dod_class_getBoundingPts


end module f_Dod_class
