!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_octahedron_class.f90
!!
!! This file contains the Fortran interface for the
!! Octahedron class.

!> \brief A fortran type class that allows the creation of
!! IRL's Octahedron class along with enabling
!! some of its methods.
module f_Octa_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_Octa
    private
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_Octa

  type, public :: Octa_type
    type(c_Octa) :: c_object
  contains
    final :: Octa_class_delete
  end type Octa_type

  interface new
    module procedure Octa_class_new
  end interface
  interface construct
    module procedure Octa_class_construct
  end interface
  interface calculateVolume
    module procedure Octa_class_calculateVolume
  end interface
  interface printToScreen
    module procedure Octa_class_printToScreen
  end interface
  interface getBoundingPts
    module procedure Octa_class_getBoundingPts
  end interface

  interface

    subroutine F_Octa_new(this) &
      bind(C, name="c_Octa_new")
      import
      implicit none
      type(c_Octa) :: this
    end subroutine F_Octa_new

    subroutine F_Octa_delete(this) &
      bind(C, name="c_Octa_delete")
      import
      implicit none
      type(c_Octa) :: this
    end subroutine F_Octa_delete

    subroutine F_Octa_construct(this, a_transported_cell) &
      bind(C, name="c_Octa_construct")
      import
      implicit none
      type(c_Octa) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_transported_cell ! dimension(1:3,1:6)
    end subroutine F_Octa_construct

    function F_Octa_calculateVolume(this) result(a_octahedron_volume) &
      bind(C, name="c_Octa_calculateVolume")
      import
      implicit none
      type(c_Octa) :: this
      real(C_DOUBLE) :: a_octahedron_volume
    end function F_Octa_calculateVolume

    subroutine F_Octa_printToScreen(this) &
      bind(C, name="c_Octa_printToScreen")
      import
      implicit none
      type(c_Octa) :: this
    end subroutine F_Octa_printToScreen

    subroutine F_Octa_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_Octa_getBoundingPts")
      import
      implicit none
      type(c_Octa) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_Octa_getBoundingPts

  end interface

  contains

    subroutine Octa_class_new(this)
      implicit none
      type(Octa_type), intent(inout) :: this
      call F_Octa_new(this%c_object)
    end subroutine Octa_class_new

    impure elemental subroutine Octa_class_delete(this)
      implicit none
      type(Octa_type), intent(in) :: this
      call F_Octa_delete(this%c_object)
    end subroutine Octa_class_delete

    subroutine Octa_class_construct(this, a_transported_cell)
      implicit none
      type(Octa_type), intent(inout) :: this
      real(IRL_double), dimension(1:3,1:6), intent(in) :: a_transported_cell
      call F_Octa_construct(this%c_object, a_transported_cell)
    end subroutine Octa_class_construct

    function Octa_class_calculateVolume(this) result(a_octahedron_volume)
      implicit none
      type(Octa_type), intent(inout) :: this
      real(IRL_double) :: a_octahedron_volume
      a_octahedron_volume = F_Octa_calculateVolume(this%c_object)
    end function Octa_class_calculateVolume

    subroutine Octa_class_printToScreen(this)
      implicit none
      type(Octa_type), intent(in) :: this
      call F_Octa_printToScreen(this%c_object)
    end subroutine Octa_class_printToScreen

    subroutine Octa_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(Octa_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_Octa_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine Octa_class_getBoundingPts


end module f_Octa_class
