!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_SymmetricPyramid_class.f90
!!
!! This file contains the Fortran interface for the
!! SymmetricPyramid class.

!> \brief A fortran type class that allows the creation of
!! IRL's SymmetricPyramid class along with enabling
!! some of its methods.
module f_SymPyrmd_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_SymPyrmd
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_SymPyrmd

  type, public :: SymPyrmd_type
    type(c_SymPyrmd) :: c_object
  contains
    final :: SymPyrmd_class_delete
  end type SymPyrmd_type

  interface new
    module procedure SymPyrmd_class_new
  end interface

  interface construct
    module procedure SymPyrmd_class_construct
  end interface

  interface calculateVolume
    module procedure SymPyrmd_class_calculateVolume
  end interface 
 
  interface printToScreen
    module procedure SymPyrmd_class_printToScreen
  end interface

  interface getBoundingPts
    module procedure SymPyrmd_class_getBoundingPts
  end interface

  interface getPt
    module procedure SymPyrmd_class_getPt
  end interface

  interface

    subroutine F_SymPyrmd_new(this) &
      bind(C, name="c_SymPyrmd_new")
      import
      implicit none
      type(c_SymPyrmd) :: this
    end subroutine F_SymPyrmd_new

    subroutine F_SymPyrmd_delete(this) &
      bind(C, name="c_SymPyrmd_delete")
      import
      implicit none
      type(c_SymPyrmd) :: this
    end subroutine F_SymPyrmd_delete

    subroutine F_SymPyrmd_construct(this, a_sym_pyramid) &
      bind(C, name="c_SymPyrmd_construct")
      import
      implicit none
      type(c_SymPyrmd) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_sym_pyramid ! dimension(1:3,1:10)
    end subroutine F_SymPyrmd_construct

    function F_SymPyrmd_calculateVolume(this) result(a_volume) &
      bind(C, name="c_SymPyrmd_calculateVolume")
      import
      implicit none
      type(c_SymPyrmd) :: this
      real(C_DOUBLE) :: a_volume
    end function F_SymPyrmd_calculateVolume
    
    subroutine F_SymPyrmd_printToScreen(this) &
      bind(C, name="c_SymPyrmd_printToScreen")
      import
      implicit none
      type(c_SymPyrmd) :: this
    end subroutine F_SymPyrmd_printToScreen
    
    subroutine F_SymPyrmd_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_SymPyrmd_getBoundingPts")
      import
      implicit none
      type(c_SymPyrmd) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_SymPyrmd_getBoundingPts

    subroutine F_SymPyrmd_getPt(this, a_index, a_pt) &
      bind(C, name="c_SymPyrmd_getPt")
      import
      implicit none
      type(c_SymPyrmd) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(out) :: a_pt ! dimension(1:3)
    end subroutine F_SymPyrmd_getPt

  end interface


  contains

    impure elemental subroutine SymPyrmd_class_delete(this)
      implicit none
      type(SymPyrmd_type), intent(in) :: this
      call F_SymPyrmd_delete(this%c_object)
    end subroutine SymPyrmd_class_delete

    subroutine SymPyrmd_class_new(this)
      implicit none
      type(SymPyrmd_type), intent(inout) :: this
      call F_SymPyrmd_new(this%c_object)
    end subroutine SymPyrmd_class_new

    subroutine SymPyrmd_class_construct(this, a_sym_pyramid)
      implicit none
      type(SymPyrmd_type), intent(inout) :: this
      real(IRL_double), dimension(1:3,1:10), intent(in) :: a_sym_pyramid
      call F_SymPyrmd_construct(this%c_object, a_sym_pyramid)
    end subroutine SymPyrmd_class_construct
    
    function SymPyrmd_class_calculateVolume(this) result(a_volume)
      implicit none
      type(SymPyrmd_type), intent(inout) :: this
      real(IRL_double) :: a_volume
      a_volume = F_SymPyrmd_calculateVolume(this%c_object)
    end function SymPyrmd_class_calculateVolume
    
    subroutine SymPyrmd_class_printToScreen(this)
      implicit none
      type(SymPyrmd_type), intent(in) :: this
      call F_SymPyrmd_printToScreen(this%c_object)
    end subroutine SymPyrmd_class_printToScreen
    
    subroutine SymPyrmd_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(SymPyrmd_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_SymPyrmd_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine SymPyrmd_class_getBoundingPts

    function SymPyrmd_class_getPt(this, a_index) result(a_pt)
      implicit none
      type(SymPyrmd_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(3) :: a_pt
      call F_SymPyrmd_getPt(this%c_object, a_index, a_pt)
      return
    end function SymPyrmd_class_getPt

end module f_SymPyrmd_class
