!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_SymmetricHexahedron_class.f90
!!
!! This file contains the Fortran interface for the
!! SymmetricHexahedron class.

!> \brief A fortran type class that allows the creation of
!! IRL's SymmetricHexahedron class along with enabling
!! some of its methods.
module f_SymHex_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_SymHex
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_SymHex

  type, public :: SymHex_type
    type(c_SymHex) :: c_object
  contains
    final :: SymHex_class_delete
  end type SymHex_type

  interface new
    module procedure SymHex_class_new
  end interface

  interface construct
    module procedure SymHex_class_construct
  end interface

  interface calculateVolume
    module procedure SymHex_class_calculateVolume
  end interface 

  interface adjustCapToMatchVolume
    module procedure SymHex_class_adjustCapToMatchVolume
  end interface
  
  interface printToScreen
    module procedure SymHex_class_printToScreen
  end interface

  interface getBoundingPts
    module procedure SymHex_class_getBoundingPts
  end interface

  interface getPt
    module procedure SymHex_class_getPt
  end interface

  interface

    subroutine F_SymHex_new(this) &
      bind(C, name="c_SymHex_new")
      import
      implicit none
      type(c_SymHex) :: this
    end subroutine F_SymHex_new

    subroutine F_SymHex_delete(this) &
      bind(C, name="c_SymHex_delete")
      import
      implicit none
      type(c_SymHex) :: this
    end subroutine F_SymHex_delete

    subroutine F_SymHex_construct(this, a_sym_hex) &
      bind(C, name="c_SymHex_construct")
      import
      implicit none
      type(c_SymHex) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_sym_hex ! dimension(1:3,1:14)
    end subroutine F_SymHex_construct

    function F_SymHex_calculateVolume(this) result(a_volume) &
      bind(C, name="c_SymHex_calculateVolume")
      import
      implicit none
      type(c_SymHex) :: this
      real(C_DOUBLE) :: a_volume
    end function F_SymHex_calculateVolume

    subroutine F_SymHex_adjustCapToMatchVolume(this, a_correct_signed_volume) &
      bind(C, name="c_SymHex_adjustCapToMatchVolume")
      import
      implicit none
      type(c_SymHex) :: this
      real(C_DOUBLE), intent(in) :: a_correct_signed_volume ! scalar
    end subroutine F_SymHex_adjustCapToMatchVolume
    
    subroutine F_SymHex_printToScreen(this) &
      bind(C, name="c_SymHex_printToScreen")
      import
      implicit none
      type(c_SymHex) :: this
    end subroutine F_SymHex_printToScreen
    
    subroutine F_SymHex_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_SymHex_getBoundingPts")
      import
      implicit none
      type(c_SymHex) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_SymHex_getBoundingPts

    subroutine F_SymHex_getPt(this, a_index, a_pt) &
      bind(C, name="c_SymHex_getPt")
      import
      implicit none
      type(c_SymHex) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(out) :: a_pt ! dimension(1:3)
    end subroutine F_SymHex_getPt

  end interface


  contains

    impure elemental subroutine SymHex_class_delete(this)
      implicit none
      type(SymHex_type), intent(in) :: this
      call F_SymHex_delete(this%c_object)
    end subroutine SymHex_class_delete

    subroutine SymHex_class_new(this)
      implicit none
      type(SymHex_type), intent(inout) :: this
      call F_SymHex_new(this%c_object)
    end subroutine SymHex_class_new

    subroutine SymHex_class_construct(this, a_sym_hex)
      implicit none
      type(SymHex_type), intent(inout) :: this
      real(IRL_double), dimension(1:3,1:14), intent(in) :: a_sym_hex
      call F_SymHex_construct(this%c_object, a_sym_hex)
    end subroutine SymHex_class_construct
    
    function SymHex_class_calculateVolume(this) result(a_volume)
      implicit none
      type(SymHex_type), intent(inout) :: this
      real(IRL_double) :: a_volume
      a_volume = F_SymHex_calculateVolume(this%c_object)
    end function SymHex_class_calculateVolume

    subroutine SymHex_class_adjustCapToMatchVolume(this, a_correct_signed_volume)
      implicit none
      type(SymHex_type), intent(inout) :: this
      real(IRL_double), intent(in) :: a_correct_signed_volume
      call F_SymHex_adjustCapToMatchVolume(this%c_object, a_correct_signed_volume)
    end subroutine SymHex_class_adjustCapToMatchVolume
    
    subroutine SymHex_class_printToScreen(this)
      implicit none
      type(SymHex_type), intent(in) :: this
      call F_SymHex_printToScreen(this%c_object)
    end subroutine SymHex_class_printToScreen
    
    subroutine SymHex_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(SymHex_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_SymHex_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine SymHex_class_getBoundingPts

    function SymHex_class_getPt(this, a_index) result(a_pt)
      implicit none
      type(SymHex_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(3) :: a_pt
      call F_SymHex_getPt(this%c_object, a_index, a_pt)
      return
    end function SymHex_class_getPt

end module f_SymHex_class
