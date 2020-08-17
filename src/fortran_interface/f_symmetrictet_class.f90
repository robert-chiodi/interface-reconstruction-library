!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_SymmetricTet_class.f90
!!
!! This file contains the Fortran interface for the
!! SymmetricTet class.

!> \brief A fortran type class that allows the creation of
!! IRL's SymmetricTet class along with enabling
!! some of its methods.
module f_SymTet_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_SymTet
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_SymTet

  type, public :: SymTet_type
    type(c_SymTet) :: c_object
  contains
    final :: SymTet_class_delete
  end type SymTet_type

  interface new
    module procedure SymTet_class_new
  end interface

  interface construct
    module procedure SymTet_class_construct
  end interface

  interface calculateVolume
    module procedure SymTet_class_calculateVolume
  end interface 
 
  interface printToScreen
    module procedure SymTet_class_printToScreen
  end interface

  interface getBoundingPts
    module procedure SymTet_class_getBoundingPts
  end interface

  interface getPt
    module procedure SymTet_class_getPt
  end interface

  interface

    subroutine F_SymTet_new(this) &
      bind(C, name="c_SymTet_new")
      import
      implicit none
      type(c_SymTet) :: this
    end subroutine F_SymTet_new

    subroutine F_SymTet_delete(this) &
      bind(C, name="c_SymTet_delete")
      import
      implicit none
      type(c_SymTet) :: this
    end subroutine F_SymTet_delete

    subroutine F_SymTet_construct(this, a_sym_tet) &
      bind(C, name="c_SymTet_construct")
      import
      implicit none
      type(c_SymTet) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_sym_tet ! dimension(1:3,1:8)
    end subroutine F_SymTet_construct

    function F_SymTet_calculateVolume(this) result(a_volume) &
      bind(C, name="c_SymTet_calculateVolume")
      import
      implicit none
      type(c_SymTet) :: this
      real(C_DOUBLE) :: a_volume
    end function F_SymTet_calculateVolume
    
    subroutine F_SymTet_printToScreen(this) &
      bind(C, name="c_SymTet_printToScreen")
      import
      implicit none
      type(c_SymTet) :: this
    end subroutine F_SymTet_printToScreen
    
    subroutine F_SymTet_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_SymTet_getBoundingPts")
      import
      implicit none
      type(c_SymTet) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_SymTet_getBoundingPts

    subroutine F_SymTet_getPt(this, a_index, a_pt) &
      bind(C, name="c_SymTet_getPt")
      import
      implicit none
      type(c_SymTet) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(out) :: a_pt ! dimension(1:3)
    end subroutine F_SymTet_getPt

  end interface


  contains

    impure elemental subroutine SymTet_class_delete(this)
      implicit none
      type(SymTet_type), intent(in) :: this
      call F_SymTet_delete(this%c_object)
    end subroutine SymTet_class_delete

    subroutine SymTet_class_new(this)
      implicit none
      type(SymTet_type), intent(inout) :: this
      call F_SymTet_new(this%c_object)
    end subroutine SymTet_class_new

    subroutine SymTet_class_construct(this, a_sym_tet)
      implicit none
      type(SymTet_type), intent(inout) :: this
      real(IRL_double), dimension(1:3,1:8), intent(in) :: a_sym_tet
      call F_SymTet_construct(this%c_object, a_sym_tet)
    end subroutine SymTet_class_construct
    
    function SymTet_class_calculateVolume(this) result(a_volume)
      implicit none
      type(SymTet_type), intent(inout) :: this
      real(IRL_double) :: a_volume
      a_volume = F_SymTet_calculateVolume(this%c_object)
    end function SymTet_class_calculateVolume
    
    subroutine SymTet_class_printToScreen(this)
      implicit none
      type(SymTet_type), intent(in) :: this
      call F_SymTet_printToScreen(this%c_object)
    end subroutine SymTet_class_printToScreen
    
    subroutine SymTet_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(SymTet_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_SymTet_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine SymTet_class_getBoundingPts

    function SymTet_class_getPt(this, a_index) result(a_pt)
      implicit none
      type(SymTet_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(3) :: a_pt
      call F_SymTet_getPt(this%c_object, a_index, a_pt)
      return
    end function SymTet_class_getPt

end module f_SymTet_class
