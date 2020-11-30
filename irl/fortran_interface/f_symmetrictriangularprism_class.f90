!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_SymmetricTriangularPrism_class.f90
!!
!! This file contains the Fortran interface for the
!! SymmetricTriangularPrism class.

!> \brief A fortran type class that allows the creation of
!! IRL's SymmetricTriangularPrism class along with enabling
!! some of its methods.
module f_SymTriPrism_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_SymTriPrism
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_SymTriPrism

  type, public :: SymTriPrism_type
    type(c_SymTriPrism) :: c_object
  contains
    final :: SymTriPrism_class_delete
  end type SymTriPrism_type

  interface new
    module procedure SymTriPrism_class_new
  end interface

  interface construct
    module procedure SymTriPrism_class_construct
  end interface

  interface calculateVolume
    module procedure SymTriPrism_class_calculateVolume
   end interface

  interface adjustCapToMatchVolume
    module procedure SymTriPrism_class_adjustCapToMatchVolume
  end interface 
 
  interface printToScreen
    module procedure SymTriPrism_class_printToScreen
  end interface

  interface getBoundingPts
    module procedure SymTriPrism_class_getBoundingPts
  end interface

  interface getPt
    module procedure SymTriPrism_class_getPt
  end interface

  interface

    subroutine F_SymTriPrism_new(this) &
      bind(C, name="c_SymTriPrism_new")
      import
      implicit none
      type(c_SymTriPrism) :: this
    end subroutine F_SymTriPrism_new

    subroutine F_SymTriPrism_delete(this) &
      bind(C, name="c_SymTriPrism_delete")
      import
      implicit none
      type(c_SymTriPrism) :: this
    end subroutine F_SymTriPrism_delete

    subroutine F_SymTriPrism_construct(this, a_sym_tri_prism) &
      bind(C, name="c_SymTriPrism_construct")
      import
      implicit none
      type(c_SymTriPrism) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_sym_tri_prism ! dimension(1:3,1:11)
    end subroutine F_SymTriPrism_construct

    function F_SymTriPrism_calculateVolume(this) result(a_volume) &
      bind(C, name="c_SymTriPrism_calculateVolume")
      import
      implicit none
      type(c_SymTriPrism) :: this
      real(C_DOUBLE) :: a_volume
    end function F_SymTriPrism_calculateVolume

    subroutine F_SymTriPrism_adjustCapToMatchVolume(this, a_correct_signed_volume) &
      bind(C, name="c_SymTriPrism_adjustCapToMatchVolume")
      import
      implicit none
      type(c_SymTriPrism) :: this
      real(C_DOUBLE), intent(in) :: a_correct_signed_volume ! scalar
    end subroutine F_SymTriPrism_adjustCapToMatchVolume
    
    subroutine F_SymTriPrism_printToScreen(this) &
      bind(C, name="c_SymTriPrism_printToScreen")
      import
      implicit none
      type(c_SymTriPrism) :: this
    end subroutine F_SymTriPrism_printToScreen
    
    subroutine F_SymTriPrism_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_SymTriPrism_getBoundingPts")
      import
      implicit none
      type(c_SymTriPrism) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_SymTriPrism_getBoundingPts

    subroutine F_SymTriPrism_getPt(this, a_index, a_pt) &
      bind(C, name="c_SymTriPrism_getPt")
      import
      implicit none
      type(c_SymTriPrism) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(out) :: a_pt ! dimension(1:3)
    end subroutine F_SymTriPrism_getPt

  end interface


  contains

    impure elemental subroutine SymTriPrism_class_delete(this)
      implicit none
      type(SymTriPrism_type), intent(in) :: this
      call F_SymTriPrism_delete(this%c_object)
    end subroutine SymTriPrism_class_delete

    subroutine SymTriPrism_class_new(this)
      implicit none
      type(SymTriPrism_type), intent(inout) :: this
      call F_SymTriPrism_new(this%c_object)
    end subroutine SymTriPrism_class_new

    subroutine SymTriPrism_class_construct(this, a_sym_tri_prism)
      implicit none
      type(SymTriPrism_type), intent(inout) :: this
      real(IRL_double), dimension(1:3,1:11), intent(in) :: a_sym_tri_prism
      call F_SymTriPrism_construct(this%c_object, a_sym_tri_prism)
    end subroutine SymTriPrism_class_construct
    
    function SymTriPrism_class_calculateVolume(this) result(a_volume)
      implicit none
      type(SymTriPrism_type), intent(inout) :: this
      real(IRL_double) :: a_volume
      a_volume = F_SymTriPrism_calculateVolume(this%c_object)
    end function SymTriPrism_class_calculateVolume

    subroutine SymTriPrism_class_adjustCapToMatchVolume(this, a_correct_signed_volume)
      implicit none
      type(SymTriPrism_type), intent(inout) :: this
      real(IRL_double), intent(in) :: a_correct_signed_volume
      call F_SymTriPrism_adjustCapToMatchVolume(this%c_object, a_correct_signed_volume)
    end subroutine SymTriPrism_class_adjustCapToMatchVolume
    
    subroutine SymTriPrism_class_printToScreen(this)
      implicit none
      type(SymTriPrism_type), intent(in) :: this
      call F_SymTriPrism_printToScreen(this%c_object)
    end subroutine SymTriPrism_class_printToScreen
    
    subroutine SymTriPrism_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(SymTriPrism_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_SymTriPrism_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine SymTriPrism_class_getBoundingPts

    function SymTriPrism_class_getPt(this, a_index) result(a_pt)
      implicit none
      type(SymTriPrism_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(3) :: a_pt
      call F_SymTriPrism_getPt(this%c_object, a_index, a_pt)
      return
    end function SymTriPrism_class_getPt

end module f_SymTriPrism_class
