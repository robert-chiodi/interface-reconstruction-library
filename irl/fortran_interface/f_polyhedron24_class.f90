!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_Polyhedron24_class.f90
!!
!! This file contains the Fortran interface for the
!! Polyhedron24 class.

!> \brief A fortran type class that allows the creation of
!! IRL's Polyhedron24 class along with enabling
!! some of its methods.
module f_Poly24_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_Poly24
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_Poly24

  type, public :: Poly24_type
    type(c_Poly24) :: c_object
  contains
    final :: Poly24_class_delete
  end type Poly24_type

  interface new
    module procedure Poly24_class_new
  end interface
  interface construct
    module procedure Poly24_class_construct
  end interface
  interface adjustCapToMatchVolume
    module procedure Poly24_class_adjustCapToMatchVolume
  end interface
  interface getBoundingPts
    module procedure Poly24_class_getBoundingPts
  end interface
  interface getPt
    module procedure Poly24_class_getPt
  end interface
  interface setPt
    module procedure Poly24_class_setPt
  end interface


  interface

    subroutine F_Poly24_new(this) &
      bind(C, name="c_Poly24_new")
      import
      implicit none
      type(c_Poly24) :: this
    end subroutine F_Poly24_new

    subroutine F_Poly24_delete(this) &
      bind(C, name="c_Poly24_delete")
      import
      implicit none
      type(c_Poly24) :: this
    end subroutine F_Poly24_delete

    subroutine F_Poly24_construct(this, a_polyhedron24) &
      bind(C, name="c_Poly24_construct")
      import
      implicit none
      type(c_Poly24) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_polyhedron24 ! dimension(1:3,1:14)
    end subroutine F_Poly24_construct

    subroutine F_Poly24_adjustCapToMatchVolume(this, a_correct_signed_volume) &
      bind(C, name="c_Poly24_adjustCapToMatchVolume")
      import
      implicit none
      type(c_Poly24) :: this
      real(C_DOUBLE), intent(in) :: a_correct_signed_volume ! scalar
    end subroutine F_Poly24_adjustCapToMatchVolume

    subroutine F_Poly24_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_Poly24_getBoundingPts")
      import
      implicit none
      type(c_Poly24) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_Poly24_getBoundingPts

    subroutine F_Poly24_getPt(this, a_index, a_pt) &
      bind(C, name="c_Poly24_getPt")
      import
      implicit none
      type(c_Poly24) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(out) :: a_pt ! dimension(1:3)
    end subroutine F_Poly24_getPt

    subroutine F_Poly24_setPt(this, a_index, a_pt) &
      bind(C, name="c_Poly24_setPt")
      import
      implicit none
      type(c_Poly24) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(in) :: a_pt ! dimension(1:3)
    end subroutine F_Poly24_setPt
  end interface


  contains

    subroutine Poly24_class_new(this)
      implicit none
      type(Poly24_type), intent(inout) :: this
      call F_Poly24_new(this%c_object)
    end subroutine Poly24_class_new

    impure elemental subroutine Poly24_class_delete(this)
      implicit none
      type(Poly24_type), intent(in) :: this
      call F_Poly24_delete(this%c_object)
    end subroutine Poly24_class_delete

    subroutine Poly24_class_construct(this, a_polyhedron24)
      implicit none
      type(Poly24_type), intent(inout) :: this
      real(IRL_double), dimension(1:3,1:14), intent(in) :: a_polyhedron24
      call F_Poly24_construct(this%c_object, a_polyhedron24)
    end subroutine Poly24_class_construct

    subroutine Poly24_class_adjustCapToMatchVolume(this, a_correct_signed_volume)
      implicit none
      type(Poly24_type), intent(inout) :: this
      real(IRL_double), intent(in) :: a_correct_signed_volume
      call F_Poly24_adjustCapToMatchVolume(this%c_object, a_correct_signed_volume)
    end subroutine Poly24_class_adjustCapToMatchVolume

    subroutine Poly24_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(Poly24_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_Poly24_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine Poly24_class_getBoundingPts

    function Poly24_class_getPt(this, a_index) result(a_pt)
      implicit none
      type(Poly24_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(3) :: a_pt
      call F_Poly24_getPt(this%c_object, a_index, a_pt)
      return
    end function Poly24_class_getPt

    subroutine Poly24_class_setPt(this, a_index, a_pt)
      implicit none
      type(Poly24_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(3) :: a_pt
      call F_Poly24_setPt(this%c_object, a_index, a_pt)
    end subroutine Poly24_class_setPt

end module f_Poly24_class
