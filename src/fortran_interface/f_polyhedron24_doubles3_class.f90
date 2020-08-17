!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_Polyhedron24_doubles3_class.f90
!!
!! This file contains the Fortran interface for the
!! Polyhedron24_doubles3 class.

!> \brief A fortran type class that allows the creation of
!! IRL's Polyhedron24_doubles3 class along with enabling
!! some of its methods.
module f_Poly24_d3_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_Poly24_d3
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_Poly24_d3

  type, public :: Poly24_d3_type
    type(c_Poly24_d3) :: c_object
  contains
    final :: Poly24_d3_class_delete
  end type Poly24_d3_type

  interface new
    module procedure Poly24_d3_class_new
  end interface
  interface construct
    module procedure Poly24_d3_class_construct
  end interface
  interface adjustCapToMatchVolume
    module procedure Poly24_d3_class_adjustCapToMatchVolume
  end interface
  interface getBoundingPts
    module procedure Poly24_d3_class_getBoundingPts
  end interface
  interface getPt
    module procedure Poly24_d3_class_getPt
  end interface
  interface setPt
    module procedure Poly24_d3_class_setPt
  end interface
  interface getData
    module procedure Poly24_d3_class_getData
  end interface
  interface setData
    module procedure Poly24_d3_class_setData
  end interface

  interface

    subroutine F_Poly24_d3_new(this) &
      bind(C, name="c_Poly24_d3_new")
      import
      implicit none
      type(c_Poly24_d3) :: this
    end subroutine F_Poly24_d3_new

    subroutine F_Poly24_d3_delete(this) &
      bind(C, name="c_Poly24_d3_delete")
      import
      implicit none
      type(c_Poly24_d3) :: this
    end subroutine F_Poly24_d3_delete

    subroutine F_Poly24_d3_construct(this, a_polyhedron24, a_data) &
      bind(C, name="c_Poly24_d3_construct")
      import
      implicit none
      type(c_Poly24_d3) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_polyhedron24 ! dimension(1:3,1:14)
      real(C_DOUBLE), dimension(*), intent(in) :: a_data ! dimension(1:3,1:14)
    end subroutine F_Poly24_d3_construct

    subroutine F_Poly24_d3_adjustCapToMatchVolume(this, a_correct_signed_volume) &
      bind(C, name="c_Poly24_d3_adjustCapToMatchVolume")
      import
      implicit none
      type(c_Poly24_d3) :: this
      real(C_DOUBLE), intent(in) :: a_correct_signed_volume ! scalar
    end subroutine F_Poly24_d3_adjustCapToMatchVolume

    subroutine F_Poly24_d3_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_Poly24_d3_getBoundingPts")
      import
      implicit none
      type(c_Poly24_d3) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_Poly24_d3_getBoundingPts

    subroutine F_Poly24_d3_getPt(this, a_index, a_pt) &
      bind(C, name="c_Poly24_d3_getPt")
      import
      implicit none
      type(c_Poly24_d3) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(out) :: a_pt ! dimension(1:3)
    end subroutine F_Poly24_d3_getPt

    subroutine F_Poly24_d3_setPt(this, a_index, a_pt) &
      bind(C, name="c_Poly24_d3_setPt")
      import
      implicit none
      type(c_Poly24_d3) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(in) :: a_pt ! dimension(1:3)
    end subroutine F_Poly24_d3_setPt

    subroutine F_Poly24_d3_getData(this, a_index, a_data) &
      bind(C, name="c_Poly24_d3_getData")
      import
      implicit none
      type(c_Poly24_d3) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(out) :: a_data ! dimension(1:3)
    end subroutine F_Poly24_d3_getData

    subroutine F_Poly24_d3_setData(this, a_index, a_data) &
      bind(C, name="c_Poly24_d3_setData")
      import
      implicit none
      type(c_Poly24_d3) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(in) :: a_data ! dimension(1:3)
    end subroutine F_Poly24_d3_setData

  end interface


  contains

    subroutine Poly24_d3_class_new(this)
      implicit none
      type(Poly24_d3_type), intent(inout) :: this
      call F_Poly24_d3_new(this%c_object)
    end subroutine Poly24_d3_class_new

    impure elemental subroutine Poly24_d3_class_delete(this)
      implicit none
      type(Poly24_d3_type), intent(in) :: this
      call F_Poly24_d3_delete(this%c_object)
    end subroutine Poly24_d3_class_delete

    subroutine Poly24_d3_class_construct(this, a_polyhedron24, a_data)
      implicit none
      type(Poly24_d3_type), intent(inout) :: this
      real(IRL_double), dimension(1:3,1:14), intent(in) :: a_polyhedron24
      real(IRL_double), dimension(1:3,1:14), intent(in) :: a_data
      call F_Poly24_d3_construct(this%c_object, a_polyhedron24, a_data)
    end subroutine Poly24_d3_class_construct

    subroutine Poly24_d3_class_adjustCapToMatchVolume(this, a_correct_signed_volume)
      implicit none
      type(Poly24_d3_type), intent(inout) :: this
      real(IRL_double), intent(in) :: a_correct_signed_volume
      call F_Poly24_d3_adjustCapToMatchVolume(this%c_object, a_correct_signed_volume)
    end subroutine Poly24_d3_class_adjustCapToMatchVolume

    subroutine Poly24_d3_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(Poly24_d3_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_Poly24_d3_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine Poly24_d3_class_getBoundingPts

    function Poly24_d3_class_getPt(this, a_index) result(a_pt)
      implicit none
      type(Poly24_d3_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(3) :: a_pt
      call F_Poly24_d3_getPt(this%c_object, a_index, a_pt)
      return
    end function Poly24_d3_class_getPt

    subroutine Poly24_d3_class_setPt(this, a_index, a_pt)
      implicit none
      type(Poly24_d3_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(3) :: a_pt
      call F_Poly24_d3_setPt(this%c_object, a_index, a_pt)
    end subroutine Poly24_d3_class_setPt

    function Poly24_d3_class_getData(this, a_index) result(a_data)
      implicit none
      type(Poly24_d3_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(3) :: a_data
      call F_Poly24_d3_getData(this%c_object, a_index, a_data)
      return
    end function Poly24_d3_class_getData

    subroutine Poly24_d3_class_setData(this, a_index, a_data)
      implicit none
      type(Poly24_d3_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(3) :: a_data
      call F_Poly24_d3_setData(this%c_object, a_index, a_data)
    end subroutine Poly24_d3_class_setData

end module f_Poly24_d3_class
