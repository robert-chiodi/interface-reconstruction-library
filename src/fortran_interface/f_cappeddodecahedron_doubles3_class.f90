!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_CappedDodecahedron_doubles3_class.f90
!!
!! This file contains the Fortran interface for the
!! CappedDodecahedron_doubles3 class.

!> \brief A fortran type class that allows the creation of
!! IRL's CappedDodecahedron_doubles3 class along with enabling
!! some of its methods.
module f_CapDod_d3_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_CapDod_d3
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_CapDod_d3

  type, public :: CapDod_d3_type
    type(c_CapDod_d3) :: c_object
  contains
    final :: CapDod_d3_class_delete
  end type CapDod_d3_type

  interface new
    module procedure CapDod_d3_class_new
  end interface

  interface construct
    module procedure CapDod_d3_class_construct
  end interface

  interface adjustCapToMatchVolume
    module procedure CapDod_d3_class_adjustCapToMatchVolume
  end interface

  interface getBoundingPts
    module procedure CapDod_d3_class_getBoundingPts
  end interface

  interface getPt
    module procedure CapDod_d3_class_getPt
  end interface

  interface setPt
    module procedure CapDod_d3_class_setPt
  end interface

  interface getData
    module procedure CapDod_d3_class_getData
  end interface

  interface setData
    module procedure CapDod_d3_class_setData
  end interface

  interface

    subroutine F_CapDod_d3_new(this) &
      bind(C, name="c_CapDod_d3_new")
      import
      implicit none
      type(c_CapDod_d3) :: this
    end subroutine F_CapDod_d3_new

    subroutine F_CapDod_d3_delete(this) &
      bind(C, name="c_CapDod_d3_delete")
      import
      implicit none
      type(c_CapDod_d3) :: this
    end subroutine F_CapDod_d3_delete

    subroutine F_CapDod_d3_construct(this, a_dodecahedron,a_attached_data) &
      bind(C, name="c_CapDod_d3_construct")
      import
      implicit none
      type(c_CapDod_d3) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_dodecahedron ! dimension(1:3,1:9)
      real(C_DOUBLE), dimension(*), intent(in) :: a_attached_data ! dimension(1:3,1:9)
    end subroutine F_CapDod_d3_construct

    subroutine F_CapDod_d3_adjustCapToMatchVolume(this, a_correct_signed_volume) &
      bind(C, name="c_CapDod_d3_adjustCapToMatchVolume")
      import
      implicit none
      type(c_CapDod_d3) :: this
      real(C_DOUBLE), intent(in) :: a_correct_signed_volume ! scalar
    end subroutine F_CapDod_d3_adjustCapToMatchVolume

    subroutine F_CapDod_d3_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_CapDod_d3_getBoundingPts")
      import
      implicit none
      type(c_CapDod_d3) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_CapDod_d3_getBoundingPts

    subroutine F_CapDod_d3_getPt(this, a_index, a_pt) &
      bind(C, name="c_CapDod_d3_getPt")
      import
      implicit none
      type(c_CapDod_d3) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(out) :: a_pt ! dimension(1:3)
    end subroutine F_CapDod_d3_getPt

    subroutine F_CapDod_d3_setPt(this, a_index, a_pt) &
      bind(C, name="c_CapDod_d3_setPt")
      import
      implicit none
      type(c_CapDod_d3) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(in) :: a_pt ! dimension(1:3)
    end subroutine F_CapDod_d3_setPt

    subroutine F_CapDod_d3_getData(this, a_index, a_data) &
      bind(C, name="c_CapDod_d3_getData")
      import
      implicit none
      type(c_CapDod_d3) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(out) :: a_data ! dimension(1:3)
    end subroutine F_CapDod_d3_getData

    subroutine F_CapDod_d3_setData(this, a_index, a_data) &
      bind(C, name="c_CapDod_d3_setData")
      import
      implicit none
      type(c_CapDod_d3) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(in) :: a_data ! dimension(1:3)
    end subroutine F_CapDod_d3_setData

  end interface


  contains

    impure elemental subroutine CapDod_d3_class_delete(this)
      implicit none
      type(CapDod_d3_type), intent(in) :: this
      call F_CapDod_d3_delete(this%c_object)
    end subroutine CapDod_d3_class_delete

    subroutine CapDod_d3_class_new(this)
      implicit none
      type(CapDod_d3_type), intent(inout) :: this
      call F_CapDod_d3_new(this%c_object)
    end subroutine CapDod_d3_class_new

    subroutine CapDod_d3_class_construct(this, a_dodecahedron, a_attached_data)
      implicit none
      type(CapDod_d3_type), intent(inout) :: this
      real(IRL_double), dimension(1:3,1:9), intent(in) :: a_dodecahedron
      real(IRL_double), dimension(1:3,1:9), intent(in) :: a_attached_data
      call F_CapDod_d3_construct(this%c_object, a_dodecahedron, a_attached_data)
    end subroutine CapDod_d3_class_construct

    subroutine CapDod_d3_class_adjustCapToMatchVolume(this, a_correct_signed_volume)
      implicit none
      type(CapDod_d3_type), intent(inout) :: this
      real(IRL_double), intent(in) :: a_correct_signed_volume
      call F_CapDod_d3_adjustCapToMatchVolume(this%c_object, a_correct_signed_volume)
    end subroutine CapDod_d3_class_adjustCapToMatchVolume

    subroutine CapDod_d3_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(CapDod_d3_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_CapDod_d3_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine CapDod_d3_class_getBoundingPts

    function CapDod_d3_class_getPt(this, a_index) result(a_pt)
      implicit none
      type(CapDod_d3_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(3) :: a_pt
      call F_CapDod_d3_getPt(this%c_object, a_index, a_pt)
      return
    end function CapDod_d3_class_getPt

    subroutine CapDod_d3_class_setPt(this, a_index, a_pt)
      implicit none
      type(CapDod_d3_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(3), intent(in) :: a_pt
      call F_CapDod_d3_setPt(this%c_object, a_index, a_pt)
      return
    end subroutine CapDod_d3_class_setPt

    function CapDod_d3_class_getData(this, a_index) result(a_data)
      implicit none
      type(CapDod_d3_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(3) :: a_data
      call F_CapDod_d3_getData(this%c_object, a_index, a_data)
      return
    end function CapDod_d3_class_getData

    subroutine CapDod_d3_class_setData(this, a_index, a_data)
      implicit none
      type(CapDod_d3_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(3), intent(in) :: a_data
      call F_CapDod_d3_setData(this%c_object, a_index, a_data)
      return
    end subroutine CapDod_d3_class_setData



end module f_CapDod_d3_class
