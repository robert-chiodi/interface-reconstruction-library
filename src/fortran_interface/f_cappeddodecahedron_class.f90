!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_CappedDodecahedron_class.f90
!!
!! This file contains the Fortran interface for the
!! CappedDodecahedron class.

!> \brief A fortran type class that allows the creation of
!! IRL's CappedDodecahedron class along with enabling
!! some of its methods.
module f_CapDod_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_CapDod
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_CapDod

  type, public :: CapDod_type
    type(c_CapDod) :: c_object
  contains
    final :: CapDod_class_delete
  end type CapDod_type

  interface new
    module procedure CapDod_class_new
  end interface

  interface construct
    module procedure CapDod_class_construct
  end interface

  interface adjustCapToMatchVolume
    module procedure CapDod_class_adjustCapToMatchVolume
  end interface

  interface getBoundingPts
    module procedure CapDod_class_getBoundingPts
  end interface

  interface getPt
    module procedure CapDod_class_getPt
  end interface

  interface

    subroutine F_CapDod_new(this) &
      bind(C, name="c_CapDod_new")
      import
      implicit none
      type(c_CapDod) :: this
    end subroutine F_CapDod_new

    subroutine F_CapDod_delete(this) &
      bind(C, name="c_CapDod_delete")
      import
      implicit none
      type(c_CapDod) :: this
    end subroutine F_CapDod_delete

    subroutine F_CapDod_construct(this, a_dodecahedron) &
      bind(C, name="c_CapDod_construct")
      import
      implicit none
      type(c_CapDod) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_dodecahedron ! dimension(1:3,1:9)
    end subroutine F_CapDod_construct

    subroutine F_CapDod_adjustCapToMatchVolume(this, a_correct_signed_volume) &
      bind(C, name="c_CapDod_adjustCapToMatchVolume")
      import
      implicit none
      type(c_CapDod) :: this
      real(C_DOUBLE), intent(in) :: a_correct_signed_volume ! scalar
    end subroutine F_CapDod_adjustCapToMatchVolume

    subroutine F_CapDod_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_CapDod_getBoundingPts")
      import
      implicit none
      type(c_CapDod) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_CapDod_getBoundingPts

    subroutine F_CapDod_getPt(this, a_index, a_pt) &
      bind(C, name="c_CapDod_getPt")
      import
      implicit none
      type(c_CapDod) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(out) :: a_pt ! dimension(1:3)
    end subroutine F_CapDod_getPt

  end interface


  contains

    impure elemental subroutine CapDod_class_delete(this)
      implicit none
      type(CapDod_type), intent(in) :: this
      call F_CapDod_delete(this%c_object)
    end subroutine CapDod_class_delete

    subroutine CapDod_class_new(this)
      implicit none
      type(CapDod_type), intent(inout) :: this
      call F_CapDod_new(this%c_object)
    end subroutine CapDod_class_new

    subroutine CapDod_class_construct(this, a_dodecahedron)
      implicit none
      type(CapDod_type), intent(inout) :: this
      real(IRL_double), dimension(1:3,1:9), intent(in) :: a_dodecahedron
      call F_CapDod_construct(this%c_object, a_dodecahedron)
    end subroutine CapDod_class_construct

    subroutine CapDod_class_adjustCapToMatchVolume(this, a_correct_signed_volume)
      implicit none
      type(CapDod_type), intent(inout) :: this
      real(IRL_double), intent(in) :: a_correct_signed_volume
      call F_CapDod_adjustCapToMatchVolume(this%c_object, a_correct_signed_volume)
    end subroutine CapDod_class_adjustCapToMatchVolume

    subroutine CapDod_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(CapDod_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_CapDod_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine CapDod_class_getBoundingPts

    function CapDod_class_getPt(this, a_index) result(a_pt)
      implicit none
      type(CapDod_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(3) :: a_pt
      call F_CapDod_getPt(this%c_object, a_index, a_pt)
      return
    end function CapDod_class_getPt

end module f_CapDod_class
