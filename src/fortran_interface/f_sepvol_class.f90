!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_SepVol_class.f90
!!
!! This file contains the Fortran interface for volume moments classes.

!> \brief A fortran type class that allows the creation of
!! IRL's SeparatedMoments<Volume>
!! class along with enabling some of its methods.
module f_SepVol_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_SepVol
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_SepVol

  type, public :: SepVol_type
    type(c_SepVol) :: c_object
  contains
    final :: SepVol_class_delete
  end type SepVol_type

  interface new
    module procedure SepVol_class_new
  end interface
  interface construct
    module procedure SepVol_class_construct
  end interface
  interface getVolume
    module procedure SepVol_class_getVolume
  end interface

  interface

    subroutine F_SepVol_new(this) &
      bind(C, name="c_SepVol_new")
      import
      implicit none
      type(c_SepVol) :: this
    end subroutine F_SepVol_new

    subroutine F_SepVol_delete(this) &
      bind(C, name="c_SepVol_delete")
      import
      implicit none
      type(c_SepVol) :: this
    end subroutine F_SepVol_delete

    subroutine F_SepVol_construct(this, a_moments_list) &
      bind(C, name="c_SepVol_construct")
      import
      implicit none
      type(c_SepVol) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_moments_list
    end subroutine F_SepVol_construct

    function F_SepVol_getVolume(this, a_index) result(a_volume) &
      bind(C, name="c_SepVol_getVolume")
      import
      implicit none
      type(c_SepVol) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE) :: a_volume
    end function F_SepVol_getVolume

  end interface


  contains

    subroutine SepVol_class_new(this)
      implicit none
      type(SepVol_type), intent(inout) :: this
      call F_SepVol_new(this%c_object)
    end subroutine SepVol_class_new

    impure elemental subroutine SepVol_class_delete(this)
      implicit none
      type(SepVol_type), intent(in) :: this
      call F_SepVol_delete(this%c_object)
    end subroutine SepVol_class_delete

    subroutine SepVol_class_construct(this, a_moments_list)
      implicit none
      type(SepVol_type), intent(inout) :: this
      real(IRL_double), dimension(2), intent(in) :: a_moments_list
      call F_SepVol_construct(this%c_object, a_moments_list)
    end subroutine SepVol_class_construct

    function SepVol_class_getVolume(this, a_index) result(a_volume)
      implicit none
      type(SepVol_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_index
      real(IRL_double) :: a_volume
      a_volume = F_SepVol_getVolume(this%c_object, a_index)
      return
    end function SepVol_class_getVolume

end module f_SepVol_class
