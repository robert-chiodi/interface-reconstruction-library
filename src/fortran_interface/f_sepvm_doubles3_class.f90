!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_SepVM_doubles3_class.f90
!!
!! This file contains the Fortran interface for volume moments classes.

!> \brief A fortran type class that allows the creation of
!! IRL's SeparatedMoments<VolumeMoments>
!! class along with enabling some of its methods.
module f_SepVM_d3_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_SepVM_d3
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_SepVM_d3

  type, public :: SepVM_d3_type
    type(c_SepVM_d3) :: c_object
  contains
    final :: SepVM_d3_class_delete
  end type SepVM_d3_type

  interface new
    module procedure SepVM_d3_class_new
  end interface
  interface normalizeByVolume
    module procedure SepVM_d3_class_normalizeByVolume
  end interface
  interface multiplyByVolume
    module procedure SepVM_d3_class_multiplyByVolume
  end interface
  interface getVolume
    module procedure SepVM_d3_class_getVolume
  end interface
  interface getCentroid
    module procedure SepVM_d3_class_getCentroid
  end interface
  interface getData
    module procedure SepVM_d3_class_getData
  end interface
  interface getVolumePtr
    module procedure SepVM_d3_class_getVolumePtr
  end interface
  interface getCentroidPtr
    module procedure SepVM_d3_class_getCentroidPtr
  end interface

  interface

    subroutine F_SepVM_d3_new(this) &
      bind(C, name="c_SepVM_d3_new")
      import
      implicit none
      type(c_SepVM_d3) :: this
    end subroutine F_SepVM_d3_new

    subroutine F_SepVM_d3_delete(this) &
      bind(C, name="c_SepVM_d3_delete")
      import
      implicit none
      type(c_SepVM_d3) :: this
    end subroutine F_SepVM_d3_delete

    subroutine F_SepVM_d3_normalizeByVolume(this) &
      bind(C, name="c_SepVM_d3_normalizeByVolume")
      import
      implicit none
      type(c_SepVM_d3) :: this
    end subroutine F_SepVM_d3_normalizeByVolume

    subroutine F_SepVM_d3_multiplyByVolume(this) &
      bind(C, name="c_SepVM_d3_multiplyByVolume")
      import
      implicit none
      type(c_SepVM_d3) :: this
    end subroutine F_SepVM_d3_multiplyByVolume

    function F_SepVM_d3_getVolume(this, a_index) result(a_volume) &
      bind(C, name="c_SepVM_d3_getVolume")
      import
      implicit none
      type(c_SepVM_d3) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE) :: a_volume
    end function F_SepVM_d3_getVolume

    subroutine F_SepVM_d3_getCentroid(this, a_index, a_centroid) &
      bind(C, name="c_SepVM_d3_getCentroid")
      import
      implicit none
      type(c_SepVM_d3) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*) :: a_centroid ! dimension(1:3)
    end subroutine F_SepVM_d3_getCentroid

    subroutine F_SepVM_d3_getData(this, a_index, a_data) &
      bind(C, name="c_SepVM_d3_getData")
      import
      implicit none
      type(c_SepVM_d3) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*) :: a_data ! dimension(1:3)
    end subroutine F_SepVM_d3_getData

    function F_SepVM_d3_getVolumePtr(this, a_index) result(a_volume_ptr) &
      bind(C, name="c_SepVM_d3_getVolumePtr")
      import
      implicit none
      type(c_SepVM_d3) :: this
      integer(C_INT) :: a_index
      type(C_PTR) :: a_volume_ptr
    end function F_SepVM_d3_getVolumePtr

    function F_SepVM_d3_getCentroidPtr(this, a_index) result(a_centroid_ptr) &
      bind(C, name="c_SepVM_d3_getCentroidPtr")
      import
      implicit none
      type(c_SepVM_d3) :: this
      integer(C_INT) :: a_index
      type(C_PTR) :: a_centroid_ptr ! Pointer to double[3]
    end function F_SepVM_d3_getCentroidPtr

  end interface


  contains

    subroutine SepVM_d3_class_new(this)
      implicit none
      type(SepVM_d3_type), intent(out) :: this
      call F_SepVM_d3_new(this%c_object)
    end subroutine SepVM_d3_class_new

    impure elemental subroutine SepVM_d3_class_delete(this)
      implicit none
      type(SepVM_d3_type), intent(in) :: this
      call F_SepVM_d3_delete(this%c_object)
    end subroutine SepVM_d3_class_delete

    subroutine SepVM_d3_class_normalizeByVolume(this)
      implicit none
      type(SepVM_d3_type), intent(inout) :: this
      call F_SepVM_d3_normalizeByVolume(this%c_object)
    end subroutine SepVM_d3_class_normalizeByVolume

    subroutine SepVM_d3_class_multiplyByVolume(this)
      implicit none
      type(SepVM_d3_type), intent(inout) :: this
      call F_SepVM_d3_multiplyByVolume(this%c_object)
    end subroutine SepVM_d3_class_multiplyByVolume

    function SepVM_d3_class_getVolume(this, a_index) result(a_volume)
      implicit none
      type(SepVM_d3_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_index
      real(IRL_double) :: a_volume
      a_volume = F_SepVM_d3_getVolume(this%c_object, a_index)
      return
    end function SepVM_d3_class_getVolume

    function SepVM_d3_class_getCentroid(this, a_index) result(a_centroid)
      implicit none
      type(SepVM_d3_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_index
      real(IRL_double), dimension(3) :: a_centroid
      call F_SepVM_d3_getCentroid(this%c_object, a_index, a_centroid)
      return
    end function SepVM_d3_class_getCentroid

    function SepVM_d3_class_getData(this, a_index) result(a_data)
      implicit none
      type(SepVM_d3_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_index
      real(IRL_double), dimension(3) :: a_data
      call F_SepVM_d3_getData(this%c_object, a_index, a_data)
      return
    end function SepVM_d3_class_getData

    function SepVM_d3_class_getVolumePtr(this, a_index) result(a_volume_ptr)
      implicit none
      type(SepVM_d3_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_index
      real(IRL_double), pointer :: a_volume_ptr
      call c_f_pointer(F_SepVM_d3_getVolumePtr(this%c_object, a_index), a_volume_ptr)
      return
    end function SepVM_d3_class_getVolumePtr

    function SepVM_d3_class_getCentroidPtr(this, a_index) result(a_centroid_ptr)
      implicit none
      type(SepVM_d3_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_index
      real(IRL_double), pointer, dimension(:) :: a_centroid_ptr
      call c_f_pointer(F_SepVM_d3_getCentroidPtr(this%c_object, a_index), a_centroid_ptr, [3])
      return
    end function SepVM_d3_class_getCentroidPtr

end module f_SepVM_d3_class
