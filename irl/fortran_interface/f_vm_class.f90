!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_VM_class.f90
!!
!! This file contains the Fortran interface for volume moments classes.

!> \brief A fortran type class that allows the creation of
!! IRL's VolumeMoments
!! class along with enabling some of its methods.
module f_VM_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_VM
    type(C_PTR), private :: object = C_NULL_PTR
    logical(C_BOOL), private :: is_owning = .false.
  end type c_VM

  type, public :: VM_type
    type(c_VM) :: c_object
  contains
    final :: VM_class_delete
  end type VM_type

  interface new
    module procedure VM_class_new
  end interface
  interface construct
    module procedure VM_class_construct
  end interface
  interface normalizeByVolume
    module procedure VM_class_normalizeByVolume
  end interface
  interface multiplyByVolume
    module procedure VM_class_multiplyByVolume
  end interface
  interface getVolume
    module procedure VM_class_getVolume
  end interface
  interface getCentroid
    module procedure VM_class_getCentroid
  end interface
  interface getVolumePtr
    module procedure VM_class_getVolumePtr
  end interface
  interface getCentroidPtr
    module procedure VM_class_getCentroidPtr
  end interface


  interface

    subroutine F_VM_new(this) &
      bind(C, name="c_VM_new")
      import
      implicit none
      type(c_VM) :: this
    end subroutine F_VM_new

    subroutine F_VM_delete(this) &
      bind(C, name="c_VM_delete")
      import
      implicit none
      type(c_VM) :: this
    end subroutine F_VM_delete

    subroutine F_VM_construct(this, a_moments_list) &
      bind(C, name="c_VM_construct")
      import
      implicit none
      type(c_VM) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_moments_list
    end subroutine F_VM_construct

    subroutine F_VM_normalizeByVolume(this) &
      bind(C, name="c_VM_normalizeByVolume")
      import
      implicit none
      type(c_VM) :: this
    end subroutine F_VM_normalizeByVolume

    subroutine F_VM_multiplyByVolume(this) &
      bind(C, name="c_VM_multiplyByVolume")
      import
      implicit none
      type(c_VM) :: this
    end subroutine F_VM_multiplyByVolume

    function F_VM_getVolume(this) result(a_volume) &
      bind(C, name="c_VM_getVolume")
      import
      implicit none
      type(c_VM) :: this
      real(C_DOUBLE) :: a_volume
    end function F_VM_getVolume

    subroutine F_VM_getCentroid(this, a_centroid) &
      bind(C, name="c_VM_getCentroid")
      import
      implicit none
      type(c_VM) :: this
      real(C_DOUBLE), dimension(*) :: a_centroid ! dimension(1:3)
    end subroutine F_VM_getCentroid

    function F_VM_getVolumePtr(this) result(a_volume_ptr) &
      bind(C, name="c_VM_getVolumePtr")
      import
      implicit none
      type(c_VM) :: this
      type(C_PTR) :: a_volume_ptr
    end function F_VM_getVolumePtr

    function F_VM_getCentroidPtr(this) result(a_centroid_ptr) &
      bind(C, name="c_VM_getCentroidPtr")
      import
      implicit none
      type(c_VM) :: this
      type(C_PTR) :: a_centroid_ptr ! Pointer to double[3]
    end function F_VM_getCentroidPtr

  end interface


  contains

    subroutine VM_class_new(this)
      implicit none
      type(VM_type), intent(inout) :: this
      call F_VM_new(this%c_object)
    end subroutine VM_class_new

    impure elemental subroutine VM_class_delete(this)
      implicit none
      type(VM_type), intent(in) :: this
      call F_VM_delete(this%c_object)
    end subroutine VM_class_delete

    subroutine VM_class_construct(this, a_moments_list)
      implicit none
      type(VM_type), intent(inout) :: this
      real(IRL_double), dimension(4), intent(in) :: a_moments_list
      call F_VM_construct(this%c_object, a_moments_list)
    end subroutine VM_class_construct

    subroutine VM_class_normalizeByVolume(this)
      implicit none
      type(VM_type), intent(inout) :: this
      call F_VM_normalizeByVolume(this%c_object)
    end subroutine VM_class_normalizeByVolume

    subroutine VM_class_multiplyByVolume(this)
      implicit none
      type(VM_type), intent(inout) :: this
      call F_VM_multiplyByVolume(this%c_object)
    end subroutine VM_class_multiplyByVolume

    function VM_class_getVolume(this) result(a_volume)
      implicit none
      type(VM_type), intent(inout) :: this
      real(IRL_double) :: a_volume
      a_volume = F_VM_getVolume(this%c_object)
      return
    end function VM_class_getVolume

    function VM_class_getCentroid(this) result(a_centroid)
      implicit none
      type(VM_type), intent(inout) :: this
      real(IRL_double), dimension(3) :: a_centroid
      call F_VM_getCentroid(this%c_object, a_centroid)
      return
    end function VM_class_getCentroid

    function VM_class_getVolumePtr(this) result(a_volume_ptr)
      implicit none
      type(VM_type), intent(inout) :: this
      real(IRL_double), pointer :: a_volume_ptr
      call c_f_pointer(F_VM_getVolumePtr(this%c_object), a_volume_ptr)
      return
    end function VM_class_getVolumePtr

    function VM_class_getCentroidPtr(this) result(a_centroid_ptr)
      implicit none
      type(VM_type), intent(inout) :: this
      real(IRL_double), pointer, dimension(:) :: a_centroid_ptr
      call c_f_pointer(F_VM_getCentroidPtr(this%c_object), a_centroid_ptr, [3])
      return
    end function VM_class_getCentroidPtr

end module f_VM_class
