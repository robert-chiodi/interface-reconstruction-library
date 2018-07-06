!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_Tagged_AccumVM_VM_class.f90
!!
!! This file contains the Fortran interface for volume moments classes.

!> \brief A fortran type class that allows the creation of
!! IRL's AccumulatedVolumeMomentsM<VolumeMoments>
!! class along with enabling some of its methods.
module f_TagAccVM_VM_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_TagAccVM_VM
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_TagAccVM_VM

  type, public :: TagAccVM_VM_type
    type(c_TagAccVM_VM) :: c_object
  contains
    final :: TagAccVM_VM_class_delete
  end type TagAccVM_VM_type

  interface new
    module procedure TagAccVM_VM_class_new
  end interface
  interface normalizeByVolume
    module procedure TagAccVM_VM_class_normalizeByVolume
  end interface
  interface multiplyByVolume
    module procedure TagAccVM_VM_class_multiplyByVolume
  end interface
  interface getVolumeAtIndex
    module procedure TagAccVM_VM_class_getVolumeAtIndex
  end interface
  interface getCentroidAtIndex
    module procedure TagAccVM_VM_class_getCentroidAtIndex
  end interface
  interface getVolumePtrAtIndex
    module procedure TagAccVM_VM_class_getVolumePtrAtIndex
  end interface
  interface getCentroidPtrAtIndex
    module procedure TagAccVM_VM_class_getCentroidPtrAtIndex
  end interface
  interface getSize
    module procedure TagAccVM_VM_class_getSize
  end interface
  interface getTagForIndex
    module procedure TagAccVM_VM_class_getTagForIndex
  end interface

  interface

    subroutine F_TagAccVM_VM_new(this) &
      bind(C, name="c_TagAccVM_VM_new")
      import
      implicit none
      type(c_TagAccVM_VM) :: this
    end subroutine F_TagAccVM_VM_new

    subroutine F_TagAccVM_VM_delete(this) &
      bind(C, name="c_TagAccVM_VM_delete")
      import
      implicit none
      type(c_TagAccVM_VM) :: this
    end subroutine F_TagAccVM_VM_delete

    subroutine F_TagAccVM_VM_normalizeByVolume(this) &
      bind(C, name="c_TagAccVM_VM_normalizeByVolume")
      import
      implicit none
      type(c_TagAccVM_VM) :: this
    end subroutine F_TagAccVM_VM_normalizeByVolume

    subroutine F_TagAccVM_VM_multiplyByVolume(this) &
      bind(C, name="c_TagAccVM_VM_multiplyByVolume")
      import
      implicit none
      type(c_TagAccVM_VM) :: this
    end subroutine F_TagAccVM_VM_multiplyByVolume

    function F_TagAccVM_VM_getVolumeAtIndex(this, a_list_index) result(a_volume) &
      bind(C, name="c_TagAccVM_VM_getVolumeAtIndex")
      import
      implicit none
      type(c_TagAccVM_VM) :: this
      integer(C_INT) :: a_list_index
      real(C_DOUBLE) :: a_volume
    end function F_TagAccVM_VM_getVolumeAtIndex

    subroutine F_TagAccVM_VM_getCentroidAtIndex(this, a_list_index, a_centroid) &
      bind(C, name="c_TagAccVM_VM_getCentroidAtIndex")
      import
      implicit none
      type(c_TagAccVM_VM) :: this
      integer(C_INT) :: a_list_index
      real(C_DOUBLE), dimension(*) :: a_centroid ! dimension(1:3)
    end subroutine F_TagAccVM_VM_getCentroidAtIndex

    function F_TagAccVM_VM_getVolumePtrAtIndex(this, a_list_index) result(a_volume_ptr) &
      bind(C, name="c_TagAccVM_VM_getVolumePtrAtIndex")
      import
      implicit none
      type(c_TagAccVM_VM) :: this
      integer(C_INT) :: a_list_index
      type(C_PTR) :: a_volume_ptr
    end function F_TagAccVM_VM_getVolumePtrAtIndex

    function F_TagAccVM_VM_getCentroidPtrAtIndex(this, a_list_index) result(a_centroid_ptr) &
      bind(C, name="c_TagAccVM_VM_getCentroidPtrAtIndex")
      import
      implicit none
      type(c_TagAccVM_VM) :: this
      integer(C_INT) :: a_list_index
      type(C_PTR) :: a_centroid_ptr ! Pointer to double[3]
    end function F_TagAccVM_VM_getCentroidPtrAtIndex

    function F_TagAccVM_VM_getSize(this) result(a_size) &
      bind(C, name="c_TagAccVM_VM_getSize")
      import
      implicit none
      type(c_TagAccVM_VM) :: this
      integer(C_INT) :: a_size
    end function F_TagAccVM_VM_getSize

    function F_TagAccVM_VM_getTagForIndex(this, a_index) result(a_tag) &
      bind(C, name="c_TagAccVM_VM_getTagForIndex")
      import
      implicit none
      type(c_TagAccVM_VM) :: this
      integer(C_INT) :: a_index
      integer(C_INT) :: a_tag
    end function F_TagAccVM_VM_getTagForIndex

  end interface

  contains

    subroutine TagAccVM_VM_class_new(this)
      implicit none
      type(TagAccVM_VM_type), intent(inout) :: this
      call F_TagAccVM_VM_new(this%c_object)
    end subroutine TagAccVM_VM_class_new

    impure elemental subroutine TagAccVM_VM_class_delete(this)
      implicit none
      type(TagAccVM_VM_type), intent(in) :: this
      call F_TagAccVM_VM_delete(this%c_object)
    end subroutine TagAccVM_VM_class_delete

    subroutine TagAccVM_VM_class_normalizeByVolume(this)
      implicit none
      type(TagAccVM_VM_type), intent(inout) :: this
      call F_TagAccVM_VM_normalizeByVolume(this%c_object)
    end subroutine TagAccVM_VM_class_normalizeByVolume

    subroutine TagAccVM_VM_class_multiplyByVolume(this)
      implicit none
      type(TagAccVM_VM_type), intent(inout) :: this
      call F_TagAccVM_VM_multiplyByVolume(this%c_object)
    end subroutine TagAccVM_VM_class_multiplyByVolume

    function TagAccVM_VM_class_getVolumeAtIndex(this, a_list_index) result(a_volume)
      implicit none
      type(TagAccVM_VM_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_list_index
      real(IRL_double) :: a_volume
      a_volume = F_TagAccVM_VM_getVolumeAtIndex(this%c_object, a_list_index)
      return
    end function TagAccVM_VM_class_getVolumeAtIndex

    function TagAccVM_VM_class_getCentroidAtIndex(this, a_list_index) result(a_centroid)
      implicit none
      type(TagAccVM_VM_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_list_index
      real(IRL_double), dimension(3) :: a_centroid
      call F_TagAccVM_VM_getCentroidAtIndex(this%c_object, a_list_index, a_centroid)
      return
    end function TagAccVM_VM_class_getCentroidAtIndex

    function TagAccVM_VM_class_getVolumePtrAtIndex(this, a_list_index) result(a_volume_ptr)
      implicit none
      type(TagAccVM_VM_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_list_index
      real(IRL_double), pointer :: a_volume_ptr
      call c_f_pointer(F_TagAccVM_VM_getVolumePtrAtIndex(this%c_object, a_list_index), a_volume_ptr)
      return
    end function TagAccVM_VM_class_getVolumePtrAtIndex

    function TagAccVM_VM_class_getCentroidPtrAtIndex(this, a_list_index) result(a_centroid_ptr)
      implicit none
      type(TagAccVM_VM_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_list_index
      real(IRL_double), pointer, dimension(:) :: a_centroid_ptr
      call c_f_pointer(F_TagAccVM_VM_getCentroidPtrAtIndex(this%c_object, a_list_index), a_centroid_ptr, [3])
      return
    end function TagAccVM_VM_class_getCentroidPtrAtIndex

    function TagAccVM_VM_class_getSize(this) result(a_size)
      implicit none
      type(TagAccVM_VM_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_size
      a_size = F_TagAccVM_VM_getSize(this%c_object)
      return
    end function TagAccVM_VM_class_getSize

    function TagAccVM_VM_class_getTagForIndex(this, a_index) result(a_tag)
      implicit none
      type(TagAccVM_VM_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_index
      integer(IRL_UnsignedIndex_t) :: a_tag
      a_tag = F_TagAccVM_VM_getTagForIndex(this%c_object, a_index)
      return
    end function TagAccVM_VM_class_getTagForIndex

end module f_TagAccVM_VM_class
