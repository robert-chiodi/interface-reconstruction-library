!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_Tagged_AccumVM_SepVM_class.f90
!!
!! This file contains the Fortran interface for volume moments classes.

!> \brief A fortran type class that allows the creation of
!! IRL's AccumulatedVolumeMomentsM<SeparatedMoments<VolumeMoments>>
!! class along with enabling some of its methods.
module f_TagAccVM_SepVM_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  use f_SepVM_class
  implicit none

  type, public, bind(C) :: c_TagAccVM_SepVM
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_TagAccVM_SepVM

  type, public :: TagAccVM_SepVM_type
    type(c_TagAccVM_SepVM) :: c_object
  contains
    final :: TagAccVM_SepVM_class_delete
  end type TagAccVM_SepVM_type

  interface new
    module procedure TagAccVM_SepVM_class_new
  end interface
  interface normalizeByVolume
    module procedure TagAccVM_SepVM_class_normalizeByVolume
  end interface
  interface multiplyByVolume
    module procedure TagAccVM_SepVM_class_multiplyByVolume
  end interface
  interface getSepVMAtIndex
    module procedure TagAccVM_SepVM_class_getSepVMAtIndex
  end interface
  interface getSepVMAtTag
    module procedure TagAccVM_SepVM_class_getSepVMAtTag
  end interface
  interface getSize
    module procedure TagAccVM_SepVM_class_getSize
  end interface
  interface getTagForIndex
    module procedure TagAccVM_SepVM_class_getTagForIndex
  end interface

  interface

    subroutine F_TagAccVM_SepVM_new(this) &
      bind(C, name="c_TagAccVM_SepVM_new")
      import
      implicit none
      type(c_TagAccVM_SepVM) :: this
    end subroutine F_TagAccVM_SepVM_new

    subroutine F_TagAccVM_SepVM_delete(this) &
      bind(C, name="c_TagAccVM_SepVM_delete")
      import
      implicit none
      type(c_TagAccVM_SepVM) :: this
    end subroutine F_TagAccVM_SepVM_delete

    subroutine F_TagAccVM_SepVM_normalizeByVolume(this) &
      bind(C, name="c_TagAccVM_SepVM_normalizeByVolume")
      import
      implicit none
      type(c_TagAccVM_SepVM) :: this
    end subroutine F_TagAccVM_SepVM_normalizeByVolume

    subroutine F_TagAccVM_SepVM_multiplyByVolume(this) &
      bind(C, name="c_TagAccVM_SepVM_multiplyByVolume")
      import
      implicit none
      type(c_TagAccVM_SepVM) :: this
    end subroutine F_TagAccVM_SepVM_multiplyByVolume

    subroutine F_TagAccVM_SepVM_getSepVMAtIndex(this, a_index, a_sepvm) &
      bind(C, name="c_TagAccVM_SepVM_getSepVMAtIndex")
      import
      implicit none
      type(c_TagAccVM_SepVM) :: this
      integer(C_INT) :: a_index
      type(c_SepVM) :: a_sepvm
    end subroutine F_TagAccVM_SepVM_getSepVMAtIndex

    subroutine F_TagAccVM_SepVM_getSepVMAtTag(this, a_tag, a_sepvm) &
      bind(C, name="c_TagAccVM_SepVM_getSepVMAtTag")
      import
      implicit none
      type(c_TagAccVM_SepVM) :: this
      integer(C_INT) :: a_tag
      type(c_SepVM) :: a_sepvm
    end subroutine F_TagAccVM_SepVM_getSepVMAtTag

    function F_TagAccVM_SepVM_getSize(this) result(a_size) &
      bind(C, name="c_TagAccVM_SepVM_getSize")
      import
      implicit none
      type(c_TagAccVM_SepVM) :: this
      integer(C_INT) :: a_size
    end function F_TagAccVM_SepVM_getSize

    function F_TagAccVM_SepVM_getTagForIndex(this, a_index) result(a_tag) &
      bind(C, name="c_TagAccVM_SepVM_getTagForIndex")
      import
      implicit none
      type(c_TagAccVM_SepVM) :: this
      integer(C_INT) :: a_index
      integer(C_INT) :: a_tag
    end function F_TagAccVM_SepVM_getTagForIndex


  end interface


  contains

    subroutine TagAccVM_SepVM_class_new(this)
      implicit none
      type(TagAccVM_SepVM_type), intent(inout) :: this
      call F_TagAccVM_SepVM_new(this%c_object)
    end subroutine TagAccVM_SepVM_class_new

    impure elemental subroutine TagAccVM_SepVM_class_delete(this)
      implicit none
      type(TagAccVM_SepVM_type), intent(in) :: this
      call F_TagAccVM_SepVM_delete(this%c_object)
    end subroutine TagAccVM_SepVM_class_delete

    subroutine TagAccVM_SepVM_class_normalizeByVolume(this)
      implicit none
      type(TagAccVM_SepVM_type), intent(inout) :: this
      call F_TagAccVM_SepVM_normalizeByVolume(this%c_object)
    end subroutine TagAccVM_SepVM_class_normalizeByVolume

    subroutine TagAccVM_SepVM_class_multiplyByVolume(this)
      implicit none
      type(TagAccVM_SepVM_type), intent(inout) :: this
      call F_TagAccVM_SepVM_multiplyByVolume(this%c_object)
    end subroutine TagAccVM_SepVM_class_multiplyByVolume

    subroutine TagAccVM_SepVM_class_getSepVMAtIndex(this, a_index, a_sepvm)
      implicit none
      type(TagAccVM_SepVM_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_index
      type(SepVM_type), intent(inout) :: a_sepvm
      call F_TagAccVM_SepVM_getSepVMAtIndex(this%c_object, a_index, a_sepvm%c_object)
      return
    end subroutine TagAccVM_SepVM_class_getSepVMAtIndex

    subroutine TagAccVM_SepVM_class_getSepVMAtTag(this, a_tag, a_sepvm)
      implicit none
      type(TagAccVM_SepVM_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_tag
      type(SepVM_type), intent(inout) :: a_sepvm
      call F_TagAccVM_SepVM_getSepVMAtTag(this%c_object, a_tag, a_sepvm%c_object)
      return
    end subroutine TagAccVM_SepVM_class_getSepVMAtTag

    function TagAccVM_SepVM_class_getSize(this) result(a_size)
      implicit none
      type(TagAccVM_SepVM_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_size
      a_size = F_TagAccVM_SepVM_getSize(this%c_object)
      return
    end function TagAccVM_SepVM_class_getSize

    function TagAccVM_SepVM_class_getTagForIndex(this, a_index) result(a_tag)
      implicit none
      type(TagAccVM_SepVM_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_index
      integer(IRL_UnsignedIndex_t) :: a_tag
      a_tag = F_TagAccVM_SepVM_getTagForIndex(this%c_object, a_index)
      return
    end function TagAccVM_SepVM_class_getTagForIndex

end module f_TagAccVM_SepVM_class
