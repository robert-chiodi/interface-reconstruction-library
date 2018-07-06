!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_Tagged_AccumVM_SepVol_class.f90
!!
!! This file contains the Fortran interface for volume moments classes.

!> \brief A fortran type class that allows the creation of
!! IRL's AccumulatedVolumeMoments<SeparatedMoments<Volume>>
!! class along with enabling some of its methods.
module f_TagAccVM_Vol_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_TagAccVM_Vol
    type(C_PTR), private :: object = C_NULL_PTR
    logical(C_BOOL), private :: is_owning = .false.
  end type c_TagAccVM_Vol

  type, public :: TagAccVM_Vol_type
    type(c_TagAccVM_Vol) :: c_object
  contains
    final :: TagAccVM_Vol_class_delete
  end type TagAccVM_Vol_type

  interface new
    module procedure TagAccVM_Vol_class_new
  end interface
  interface getVolumeAtIndex
    module procedure TagAccVM_Vol_class_getVolumeAtIndex
  end interface
  interface getVolumeAtTag
    module procedure TagAccVM_Vol_class_getVolumeAtTag
  end interface
  interface getVolumePtrAtIndex
    module procedure TagAccVM_Vol_class_getVolumePtrAtIndex
  end interface
  interface getSize
    module procedure TagAccVM_Vol_class_getSize
  end interface
  interface getTagForIndex
    module procedure TagAccVM_Vol_class_getTagForIndex
  end interface getTagForIndex

  interface

    subroutine F_TagAccVM_Vol_new(this) &
      bind(C, name="c_TagAccVM_Vol_new")
      import
      implicit none
      type(c_TagAccVM_Vol) :: this
    end subroutine F_TagAccVM_Vol_new

    subroutine F_TagAccVM_Vol_delete(this) &
      bind(C, name="c_TagAccVM_Vol_delete")
      import
      implicit none
      type(c_TagAccVM_Vol) :: this
    end subroutine F_TagAccVM_Vol_delete

    function F_TagAccVM_Vol_getVolumeAtIndex(this, a_list_index) result(a_volume) &
      bind(C, name="c_TagAccVM_Vol_getVolumeAtIndex")
      import
      implicit none
      type(c_TagAccVM_Vol) :: this
      integer(C_INT) :: a_list_index
      real(C_DOUBLE) :: a_volume
    end function F_TagAccVM_Vol_getVolumeAtIndex

    function F_TagAccVM_Vol_getVolumeAtTag(this, a_tag) result(a_volume) &
      bind(C, name="c_TagAccVM_Vol_getVolumeAtTag")
      import
      implicit none
      type(c_TagAccVM_Vol) :: this
      integer(C_INT) :: a_tag
      real(C_DOUBLE) :: a_volume
    end function F_TagAccVM_Vol_getVolumeAtTag

    function F_TagAccVM_Vol_getVolumePtrAtIndex(this, a_list_index) result(a_volume_ptr) &
      bind(C, name="c_TagAccVM_Vol_getVolumePtrAtIndex")
      import
      implicit none
      type(c_TagAccVM_Vol) :: this
      integer(C_INT) :: a_list_index
      type(C_PTR) :: a_volume_ptr
    end function F_TagAccVM_Vol_getVolumePtrAtIndex

    function F_TagAccVM_Vol_getSize(this) result(a_size) &
      bind(C, name="c_TagAccVM_Vol_getSize")
      import
      implicit none
      type(c_TagAccVM_Vol) :: this
      integer(C_INT) :: a_size
    end function F_TagAccVM_Vol_getSize

    function F_TagAccVM_Vol_getTagForIndex(this, a_index) result(a_tag) &
      bind(C, name="c_TagAccVM_Vol_getTagForIndex")
      import
      implicit none
      type(c_TagAccVM_Vol) :: this
      integer(C_INT) :: a_index
      integer(C_INT) :: a_tag
    end function F_TagAccVM_Vol_getTagForIndex

    subroutine F_TagAccVM_Vol_printCLoc(this) &
      bind(C, name="c_TagAccVM_Vol_printCLoc")
      import
      implicit none
      type(c_TagAccVM_Vol) :: this
    end subroutine F_TagAccVM_Vol_printCLoc

  end interface

  contains

    subroutine TagAccVM_Vol_class_new(this)
      implicit none
      type(TagAccVM_Vol_type), intent(inout) :: this
      call F_TagAccVM_Vol_new(this%c_object)
    end subroutine TagAccVM_Vol_class_new

    impure elemental subroutine TagAccVM_Vol_class_delete(this)
      implicit none
      type(TagAccVM_Vol_type), intent(in) :: this
      call F_TagAccVM_Vol_delete(this%c_object)
    end subroutine TagAccVM_Vol_class_delete

    function TagAccVM_Vol_class_getVolumeAtIndex(this, a_list_index) result(a_volume)
      implicit none
      type(TagAccVM_Vol_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_list_index
      real(IRL_double) :: a_volume
      a_volume = F_TagAccVM_Vol_getVolumeAtIndex(this%c_object, a_list_index)
      return
    end function TagAccVM_Vol_class_getVolumeAtIndex

    function TagAccVM_Vol_class_getVolumeAtTag(this, a_tag) result(a_volume)
      implicit none
      type(TagAccVM_Vol_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_tag
      real(IRL_double) :: a_volume
      a_volume = F_TagAccVM_Vol_getVolumeAtTag(this%c_object, a_tag)
      return
    end function TagAccVM_Vol_class_getVolumeAtTag

    function TagAccVM_Vol_class_getVolumePtrAtIndex(this, a_list_index) result(a_volume_ptr)
      implicit none
      type(TagAccVM_Vol_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_list_index
      real(IRL_double), pointer :: a_volume_ptr
      call c_f_pointer(F_TagAccVM_Vol_getVolumePtrAtIndex(this%c_object, a_list_index), a_volume_ptr)
      return
    end function TagAccVM_Vol_class_getVolumePtrAtIndex

    function TagAccVM_Vol_class_getSize(this) result(a_size)
      implicit none
      type(TagAccVM_Vol_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_size
      a_size = F_TagAccVM_Vol_getSize(this%c_object)
      return
    end function TagAccVM_Vol_class_getSize

    function TagAccVM_Vol_class_getTagForIndex(this, a_index) result(a_tag)
      implicit none
      type(TagAccVM_Vol_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_index
      integer(IRL_UnsignedIndex_t) :: a_tag
      a_tag = F_TagAccVM_Vol_getTagForIndex(this%c_object, a_index)
      return
    end function TagAccVM_Vol_class_getTagForIndex
    
end module f_TagAccVM_Vol_class
