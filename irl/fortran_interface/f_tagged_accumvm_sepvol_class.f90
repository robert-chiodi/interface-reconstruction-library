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
module f_TagAccVM_SepVol_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_TagAccVM_SepVol
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_TagAccVM_SepVol

  type, public :: TagAccVM_SepVol_type
    type(c_TagAccVM_SepVol) :: c_object
  contains
    final :: TagAccVM_SepVol_class_delete
  end type TagAccVM_SepVol_type

  interface new
    module procedure TagAccVM_SepVol_class_new
  end interface
  interface getVolumeAtIndex
    module procedure TagAccVM_SepVol_class_getVolumeAtIndex
  end interface
  interface getVolumeAtTag
    module procedure TagAccVM_SepVol_class_getVolumeAtTag
  end interface
  interface getVolumePtrAtIndex
    module procedure TagAccVM_SepVol_class_getVolumePtrAtIndex
  end interface
  interface getSize
    module procedure TagAccVM_SepVol_class_getSize
  end interface
  interface getTagForIndex
    module procedure TagAccVM_SepVol_class_getTagForIndex
  end interface

  interface

    subroutine F_TagAccVM_SepVol_new(this) &
      bind(C, name="c_TagAccVM_SepVol_new")
      import
      implicit none
      type(c_TagAccVM_SepVol) :: this
    end subroutine F_TagAccVM_SepVol_new

    subroutine F_TagAccVM_SepVol_delete(this) &
      bind(C, name="c_TagAccVM_SepVol_delete")
      import
      implicit none
      type(c_TagAccVM_SepVol) :: this
    end subroutine F_TagAccVM_SepVol_delete

    function F_TagAccVM_SepVol_getVolumeAtIndex(this, a_list_index, a_index) result(a_volume) &
      bind(C, name="c_TagAccVM_SepVol_getVolumeAtIndex")
      import
      implicit none
      type(c_TagAccVM_SepVol) :: this
      integer(C_INT) :: a_list_index
      integer(C_INT) :: a_index
      real(C_DOUBLE) :: a_volume
    end function F_TagAccVM_SepVol_getVolumeAtIndex

    function F_TagAccVM_SepVol_getVolumeAtTag(this, a_tag, a_index) result(a_volume) &
      bind(C, name="c_TagAccVM_SepVol_getVolumeAtTag")
      import
      implicit none
      type(c_TagAccVM_SepVol) :: this
      integer(C_INT) :: a_tag
      integer(C_INT) :: a_index
      real(C_DOUBLE) :: a_volume
    end function F_TagAccVM_SepVol_getVolumeAtTag

    function F_TagAccVM_SepVol_getVolumePtrAtIndex(this, a_list_index, a_index) result(a_volume_ptr) &
      bind(C, name="c_TagAccVM_SepVol_getVolumePtrAtIndex")
      import
      implicit none
      type(c_TagAccVM_SepVol) :: this
      integer(C_INT) :: a_list_index
      integer(C_INT) :: a_index
      type(C_PTR) :: a_volume_ptr
    end function F_TagAccVM_SepVol_getVolumePtrAtIndex

    function F_TagAccVM_SepVol_getSize(this) result(a_size) &
      bind(C, name="c_TagAccVM_SepVol_getSize")
      import
      implicit none
      type(c_TagAccVM_SepVol) :: this
      integer(C_INT) :: a_size
    end function F_TagAccVM_SepVol_getSize

    function F_TagAccVM_SepVol_getTagForIndex(this, a_index) result(a_tag) &
      bind(C, name="c_TagAccVM_SepVol_getTagForIndex")
      import
      implicit none
      type(c_TagAccVM_SepVol) :: this
      integer(C_INT) :: a_index
      integer(C_INT) :: a_tag
    end function F_TagAccVM_SepVol_getTagForIndex


  end interface


  contains

    subroutine TagAccVM_SepVol_class_new(this)
      implicit none
      type(TagAccVM_SepVol_type), intent(inout) :: this
      call F_TagAccVM_SepVol_new(this%c_object)
    end subroutine TagAccVM_SepVol_class_new

    impure elemental subroutine TagAccVM_SepVol_class_delete(this)
      implicit none
      type(TagAccVM_SepVol_type), intent(in) :: this
      call F_TagAccVM_SepVol_delete(this%c_object)
    end subroutine TagAccVM_SepVol_class_delete

    function TagAccVM_SepVol_class_getVolumeAtIndex(this, a_list_index, a_index) result(a_volume)
      implicit none
      type(TagAccVM_SepVol_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_list_index
      integer(IRL_UnsignedIndex_t) :: a_index
      real(IRL_double) :: a_volume
      a_volume = F_TagAccVM_SepVol_getVolumeAtIndex(this%c_object, a_list_index, a_index)
      return
    end function TagAccVM_SepVol_class_getVolumeAtIndex

    function TagAccVM_SepVol_class_getVolumeAtTag(this, a_tag, a_index) result(a_volume)
      implicit none
      type(TagAccVM_SepVol_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_tag
      integer(IRL_UnsignedIndex_t) :: a_index
      real(IRL_double) :: a_volume
      a_volume = F_TagAccVM_SepVol_getVolumeAtTag(this%c_object, a_tag, a_index)
      return
    end function TagAccVM_SepVol_class_getVolumeAtTag

    function TagAccVM_SepVol_class_getVolumePtrAtIndex(this, a_list_index, a_index) result(a_volume_ptr)
      implicit none
      type(TagAccVM_SepVol_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_list_index
      integer(IRL_UnsignedIndex_t) :: a_index
      real(IRL_double), pointer :: a_volume_ptr
      call c_f_pointer(F_TagAccVM_SepVol_getVolumePtrAtIndex(this%c_object, a_list_index, a_index), a_volume_ptr)
      return
    end function TagAccVM_SepVol_class_getVolumePtrAtIndex

    function TagAccVM_SepVol_class_getSize(this) result(a_size)
      implicit none
      type(TagAccVM_SepVol_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_size
      a_size = F_TagAccVM_SepVol_getSize(this%c_object)
      return
    end function TagAccVM_SepVol_class_getSize

    function TagAccVM_SepVol_class_getTagForIndex(this, a_index) result(a_tag)
      implicit none
      type(TagAccVM_SepVol_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_index
      integer(IRL_UnsignedIndex_t) :: a_tag
      a_tag = F_TagAccVM_SepVol_getTagForIndex(this%c_object, a_index)
      return
    end function TagAccVM_SepVol_class_getTagForIndex

end module f_TagAccVM_SepVol_class
