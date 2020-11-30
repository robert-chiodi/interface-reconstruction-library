!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_Tagged_AccumVM_Tagged_AccumVM_Vol_class.f90
!!
!! This file contains the Fortran interface for volume moments classes.

!> \brief A fortran type class that allows the creation of
!! IRL's AccumulatedVolumeMoments<AccumulatedVolumeMoments<Volume>>
!! class along with enabling some of its methods.
module f_TagAccVM2_Vol_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  use f_TagAccVM_Vol_class
  implicit none

  type, public, bind(C) :: c_TagAccVM2_Vol
    type(C_PTR), private :: object = C_NULL_PTR
    logical(C_BOOL), private :: is_owning = .false.
  end type c_TagAccVM2_Vol

  type, public :: TagAccVM2_Vol_type
    type(c_TagAccVM2_Vol) :: c_object
  contains
    final :: TagAccVM2_Vol_class_delete
  end type TagAccVM2_Vol_type

  interface new
    module procedure TagAccVM2_Vol_class_new
  end interface
  interface getAtIndex
    module procedure TagAccVM2_Vol_class_getAtIndex
  end interface
  interface getAtTag
    module procedure TagAccVM2_Vol_class_getAtTag
  end interface
  interface getSize
    module procedure TagAccVM2_Vol_class_getSize
  end interface
  interface getTagForIndex
    module procedure TagAccVM2_Vol_class_getTagForIndex
  end interface

  interface

    subroutine F_TagAccVM2_Vol_new(this) &
      bind(C, name="c_TagAccVM2_Vol_new")
      import
      implicit none
      type(c_TagAccVM2_Vol) :: this
    end subroutine F_TagAccVM2_Vol_new

    subroutine F_TagAccVM2_Vol_delete(this) &
      bind(C, name="c_TagAccVM2_Vol_delete")
      import
      implicit none
      type(c_TagAccVM2_Vol) :: this
    end subroutine F_TagAccVM2_Vol_delete

    subroutine F_TagAccVM2_Vol_getAtIndex(this, a_list_index, a_tagged_list) &
      bind(C, name="c_TagAccVM2_Vol_getAtIndex")
      import
      implicit none
      type(c_TagAccVM2_Vol) :: this
      integer(C_INT) :: a_list_index
      type(c_TagAccVM_Vol) :: a_tagged_list
    end subroutine F_TagAccVM2_Vol_getAtIndex

    subroutine F_TagAccVM2_Vol_getAtTag(this, a_tag, a_tagged_list) &
      bind(C, name="c_TagAccVM2_Vol_getAtTag")
      import
      implicit none
      type(c_TagAccVM2_Vol) :: this
      integer(C_INT) :: a_tag
      type(c_TagAccVM_Vol) :: a_tagged_list
    end subroutine F_TagAccVM2_Vol_getAtTag

    function F_TagAccVM2_Vol_getSize(this) result(a_size) &
      bind(C, name="c_TagAccVM2_Vol_getSize")
      import
      implicit none
      type(c_TagAccVM2_Vol) :: this
      integer(C_INT) :: a_size
    end function F_TagAccVM2_Vol_getSize

    function F_TagAccVM2_Vol_getTagForIndex(this, a_index) result(a_tag) &
      bind(C, name="c_TagAccVM2_Vol_getTagForIndex")
      import
      implicit none
      type(c_TagAccVM2_Vol) :: this
      integer(C_INT) :: a_index
      integer(C_INT) :: a_tag
    end function F_TagAccVM2_Vol_getTagForIndex


  end interface


  contains

    subroutine TagAccVM2_Vol_class_new(this)
      implicit none
      type(TagAccVM2_Vol_type), intent(inout) :: this
      call F_TagAccVM2_Vol_new(this%c_object)
    end subroutine TagAccVM2_Vol_class_new

    impure elemental subroutine TagAccVM2_Vol_class_delete(this)
      implicit none
      type(TagAccVM2_Vol_type), intent(in) :: this
      call F_TagAccVM2_Vol_delete(this%c_object)
    end subroutine TagAccVM2_Vol_class_delete

    subroutine TagAccVM2_Vol_class_getAtIndex(this, a_list_index, a_tagged_list)
      implicit none
      type(TagAccVM2_Vol_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_list_index
      type(TagAccVM_Vol_type), intent(inout) :: a_tagged_list
      call F_TagAccVM2_Vol_getAtIndex(this%c_object, a_list_index, a_tagged_list%c_object)
      return
    end subroutine TagAccVM2_Vol_class_getAtIndex

    subroutine TagAccVM2_Vol_class_getAtTag(this, a_tag, a_tagged_list)
      implicit none
      type(TagAccVM2_Vol_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_tag
      type(TagAccVM_Vol_type), intent(inout) :: a_tagged_list

      call F_TagAccVM2_Vol_getAtTag(this%c_object, a_tag, a_tagged_list%c_object)
      return
    end subroutine TagAccVM2_Vol_class_getAtTag

    function TagAccVM2_Vol_class_getSize(this) result(a_size)
      implicit none
      type(TagAccVM2_Vol_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_size
      a_size = F_TagAccVM2_Vol_getSize(this%c_object)
      return
    end function TagAccVM2_Vol_class_getSize

    function TagAccVM2_Vol_class_getTagForIndex(this, a_index) result(a_tag)
      implicit none
      type(TagAccVM2_Vol_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_index
      integer(IRL_UnsignedIndex_t) :: a_tag
      a_tag = F_TagAccVM2_Vol_getTagForIndex(this%c_object, a_index)
      return
    end function TagAccVM2_Vol_class_getTagForIndex

end module f_TagAccVM2_Vol_class
