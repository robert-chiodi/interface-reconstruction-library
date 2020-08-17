!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_Tagged_AccumListVM_VMAN_class.f90
!!
!! This file contains the Fortran interface for volume moments classes.

!> \brief A fortran type class that allows the creation of
!! IRL's TaggedAccumulatedListedVolumeMomentsM<VolumeMomentsAndNormal>
!! class along with enabling some of its methods.
module f_TagAccListVM_VMAN_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  use f_ListVM_VMAN_class
  implicit none

  type, public, bind(C) :: c_TagAccListVM_VMAN
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_TagAccListVM_VMAN

  type, public :: TagAccListVM_VMAN_type
    type(c_TagAccListVM_VMAN) :: c_object
  contains
    final :: TagAccListVM_VMAN_class_delete
  end type TagAccListVM_VMAN_type

  interface new
    module procedure TagAccListVM_VMAN_class_new
  end interface
  interface getListAtIndex
    module procedure TagAccListVM_VMAN_class_getListAtIndex
  end interface
  interface append
    module procedure TagAccListVM_VMAN_class_append
  end interface
  interface clear
    module procedure TagAccListVM_VMAN_class_clear
  end interface
  interface getSize
    module procedure TagAccListVM_VMAN_class_getSize
  end interface
  interface getTagForIndex
    module procedure TagAccListVM_VMAN_class_getTagForIndex
  end interface


  interface

    subroutine F_TagAccListVM_VMAN_new(this) &
      bind(C, name="c_TagAccListVM_VMAN_new")
      import
      implicit none
      type(c_TagAccListVM_VMAN) :: this
    end subroutine F_TagAccListVM_VMAN_new

    subroutine F_TagAccListVM_VMAN_delete(this) &
      bind(C, name="c_TagAccListVM_VMAN_delete")
      import
      implicit none
      type(c_TagAccListVM_VMAN) :: this
    end subroutine F_TagAccListVM_VMAN_delete

    subroutine F_TagAccListVM_VMAN_getListAtIndex(this, a_index, a_other_list) &
      bind(C, name="c_TagAccListVM_VMAN_getListAtIndex")
      import
      implicit none
      type(c_TagAccListVM_VMAN) :: this
      integer(C_INT), intent(in) :: a_index
      type(c_ListVM_VMAN) :: a_other_list
    end subroutine F_TagAccListVM_VMAN_getListAtIndex

    subroutine F_TagAccListVM_VMAN_append(this, a_other_list) &
      bind(C, name="c_TagAccListVM_VMAN_append")
      import
      implicit none
      type(c_TagAccListVM_VMAN) :: this
      type(c_TagAccListVM_VMAN) :: a_other_list
    end subroutine F_TagAccListVM_VMAN_append

    subroutine F_TagAccListVM_VMAN_clear(this) &
      bind(C, name="c_TagAccListVM_VMAN_clear")
      import
      implicit none
      type(c_TagAccListVM_VMAN) :: this
    end subroutine F_TagAccListVM_VMAN_clear

    function F_TagAccListVM_VMAN_getSize(this) result(a_size) &
      bind(C, name="c_TagAccListVM_VMAN_getSize")
      import
      implicit none
      type(c_TagAccListVM_VMAN) :: this
      integer(C_INT) :: a_size
    end function F_TagAccListVM_VMAN_getSize

    function F_TagAccListVM_VMAN_getTagForIndex(this, a_index) result(a_tag) &
      bind(C, name="c_TagAccListVM_VMAN_getTagForIndex")
      import
      implicit none
      type(c_TagAccListVM_VMAN) :: this
      integer(C_INT) :: a_index
      integer(C_INT) :: a_tag
    end function F_TagAccListVM_VMAN_getTagForIndex

  end interface

  contains

    subroutine TagAccListVM_VMAN_class_new(this)
      implicit none
      type(TagAccListVM_VMAN_type), intent(inout) :: this
      call F_TagAccListVM_VMAN_new(this%c_object)
    end subroutine TagAccListVM_VMAN_class_new

    impure elemental subroutine TagAccListVM_VMAN_class_delete(this)
      implicit none
      type(TagAccListVM_VMAN_type), intent(in) :: this
      call F_TagAccListVM_VMAN_delete(this%c_object)
    end subroutine TagAccListVM_VMAN_class_delete

    subroutine TagAccListVM_VMAN_class_getListAtIndex(this, a_index, a_other_list)
      implicit none
      type(TagAccListVM_VMAN_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      type(ListVM_VMAN_type), intent(inout) :: a_other_list
      call F_TagAccListVM_VMAN_getListAtIndex(this%c_object, a_index, a_other_list%c_object)
    end subroutine TagAccListVM_VMAN_class_getListAtIndex

    subroutine TagAccListVM_VMAN_class_append(this, a_other_list)
      implicit none
      type(TagAccListVM_VMAN_type), intent(inout) :: this
      type(TagAccListVM_VMAN_type), intent(in) :: a_other_list
      call F_TagAccListVM_VMAN_append(this%c_object, a_other_list%c_object)
    end subroutine TagAccListVM_VMAN_class_append

    subroutine TagAccListVM_VMAN_class_clear(this)
      implicit none
      type(TagAccListVM_VMAN_type), intent(inout) :: this
      call F_TagAccListVM_VMAN_clear(this%c_object)
    end subroutine TagAccListVM_VMAN_class_clear

    function TagAccListVM_VMAN_class_getSize(this) result(a_size)
      implicit none
      type(TagAccListVM_VMAN_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_size
      a_size = F_TagAccListVM_VMAN_getSize(this%c_object)
      return
    end function TagAccListVM_VMAN_class_getSize

    function TagAccListVM_VMAN_class_getTagForIndex(this, a_index) result(a_tag)
      implicit none
      type(TagAccListVM_VMAN_type), intent(inout) :: this
      integer(IRL_UnsignedIndex_t) :: a_index
      integer(IRL_UnsignedIndex_t) :: a_tag
      a_tag = F_TagAccListVM_VMAN_getTagForIndex(this%c_object, a_index)
      return
    end function TagAccListVM_VMAN_class_getTagForIndex

end module f_TagAccListVM_VMAN_class
