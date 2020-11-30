!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_volumemoments_class.f90
!!
!! This file contains the Fortran interface for volume moments classes.

!> \brief A fortran type class that allows the creation of
!! IRL's ListedVolumeMomentsM<VolumeMomentsAndNormal>
!! class along with enabling some of its methods.
module f_ListVM_VMAN_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  use f_VMAN_class
  implicit none

  type, public, bind(C) :: c_ListVM_VMAN
    private
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_ListVM_VMAN

  type, public :: ListVM_VMAN_type
    type(c_ListVM_VMAN) :: c_object
  contains
    final :: ListVM_VMAN_class_delete
  end type ListVM_VMAN_type

  interface new
    module procedure ListVM_VMAN_class_new
  end interface
  interface append
    module procedure ListVM_VMAN_class_append
  end interface
  interface clear
    module procedure ListVM_VMAN_class_clear
  end interface
  interface getSize
    module procedure ListVM_VMAN_class_getSize
  end interface
  interface getMoments
    module procedure ListVM_VMAN_class_getMoments
  end interface
  interface zeroNormalComponent
    module procedure ListVM_VMAN_class_zeroNormalComponent
  end interface
  interface erase
    module procedure ListVM_VMAN_class_erase
  end interface

  interface

    subroutine F_ListVM_VMAN_new(this) &
      bind(C, name="c_ListVM_VMAN_new")
      import
      implicit none
      type(c_ListVM_VMAN) :: this
    end subroutine F_ListVM_VMAN_new

    subroutine F_ListVM_VMAN_delete(this) &
      bind(C, name="c_ListVM_VMAN_delete")
      import
      implicit none
      type(c_ListVM_VMAN) :: this
    end subroutine F_ListVM_VMAN_delete

    subroutine F_ListVM_VMAN_append(this, a_other_list) &
      bind(C, name="c_ListVM_VMAN_append")
      import
      implicit none
      type(c_ListVM_VMAN) :: this
      type(c_ListVM_VMAN) :: a_other_list
    end subroutine F_ListVM_VMAN_append

    subroutine F_ListVM_VMAN_clear(this) &
      bind(C, name="c_ListVM_VMAN_clear")
      import
      implicit none
      type(c_ListVM_VMAN) :: this
    end subroutine F_ListVM_VMAN_clear

    function F_ListVM_VMAN_getSize(this) result(a_size) &
      bind(C, name="c_ListVM_VMAN_getSize")
      import
      implicit none
      type(c_ListVM_VMAN) :: this
      integer(C_INT) :: a_size
    end function F_ListVM_VMAN_getSize

    subroutine F_ListVM_VMAN_getMoments(this, a_index, a_moments) &
      bind(C, name="c_ListVM_VMAN_getMoments")
      import
      implicit none
      type(c_ListVM_VMAN) :: this
      integer(C_INT) :: a_index
      type(c_VMAN) :: a_moments ! Pointer to VolumeMomentsAndNormal
    end subroutine F_ListVM_VMAN_getMoments

    subroutine F_ListVM_VMAN_zeroNormalComponent(this, a_index) &
      bind(C, name="c_ListVM_VMAN_zeroNormalComponent")
      import
      implicit none
      type(c_ListVM_VMAN) :: this
      integer(C_INT) :: a_index
    end subroutine F_ListVM_VMAN_zeroNormalComponent

    subroutine F_ListVM_VMAN_erase(this, a_index) &
      bind(C, name="c_ListVM_VMAN_erase")
      import
      implicit none
      type(c_ListVM_VMAN) :: this
      integer(C_INT) :: a_index
    end subroutine F_ListVM_VMAN_erase

  end interface


  contains

    subroutine ListVM_VMAN_class_new(this)
      implicit none
      type(ListVM_VMAN_type), intent(inout) :: this
      call F_ListVM_VMAN_new(this%c_object)
    end subroutine ListVM_VMAN_class_new

    impure elemental subroutine ListVM_VMAN_class_delete(this)
      implicit none
      type(ListVM_VMAN_type), intent(in) :: this
      call F_ListVM_VMAN_delete(this%c_object)
    end subroutine ListVM_VMAN_class_delete

    subroutine ListVM_VMAN_class_append(this, a_other_list)
      implicit none
      type(ListVM_VMAN_type) :: this
      type(ListVM_VMAN_type), intent(in) :: a_other_list
      call F_ListVM_VMAN_append(this%c_object, a_other_list%c_object)
    end subroutine ListVM_VMAN_class_append

    subroutine ListVM_VMAN_class_clear(this)
      implicit none
      type(ListVM_VMAN_type), intent(inout) :: this
      call F_ListVM_VMAN_clear(this%c_object)
    end subroutine ListVM_VMAN_class_clear

    function ListVM_VMAN_class_getSize(this) result(a_size)
      implicit none
      type(ListVM_VMAN_type) :: this
      integer(IRL_UnsignedIndex_t) :: a_size
      a_size = F_ListVM_VMAN_getSize(this%c_object)
      return
    end function ListVM_VMAN_class_getSize

    subroutine ListVM_VMAN_class_getMoments(this, a_index, a_moments)
      implicit none
      type(ListVM_VMAN_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      type(VMAN_type), intent(in) :: a_moments
      call F_ListVM_VMAN_getMoments(this%c_object, a_index, a_moments%c_object)
      return
    end subroutine ListVM_VMAN_class_getMoments

    subroutine ListVM_VMAN_class_zeroNormalComponent(this, a_index)
      implicit none
      type(ListVM_VMAN_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      call F_ListVM_VMAN_zeroNormalComponent(this%c_object, a_index)
      return
    end subroutine ListVM_VMAN_class_zeroNormalComponent

    subroutine ListVM_VMAN_class_erase(this, a_index)
      implicit none
      type(ListVM_VMAN_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      call F_ListVM_VMAN_erase(this%c_object, a_index)
      return
    end subroutine ListVM_VMAN_class_erase


end module f_ListVM_VMAN_class
