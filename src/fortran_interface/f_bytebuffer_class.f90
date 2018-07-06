!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_ByteBuffer_class.f90
!!
!! This file contains the Fortran interface for the
!! ByteBuffer class.

!> \brief A fortran type class that allows the creation of
!! IRL's ByteBuffer class along with enabling
!! some of its methods.
module f_ByteBuffer_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_ByteBuffer
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_ByteBuffer

  type, public :: ByteBuffer_type
    type(c_ByteBuffer) :: c_object
  contains
    final :: ByteBuffer_class_delete
  end type ByteBuffer_type

  interface new
    module procedure ByteBuffer_class_new
  end interface

  interface getSize
    module procedure ByteBuffer_class_getSize
  end interface

  interface setSize
    module procedure ByteBuffer_class_setSize
  end interface

  interface resetBufferPointer
    module procedure ByteBuffer_class_resetBufferPointer
  end interface

  interface dataPtr
    module procedure ByteBuffer_class_dataPtr
  end interface

  interface

    subroutine F_ByteBuffer_new(this) &
      bind(C, name="c_ByteBuffer_new")
      import
      implicit none
      type(c_ByteBuffer) :: this
    end subroutine F_ByteBuffer_new

    subroutine F_ByteBuffer_delete(this) &
      bind(C, name="c_ByteBuffer_delete")
      import
      implicit none
      type(c_ByteBuffer) :: this
    end subroutine F_ByteBuffer_delete

    function F_ByteBuffer_getSize(this) result(a_size) &
      bind(C, name="c_ByteBuffer_getSize")
      import
      implicit none
      type(c_ByteBuffer) :: this
      integer(C_SIZE_T) :: a_size
    end function F_ByteBuffer_getSize

    subroutine F_ByteBuffer_setSize(this,a_size) &
      bind(C, name="c_ByteBuffer_setSize")
      import
      implicit none
      type(c_ByteBuffer) :: this
      integer(C_SIZE_T) :: a_size
    end subroutine F_ByteBuffer_setSize

    subroutine F_ByteBuffer_resetBufferPointer(this) &
      bind(C, name="c_ByteBuffer_resetBufferPointer")
      import
      implicit none
      type(c_ByteBuffer) :: this
    end subroutine F_ByteBuffer_resetBufferPointer

    function F_ByteBuffer_dataPtr(this) result(a_data_ptr) &
      bind(C, name="c_ByteBuffer_dataPtr")
      import
      implicit none
      type(c_ByteBuffer) :: this
      type(C_PTR) :: a_data_ptr
    end function F_ByteBuffer_dataPtr

  end interface

  contains

    impure elemental subroutine ByteBuffer_class_delete(this)
      implicit none
      type(ByteBuffer_type), intent(in) :: this
      call F_ByteBuffer_delete(this%c_object)
    end subroutine ByteBuffer_class_delete

    subroutine ByteBuffer_class_new(this)
      implicit none
      type(ByteBuffer_type), intent(inout) :: this
      call F_ByteBuffer_new(this%c_object)
    end subroutine ByteBuffer_class_new

    function ByteBuffer_class_getSize(this) result(a_size)
      implicit none
      type(ByteBuffer_type), intent(in) :: this
      integer(IRL_LargeOffsetIndex_t) :: a_size
      a_size = F_ByteBuffer_getSize(this%c_object)
      return
    end function ByteBuffer_class_getSize

    subroutine ByteBuffer_class_setsize(this, a_size)
      implicit none
      type(ByteBuffer_type), intent(inout) :: this
      integer(IRL_LargeOffsetIndex_t), intent(in) :: a_size
      call F_ByteBuffer_setSize(this%c_object, a_size)
    end subroutine ByteBuffer_class_setsize

    subroutine ByteBuffer_class_resetBufferPointer(this)
      implicit none
      type(ByteBuffer_type), intent(inout) :: this
      call F_ByteBuffer_resetBufferPointer(this%c_object)
    end subroutine ByteBuffer_class_resetBufferPointer

    function ByteBuffer_class_dataPtr(this) result(a_ptr_to_byte_array)
      use iso_c_binding
      implicit none
      type(ByteBuffer_type), intent(in) :: this
      integer(IRL_Byte_t), pointer, dimension(:) :: a_ptr_to_byte_array ! integer(IRL_Byte_t) is a byte
      call c_f_pointer(F_ByteBuffer_dataPtr(this%c_object), a_ptr_to_byte_array, [getSize(this)])
      return
    end function ByteBuffer_class_dataPtr

end module f_ByteBuffer_class
