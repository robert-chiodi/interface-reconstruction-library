!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_serializer.f90
!!
!! This file deals with serializing IRL
!! class objects into byte buffers. This is
!! usually done before parallel communication via
!! MPI using MPI_BYTE.

!> This module contains mappings to the
!! IRL C interface that deal with serializing
!! IRL class objects into an array of bytes and 
!! packing them into a byte buffer.
module f_Serializer
  use f_DefinedTypes
  use f_PlanarSep_class
  use f_ByteBuffer_class
  implicit none

  interface serializeAndPack
    ! Pack a PlanarSep object into a ByteBuffer
    module procedure serializeAndPack_PlanarSep_ByteBuffer
  end interface serializeAndPack

  interface unpackAndStore
    ! Unpack a ByteBuffer to store into a PlanarSep
    module procedure unpackAndStore_PlanarSep_ByteBuffer
  end interface unpackAndStore

  interface
    subroutine F_serializeAndPack_PlanarSep_ByteBuffer(a_separator, a_byte_buffer) &
    bind(C, name="c_serializeAndPack_PlanarSep_ByteBuffer")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_PlanarSep) :: a_separator ! Pointer to PlanarSep object
      type(c_ByteBuffer) :: a_byte_buffer ! Pointer to ByteBuffer object
    end subroutine F_serializeAndPack_PlanarSep_ByteBuffer
  end interface

  interface
    subroutine F_unpackAndStore_PlanarSep_ByteBuffer(a_separator, a_byte_buffer) &
    bind(C, name="c_unpackAndStore_PlanarSep_ByteBuffer")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_PlanarSep) :: a_separator ! Pointer to PlanarSep object
      type(c_ByteBuffer) :: a_byte_buffer ! Pointer to ByteBuffer object
    end subroutine F_unpackAndStore_PlanarSep_ByteBuffer
  end interface

contains

  subroutine serializeAndPack_PlanarSep_ByteBuffer(a_separator, a_byte_buffer)
    implicit none
      type(PlanarSep_type) :: a_separator
      type(ByteBuffer_type) :: a_byte_buffer

      call F_serializeAndPack_PlanarSep_ByteBuffer &
          (a_separator%c_object, a_byte_buffer%c_object)
  end subroutine serializeAndPack_PlanarSep_ByteBuffer

  subroutine unpackAndStore_PlanarSep_ByteBuffer(a_separator, a_byte_buffer)
    implicit none
      type(PlanarSep_type) :: a_separator
      type(ByteBuffer_type) :: a_byte_buffer

      call F_unpackAndStore_PlanarSep_ByteBuffer &
          (a_separator%c_object, a_byte_buffer%c_object)
  end subroutine unpackAndStore_PlanarSep_ByteBuffer


end module f_Serializer
