!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_defined_types.f90
!!
!! This file sets precision values for the IRL fortran interface.
module f_DefinedTypes
  ! 4-Byte integer to represent indices.
  ! Can support only half of what C++ IRL can because Fortran does
  ! not use unsigned integers.
  use, intrinsic :: iso_fortran_env, only: IRL_UnsignedIndex_t => INT32
  ! 4-Byte integer used for signed indices.
  use, intrinsic :: iso_fortran_env, only: IRL_SignedIndex_t => INT32
  ! Integer to represent large offsets.
  use, intrinsic :: iso_fortran_env, only: IRL_LargeOffsetIndex_t => INT64
  ! Integer that is really a byte.
  use, intrinsic :: iso_fortran_env, only: IRL_Byte_t => INT8
  ! Precision of doubles used for IRL
  use, intrinsic :: iso_fortran_env, only: IRL_double => REAL64
  public

end module f_DefinedTypes
