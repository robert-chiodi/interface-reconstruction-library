!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_elviraneighborhood.f90
!!
!! This file contains functions reproducing
!! the functionality of the IRL class
!! ELVIRANeighborhood. The purpose of this
!! is to allow building the stencil
!! through references to then supply
!! to obtain a PlanarSeparator using
!! the ELVIRA method.

!> \brief A fortran type class to 
!! provide the functionality of 
!! ELVIRANeighborhood.
module f_ELVIRANeigh_class
  use f_RectCub_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_ELVIRANeigh
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_ELVIRANeigh

  type, public :: ELVIRANeigh_type
    type(c_ELVIRANeigh) :: c_object
  contains
    final :: ELVIRANeigh_class_delete
  end type ELVIRANeigh_type

  interface new
    module procedure ELVIRANeigh_class_new
  end interface
  interface setSize
    module procedure ELVIRANeigh_class_setSize
  end interface
  interface setMember
    module procedure ELVIRANeigh_class_setMember
  end interface

  interface

    subroutine F_ELVIRANeigh_new(this) &
      bind(C, name="c_ELVIRANeigh_new")
      import
      implicit none
      type(c_ELVIRANeigh) :: this
    end subroutine F_ELVIRANeigh_new

    subroutine F_ELVIRANeigh_delete(this) &
      bind(C, name="c_ELVIRANeigh_delete")
      import
      implicit none
      type(c_ELVIRANeigh) :: this
    end subroutine F_ELVIRANeigh_delete

    subroutine F_ELVIRANeigh_setSize(this, a_size) &
      bind(C, name="c_ELVIRANeigh_setSize")
      import
      implicit none
      type(c_ELVIRANeigh) :: this
      integer(C_INT) :: a_size
    end subroutine F_ELVIRANeigh_setSize

    subroutine F_ELVIRANeigh_setMember(this, a_rectangular_cuboid, &
        a_liquid_volume_fraction, i, j, k) &
      bind(C, name="c_ELVIRANeigh_setMember")
      import
      implicit none
      type(c_ELVIRANeigh) :: this
      type(c_RectCub) :: a_rectangular_cuboid ! Pointer to RectCub
      real(C_DOUBLE), intent(in) :: a_liquid_volume_fraction ! scalar
      integer(C_INT), intent(in) :: i
      integer(C_INT), intent(in) :: j
      integer(C_INT), intent(in) :: k
    end subroutine F_ELVIRANeigh_setMember

  end interface


  contains

    subroutine ELVIRANeigh_class_new(this)
      implicit none
      type(ELVIRANeigh_type), intent(inout) :: this
      call F_ELVIRANeigh_new(this%c_object)
    end subroutine ELVIRANeigh_class_new

    impure elemental subroutine ELVIRANeigh_class_delete(this)
      implicit none
      type(ELVIRANeigh_type), intent(in) :: this
      call F_ELVIRANeigh_delete(this%c_object)
    end subroutine ELVIRANeigh_class_delete

    subroutine ELVIRANeigh_class_setSize(this, a_size)
      implicit none
      type(ELVIRANeigh_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_size
      call F_ELVIRANeigh_setSize(this%c_object,a_size)
    end subroutine ELVIRANeigh_class_setSize

    subroutine ELVIRANeigh_class_setMember(this, a_rectangular_cuboid, &
          a_liquid_volume_fraction, i, j, k)
      implicit none
      type(ELVIRANeigh_type), intent(in) :: this
      type(RectCub_type), intent(in) :: a_rectangular_cuboid
      real(IRL_double), intent(in) :: a_liquid_volume_fraction
      integer(IRL_SignedIndex_t), intent(in) :: i
      integer(IRL_SignedIndex_t), intent(in) :: j
      integer(IRL_SignedIndex_t), intent(in) :: k
      call F_ELVIRANeigh_setMember(this%c_object,a_rectangular_cuboid%c_object, &
          a_liquid_volume_fraction,i,j,k)
    end subroutine ELVIRANeigh_class_setMember

end module f_ELVIRANeigh_class
