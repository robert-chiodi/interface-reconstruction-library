!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_lviraneighborhood.f90
!!
!! This file contains functions reproducing
!! the functionality of the IRL class
!! LVIRANeighborhood<RectangularCuboid>. The purpose of this
!! is to allow building the stencil
!! through references to then be supplied
!! to obtain a PlanarSeparator using
!! the LVIRA method.

!> \brief A fortran type class to 
!! provide the functionality of 
!! LVIRANeighborhood.
module f_LVIRANeigh_RectCub_class
  use f_RectCub_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_LVIRANeigh_RectCub
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_LVIRANeigh_RectCub

  type, public :: LVIRANeigh_RectCub_type
    type(c_LVIRANeigh_RectCub) :: c_object
  contains
    final :: LVIRANeigh_RectCub_class_delete
  end type LVIRANeigh_RectCub_type

  interface new
    module procedure LVIRANeigh_RectCub_class_new
  end interface
  interface setSize
    module procedure LVIRANeigh_RectCub_class_setSize
  end interface
  interface setMember
    module procedure LVIRANeigh_RectCub_class_setMember
  end interface
  interface addMember
    module procedure LVIRANeigh_RectCub_class_addMember
  end interface
  interface emptyNeighborhood
    module procedure LVIRANeigh_RectCub_class_emptyNeighborhood
  end interface
  interface setCenterOfStencil
    module procedure LVIRANeigh_RectCub_class_setCenterOfStencil
  end interface

  interface

    subroutine F_LVIRANeigh_RectCub_new(this) &
      bind(C, name="c_LVIRANeigh_RectCub_new")
      import
      implicit none
      type(c_LVIRANeigh_RectCub) :: this
    end subroutine F_LVIRANeigh_RectCub_new

    subroutine F_LVIRANeigh_RectCub_delete(this) &
      bind(C, name="c_LVIRANeigh_RectCub_delete")
      import
      implicit none
      type(c_LVIRANeigh_RectCub) :: this
    end subroutine F_LVIRANeigh_RectCub_delete

    subroutine F_LVIRANeigh_RectCub_setSize(this, a_size) &
      bind(C, name="c_LVIRANeigh_RectCub_setSize")
      import
      implicit none
      type(c_LVIRANeigh_RectCub) :: this
      integer(C_INT) :: a_size
    end subroutine F_LVIRANeigh_RectCub_setSize

    subroutine F_LVIRANeigh_RectCub_setMember(this, a_index, a_rectangular_cuboid, &
        a_liquid_volume_fraction) &
      bind(C, name="c_LVIRANeigh_RectCub_setMember")
      import
      implicit none
      integer(C_INT), intent(in) :: a_index
      type(c_LVIRANeigh_RectCub) :: this
      type(c_RectCub) :: a_rectangular_cuboid ! Pointer to RectCub
      real(C_DOUBLE), intent(in) :: a_liquid_volume_fraction ! scalar
    end subroutine F_LVIRANeigh_RectCub_setMember

    subroutine F_LVIRANeigh_RectCub_addMember(this, a_rectangular_cuboid, &
        a_volume_fraction) &
      bind(C, name="c_LVIRANeigh_RectCub_addMember")
      import
      implicit none
      type(c_LVIRANeigh_RectCub) :: this
      type(c_RectCub) :: a_rectangular_cuboid ! Pointer to RectCub
      real(C_DOUBLE) :: a_volume_fraction ! Pointer to double with volume fraction
    end subroutine F_LVIRANeigh_RectCub_addMember

    subroutine F_LVIRANeigh_RectCub_emptyNeighborhood(this) &
      bind(C, name="c_LVIRANeigh_RectCub_emptyNeighborhood")
      import
      implicit none
      type(c_LVIRANeigh_RectCub) :: this
    end subroutine F_LVIRANeigh_RectCub_emptyNeighborhood

    subroutine F_LVIRANeigh_RectCub_setCenterOfStencil(this, a_center_cell_index) &
      bind(C, name="c_LVIRANeigh_RectCub_setCenterOfStencil")
      import
      implicit none
      type(c_LVIRANeigh_RectCub) :: this
      integer(C_INT) :: a_center_cell_index
    end subroutine F_LVIRANeigh_RectCub_setCenterOfStencil

  end interface


  contains

    subroutine LVIRANeigh_RectCub_class_new(this)
      implicit none
      type(LVIRANeigh_RectCub_type), intent(inout) :: this
      call F_LVIRANeigh_RectCub_new(this%c_object)
    end subroutine LVIRANeigh_RectCub_class_new

    impure elemental subroutine LVIRANeigh_RectCub_class_delete(this)
      implicit none
      type(LVIRANeigh_RectCub_type), intent(in) :: this
      call F_LVIRANeigh_RectCub_delete(this%c_object)
    end subroutine LVIRANeigh_RectCub_class_delete

    subroutine LVIRANeigh_RectCub_class_setSize(this, a_size)
      implicit none
      type(LVIRANeigh_RectCub_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_size
      call F_LVIRANeigh_RectCub_setSize(this%c_object,a_size)
    end subroutine LVIRANeigh_RectCub_class_setSize

    subroutine LVIRANeigh_RectCub_class_setMember(this, a_index, a_rectangular_cuboid, &
          a_liquid_volume_fraction)
      implicit none
      type(LVIRANeigh_RectCub_type), intent(in) :: this
      integer(IRL_SignedIndex_t), intent(in) :: a_index
      type(RectCub_type), intent(in) :: a_rectangular_cuboid
      real(IRL_double), intent(in) :: a_liquid_volume_fraction
      call F_LVIRANeigh_RectCub_setMember(this%c_object,a_index, a_rectangular_cuboid%c_object, &
          a_liquid_volume_fraction)
    end subroutine LVIRANeigh_RectCub_class_setMember

    subroutine LVIRANeigh_RectCub_class_addMember(this, a_rectangular_cuboid, &
          a_volume_fraction)
      implicit none
      type(LVIRANeigh_RectCub_type), intent(in) :: this
      type(RectCub_type), intent(in) :: a_rectangular_cuboid
      real(IRL_double), intent(in) :: a_volume_fraction
      call F_LVIRANeigh_RectCub_addMember(this%c_object, a_rectangular_cuboid%c_object, &
        a_volume_fraction)
    end subroutine LVIRANeigh_RectCub_class_addMember

    subroutine LVIRANeigh_RectCub_class_emptyNeighborhood(this)
      implicit none
      type(LVIRANeigh_RectCub_type), intent(in) :: this
      call F_LVIRANeigh_RectCub_emptyNeighborhood(this%c_object)
    end subroutine LVIRANeigh_RectCub_class_emptyNeighborhood

    subroutine LVIRANeigh_RectCub_class_setCenterOfStencil(this, a_center_cell_index)
      implicit none
      type(LVIRANeigh_RectCub_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_center_cell_index
      call F_LVIRANeigh_RectCub_setCenterOfStencil(this%c_object, a_center_cell_index)
    end subroutine LVIRANeigh_RectCub_class_setCenterOfStencil

end module f_LVIRANeigh_RectCub_class
