!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_r2pneighborhood_rectangularcuboid_class.f90
!!
!! This file contains functions reproducing
!! the functionality of the IRL class
!! R2PNeighborhood_RectangularCuboid. The purpose of this
!! is to allow building the stencil
!! through references to then supply
!! to obtain a PlanarSeparator using
!! the R2P method.

!> \brief A fortran type class to 
!! provide the functionality of 
!! R2PNeighborhood_RectangularCuboid.
module f_R2PNeigh_RectCub_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  use f_RectCub_class
  use f_SepVM_class
  implicit none

  type, public, bind(C) :: c_R2PNeigh_RectCub
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_R2PNeigh_RectCub

  type, public :: R2PNeigh_RectCub_type
    type(c_R2PNeigh_RectCub) :: c_object
  contains
    final :: R2PNeigh_RectCub_class_delete
  end type R2PNeigh_RectCub_type

  interface new
    module procedure R2PNeigh_RectCub_class_new
  end interface
  interface setSize
    module procedure R2PNeigh_RectCub_class_setSize
  end interface
  interface setMember
    module procedure R2PNeigh_RectCub_class_setMember
  end interface
  interface addMember
    module procedure R2PNeigh_RectCub_class_addMember
  end interface
  interface emptyNeighborhood
    module procedure R2PNeigh_RectCub_class_emptyNeighborhood
  end interface
  interface setCenterOfStencil
    module procedure R2PNeigh_RectCub_class_setCenterOfStencil
  end interface
  interface setSurfaceArea
    module procedure R2PNeigh_RectCub_class_setSurfaceArea
  end interface

  interface

    subroutine F_R2PNeigh_RectCub_new(this) &
      bind(C, name="c_R2PNeigh_RectCub_new")
      import
      implicit none
      type(c_R2PNeigh_RectCub) :: this
    end subroutine F_R2PNeigh_RectCub_new

    subroutine F_R2PNeigh_RectCub_delete(this) &
      bind(C, name="c_R2PNeigh_RectCub_delete")
      import
      implicit none
      type(c_R2PNeigh_RectCub) :: this
    end subroutine F_R2PNeigh_RectCub_delete

    subroutine F_R2PNeigh_RectCub_setSize(this, a_size) &
      bind(C, name="c_R2PNeigh_RectCub_setSize")
      import
      implicit none
      type(c_R2PNeigh_RectCub) :: this
      integer(C_INT) :: a_size
    end subroutine F_R2PNeigh_RectCub_setSize

    subroutine F_R2PNeigh_RectCub_setMember(this, a_rectangular_cuboid, &
        a_separated_volume_moments, a_index) &
      bind(C, name="c_R2PNeigh_RectCub_setMember")
      import
      implicit none
      type(c_R2PNeigh_RectCub) :: this
      type(c_RectCub) :: a_rectangular_cuboid ! Pointer to RectCub
      type(c_SepVM) :: a_separated_volume_moments ! Pointer to SeparatedMoments<VolumeMoments>
      integer(C_INT), intent(in) :: a_index
    end subroutine F_R2PNeigh_RectCub_setMember

    subroutine F_R2PNeigh_RectCub_addMember(this, a_rectangular_cuboid, &
        a_separated_volume_moments) &
      bind(C, name="c_R2PNeigh_RectCub_addMember")
      import
      implicit none
      type(c_R2PNeigh_RectCub) :: this
      type(c_RectCub) :: a_rectangular_cuboid ! Pointer to RectCub
      type(c_SepVM) :: a_separated_volume_moments ! Pointer to SeparatedMoments<VolumeMoments>
    end subroutine F_R2PNeigh_RectCub_addMember

    subroutine F_R2PNeigh_RectCub_emptyNeighborhood(this) &
      bind(C, name="c_R2PNeigh_RectCub_emptyNeighborhood")
      import
      implicit none
      type(c_R2PNeigh_RectCub) :: this
    end subroutine F_R2PNeigh_RectCub_emptyNeighborhood

    subroutine F_R2PNeigh_RectCub_setCenterOfStencil(this, a_center_cell_index) &
      bind(C, name="c_R2PNeigh_RectCub_setCenterOfStencil")
      import
      implicit none
      type(c_R2PNeigh_RectCub) :: this
      integer(C_INT) :: a_center_cell_index
    end subroutine F_R2PNeigh_RectCub_setCenterOfStencil

    subroutine F_R2PNeigh_RectCub_setSurfaceArea(this, a_surface_area) &
      bind(C, name="c_R2PNeigh_RectCub_setSurfaceArea")
      import
      implicit none
      type(c_R2PNeigh_RectCub) :: this
      real(C_DOUBLE) :: a_surface_area
    end subroutine F_R2PNeigh_RectCub_setSurfaceArea

  end interface

  contains

    subroutine R2PNeigh_RectCub_class_new(this)
      implicit none
      type(R2PNeigh_RectCub_type), intent(inout) :: this
      call F_R2PNeigh_RectCub_new(this%c_object)
    end subroutine R2PNeigh_RectCub_class_new

    impure elemental subroutine R2PNeigh_RectCub_class_delete(this)
      implicit none
      type(R2PNeigh_RectCub_type), intent(in) :: this
      call F_R2PNeigh_RectCub_delete(this%c_object)
    end subroutine R2PNeigh_RectCub_class_delete

    subroutine R2PNeigh_RectCub_class_setSize(this, a_size)
      implicit none
      type(R2PNeigh_RectCub_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_size
      call F_R2PNeigh_RectCub_setSize(this%c_object,a_size)
    end subroutine R2PNeigh_RectCub_class_setSize

    subroutine R2PNeigh_RectCub_class_setMember(this, a_rectangular_cuboid, &
          a_separated_volume_moments, a_index)
      implicit none
      type(R2PNeigh_RectCub_type), intent(in) :: this
      type(RectCub_type), intent(in) :: a_rectangular_cuboid
      type(SepVM_type), intent(in) :: a_separated_volume_moments
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      call F_R2PNeigh_RectCub_setMember(this%c_object, a_rectangular_cuboid%c_object, &
        a_separated_volume_moments%c_object,a_index)
    end subroutine R2PNeigh_RectCub_class_setMember

    subroutine R2PNeigh_RectCub_class_addMember(this, a_rectangular_cuboid, &
          a_separated_volume_moments)
      implicit none
      type(R2PNeigh_RectCub_type), intent(in) :: this
      type(RectCub_type), intent(in) :: a_rectangular_cuboid
      type(SepVM_type), intent(in) :: a_separated_volume_moments
      call F_R2PNeigh_RectCub_addMember(this%c_object, a_rectangular_cuboid%c_object, &
        a_separated_volume_moments%c_object)
    end subroutine R2PNeigh_RectCub_class_addMember

    subroutine R2PNeigh_RectCub_class_emptyNeighborhood(this)
      implicit none
      type(R2PNeigh_RectCub_type), intent(in) :: this
      call F_R2PNeigh_RectCub_emptyNeighborhood(this%c_object)
    end subroutine R2PNeigh_RectCub_class_emptyNeighborhood

    subroutine R2PNeigh_RectCub_class_setCenterOfStencil(this, a_center_cell_index)
      implicit none
      type(R2PNeigh_RectCub_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_center_cell_index
      call F_R2PNeigh_RectCub_setCenterOfStencil(this%c_object, a_center_cell_index)
    end subroutine R2PNeigh_RectCub_class_setCenterOfStencil

    subroutine R2PNeigh_RectCub_class_setSurfaceArea(this, a_surface_area)
      implicit none
      type(R2PNeigh_RectCub_type), intent(in) :: this
      real(IRL_double), intent(in) :: a_surface_area
      call F_R2PNeigh_RectCub_setSurfaceArea(this%c_object, a_surface_area)
    end subroutine R2PNeigh_RectCub_class_setSurfaceArea


end module f_R2PNeigh_RectCub_class
