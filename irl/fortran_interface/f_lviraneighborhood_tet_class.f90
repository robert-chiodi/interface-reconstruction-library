!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_lviraneighborhood_hexahedron_class.f90
!!
!! This file contains functions reproducing
!! the functionality of the IRL class
!! LVIRANeighborhood<Tetahedron>. The purpose of this
!! is to allow building the stencil
!! through references to then be supplied
!! to obtain a PlanarSeparator using
!! the LVIRA method.

!> \brief A fortran type class to 
!! provide the functionality of 
!! LVIRANeighborhood.
module f_LVIRANeigh_Tet_class
  use f_Tet_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_LVIRANeigh_Tet
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_LVIRANeigh_Tet

  type, public :: LVIRANeigh_Tet_type
    type(c_LVIRANeigh_Tet) :: c_object
  contains
    final :: LVIRANeigh_Tet_class_delete
  end type LVIRANeigh_Tet_type

  interface new
    module procedure LVIRANeigh_Tet_class_new
  end interface
  interface setSize
    module procedure LVIRANeigh_Tet_class_setSize
  end interface
  interface setMember
    module procedure LVIRANeigh_Tet_class_setMember
  end interface
  interface addMember
    module procedure LVIRANeigh_Tet_class_addMember
  end interface
  interface emptyNeighborhood
    module procedure LVIRANeigh_Tet_class_emptyNeighborhood
  end interface
  interface setCenterOfStencil
    module procedure LVIRANeigh_Tet_class_setCenterOfStencil
  end interface

  interface

    subroutine F_LVIRANeigh_Tet_new(this) &
      bind(C, name="c_LVIRANeigh_Tet_new")
      import
      implicit none
      type(c_LVIRANeigh_Tet) :: this
    end subroutine F_LVIRANeigh_Tet_new

    subroutine F_LVIRANeigh_Tet_delete(this) &
      bind(C, name="c_LVIRANeigh_Tet_delete")
      import
      implicit none
      type(c_LVIRANeigh_Tet) :: this
    end subroutine F_LVIRANeigh_Tet_delete

    subroutine F_LVIRANeigh_Tet_setSize(this, a_size) &
      bind(C, name="c_LVIRANeigh_Tet_setSize")
      import
      implicit none
      type(c_LVIRANeigh_Tet) :: this
      integer(C_INT) :: a_size
    end subroutine F_LVIRANeigh_Tet_setSize

    subroutine F_LVIRANeigh_Tet_setMember(this, a_index, a_rectangular_cuboid, &
        a_liquid_volume_fraction) &
      bind(C, name="c_LVIRANeigh_Tet_setMember")
      import
      implicit none
      integer(C_INT), intent(in) :: a_index
      type(c_LVIRANeigh_Tet) :: this
      type(c_Tet) :: a_rectangular_cuboid ! Pointer to Tet
      real(C_DOUBLE), intent(in) :: a_liquid_volume_fraction ! scalar
    end subroutine F_LVIRANeigh_Tet_setMember

    subroutine F_LVIRANeigh_Tet_addMember(this, a_rectangular_cuboid, &
        a_volume_fraction) &
      bind(C, name="c_LVIRANeigh_Tet_addMember")
      import
      implicit none
      type(c_LVIRANeigh_Tet) :: this
      type(c_Tet) :: a_rectangular_cuboid ! Pointer to Tet
      real(C_DOUBLE) :: a_volume_fraction ! Pointer to double with volume fraction
    end subroutine F_LVIRANeigh_Tet_addMember

    subroutine F_LVIRANeigh_Tet_emptyNeighborhood(this) &
      bind(C, name="c_LVIRANeigh_Tet_emptyNeighborhood")
      import
      implicit none
      type(c_LVIRANeigh_Tet) :: this
    end subroutine F_LVIRANeigh_Tet_emptyNeighborhood

    subroutine F_LVIRANeigh_Tet_setCenterOfStencil(this, a_center_cell_index) &
      bind(C, name="c_LVIRANeigh_Tet_setCenterOfStencil")
      import
      implicit none
      type(c_LVIRANeigh_Tet) :: this
      integer(C_INT) :: a_center_cell_index
    end subroutine F_LVIRANeigh_Tet_setCenterOfStencil

  end interface


  contains

    subroutine LVIRANeigh_Tet_class_new(this)
      implicit none
      type(LVIRANeigh_Tet_type), intent(inout) :: this
      call F_LVIRANeigh_Tet_new(this%c_object)
    end subroutine LVIRANeigh_Tet_class_new

    impure elemental subroutine LVIRANeigh_Tet_class_delete(this)
      implicit none
      type(LVIRANeigh_Tet_type), intent(in) :: this
      call F_LVIRANeigh_Tet_delete(this%c_object)
    end subroutine LVIRANeigh_Tet_class_delete

    subroutine LVIRANeigh_Tet_class_setSize(this, a_size)
      implicit none
      type(LVIRANeigh_Tet_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_size
      call F_LVIRANeigh_Tet_setSize(this%c_object,a_size)
    end subroutine LVIRANeigh_Tet_class_setSize

    subroutine LVIRANeigh_Tet_class_setMember(this, a_index, a_tet, a_liquid_volume_fraction)
      implicit none
      type(LVIRANeigh_Tet_type), intent(in) :: this
      integer(IRL_SignedIndex_t), intent(in) :: a_index
      type(Tet_type), intent(in) :: a_tet
      real(IRL_double), intent(in) :: a_liquid_volume_fraction
      call F_LVIRANeigh_Tet_setMember(this%c_object,a_index, a_tet%c_object, &
          a_liquid_volume_fraction)
    end subroutine LVIRANeigh_Tet_class_setMember

    subroutine LVIRANeigh_Tet_class_addMember(this, a_tet, a_volume_fraction)
      implicit none
      type(LVIRANeigh_Tet_type), intent(in) :: this
      type(Tet_type), intent(in) :: a_tet
      real(IRL_double), intent(in) :: a_volume_fraction
      call F_LVIRANeigh_Tet_addMember(this%c_object, a_tet%c_object, &
        a_volume_fraction)
    end subroutine LVIRANeigh_Tet_class_addMember

    subroutine LVIRANeigh_Tet_class_emptyNeighborhood(this)
      implicit none
      type(LVIRANeigh_Tet_type), intent(in) :: this
      call F_LVIRANeigh_Tet_emptyNeighborhood(this%c_object)
    end subroutine LVIRANeigh_Tet_class_emptyNeighborhood

    subroutine LVIRANeigh_Tet_class_setCenterOfStencil(this, a_center_cell_index)
      implicit none
      type(LVIRANeigh_Tet_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_center_cell_index
      call F_LVIRANeigh_Tet_setCenterOfStencil(this%c_object, a_center_cell_index)
    end subroutine LVIRANeigh_Tet_class_setCenterOfStencil

end module f_LVIRANeigh_Tet_class
