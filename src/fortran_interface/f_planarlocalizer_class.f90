!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_PlanarLocalizer_class.f90
!!
!! This file allows use of the IRL PlanarLocalizer
!! class through a fortran interface.

!> \brief A fortran type class that allows the creation of
!! IRL's PlanarLocalizer class along with enabling
!! some of its methods.
module f_PlanarLoc_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  use f_ObjServer_PlanarLoc_class
  implicit none

  type, public, bind(C) :: c_PlanarLoc
    type(C_PTR), private :: object = C_NULL_PTR
    logical(C_BOOL), private :: is_owning = .false.
  end type c_PlanarLoc

  type, public :: PlanarLoc_type
    type(c_PlanarLoc) :: c_object
  contains
    final :: PlanarLoc_class_delete
  end type PlanarLoc_type

  interface new
    module procedure PlanarLoc_class_new
    module procedure PlanarLoc_class_newFromObjectAllocationServer
  end interface
  interface addPlane
    module procedure PlanarLoc_class_addPlane
  end interface
  interface setNumberOfPlanes
    module procedure PlanarLoc_class_setNumberOfPlanes
  end interface
  interface setPlane
    module procedure PlanarLoc_class_setPlane
  end interface
  interface setFromRectangularCuboid
    module procedure PlanarLoc_class_setFromRectangularCuboid
  end interface
  interface printToScreen
    module procedure PlanarLoc_class_printToScreen
  end interface

  interface

    subroutine F_PlanarLoc_new(this) &
      bind(C, name="c_PlanarLoc_new")
      import
      implicit none
      type(c_PlanarLoc) :: this
    end subroutine F_PlanarLoc_new

    subroutine F_PlanarLoc_newFromObjectAllocationServer(this, a_object_allocation_server) &
      bind(C, name="c_PlanarLoc_newFromObjectAllocationServer")
      import
      implicit none
      type(c_PlanarLoc) :: this
      type(c_ObjServer_PlanarLoc) :: a_object_allocation_server ! Pointer to ObjectAllocationServer<PlanarLoc>
    end subroutine F_PlanarLoc_newFromObjectAllocationServer

    subroutine F_PlanarLoc_delete(this) &
      bind(C, name="c_PlanarLoc_delete")
      import
      implicit none
      type(c_PlanarLoc) :: this
    end subroutine F_PlanarLoc_delete

    subroutine F_PlanarLoc_addPlane(this, a_normal, a_distance) &
      bind(C, name="c_PlanarLoc_addPlane")
      import
      implicit none
      type(c_PlanarLoc) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_normal !  dimension(1:3)
      real(C_DOUBLE), intent(in) :: a_distance ! scalar
    end subroutine F_PlanarLoc_addPlane

    subroutine F_PlanarLoc_setNumberOfPlanes(this, a_number_to_set) &
      bind(C, name="c_PlanarLoc_setNumberOfPlanes")
      import
      implicit none
      type(c_PlanarLoc) :: this
      integer(C_INT), intent(in) :: a_number_to_set ! scalar
    end subroutine F_PlanarLoc_setNumberOfPlanes

    subroutine F_PlanarLoc_setPlane(this, a_plane_index_to_set,a_normal, a_distance) &
      bind(C, name="c_PlanarLoc_setPlane")
      import
      implicit none
      type(c_PlanarLoc) :: this
      integer(C_INT), intent(in) :: a_plane_index_to_set ! scalar
      real(C_DOUBLE), dimension(*), intent(in) :: a_normal !  dimension(1:3)
      real(C_DOUBLE), intent(in) :: a_distance ! scalar
    end subroutine F_PlanarLoc_setPlane

    subroutine F_PlanarLoc_setFromRectangularCuboid(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_PlanarLoc_setFromRectangularCuboid")
      import
      implicit none
      type(c_PlanarLoc) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_lower_pt !  dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(in) :: a_upper_pt !  dimension(1:3)
    end subroutine F_PlanarLoc_setFromRectangularCuboid

    subroutine F_PlanarLoc_printToScreen(this) &
      bind(C, name="c_PlanarLoc_printToScreen")
      import
      implicit none
      type(c_PlanarLoc) :: this
    end subroutine F_PlanarLoc_printToScreen

  end interface

  contains

    subroutine PlanarLoc_class_new(this)
      implicit none
      type(PlanarLoc_type), intent(inout) :: this
      call F_PlanarLoc_new(this%c_object)
    end subroutine PlanarLoc_class_new

    subroutine PlanarLoc_class_newFromObjectAllocationServer(this, a_object_allocation_server)
      implicit none
      type(PlanarLoc_type), intent(inout) :: this
      type(ObjServer_PlanarLoc_type), intent(in) :: a_object_allocation_server
      call F_PlanarLoc_newFromObjectAllocationServer(this%c_object, a_object_allocation_server%c_object)
    end subroutine PlanarLoc_class_newFromObjectAllocationServer

    impure elemental subroutine PlanarLoc_class_delete(this)
      implicit none
      type(PlanarLoc_type), intent(in) :: this
      call F_PlanarLoc_delete(this%c_object)
    end subroutine PlanarLoc_class_delete

    subroutine PlanarLoc_class_addPlane(this, a_normal, a_distance)
      implicit none
      type(PlanarLoc_type), intent(in) :: this
      real(IRL_double), dimension(1:3), intent(in) :: a_normal
      real(IRL_double), intent(in) :: a_distance
      call F_PlanarLoc_addPlane(this%c_object, a_normal, a_distance)
    end subroutine PlanarLoc_class_addPlane

    subroutine PlanarLoc_class_setNumberOfPlanes(this, a_number_to_set)
      implicit none
      type(PlanarLoc_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_number_to_set
      call F_PlanarLoc_setNumberOfPlanes(this%c_object, a_number_to_set)
    end subroutine PlanarLoc_class_setNumberOfPlanes

    subroutine PlanarLoc_class_setPlane(this, a_plane_index_to_set,a_normal, a_distance)
      implicit none
      type(PlanarLoc_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_plane_index_to_set
      real(IRL_double), dimension(1:3), intent(in) :: a_normal
      real(IRL_double), intent(in) :: a_distance
      call F_PlanarLoc_setPlane(this%c_object, a_plane_index_to_set, a_normal, a_distance)
    end subroutine PlanarLoc_class_setPlane

    subroutine PlanarLoc_class_setFromRectangularcuboid(this, a_lower_pt, a_upper_pt)
      implicit none
      type(PlanarLoc_type), intent(in) :: this
      real(IRL_double), dimension(1:3), intent(in) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(in) :: a_upper_pt
      call F_PlanarLoc_setFromRectangularCuboid(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine PlanarLoc_class_setFromRectangularcuboid

    subroutine PlanarLoc_class_printToScreen(this)
      implicit none
      type(PlanarLoc_type), intent(in) :: this
      call F_PlanarLoc_printToScreen(this%c_object)
    end subroutine PlanarLoc_class_printToScreen

end module f_PlanarLoc_class
