!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_PlanarSeparator_class.f90
!!
!! This file allows use of the IRL PlanarSeparator
!! class through a fortran interface.

!> \brief A fortran type class that allows the creation of
!! IRL's PlanarSeparator class along with enabling
!! some of its methods.
module f_PlanarSep_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  use f_ObjServer_PlanarSep_class
  implicit none

  type, public, bind(C) :: c_PlanarSep
    type(C_PTR), private :: object = C_NULL_PTR
    logical(C_BOOL), private :: is_owning  = .false.
  end type c_PlanarSep

  type, public :: PlanarSep_type
    type(c_PlanarSep) :: c_object
  contains
    final :: PlanarSep_class_delete
  end type PlanarSep_type

  interface new
    module procedure PlanarSep_class_new
    module procedure PlanarSep_class_newFromObjectAllocationServer
  end interface
  interface addPlane
    module procedure PlanarSep_class_addPlane
  end interface
  interface setNumberOfPlanes
    module procedure PlanarSep_class_setNumberOfPlanes
  end interface
  interface setPlane
    module procedure PlanarSep_class_setPlane
  end interface
  interface copy
    module procedure PlanarSep_class_copy
  end interface
  interface getNumberOfPlanes
    module procedure PlanarSep_class_getNumberOfPlanes
  end interface
  interface getPlane
    module procedure PlanarSep_class_getPlane
  end interface
  interface isFlipped
    module procedure PlanarSep_class_isFlipped
  end interface
  interface printToScreen
    module procedure PlanarSep_class_printToScreen
  end interface


  interface

    subroutine F_PlanarSep_new(this) &
      bind(C, name="c_PlanarSep_new")
      import
      implicit none
      type(c_PlanarSep) :: this
    end subroutine F_PlanarSep_new

    subroutine F_PlanarSep_newFromObjectAllocationServer(this, a_object_allocation_server) &
      bind(C, name="c_PlanarSep_newFromObjectAllocationServer")
      import
      implicit none
      type(c_PlanarSep) :: this
      type(c_ObjServer_PlanarSep) :: a_object_allocation_server ! Pointer to ObjectAllocationServer<PlanarSep>
    end subroutine F_PlanarSep_newFromObjectAllocationServer

    subroutine F_PlanarSep_delete(this) &
      bind(C, name="c_PlanarSep_delete")
      import
      implicit none
      type(c_PlanarSep) :: this
    end subroutine F_PlanarSep_delete

    subroutine F_PlanarSep_addPlane(this, a_normal, a_distance) &
      bind(C, name="c_PlanarSep_addPlane")
      import
      implicit none
      type(c_PlanarSep) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_normal !  dimension(1:3)
      real(C_DOUBLE), intent(in) :: a_distance ! scalar
    end subroutine F_PlanarSep_addPlane

    subroutine F_PlanarSep_setNumberOfPlanes(this, a_number_to_set) &
      bind(C, name="c_PlanarSep_setNumberOfPlanes")
      import
      implicit none
      type(c_PlanarSep) :: this
      integer(C_INT), intent(in) :: a_number_to_set ! scalar
    end subroutine F_PlanarSep_setNumberOfPlanes

    subroutine F_PlanarSep_setPlane(this, a_plane_index_to_set,a_normal, a_distance) &
      bind(C, name="c_PlanarSep_setPlane")
      import
      implicit none
      type(c_PlanarSep) :: this
      integer(C_INT), intent(in) :: a_plane_index_to_set ! scalar
      real(C_DOUBLE), dimension(*), intent(in) :: a_normal !  dimension(1:3)
      real(C_DOUBLE), intent(in) :: a_distance ! scalar
    end subroutine F_PlanarSep_setPlane

    subroutine F_PlanarSep_copy(this, a_other_PlanarSep) &
      bind(C, name="c_PlanarSep_copy")
      import
      implicit none
      type(c_PlanarSep) :: this
      type(c_PlanarSep) :: a_other_PlanarSep
    end subroutine F_PlanarSep_copy

    function F_PlanarSep_getNumberOfPlanes(this) result(a_number_of_planes) &
      bind(C, name="c_PlanarSep_getNumberOfPlanes")
      import
      implicit none
      type(c_PlanarSep) :: this
      integer(C_INT) :: a_number_of_planes
    end function F_PlanarSep_getNumberOfPlanes

    subroutine F_PlanarSep_getPlane(this, a_index, a_plane_listed) &
      bind(C, name="c_PlanarSep_getPlane")
      import
      implicit none
      type(c_PlanarSep) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(out) :: a_plane_listed
    end subroutine F_PlanarSep_getPlane

    function F_PlanarSep_isFlipped(this) result(a_flipped) &
      bind(C, name="c_PlanarSep_isFlipped")
      import
      implicit none
      type(c_PlanarSep) :: this
      logical(C_BOOL) :: a_flipped
    end function F_PlanarSep_isFlipped

    subroutine F_PlanarSep_printToScreen(this) &
      bind(C, name="c_PlanarSep_printToScreen")
      import
      implicit none
      type(c_PlanarSep) :: this
    end subroutine F_PlanarSep_printToScreen

  end interface


  contains

    subroutine PlanarSep_class_new(this)
      implicit none
      type(PlanarSep_type), intent(inout) :: this
      call F_PlanarSep_new(this%c_object)
    end subroutine PlanarSep_class_new

    subroutine PlanarSep_class_newFromObjectAllocationServer(this, a_object_allocation_server)
      implicit none
      type(PlanarSep_type), intent(inout) :: this
      type(ObjServer_PlanarSep_type), intent(in) :: a_object_allocation_server
      call F_PlanarSep_newFromObjectAllocationServer(this%c_object, a_object_allocation_server%c_object)
    end subroutine PlanarSep_class_newFromObjectAllocationServer

    impure elemental subroutine PlanarSep_class_delete(this)
      implicit none
      type(PlanarSep_type), intent(in) :: this
      call F_PlanarSep_delete(this%c_object)
    end subroutine PlanarSep_class_delete

    subroutine PlanarSep_class_addPlane(this, a_normal, a_distance)
      implicit none
      type(PlanarSep_type), intent(in) :: this
      real(IRL_double), dimension(1:3), intent(in) :: a_normal
      real(IRL_double), intent(in) :: a_distance
      call F_PlanarSep_addPlane(this%c_object, a_normal, a_distance)
    end subroutine PlanarSep_class_addPlane

    subroutine PlanarSep_class_setNumberOfPlanes(this, a_number_to_set)
      implicit none
      type(PlanarSep_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_number_to_set
      call F_PlanarSep_setNumberOfPlanes(this%c_object, a_number_to_set)
    end subroutine PlanarSep_class_setNumberOfPlanes

    subroutine PlanarSep_class_setPlane(this, a_plane_index_to_set,a_normal, a_distance)
      implicit none
      type(PlanarSep_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_plane_index_to_set
      real(IRL_double), dimension(1:3), intent(in) :: a_normal
      real(IRL_double), intent(in) :: a_distance
      call F_PlanarSep_setPlane(this%c_object, a_plane_index_to_set, a_normal, a_distance)
    end subroutine PlanarSep_class_setPlane

    subroutine PlanarSep_class_copy(this, a_other_PlanarSep)
      implicit none
      type(PlanarSep_type), intent(inout) :: this
      type(PlanarSep_type), intent(in) :: a_other_PlanarSep
      call F_PlanarSep_copy(this%c_object, a_other_PlanarSep%c_object)
    end subroutine PlanarSep_class_copy

    function PlanarSep_class_getNumberOfPlanes(this) result(a_number_of_planes)
      implicit none
      type(PlanarSep_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t) :: a_number_of_planes
      a_number_of_planes = F_PlanarSep_getNumberOfPlanes(this%c_object)
    end function PlanarSep_class_getNumberOfPlanes

    function PlanarSep_class_getPlane(this, a_index) result(a_plane_listed)
      implicit none
      type(PlanarSep_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(4) :: a_plane_listed
      call F_PlanarSep_getPlane(this%c_object, a_index, a_plane_listed)
    end function PlanarSep_class_getPlane

    function PlanarSep_class_isFlipped(this) result(a_flipped)
      implicit none
      type(PlanarSep_type), intent(in) :: this
      logical(1) :: a_flipped
      a_flipped = F_PlanarSep_isFlipped(this%c_object)
      return
    end function PlanarSep_class_isFlipped

    subroutine PlanarSep_class_printToScreen(this)
      implicit none
      type(PlanarSep_type), intent(in) :: this
      call F_PlanarSep_printToScreen(this%c_object)
    end subroutine PlanarSep_class_printToScreen


end module f_PlanarSep_class
