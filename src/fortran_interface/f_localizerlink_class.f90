!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_LocalizerLink_class.f90
!!
!! This file allows use of the IRL LocalizerLink
!! class through a fortran interface.

!> \brief A fortran type class that allows the creation of
!! IRL's LocalizerLink class along with enabling
!! some of its methods.
module f_LocLink_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  use f_ObjServer_LocLink_class
  use f_PlanarLoc_class
  implicit none

  type, public, bind(C) :: c_LocLink
    type(C_PTR), private :: object = C_NULL_PTR
    logical(C_BOOL), private :: is_owning = .false.
  end type c_LocLink

  type, public :: LocLink_type
    type(c_LocLink) :: c_object
  contains
    final :: LocLink_class_delete
  end type LocLink_type

  interface new
    module procedure LocLink_class_new
    module procedure LocLink_class_newFromObjectAllocationServer
  end interface
  interface setId
    module procedure LocLink_class_setId
  end interface
  interface getId
    module procedure LocLink_class_getId
  end interface
  interface setEdgeConnectivity
    module procedure LocLink_class_setEdgeConnectivity
  end interface
  interface setEdgeConnectivityNull
    module procedure LocLink_class_setEdgeConnectivityNull
  end interface

  interface

    subroutine F_LocLink_new(this, a_planar_localizer) &
      bind(C, name="c_LocLink_new")
      import
      implicit none
      type(c_LocLink) :: this
      type(c_PlanarLoc) :: a_planar_localizer ! Pointer to PlanarLoc
    end subroutine F_LocLink_new

    subroutine F_LocLink_newFromObjectAllocationServer(this, a_object_allocation_server, &
                                                            a_planar_localizer) &
      bind(C, name="c_LocLink_newFromObjectAllocationServer")
      import
      implicit none
      type(c_LocLink) :: this
      type(c_ObjServer_LocLink) :: a_object_allocation_server ! Pointer to ObjectAllocationServer<LocLink>
      type(c_PlanarLoc) :: a_planar_localizer ! Pointer to PlanarLoc
    end subroutine F_LocLink_newFromObjectAllocationServer

    subroutine F_LocLink_delete(this) &
      bind(C, name="c_LocLink_delete")
      import
      implicit none
      type(c_LocLink) :: this
    end subroutine F_LocLink_delete

    subroutine F_LocLink_setId(this, a_id) &
      bind(C, name="c_LocLink_setId")
      import
      implicit none
      type(c_LocLink) :: this
      integer(C_INT), intent(in) :: a_id ! Scalar >= 0
    end subroutine F_LocLink_setId

    function F_LocLink_getId(this) result(a_id) &
      bind(C, name="c_LocLink_getId")
      import
      implicit none
      type(c_LocLink) :: this
      integer(C_INT) :: a_id ! Scalar >= 0
    end function F_LocLink_getId

    subroutine F_LocLink_setEdgeConnectivity(this, a_plane_index, a_ptr_to_neighbor) &
      bind(C, name="c_LocLink_setEdgeConnectivity")
      import
      implicit none
      type(c_LocLink) :: this
      integer(C_INT), intent(in) :: a_plane_index ! Index for plane of localizer
      type(c_LocLink) :: a_ptr_to_neighbor ! Ptr to neighboring LocLink on other side of plane.
    end subroutine F_LocLink_setEdgeConnectivity

    subroutine F_LocLink_setEdgeConnectivityNull(this, a_plane_index) &
      bind(C, name="c_LocLink_setEdgeConnectivityNull")
      import
      implicit none
      type(c_LocLink) :: this
      integer(C_INT), intent(in) :: a_plane_index ! Index for plane of localizer
    end subroutine F_LocLink_setEdgeConnectivityNull

  end interface


  contains

    subroutine LocLink_class_new(this, a_planar_localizer)
      implicit none
      type(LocLink_type), intent(inout) :: this
      type(PlanarLoc_type), intent(in) :: a_planar_localizer
      call F_LocLink_new(this%c_object, a_planar_localizer%c_object)
    end subroutine LocLink_class_new

    subroutine LocLink_class_newFromObjectAllocationServer(this, a_object_allocation_server, &
                                                                a_planar_localizer)
      implicit none
      type(LocLink_type), intent(inout) :: this
      type(ObjServer_LocLink_type), intent(in) :: a_object_allocation_server
      type(PlanarLoc_type), intent(in) :: a_planar_localizer
      call F_LocLink_newFromObjectAllocationServer(this%c_object, a_object_allocation_server%c_object, &
                                                         a_planar_localizer%c_object)
    end subroutine LocLink_class_newFromObjectAllocationServer

    impure elemental subroutine LocLink_class_delete(this)
      implicit none
      type(LocLink_type), intent(in) :: this
      call F_LocLink_delete(this%c_object)
    end subroutine LocLink_class_delete

    subroutine LocLink_class_setId(this, a_id)
      implicit none
      type(LocLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_id
      ! a_id must be >= 0
      call F_LocLink_setId(this%c_object, a_id)
    end subroutine LocLink_class_setId

    function LocLink_class_getId(this) result(a_id)
      implicit none
      type(LocLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t):: a_id
      ! a_id must be >= 0
      a_id = F_LocLink_getId(this%c_object)
    end function LocLink_class_getId

    subroutine LocLink_class_setEdgeConnectivity(this, a_plane_index, &
        a_neighboring_LocLink)
      implicit none
      type(LocLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_plane_index
      type(LocLink_type) :: a_neighboring_LocLink
      call F_LocLink_setEdgeConnectivity(this%c_object, a_plane_index, a_neighboring_LocLink%c_object)
    end subroutine LocLink_class_setEdgeConnectivity

    subroutine LocLink_class_setEdgeConnectivityNull(this, a_plane_index)
      implicit none
      type(LocLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_plane_index
      call F_LocLink_setEdgeConnectivityNull(this%c_object, a_plane_index)
    end subroutine LocLink_class_setEdgeConnectivityNull

end module f_LocLink_class
