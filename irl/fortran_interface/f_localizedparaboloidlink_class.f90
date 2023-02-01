!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_LocalizedSeparatorLink_class.f90
!!
!! This file allows use of the IRL LocalizedSeparatorLink
!! class through a fortran interface.

!> \brief A fortran type class that allows the creation of
!! IRL's LocalizedSeparatorLink class along with enabling
!! some of its methods.
module f_LocParabLink_class
  use f_DefinedTypes
  use f_ObjServer_LocParabLink_class
  use f_PlanarLoc_class
  use f_Paraboloid_class
  use, intrinsic :: iso_c_binding
  implicit none

  type, public, bind(C) :: c_LocParabLink
    type(C_PTR), private :: object = C_NULL_PTR
    logical(C_BOOL), private :: is_owning = .false.
  end type c_LocParabLink

  type, public :: LocParabLink_type
    type(c_LocParabLink) :: c_object
  contains
    final :: LocParabLink_class_delete
  end type LocParabLink_type

  interface new
    module procedure LocParabLink_class_new
    module procedure LocParabLink_class_newFromObjectAllocationServer
  end interface
  interface setId
    module procedure LocParabLink_class_setId
  end interface
  interface getId
    module procedure LocParabLink_class_getId
  end interface
  interface setEdgeConnectivity
    module procedure LocParabLink_class_setEdgeConnectivity
  end interface
  interface setEdgeConnectivityNull
    module procedure LocParabLink_class_setEdgeConnectivityNull
  end interface
  interface printToScreen
    module procedure LocParabLink_class_printToScreen
  end interface

  interface

    subroutine F_LocParabLink_new(this, a_planar_localizer, a_paraboloid) &
      bind(C, name="c_LocParabLink_new")
      import
      implicit none
      type(c_LocParabLink) :: this
      type(c_PlanarLoc) :: a_planar_localizer ! Pointer to PlanarLoc
      type(c_Paraboloid) :: a_paraboloid ! Pointer to Paraboloid
    end subroutine F_LocParabLink_new

    subroutine F_LocParabLink_newFromObjectAllocationServer(this, a_object_allocation_server, &
                                                            a_planar_localizer, a_paraboloid) &
      bind(C, name="c_LocParabLink_newFromObjectAllocationServer")
      import
      implicit none
      type(c_LocParabLink) :: this
      type(c_ObjServer_LocParabLink) :: a_object_allocation_server ! Pointer to ObjectAllocationServer<LocParabLink>
      type(c_PlanarLoc) :: a_planar_localizer ! Pointer to PlanarLoc
      type(c_Paraboloid) :: a_paraboloid ! Pointer to Paraboloid
    end subroutine F_LocParabLink_newFromObjectAllocationServer

    subroutine F_LocParabLink_delete(this) &
      bind(C, name="c_LocParabLink_delete")
      import
      implicit none
      type(c_LocParabLink) :: this
    end subroutine F_LocParabLink_delete

    subroutine F_LocParabLink_setId(this, a_id) &
      bind(C, name="c_LocParabLink_setId")
      import
      implicit none
      type(c_LocParabLink) :: this
      integer(C_INT), intent(in) :: a_id ! Scalar >= 0
    end subroutine F_LocParabLink_setId

    function F_LocParabLink_getId(this) result(a_id) &
      bind(C, name="c_LocParabLink_getId")
      import
      implicit none
      type(c_LocParabLink) :: this
      integer(C_INT) :: a_id ! Scalar >= 0
    end function F_LocParabLink_getId

    subroutine F_LocParabLink_setEdgeConnectivity(this, a_plane_index, a_ptr_to_neighbor) &
      bind(C, name="c_LocParabLink_setEdgeConnectivity")
      import
      implicit none
      type(c_LocParabLink) :: this
      integer(C_INT), intent(in) :: a_plane_index ! Index for plane of localizer
      type(c_LocParabLink) :: a_ptr_to_neighbor ! Ptr to neighboring LocParabLink on other side of plane.
    end subroutine F_LocParabLink_setEdgeConnectivity

    subroutine F_LocParabLink_setEdgeConnectivityNull(this, a_plane_index) &
      bind(C, name="c_LocParabLink_setEdgeConnectivityNull")
      import
      implicit none
      type(c_LocParabLink) :: this
      integer(C_INT), intent(in) :: a_plane_index ! Index for plane of localizer
    end subroutine F_LocParabLink_setEdgeConnectivityNull
    
    subroutine F_LocParabLink_printToScreen(this) &
      bind(C, name="c_LocParabLink_printToScreen")
      import
      implicit none
      type(c_LocParabLink) :: this
    end subroutine F_LocParabLink_printToScreen

  end interface

  contains

    subroutine LocParabLink_class_new(this, a_planar_localizer, a_paraboloid)
      implicit none
      type(LocParabLink_type), intent(inout) :: this
      type(PlanarLoc_type), intent(in) :: a_planar_localizer
      type(Paraboloid_type), intent(in) :: a_paraboloid
      call F_LocParabLink_new(this%c_object, a_planar_localizer%c_object, a_paraboloid%c_object)
    end subroutine LocParabLink_class_new

    subroutine LocParabLink_class_newFromObjectAllocationServer(this, a_object_allocation_server, &
                                                                a_planar_localizer, a_paraboloid)
      implicit none
      type(LocParabLink_type), intent(inout) :: this
      type(ObjServer_LocParabLink_type), intent(in) :: a_object_allocation_server
      type(PlanarLoc_type), intent(in) :: a_planar_localizer
      type(Paraboloid_type), intent(in) :: a_paraboloid
      call F_LocParabLink_newFromObjectAllocationServer(this%c_object, &
          a_object_allocation_server%c_object, a_planar_localizer%c_object, a_paraboloid%c_object)
    end subroutine LocParabLink_class_newFromObjectAllocationServer

    impure elemental subroutine LocParabLink_class_delete(this)
      implicit none
      type(LocParabLink_type), intent(in) :: this
      call F_LocParabLink_delete(this%c_object)
    end subroutine LocParabLink_class_delete

    subroutine LocParabLink_class_setId(this, a_id)
      implicit none
      type(LocParabLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_id
      ! a_id must be >= 0
      call F_LocParabLink_setId(this%c_object, a_id)
    end subroutine LocParabLink_class_setId

    function LocParabLink_class_getId(this) result(a_id)
      implicit none
      type(LocParabLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t) :: a_id
      ! a_id must be >= 0
      a_id = F_LocParabLink_getId(this%c_object)
    end function LocParabLink_class_getId

    subroutine LocParabLink_class_setEdgeConnectivity(this, a_plane_index, &
        a_neighboring_LocParabLink)
      implicit none
      type(LocParabLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_plane_index
      type(LocParabLink_type) :: a_neighboring_LocParabLink
      call F_LocParabLink_setEdgeConnectivity(this%c_object, a_plane_index, a_neighboring_LocParabLink%c_object)
    end subroutine LocParabLink_class_setEdgeConnectivity

    subroutine LocParabLink_class_setEdgeConnectivityNull(this, a_plane_index)
      implicit none
      type(LocParabLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_plane_index
      call F_LocParabLink_setEdgeConnectivityNull(this%c_object, a_plane_index)
    end subroutine LocParabLink_class_setEdgeConnectivityNull
    
    subroutine LocParabLink_class_printToScreen(this)
      implicit none
      type(LocParabLink_type), intent(in) :: this
      call F_LocParabLink_printToScreen(this%c_object)
    end subroutine LocParabLink_class_printToScreen

end module f_LocParabLink_class
