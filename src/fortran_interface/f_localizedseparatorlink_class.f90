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
module f_LocSepLink_class
  use f_DefinedTypes
  use f_ObjServer_LocSepLink_class
  use f_PlanarLoc_class
  use f_PlanarSep_class
  use, intrinsic :: iso_c_binding
  implicit none

  type, public, bind(C) :: c_LocSepLink
    type(C_PTR), private :: object = C_NULL_PTR
    logical(C_BOOL), private :: is_owning = .false.
  end type c_LocSepLink

  type, public :: LocSepLink_type
    type(c_LocSepLink) :: c_object
  contains
    final :: LocSepLink_class_delete
  end type LocSepLink_type

  interface new
    module procedure LocSepLink_class_new
    module procedure LocSepLink_class_newFromObjectAllocationServer
  end interface
  interface setId
    module procedure LocSepLink_class_setId
  end interface
  interface getId
    module procedure LocSepLink_class_getId
  end interface
  interface setEdgeConnectivity
    module procedure LocSepLink_class_setEdgeConnectivity
  end interface
  interface setEdgeConnectivityNull
    module procedure LocSepLink_class_setEdgeConnectivityNull
  end interface
  interface printToScreen
    module procedure LocSepLink_class_printToScreen
  end interface

  interface

    subroutine F_LocSepLink_new(this, a_planar_localizer, a_planar_separator) &
      bind(C, name="c_LocSepLink_new")
      import
      implicit none
      type(c_LocSepLink) :: this
      type(c_PlanarLoc) :: a_planar_localizer ! Pointer to PlanarLoc
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep
    end subroutine F_LocSepLink_new

    subroutine F_LocSepLink_newFromObjectAllocationServer(this, a_object_allocation_server, &
                                                            a_planar_localizer, a_planar_separator) &
      bind(C, name="c_LocSepLink_newFromObjectAllocationServer")
      import
      implicit none
      type(c_LocSepLink) :: this
      type(c_ObjServer_LocSepLink) :: a_object_allocation_server ! Pointer to ObjectAllocationServer<LocSepLink>
      type(c_PlanarLoc) :: a_planar_localizer ! Pointer to PlanarLoc
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep
    end subroutine F_LocSepLink_newFromObjectAllocationServer

    subroutine F_LocSepLink_delete(this) &
      bind(C, name="c_LocSepLink_delete")
      import
      implicit none
      type(c_LocSepLink) :: this
    end subroutine F_LocSepLink_delete

    subroutine F_LocSepLink_setId(this, a_id) &
      bind(C, name="c_LocSepLink_setId")
      import
      implicit none
      type(c_LocSepLink) :: this
      integer(C_INT), intent(in) :: a_id ! Scalar >= 0
    end subroutine F_LocSepLink_setId

    function F_LocSepLink_getId(this) result(a_id) &
      bind(C, name="c_LocSepLink_getId")
      import
      implicit none
      type(c_LocSepLink) :: this
      integer(C_INT) :: a_id ! Scalar >= 0
    end function F_LocSepLink_getId

    subroutine F_LocSepLink_setEdgeConnectivity(this, a_plane_index, a_ptr_to_neighbor) &
      bind(C, name="c_LocSepLink_setEdgeConnectivity")
      import
      implicit none
      type(c_LocSepLink) :: this
      integer(C_INT), intent(in) :: a_plane_index ! Index for plane of localizer
      type(c_LocSepLink) :: a_ptr_to_neighbor ! Ptr to neighboring LocSepLink on other side of plane.
    end subroutine F_LocSepLink_setEdgeConnectivity

    subroutine F_LocSepLink_setEdgeConnectivityNull(this, a_plane_index) &
      bind(C, name="c_LocSepLink_setEdgeConnectivityNull")
      import
      implicit none
      type(c_LocSepLink) :: this
      integer(C_INT), intent(in) :: a_plane_index ! Index for plane of localizer
    end subroutine F_LocSepLink_setEdgeConnectivityNull
    
    subroutine F_LocSepLink_printToScreen(this) &
      bind(C, name="c_LocSepLink_printToScreen")
      import
      implicit none
      type(c_LocSepLink) :: this
    end subroutine F_LocSepLink_printToScreen

  end interface

  contains

    subroutine LocSepLink_class_new(this, a_planar_localizer, a_planar_separator)
      implicit none
      type(LocSepLink_type), intent(inout) :: this
      type(PlanarLoc_type), intent(in) :: a_planar_localizer
      type(PlanarSep_type), intent(in) :: a_planar_separator
      call F_LocSepLink_new(this%c_object, a_planar_localizer%c_object, a_planar_separator%c_object)
    end subroutine LocSepLink_class_new

    subroutine LocSepLink_class_newFromObjectAllocationServer(this, a_object_allocation_server, &
                                                                a_planar_localizer, a_planar_separator)
      implicit none
      type(LocSepLink_type), intent(inout) :: this
      type(ObjServer_LocSepLink_type), intent(in) :: a_object_allocation_server
      type(PlanarLoc_type), intent(in) :: a_planar_localizer
      type(PlanarSep_type), intent(in) :: a_planar_separator
      call F_LocSepLink_newFromObjectAllocationServer(this%c_object, &
          a_object_allocation_server%c_object, a_planar_localizer%c_object, a_planar_separator%c_object)
    end subroutine LocSepLink_class_newFromObjectAllocationServer

    impure elemental subroutine LocSepLink_class_delete(this)
      implicit none
      type(LocSepLink_type), intent(in) :: this
      call F_LocSepLink_delete(this%c_object)
    end subroutine LocSepLink_class_delete

    subroutine LocSepLink_class_setId(this, a_id)
      implicit none
      type(LocSepLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_id
      ! a_id must be >= 0
      call F_LocSepLink_setId(this%c_object, a_id)
    end subroutine LocSepLink_class_setId

    function LocSepLink_class_getId(this) result(a_id)
      implicit none
      type(LocSepLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t) :: a_id
      ! a_id must be >= 0
      a_id = F_LocSepLink_getId(this%c_object)
    end function LocSepLink_class_getId

    subroutine LocSepLink_class_setEdgeConnectivity(this, a_plane_index, &
        a_neighboring_LocSepLink)
      implicit none
      type(LocSepLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_plane_index
      type(LocSepLink_type) :: a_neighboring_LocSepLink
      call F_LocSepLink_setEdgeConnectivity(this%c_object, a_plane_index, a_neighboring_LocSepLink%c_object)
    end subroutine LocSepLink_class_setEdgeConnectivity

    subroutine LocSepLink_class_setEdgeConnectivityNull(this, a_plane_index)
      implicit none
      type(LocSepLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_plane_index
      call F_LocSepLink_setEdgeConnectivityNull(this%c_object, a_plane_index)
    end subroutine LocSepLink_class_setEdgeConnectivityNull
    
    subroutine LocSepLink_class_printToScreen(this)
      implicit none
      type(LocSepLink_type), intent(in) :: this
      call F_LocSepLink_printToScreen(this%c_object)
    end subroutine LocSepLink_class_printToScreen

end module f_LocSepLink_class
