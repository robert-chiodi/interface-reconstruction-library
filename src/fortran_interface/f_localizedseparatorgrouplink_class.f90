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
module f_LocSepGroupLink_class
  use f_DefinedTypes
  use f_PlanarLoc_class
  use f_PlanarSepPathGroup_class
  use, intrinsic :: iso_c_binding
  implicit none

  type, public, bind(C) :: c_LocSepGroupLink
    type(C_PTR), private :: object = C_NULL_PTR
    logical(C_BOOL), private :: is_owning  = .false.
  end type c_LocSepGroupLink

  type, public :: LocSepGroupLink_type
    type(c_LocSepGroupLink) :: c_object
  contains
    final :: LocSepGroupLink_class_delete
  end type LocSepGroupLink_type

  interface new
    module procedure LocSepGroupLink_class_new
  end interface
  interface setId
    module procedure LocSepGroupLink_class_setId
  end interface
  interface getId
    module procedure LocSepGroupLink_class_getId
  end interface
  interface setEdgeConnectivity
    module procedure LocSepGroupLink_class_setEdgeConnectivity
  end interface
  interface setEdgeConnectivityNull
    module procedure LocSepGroupLink_class_setEdgeConnectivityNull
  end interface

  interface

    subroutine F_LocSepGroupLink_new(this, a_planar_localizer, a_planar_separator) &
      bind(C, name="c_LocSepGroupLink_new")
      import
      implicit none
      type(c_LocSepGroupLink) :: this
      type(c_PlanarLoc) :: a_planar_localizer ! Pointer to PlanarLoc
      type(c_PlanarSepPathGroup) :: a_planar_separator ! Pointer to PlanarSep
    end subroutine F_LocSepGroupLink_new

    subroutine F_LocSepGroupLink_delete(this) &
      bind(C, name="c_LocSepGroupLink_delete")
      import
      implicit none
      type(c_LocSepGroupLink) :: this
    end subroutine F_LocSepGroupLink_delete

    subroutine F_LocSepGroupLink_setId(this, a_id) &
      bind(C, name="c_LocSepGroupLink_setId")
      import
      implicit none
      type(c_LocSepGroupLink) :: this
      integer(C_INT), intent(in) :: a_id ! Scalar >= 0
    end subroutine F_LocSepGroupLink_setId

    function F_LocSepGroupLink_getId(this) result(a_id) &
      bind(C, name="c_LocSepGroupLink_getId")
      import
      implicit none
      type(c_LocSepGroupLink) :: this
      integer(C_INT) :: a_id ! Scalar >= 0
    end function F_LocSepGroupLink_getId

    subroutine F_LocSepGroupLink_setEdgeConnectivity(this, a_plane_index, a_ptr_to_neighbor) &
      bind(C, name="c_LocSepGroupLink_setEdgeConnectivity")
      import
      implicit none
      type(c_LocSepGroupLink) :: this
      integer(C_INT), intent(in) :: a_plane_index ! Index for plane of localizer
      type(c_LocSepGroupLink) :: a_ptr_to_neighbor ! Ptr to neighboring LocSepGroupLink on other side of plane.
    end subroutine F_LocSepGroupLink_setEdgeConnectivity

    subroutine F_LocSepGroupLink_setEdgeConnectivityNull(this, a_plane_index) &
      bind(C, name="c_LocSepGroupLink_setEdgeConnectivityNull")
      import
      implicit none
      type(c_LocSepGroupLink) :: this
      integer(C_INT), intent(in) :: a_plane_index ! Index for plane of localizer
    end subroutine F_LocSepGroupLink_setEdgeConnectivityNull

  end interface

  contains

    subroutine LocSepGroupLink_class_new(this, a_planar_localizer, a_planar_separator)
      implicit none
      type(LocSepGroupLink_type), intent(inout) :: this
      type(PlanarLoc_type), intent(in) :: a_planar_localizer
      type(PlanarSepPathGroup_type), intent(in) :: a_planar_separator
      call F_LocSepGroupLink_new(this%c_object, a_planar_localizer%c_object, a_planar_separator%c_object)
    end subroutine LocSepGroupLink_class_new

    impure elemental subroutine LocSepGroupLink_class_delete(this)
      implicit none
      type(LocSepGroupLink_type), intent(in) :: this
      call F_LocSepGroupLink_delete(this%c_object)
    end subroutine LocSepGroupLink_class_delete

    subroutine LocSepGroupLink_class_setId(this, a_id)
      implicit none
      type(LocSepGroupLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_id
      ! a_id must be >= 0
      call F_LocSepGroupLink_setId(this%c_object, a_id)
    end subroutine LocSepGroupLink_class_setId

    function LocSepGroupLink_class_getId(this) result(a_id)
      implicit none
      type(LocSepGroupLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t) :: a_id
      ! a_id must be >= 0
      a_id = F_LocSepGroupLink_getId(this%c_object)
    end function LocSepGroupLink_class_getId

    subroutine LocSepGroupLink_class_setEdgeConnectivity(this, a_plane_index, &
        a_neighboring_LocSepGroupLink)
      implicit none
      type(LocSepGroupLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_plane_index
      type(LocSepGroupLink_type) :: a_neighboring_LocSepGroupLink
      call F_LocSepGroupLink_setEdgeConnectivity(this%c_object, a_plane_index, a_neighboring_LocSepGroupLink%c_object)
    end subroutine LocSepGroupLink_class_setEdgeConnectivity

    subroutine LocSepGroupLink_class_setEdgeConnectivityNull(this, a_plane_index)
      implicit none
      type(LocSepGroupLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_plane_index
      call F_LocSepGroupLink_setEdgeConnectivityNull(this%c_object, a_plane_index)
    end subroutine LocSepGroupLink_class_setEdgeConnectivityNull

end module f_LocSepGroupLink_class
