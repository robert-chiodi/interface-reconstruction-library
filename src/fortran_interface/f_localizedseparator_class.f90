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
!! IRL's LocalizedSeparator class along with enabling
!! some of its methods.
module f_LocSep_class
  use f_DefinedTypes
  use f_ObjServer_LocSep_class
  use f_PlanarLoc_class
  use f_PlanarSep_class
  use, intrinsic :: iso_c_binding
  implicit none

  type, public, bind(C) :: c_LocSep
    type(C_PTR), private :: object = C_NULL_PTR
    logical(C_BOOL), private :: is_owning = .false.
  end type c_LocSep

  type, public :: LocSep_type
    type(c_LocSep) :: c_object
  contains
    final :: LocSep_class_delete
  end type LocSep_type

  interface new
    module procedure LocSep_class_new
    module procedure LocSep_class_newFromObjectAllocationServer
  end interface

  interface

    subroutine F_LocSep_new(this, a_planar_localizer, a_planar_separator) &
      bind(C, name="c_LocSep_new")
      import
      implicit none
      type(c_LocSep) :: this
      type(c_PlanarLoc) :: a_planar_localizer ! Pointer to PlanarLoc
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep
    end subroutine F_LocSep_new

    subroutine F_LocSep_newFromObjectAllocationServer(this, a_object_allocation_server, &
                                                            a_planar_localizer, a_planar_separator) &
      bind(C, name="c_LocSep_newFromObjectAllocationServer")
      import
      implicit none
      type(c_LocSep) :: this
      type(c_ObjServer_LocSep) :: a_object_allocation_server ! Pointer to ObjectAllocationServer<LocSep>
      type(c_PlanarLoc) :: a_planar_localizer ! Pointer to PlanarLoc
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep
    end subroutine F_LocSep_newFromObjectAllocationServer

    subroutine F_LocSep_delete(this) &
      bind(C, name="c_LocSep_delete")
      import
      implicit none
      type(c_LocSep) :: this
    end subroutine F_LocSep_delete

  end interface

  contains

    subroutine LocSep_class_new(this, a_planar_localizer, a_planar_separator)
      implicit none
      type(LocSep_type), intent(inout) :: this
      type(PlanarLoc_type), intent(in) :: a_planar_localizer
      type(PlanarSep_type), intent(in) :: a_planar_separator
      call F_LocSep_new(this%c_object, a_planar_localizer%c_object, a_planar_separator%c_object)
    end subroutine LocSep_class_new

    subroutine LocSep_class_newFromObjectAllocationServer(this, a_object_allocation_server, &
                                                                a_planar_localizer, a_planar_separator)
      implicit none
      type(LocSep_type), intent(inout) :: this
      type(ObjServer_LocSep_type), intent(in) :: a_object_allocation_server
      type(PlanarLoc_type), intent(in) :: a_planar_localizer
      type(PlanarSep_type), intent(in) :: a_planar_separator
      call F_LocSep_newFromObjectAllocationServer(this%c_object, &
          a_object_allocation_server%c_object, a_planar_localizer%c_object, a_planar_separator%c_object)
    end subroutine LocSep_class_newFromObjectAllocationServer

    impure elemental subroutine LocSep_class_delete(this)
      implicit none
      type(LocSep_type), intent(in) :: this
      call F_LocSep_delete(this%c_object)
    end subroutine LocSep_class_delete

end module f_LocSep_class
