!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_Paraboloidarator_class.f90
!!
!! This file allows use of the IRL Paraboloidarator
!! class through a fortran interface.

!> \brief A fortran type class that allows the creation of
!! IRL's Paraboloidarator class along with enabling
!! some of its methods.
module f_TriangulatedParaboloid_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  use f_ObjServer_TriangulatedParaboloid_class
  implicit none

  type, public, bind(C) :: c_TriangulatedParaboloid
    type(C_PTR), private :: object = C_NULL_PTR
    logical(C_BOOL), private :: is_owning  = .false.
  end type c_TriangulatedParaboloid

  type, public :: TriangulatedParaboloid_type
    type(c_TriangulatedParaboloid) :: c_object
  contains
    final :: TriangulatedParaboloid_class_delete
  end type TriangulatedParaboloid_type

  interface new
    module procedure TriangulatedParaboloid_class_new
    module procedure TriangulatedParaboloid_class_newFromObjectAllocationServer
  end interface
  interface getNumberOfTriangles
    module procedure TriangulatedParaboloid_class_getNumberOfTriangles
  end interface  
  interface getPt
    module procedure TriangulatedParaboloid_class_getPt
  end interface
  interface zeroTriangulatedParaboloid
    module procedure TriangulatedParaboloid_class_clear
  end interface  


  interface

    subroutine F_TriangulatedParaboloid_new(this) &
      bind(C, name="c_TriangulatedParaboloid_new")
      import
      implicit none
      type(c_TriangulatedParaboloid) :: this
    end subroutine F_TriangulatedParaboloid_new

    subroutine F_TriangulatedParaboloid_newFromObjectAllocationServer(this, a_object_allocation_server) &
      bind(C, name="c_TriangulatedParaboloid_newFromObjectAllocationServer")
      import
      implicit none
      type(c_TriangulatedParaboloid) :: this
      type(c_ObjServer_TriangulatedParaboloid) :: a_object_allocation_server ! Pointer to ObjectAllocationServer<TriangulatedParaboloid>
    end subroutine F_TriangulatedParaboloid_newFromObjectAllocationServer

    subroutine F_TriangulatedParaboloid_delete(this) &
      bind(C, name="c_TriangulatedParaboloid_delete")
      import
      implicit none
      type(c_TriangulatedParaboloid) :: this
    end subroutine F_TriangulatedParaboloid_delete

    function F_TriangulatedParaboloid_getNumberOfTriangles(this) result(a_number_of_tris) &
      bind(C, name="c_TriangulatedParaboloid_getNumberOfTriangles")
      import
      implicit none
      type(c_TriangulatedParaboloid) :: this
      integer(C_INT) :: a_number_of_tris
    end function F_TriangulatedParaboloid_getNumberOfTriangles

    subroutine F_TriangulatedParaboloid_getPt(this, a_index_tri, a_index_vert, a_pt) &
      bind(C, name="c_TriangulatedParaboloid_getPt")
      import
      implicit none
      type(c_TriangulatedParaboloid) :: this
      integer(C_INT) :: a_index_tri
      integer(C_INT) :: a_index_vert
      real(C_DOUBLE), dimension(*), intent(out) :: a_pt ! dimension(1:3)
    end subroutine F_TriangulatedParaboloid_getPt

    subroutine F_TriangulatedParaboloid_clear(this) &
      bind(C, name="c_TriangulatedParaboloid_clear")
      import
      implicit none
      type(c_TriangulatedParaboloid) :: this
    end subroutine F_TriangulatedParaboloid_clear

  end interface


  contains

    subroutine TriangulatedParaboloid_class_new(this)
      implicit none
      type(TriangulatedParaboloid_type), intent(inout) :: this
      call F_TriangulatedParaboloid_new(this%c_object)
    end subroutine TriangulatedParaboloid_class_new

    subroutine TriangulatedParaboloid_class_newFromObjectAllocationServer(this, a_object_allocation_server)
      implicit none
      type(TriangulatedParaboloid_type), intent(inout) :: this
      type(ObjServer_TriangulatedParaboloid_type), intent(in) :: a_object_allocation_server
      call F_TriangulatedParaboloid_newFromObjectAllocationServer(this%c_object, a_object_allocation_server%c_object)
    end subroutine TriangulatedParaboloid_class_newFromObjectAllocationServer

    impure elemental subroutine TriangulatedParaboloid_class_delete(this)
      implicit none
      type(TriangulatedParaboloid_type), intent(in) :: this
      call F_TriangulatedParaboloid_delete(this%c_object)
    end subroutine TriangulatedParaboloid_class_delete

    function TriangulatedParaboloid_class_getNumberOfTriangles(this) result(a_number_of_tris)
      implicit none
      type(TriangulatedParaboloid_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t) :: a_number_of_tris
      a_number_of_tris = F_TriangulatedParaboloid_getNumberOfTriangles(this%c_object)
      return
    end function TriangulatedParaboloid_class_getNumberOfTriangles

    function TriangulatedParaboloid_class_getPt(this, a_index_tri, a_index_vert) result(a_pt)
      implicit none
      type(TriangulatedParaboloid_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index_tri
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index_vert
      real(IRL_double), dimension(3) :: a_pt
      call F_TriangulatedParaboloid_getPt(this%c_object, a_index_tri, a_index_vert, a_pt)
      return
    end function TriangulatedParaboloid_class_getPt

    impure elemental subroutine TriangulatedParaboloid_class_clear(this)
      implicit none
      type(TriangulatedParaboloid_type), intent(in) :: this
      call F_TriangulatedParaboloid_clear(this%c_object)
    end subroutine TriangulatedParaboloid_class_clear


end module f_TriangulatedParaboloid_class
