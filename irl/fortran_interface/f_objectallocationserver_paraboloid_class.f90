!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_objectallocationserver_Paraboloidarator_class.f90
!!
!! This file allows use of the IRL ObjectAllocationServer<Paraboloidarator>
!! class through a fortran interface.

!> \brief A fortran type class that allows the creation of
!! IRL's ObjectAllocationServer<Paraboloidarator> class along with enabling
!! some of its methods.
module f_ObjServer_Paraboloid_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_ObjServer_Paraboloid
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_ObjServer_Paraboloid

  type, public :: ObjServer_Paraboloid_type
    type(c_ObjServer_Paraboloid) :: c_object
  contains
    final :: ObjServer_Paraboloid_class_delete
  end type ObjServer_Paraboloid_type

  interface new
    module procedure ObjServer_Paraboloid_class_new
  end interface

  interface

    subroutine F_ObjServer_Paraboloid_new(this, a_number_to_allocate) &
      bind(C, name="c_ObjServer_Paraboloid_new")
      import
      implicit none
      type(c_ObjServer_Paraboloid) :: this
      integer(C_SIZE_T) :: a_number_to_allocate
    end subroutine F_ObjServer_Paraboloid_new

    subroutine F_ObjServer_Paraboloid_delete(this) &
      bind(C, name="c_ObjServer_Paraboloid_delete")
      import
      implicit none
      type(c_ObjServer_Paraboloid) :: this
    end subroutine F_ObjServer_Paraboloid_delete

  end interface

  contains

    subroutine ObjServer_Paraboloid_class_new(this, a_number_to_allocate)
      implicit none
      type(ObjServer_Paraboloid_type), intent(inout) :: this
      integer(IRL_LargeOffsetIndex_t) :: a_number_to_allocate
      call F_ObjServer_Paraboloid_new(this%c_object, a_number_to_allocate)
    end subroutine ObjServer_Paraboloid_class_new

    impure elemental subroutine ObjServer_Paraboloid_class_delete(this)
      implicit none
      type(ObjServer_Paraboloid_type), intent(in) :: this
      call F_ObjServer_Paraboloid_delete(this%c_object)
    end subroutine ObjServer_Paraboloid_class_delete

end module f_ObjServer_Paraboloid_class
