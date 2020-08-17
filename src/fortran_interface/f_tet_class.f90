!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_Tet_class.f90
!!
!! This file contains the Fortran interface for the
!! Tet class.

!> \brief A fortran type class that allows the creation of
!! IRL's Tet class along with enabling
!! some of its methods.
module f_Tet_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_Tet
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_Tet

  type, public :: Tet_type
    type(c_Tet) :: c_object
  contains
    final :: Tet_class_delete
  end type Tet_type

  interface new
    module procedure Tet_class_new
  end interface
  interface construct
    module procedure Tet_class_construct
  end interface
  interface calculateVolume
    module procedure Tet_class_calculateVolume
  end interface
  interface getBoundingPts
    module procedure Tet_class_getBoundingPts
  end interface
  interface printToScreen
    module procedure Tet_class_printToScreen
  end interface printToScreen

  interface

    subroutine F_Tet_new(this) &
      bind(C, name="c_Tet_new")
      import
      implicit none
      type(c_Tet) :: this
    end subroutine F_Tet_new

    subroutine F_Tet_delete(this) &
      bind(C, name="c_Tet_delete")
      import
      implicit none
      type(c_Tet) :: this
    end subroutine F_Tet_delete

    subroutine F_Tet_construct(this, a_Tet_pts) &
      bind(C, name="c_Tet_construct")
      import
      implicit none
      type(c_Tet) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_Tet_pts ! dimension(1:3,1:4)
    end subroutine F_Tet_construct
    
    function F_Tet_calculateVolume(this) result(a_tet_volume) &
      bind(C, name="c_Tet_calculateVolume")
      import
      implicit none
      type(c_Tet) :: this
      real(C_DOUBLE) :: a_tet_volume
    end function F_Tet_calculateVolume

    subroutine F_Tet_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_Tet_getBoundingPts")
      import
      implicit none
      type(c_Tet) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_Tet_getBoundingPts
    
    subroutine F_Tet_printToScreen(this) &
      bind(C, name="c_Tet_printToScreen")
      import
      implicit none
      type(c_Tet) :: this
    end subroutine F_Tet_printToScreen

  end interface

  contains

    subroutine Tet_class_new(this)
      implicit none
      type(Tet_type), intent(inout) :: this
      call F_Tet_new(this%c_object)
    end subroutine Tet_class_new

    impure elemental subroutine Tet_class_delete(this)
      implicit none
      type(Tet_type), intent(in) :: this
      call F_Tet_delete(this%c_object)
    end subroutine Tet_class_delete

    subroutine Tet_class_construct(this, a_Tet_pts)
      implicit none
      type(Tet_type), intent(inout) :: this
      real(IRL_double), dimension(1:3,1:4), intent(in) :: a_Tet_pts
      call F_Tet_construct(this%c_object, a_Tet_pts)
    end subroutine Tet_class_construct
    
    function Tet_class_calculateVolume(this) result(a_tet_volume)
      implicit none
      type(Tet_type), intent(inout) :: this
      real(IRL_double) :: a_tet_volume
      a_tet_volume = F_Tet_calculateVolume(this%c_object)
    end function Tet_class_calculateVolume

    subroutine Tet_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(Tet_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_Tet_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine Tet_class_getBoundingPts

    subroutine Tet_class_printToScreen(this)
      implicit none
      type(Tet_type), intent(inout) :: this
      call F_Tet_printToScreen(this%c_object)
    end subroutine Tet_class_printToScreen

end module f_Tet_class
