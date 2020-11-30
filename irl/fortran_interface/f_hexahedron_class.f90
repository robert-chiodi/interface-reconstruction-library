!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_Hexahedron_class.f90
!!
!! This file contains the Fortran interface for the
!! Hexahedron class.

!> \brief A fortran type class that allows the creation of
!! IRL's Hexahedron class along with enabling
!! some of its methods.
module f_Hex_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_Hex
    type(C_PTR), public :: object = C_NULL_PTR
  end type c_Hex

  type, public :: Hex_type
    type(c_Hex) :: c_object
  contains
    final :: Hex_class_delete
  end type Hex_type

  interface new
    module procedure Hex_class_new
  end interface
  interface construct
    module procedure Hex_class_construct
  end interface
  interface calculateVolume
    module procedure Hex_class_calculateVolume
  end interface
  interface getVertices
    module procedure Hex_class_getVertices
  end interface
  interface getBoundingPts
    module procedure Hex_class_getBoundingPts
  end interface
  interface printToScreen
    module procedure Hex_class_printToScreen
  end interface 

  interface

    subroutine F_Hex_new(this) &
      bind(C, name="c_Hex_new")
      import
      implicit none
      type(c_Hex) :: this
    end subroutine F_Hex_new

    subroutine F_Hex_delete(this) &
      bind(C, name="c_Hex_delete")
      import
      implicit none
      type(c_Hex) :: this
    end subroutine F_Hex_delete

    subroutine F_Hex_construct(this, a_cell) &
      bind(C, name="c_Hex_construct")
      import
      implicit none
      type(c_Hex) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_cell ! dimension(1:3,1:8)
    end subroutine F_Hex_construct

    function F_Hex_calculateVolume(this) result(a_hexahedron_volume) &
      bind(C, name="c_Hex_calculateVolume")
      import
      implicit none
      type(c_Hex) :: this
      real(C_DOUBLE) :: a_hexahedron_volume
    end function F_Hex_calculateVolume
    
    subroutine F_Hex_getVertices(this, a_pts) &
      bind(C, name="c_Hex_getVertices")
      import
      implicit none
      type(c_Hex) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_pts ! dimension(1:3,1:8)
    end subroutine F_Hex_getVertices

    subroutine F_Hex_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_Hex_getBoundingPts")
      import
      implicit none
      type(c_Hex) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_Hex_getBoundingPts
 
    subroutine F_Hex_printToScreen(this) &
      bind(C, name="c_Hex_printToScreen")
      import
      implicit none
      type(c_Hex) :: this
    end subroutine F_Hex_printToScreen 

  end interface


  contains

    subroutine Hex_class_new(this)
      implicit none
      type(Hex_type), intent(inout) :: this
      call F_Hex_new(this%c_object)
    end subroutine Hex_class_new

    impure elemental subroutine Hex_class_delete(this)
      use, intrinsic :: iso_c_binding
      implicit none
      type(Hex_type), intent(in) :: this
      call F_Hex_delete(this%c_object)
    end subroutine Hex_class_delete

    subroutine Hex_class_construct(this, a_cell)
      implicit none
      type(Hex_type), intent(inout) :: this
      real(IRL_double), dimension(1:3,1:8), intent(in) :: a_cell
      call F_Hex_construct(this%c_object, a_cell)
    end subroutine Hex_class_construct

    function Hex_class_calculateVolume(this) result(a_hexahedron_volume)
      implicit none
      type(Hex_type), intent(inout) :: this
      real(IRL_double) :: a_hexahedron_volume
      a_hexahedron_volume = F_Hex_calculateVolume(this%c_object)
    end function Hex_class_calculateVolume
    
    function Hex_class_getVertices(this) result(a_pts)
      implicit none
      type(Hex_type), intent(in) :: this
      real(IRL_double), dimension(1:3,1:8) :: a_pts
      call F_Hex_getVertices(this%c_object, a_pts)
      return
    end function Hex_class_getVertices

    subroutine Hex_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(Hex_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_Hex_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine Hex_class_getBoundingPts

    subroutine Hex_class_printToScreen(this)
      implicit none
      type(Hex_type), intent(inout) :: this
      call F_Hex_printToScreen(this%c_object)
    end subroutine Hex_class_printToScreen 


end module f_Hex_class
