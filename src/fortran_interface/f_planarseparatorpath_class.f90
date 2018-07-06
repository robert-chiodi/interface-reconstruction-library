!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_PlanarSeparatorPath_class.f90
!!
!! This file allows use of the IRL PlanarSeparatorPath
!! class through a fortran interface.

!> \brief A fortran type class that allows the creation of
!! IRL's PlanarSeparatorPath class along with enabling
!! some of its methods.
module f_PlanarSepPath_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  use f_PlanarSep_class
  implicit none

  type, public, bind(C) :: c_PlanarSepPath
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_PlanarSepPath

  type, public :: PlanarSepPath_type
    type(c_PlanarSepPath) :: c_object
  contains
    final :: PlanarSepPath_class_delete
  end type PlanarSepPath_type

  interface new
    module procedure PlanarSepPath_class_new
    module procedure PlanarSepPath_class_new_PlanarSep    
  end interface new
  interface construct
    module procedure PlanarSepPath_class_construct
  end interface
  interface setId
    module procedure PlanarSepPath_class_setId
  end interface
  interface getId
    module procedure PlanarSepPath_class_getId
  end interface
  interface setEdgeConnectivity
    module procedure PlanarSepPath_class_setEdgeConnectivity
  end interface
  interface setEdgeConnectivityNull
    module procedure PlanarSepPath_class_setEdgeConnectivityNull
  end interface

  interface

    subroutine F_PlanarSepPath_new(this) &
      bind(C, name="c_PlanarSepPath_new")
      import
      implicit none
      type(c_PlanarSepPath) :: this
    end subroutine F_PlanarSepPath_new
    
    subroutine F_PlanarSepPath_new_PlanarSep(this, a_planar_separator) &
      bind(C, name="c_PlanarSepPath_new_PlanarSep")
      import
      implicit none
      type(c_PlanarSepPath) :: this
      type(c_PlanarSep) :: a_planar_separator
    end subroutine F_PlanarSepPath_new_PlanarSep

    subroutine F_PlanarSepPath_construct(this, a_planar_separator) &
      bind(C, name="c_PlanarSepPath_construct")
      import
      implicit none
      type(c_PlanarSepPath) :: this
      type(c_PlanarSep) :: a_planar_separator
    end subroutine F_PlanarSepPath_construct

    subroutine F_PlanarSepPath_delete(this) &
      bind(C, name="c_PlanarSepPath_delete")
      import
      implicit none
      type(c_PlanarSepPath) :: this
    end subroutine F_PlanarSepPath_delete

    subroutine F_PlanarSepPath_setId(this, a_id) &
      bind(C, name="c_PlanarSepPath_setId")
      import
      implicit none
      type(c_PlanarSepPath) :: this
      integer(C_INT), intent(in) :: a_id ! Scalar >= 0
    end subroutine F_PlanarSepPath_setId

    function F_PlanarSepPath_getId(this) result(a_id) &
      bind(C, name="c_PlanarSepPath_getId")
      import
      implicit none
      type(c_PlanarSepPath) :: this
      integer(C_INT) :: a_id ! Scalar >= 0
    end function F_PlanarSepPath_getId

    subroutine F_PlanarSepPath_setEdgeConnectivity(this, a_ptr_to_neighbor) &
      bind(C, name="c_PlanarSepPath_setEdgeConnectivity")
      import
      implicit none
      type(c_PlanarSepPath) :: this
      type(c_PlanarSepPath) :: a_ptr_to_neighbor ! Ptr to next PlanarSepPath
    end subroutine F_PlanarSepPath_setEdgeConnectivity

    subroutine F_PlanarSepPath_setEdgeConnectivityNull(this) &
      bind(C, name="c_PlanarSepPath_setEdgeConnectivityNull")
      import
      implicit none
      type(c_PlanarSepPath) :: this
    end subroutine F_PlanarSepPath_setEdgeConnectivityNull

  end interface


contains

    subroutine PlanarSepPath_class_new(this)
      implicit none
      type(PlanarSepPath_type), intent(inout) :: this
      call F_PlanarSepPath_new(this%c_object)
    end subroutine PlanarSepPath_class_new  

    subroutine PlanarSepPath_class_new_PlanarSep(this, a_planar_separator)
      implicit none
      type(PlanarSepPath_type), intent(inout) :: this
      type(PlanarSep_type), intent(in) :: a_planar_separator
      call F_PlanarSepPath_new_PlanarSep(this%c_object, a_planar_separator%c_object)
    end subroutine PlanarSepPath_class_new_PlanarSep

    subroutine PlanarSepPath_class_construct(this, a_planar_separator)
      implicit none
      type(PlanarSepPath_type), intent(inout) :: this
      type(PlanarSep_type), intent(in) :: a_planar_separator
      call F_PlanarSepPath_construct(this%c_object, a_planar_separator%c_object)
    end subroutine PlanarSepPath_class_construct
    
    impure elemental subroutine PlanarSepPath_class_delete(this)
      implicit none
      type(PlanarSepPath_type), intent(in) :: this
      call F_PlanarSepPath_delete(this%c_object)
    end subroutine PlanarSepPath_class_delete

    subroutine PlanarSepPath_class_setId(this, a_id)
      implicit none
      type(PlanarSepPath_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_id
      ! a_id must be >= 0
      call F_PlanarSepPath_setId(this%c_object, a_id)
    end subroutine PlanarSepPath_class_setId

    function PlanarSepPath_class_getId(this) result(a_id)
      implicit none
      type(PlanarSepPath_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t):: a_id
      ! a_id must be >= 0
      a_id = F_PlanarSepPath_getId(this%c_object)
    end function PlanarSepPath_class_getId

    subroutine PlanarSepPath_class_setEdgeConnectivity(this,  &
        a_neighboring_PlanarSepPath)
      implicit none
      type(PlanarSepPath_type), intent(in) :: this
      type(PlanarSepPath_type) :: a_neighboring_PlanarSepPath
      call F_PlanarSepPath_setEdgeConnectivity(this%c_object, a_neighboring_PlanarSepPath%c_object)
    end subroutine PlanarSepPath_class_setEdgeConnectivity

    subroutine PlanarSepPath_class_setEdgeConnectivityNull(this)
      implicit none
      type(PlanarSepPath_type), intent(in) :: this
      call F_PlanarSepPath_setEdgeConnectivityNull(this%c_object)
    end subroutine PlanarSepPath_class_setEdgeConnectivityNull

end module f_PlanarSepPath_class
