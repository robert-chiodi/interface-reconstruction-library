!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_PlanarSeparatorPathGroup_class.f90
!!
!! This file allows use of the IRL PlanarSeparatorPathGroup
!! class through a fortran interface.

!> \brief A fortran type class that allows the creation of
!! IRL's PlanarSeparatorPathGroup class along with enabling
!! some of its methods.
module f_PlanarSepPathGroup_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  use f_PlanarSepPath_class
  implicit none

  type, public, bind(C) :: c_PlanarSepPathGroup
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_PlanarSepPathGroup

  type, public :: PlanarSepPathGroup_type
    type(c_PlanarSepPathGroup) :: c_object
  contains
    final :: PlanarSepPathGroup_class_delete
  end type PlanarSepPathGroup_type

  interface new
    module procedure PlanarSepPathGroup_class_new
  end interface
  interface addPlanarSeparatorPath
    module procedure PlanarSepPathGroup_class_addPlanarSeparatorPath    
    module procedure PlanarSepPathGroup_class_addPlanarSeparatorPath_Id
  end interface 
  interface setPriorityOrder
    module procedure PlanarSepPathGroup_class_setPriorityOrder
  end interface 
  interface getPriorityOrderSize
    module procedure PlanarSepPathGroup_class_getPriorityOrderSize
  end interface 
  interface getPriorityOrderTag
    module procedure PlanarSepPathGroup_class_getPriorityOrderTag
  end interface 
  
  interface

    subroutine F_PlanarSepPathGroup_new(this) &
      bind(C, name="c_PlanarSepPathGroup_new")
      import
      implicit none
      type(c_PlanarSepPathGroup) :: this
    end subroutine F_PlanarSepPathGroup_new

    subroutine F_PlanarSepPathGroup_delete(this) &
      bind(C, name="c_PlanarSepPathGroup_delete")
      import
      implicit none
      type(c_PlanarSepPathGroup) :: this
    end subroutine F_PlanarSepPathGroup_delete

    subroutine F_PlanarSepPathGroup_addPlanarSeparatorPath(this, a_planar_separator_path) &
      bind(C, name="c_PlanarSepPathGroup_addPlanarSeparatorPath")
      import
      implicit none
      type(c_PlanarSepPathGroup) :: this
      type(c_PlanarSepPath) :: a_planar_separator_path
    end subroutine F_PlanarSepPathGroup_addPlanarSeparatorPath

    subroutine F_PlanarSepPathGroup_addPlanarSeparatorPath_Id(this, a_planar_separator_path, a_id) &
      bind(C, name="c_PlanarSepPathGroup_addPlanarSeparatorPath_Id")
      import
      implicit none
      type(c_PlanarSepPathGroup) :: this
      type(c_PlanarSepPath) :: a_planar_separator_path
      integer(C_INT) :: a_id
    end subroutine F_PlanarSepPathGroup_addPlanarSeparatorPath_Id

    subroutine F_PlanarSepPathGroup_setPriorityOrder(this, a_size, a_priority_order) &
      bind(C, name="c_PlanarSepPathGroup_setPriorityOrder")
      import
      implicit none
      type(c_PlanarSepPathGroup) :: this
      integer(C_INT) :: a_size
      integer(C_INT), dimension(*) :: a_priority_order ! Dimension [1:a_size]
    end subroutine F_PlanarSepPathGroup_setPriorityOrder

    function F_PlanarSepPathGroup_getPriorityOrderSize(this) result(a_size)&
      bind(C, name="c_PlanarSepPathGroup_getPriorityOrderSize")
      import
      implicit none
      type(c_PlanarSepPathGroup) :: this
      integer(C_INT) :: a_size
    end function F_PlanarSepPathGroup_getPriorityOrderSize

    function F_PlanarSepPathGroup_getPriorityOrderTag(this, a_index) result(a_tag)&
      bind(C, name="c_PlanarSepPathGroup_getPriorityOrderTag")
      import
      implicit none
      type(c_PlanarSepPathGroup) :: this
      integer(C_INT) :: a_index
      integer(C_INT) :: a_tag
    end function F_PlanarSepPathGroup_getPriorityOrderTag
    
  end interface
    
  contains

    subroutine PlanarSepPathGroup_class_new(this)
      implicit none
      type(PlanarSepPathGroup_type), intent(inout) :: this
      call F_PlanarSepPathGroup_new(this%c_object)
    end subroutine PlanarSepPathGroup_class_new

    impure elemental subroutine PlanarSepPathGroup_class_delete(this)
      implicit none
      type(PlanarSepPathGroup_type), intent(in) :: this
      call F_PlanarSepPathGroup_delete(this%c_object)
    end subroutine PlanarSepPathGroup_class_delete

    subroutine PlanarSepPathGroup_class_addPlanarSeparatorPath(this, a_planar_separator_path)
      implicit none
      type(PlanarSepPathGroup_type), intent(in) :: this
      type(PlanarSepPath_type) :: a_planar_separator_path

      call F_PlanarSepPathGroup_addPlanarSeparatorPath(this%c_object, &
                                                             a_planar_separator_path%c_object)
    end subroutine PlanarSepPathGroup_class_addPlanarSeparatorPath

    subroutine PlanarSepPathGroup_class_addPlanarSeparatorPath_Id(this, a_planar_separator_path, a_id)
      implicit none
      type(PlanarSepPathGroup_type), intent(in) :: this
      type(PlanarSepPath_type) :: a_planar_separator_path
      integer, intent(in) :: a_id

      call F_PlanarSepPathGroup_addPlanarSeparatorPath_Id(this%c_object, &
                                                             a_planar_separator_path%c_object, a_id)
    end subroutine PlanarSepPathGroup_class_addPlanarSeparatorPath_Id

    subroutine PlanarSepPathGroup_class_setPriorityOrder(this, a_size, a_priority_order)
      implicit none
      type(PlanarSepPathGroup_type), intent(in) :: this
      integer, intent(in) :: a_size
      integer, dimension(:), intent(in) :: a_priority_order

      call F_PlanarSepPathGroup_setPriorityOrder(this%c_object, a_size, a_priority_order)
      
    end subroutine PlanarSepPathGroup_class_setPriorityOrder

    function PlanarSepPathGroup_class_getPriorityOrderSize(this) result(a_size)
      implicit none
      type(PlanarSepPathGroup_type), intent(in) :: this
      integer:: a_size

      a_size = F_PlanarSepPathGroup_getPriorityOrderSize(this%c_object)
      return
    end function PlanarSepPathGroup_class_getPriorityOrderSize

    function PlanarSepPathGroup_class_getPriorityOrderTag(this, a_index) result(a_tag)
      implicit none
      type(PlanarSepPathGroup_type), intent(in) :: this
      integer, intent(in) :: a_index
      integer:: a_tag

      a_tag = F_PlanarSepPathGroup_getPriorityOrderTag(this%c_object, a_index)
      return
    end function PlanarSepPathGroup_class_getPriorityOrderTag
    
end module f_PlanarSepPathGroup_class
