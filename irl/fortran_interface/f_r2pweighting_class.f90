!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2021 Austin Han <austinhhan@outlook.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_r2pweighting_class.f90
!!
!! This file contains the Fortran interface for the
!! R2PWeighting class.

!> \brief A fortran type class that allows the creation of
!! IRL's R2PWeighting type along with enabling
!! some of its methods.
module f_R2PWeighting_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_R2PWeighting
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_R2PWeighting

  type, public :: R2PWeighting_type
    type(c_R2PWeighting) :: c_object
  contains
    final :: R2PWeighting_class_delete
  end type R2PWeighting_type

  interface new
    module procedure R2PWeighting_class_new
  end interface
  interface setImportances
    module procedure R2PWeighting_class_setImportances
  end interface
  interface setImportanceOfLiquidVolumeFraction
    module procedure R2PWeighting_class_setImpOfLiqVolFrac
  end interface
  interface setImportanceOfLiquidCentroidRelativeToGas
    module procedure R2PWeighting_class_setImpOfLiqCentRelToGas
  end interface
  interface setImportanceOfCentroid
    module procedure R2PWeighting_class_setImpOfCentroid
  end interface
  interface setImportanceOfSurfaceArea
    module procedure R2PWeighting_class_setImpOfSurfArea
  end interface      
  interface getImportances
    module procedure R2PWeighting_class_getImportances
  end interface


  interface

    subroutine F_R2PWeighting_new(this) &
      bind(C, name="c_R2PWeighting_new")
      import
      implicit none
      type(c_R2PWeighting) :: this
    end subroutine F_R2PWeighting_new

    subroutine F_R2PWeighting_delete(this) &
      bind(C, name="c_R2PWeighting_delete")
      import
      implicit none
      type(c_R2PWeighting) :: this
    end subroutine F_R2PWeighting_delete

    subroutine F_R2PWeighting_setImportances(this, a_importances) &
      bind(C, name="c_R2PWeighting_setImportances")
      import
      implicit none
      type(c_R2PWeighting) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_importances !  dimension(4)
    end subroutine F_R2PWeighting_setImportances

    subroutine F_R2PWeighting_setImpOfLiqVolFrac(this, a_importance) &
      bind(C, name="c_R2PWeighting_setImpOfLiqVolFrac")
      import
      implicit none
      type(c_R2PWeighting) :: this
      real(C_DOUBLE), intent(in) :: a_importance
    end subroutine F_R2PWeighting_setImpOfLiqVolFrac

    subroutine F_R2PWeighting_setImpOfLiqCentRelToGas(this, a_importance) &
      bind(C, name="c_R2PWeighting_setImpOfLiqCentRelToGas")
      import
      implicit none
      type(c_R2PWeighting) :: this
      real(C_DOUBLE), intent(in) :: a_importance
    end subroutine F_R2PWeighting_setImpOfLiqCentRelToGas

    subroutine F_R2PWeighting_setImpOfCentroid(this, a_importance) &
      bind(C, name="c_R2PWeighting_setImpOfCentroid")
      import
      implicit none
      type(c_R2PWeighting) :: this
      real(C_DOUBLE), intent(in) :: a_importance
    end subroutine F_R2PWeighting_setImpOfCentroid

    subroutine F_R2PWeighting_setImpOfSurfArea(this, a_importance) &
      bind(C, name="c_R2PWeighting_setImpOfSurfArea")
      import
      implicit none
      type(c_R2PWeighting) :: this
      real(C_DOUBLE), intent(in) :: a_importance
    end subroutine F_R2PWeighting_setImpOfSurfArea

    subroutine F_R2PWeighting_getImportances(this, a_importances) &
      bind(C, name="c_R2PWeighting_getImportances")
      import
      implicit none
      type(c_R2PWeighting) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_importances !  dimension(4)
    end subroutine F_R2PWeighting_getImportances

  end interface


  contains

    subroutine R2PWeighting_class_new(this)
      implicit none
      type(R2PWeighting_type), intent(inout) :: this
      call F_R2PWeighting_new(this%c_object)
    end subroutine R2PWeighting_class_new

    impure elemental subroutine R2PWeighting_class_delete(this)
      implicit none
      type(R2PWeighting_type), intent(in) :: this
      call F_R2PWeighting_delete(this%c_object)
    end subroutine R2PWeighting_class_delete

    subroutine R2PWeighting_class_setImportances(this, a_importances)
      implicit none
      type(R2PWeighting_type), intent(in) :: this
      real(IRL_double), dimension(4), intent(in) :: a_importances
      call F_R2PWeighting_setImportances(this%c_object, a_importances)
    end subroutine R2PWeighting_class_setImportances

    subroutine R2PWeighting_class_setImpOfLiqVolFrac(this, a_importance)
      implicit none
      type(R2PWeighting_type), intent(in) :: this
      real(IRL_double), intent(in) :: a_importance
      call F_R2PWeighting_setImpOfLiqVolFrac(this%c_object, a_importance)
    end subroutine R2PWeighting_class_setImpOfLiqVolFrac

    subroutine R2PWeighting_class_setImpOfLiqCentRelToGas(this, a_importance)
      implicit none
      type(R2PWeighting_type), intent(in) :: this
      real(IRL_double), intent(in) :: a_importance
      call F_R2PWeighting_setImpOfLiqCentRelToGas(this%c_object, a_importance)
    end subroutine R2PWeighting_class_setImpOfLiqCentRelToGas

    subroutine R2PWeighting_class_setImpOfCentroid(this, a_importance)
      implicit none
      type(R2PWeighting_type), intent(in) :: this
      real(IRL_double), intent(in) :: a_importance
      call F_R2PWeighting_setImpOfCentroid(this%c_object, a_importance)
    end subroutine R2PWeighting_class_setImpOfCentroid

    subroutine R2PWeighting_class_setImpOfSurfArea(this, a_importance)
      implicit none
      type(R2PWeighting_type), intent(in) :: this
      real(IRL_double), intent(in) :: a_importance
      call F_R2PWeighting_setImpOfSurfArea(this%c_object, a_importance)
    end subroutine R2PWeighting_class_setImpOfSurfArea

    function R2PWeighting_class_getImportances(this) result(a_importances)
      implicit none
      type(R2PWeighting_type), intent(in) :: this
      real(IRL_double), dimension(4) :: a_importances
      call F_R2PWeighting_getImportances(this%c_object, a_importances)
    end function R2PWeighting_class_getImportances


end module f_R2PWeighting_class
