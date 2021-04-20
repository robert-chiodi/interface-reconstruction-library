!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2021 Austin Han <austinhhan@outlook.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_OptimizationBehavior_class.f90
!!
!! This file contains the Fortran interface for the
!! OptimizationBehavior class.

!> \brief A fortran type class that allows the creation of
!! IRL's OptimizationBehavior type along with enabling
!! some of its methods.
module f_OptimizationBehavior_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_OptimizationBehavior
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_OptimizationBehavior

  type, public :: OptimizationBehavior_type
    type(c_OptimizationBehavior) :: c_object
  contains
    final :: OptimizationBehavior_class_delete
  end type OptimizationBehavior_type

  interface new
    module procedure OptimizationBehavior_class_new
  end interface
  interface setAcceptableError
    module procedure OptimizationBehavior_class_setAcceptableError
  end interface
  interface setMaxIterations
    module procedure OptimizationBehavior_class_setMaxIterations
  end interface
  interface setMinAngleChange
    module procedure OptimizationBehavior_class_setMinAngleChange
  end interface
  interface setMinDistChange
    module procedure OptimizationBehavior_class_setMinDistChange
  end interface
  interface setLambdaIncrease
    module procedure OptimizationBehavior_class_setLambdaIncrease
  end interface
  interface setLambdaDecrease
    module procedure OptimizationBehavior_class_setLambdaDecrease
  end interface   
  interface setDelayJacobianAmt
    module procedure OptimizationBehavior_class_setDelayJacobianAmt
  end interface   
  interface setInitialAngle
    module procedure OptimizationBehavior_class_setInitialAngle
  end interface   
  interface setInitialDistance
    module procedure OptimizationBehavior_class_setInitialDistance
  end interface 
  interface setFinDiffAngle
    module procedure OptimizationBehavior_class_setFinDiffAngle
  end interface          
  interface getParameters
    module procedure OptimizationBehavior_class_getParameters
  end interface


  interface

    subroutine F_OptimizationBehavior_new(this) &
      bind(C, name="c_OptimizationBehavior_new")
      import
      implicit none
      type(c_OptimizationBehavior) :: this
    end subroutine F_OptimizationBehavior_new

    subroutine F_OptimizationBehavior_delete(this) &
      bind(C, name="c_OptimizationBehavior_delete")
      import
      implicit none
      type(c_OptimizationBehavior) :: this
    end subroutine F_OptimizationBehavior_delete

    subroutine F_OptimizationBehavior_setAcceptableError(this, a_parameter) &
      bind(C, name="c_OptimizationBehavior_setAcceptableError")
      import
      implicit none
      type(c_OptimizationBehavior) :: this
      real(C_DOUBLE), intent(in) :: a_parameter
    end subroutine F_OptimizationBehavior_setAcceptableError

    subroutine F_OptimizationBehavior_setMaxIterations(this, a_parameter) &
      bind(C, name="c_OptimizationBehavior_setMaxIterations")
      import
      implicit none
      type(c_OptimizationBehavior) :: this
      integer(C_INT), intent(in) :: a_parameter
    end subroutine F_OptimizationBehavior_setMaxIterations

    subroutine F_OptimizationBehavior_setMinAngleChange(this, a_parameter) &
      bind(C, name="c_OptimizationBehavior_setMinAngleChange")
      import
      implicit none
      type(c_OptimizationBehavior) :: this
      real(C_DOUBLE), intent(in) :: a_parameter
    end subroutine F_OptimizationBehavior_setMinAngleChange

    subroutine F_OptimizationBehavior_setMinDistChange(this, a_parameter) &
      bind(C, name="c_OptimizationBehavior_setMinDistChange")
      import
      implicit none
      type(c_OptimizationBehavior) :: this
      real(C_DOUBLE), intent(in) :: a_parameter
    end subroutine F_OptimizationBehavior_setMinDistChange

    subroutine F_OptimizationBehavior_setLambdaIncrease(this, a_parameter) &
      bind(C, name="c_OptimizationBehavior_setLambdaIncrease")
      import
      implicit none
      type(c_OptimizationBehavior) :: this
      real(C_DOUBLE), intent(in) :: a_parameter
    end subroutine F_OptimizationBehavior_setLambdaIncrease

    subroutine F_OptimizationBehavior_setLambdaDecrease(this, a_parameter) &
      bind(C, name="c_OptimizationBehavior_setLambdaDecrease")
      import
      implicit none
      type(c_OptimizationBehavior) :: this
      real(C_DOUBLE), intent(in) :: a_parameter
    end subroutine F_OptimizationBehavior_setLambdaDecrease

    subroutine F_OptimizationBehavior_setDelayJacobianAmt(this, a_parameter) &
      bind(C, name="c_OptimizationBehavior_setDelayJacobianAmt")
      import
      implicit none
      type(c_OptimizationBehavior) :: this
      integer(C_INT), intent(in) :: a_parameter
    end subroutine F_OptimizationBehavior_setDelayJacobianAmt

    subroutine F_OptimizationBehavior_setInitialAngle(this, a_parameter) &
      bind(C, name="c_OptimizationBehavior_setInitialAngle")
      import
      implicit none
      type(c_OptimizationBehavior) :: this
      real(C_DOUBLE), intent(in) :: a_parameter
    end subroutine F_OptimizationBehavior_setInitialAngle

    subroutine F_OptimizationBehavior_setInitialDistance(this, a_parameter) &
      bind(C, name="c_OptimizationBehavior_setInitialDistance")
      import
      implicit none
      type(c_OptimizationBehavior) :: this
      real(C_DOUBLE), intent(in) :: a_parameter
    end subroutine F_OptimizationBehavior_setInitialDistance

    subroutine F_OptimizationBehavior_setFinDiffAngle(this, a_parameter) &
      bind(C, name="c_OptimizationBehavior_setFinDiffAngle")
      import
      implicit none
      type(c_OptimizationBehavior) :: this
      real(C_DOUBLE), intent(in) :: a_parameter
    end subroutine F_OptimizationBehavior_setFinDiffAngle

    subroutine F_OptimizationBehavior_getParameters(this, a_parameters) &
      bind(C, name="c_OptimizationBehavior_getParameters")
      import
      implicit none
      type(c_OptimizationBehavior) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_parameters !  dimension(4)
    end subroutine F_OptimizationBehavior_getParameters

  end interface


  contains

    subroutine OptimizationBehavior_class_new(this)
      implicit none
      type(OptimizationBehavior_type), intent(inout) :: this
      call F_OptimizationBehavior_new(this%c_object)
    end subroutine OptimizationBehavior_class_new

    impure elemental subroutine OptimizationBehavior_class_delete(this)
      implicit none
      type(OptimizationBehavior_type), intent(in) :: this
      call F_OptimizationBehavior_delete(this%c_object)
    end subroutine OptimizationBehavior_class_delete

    subroutine OptimizationBehavior_class_setAcceptableError(this, a_parameter)
      implicit none
      type(OptimizationBehavior_type), intent(in) :: this
      real(IRL_double), intent(in) :: a_parameter
      call F_OptimizationBehavior_setAcceptableError(this%c_object, a_parameter)
    end subroutine OptimizationBehavior_class_setAcceptableError

    subroutine OptimizationBehavior_class_setMaxIterations(this, a_parameter)
      implicit none
      type(OptimizationBehavior_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_parameter
      call F_OptimizationBehavior_setMaxIterations(this%c_object, a_parameter)
    end subroutine OptimizationBehavior_class_setMaxIterations

    subroutine OptimizationBehavior_class_setMinAngleChange(this, a_parameter)
      implicit none
      type(OptimizationBehavior_type), intent(in) :: this
      real(IRL_double), intent(in) :: a_parameter
      call F_OptimizationBehavior_setMinAngleChange(this%c_object, a_parameter)
    end subroutine OptimizationBehavior_class_setMinAngleChange

    subroutine OptimizationBehavior_class_setMinDistChange(this, a_parameter)
      implicit none
      type(OptimizationBehavior_type), intent(in) :: this
      real(IRL_double), intent(in) :: a_parameter
      call F_OptimizationBehavior_setMinDistChange(this%c_object, a_parameter)
    end subroutine OptimizationBehavior_class_setMinDistChange

    subroutine OptimizationBehavior_class_setLambdaIncrease(this, a_parameter)
      implicit none
      type(OptimizationBehavior_type), intent(in) :: this
      real(IRL_double), intent(in) :: a_parameter
      call F_OptimizationBehavior_setLambdaIncrease(this%c_object, a_parameter)
    end subroutine OptimizationBehavior_class_setLambdaIncrease

    subroutine OptimizationBehavior_class_setLambdaDecrease(this, a_parameter)
      implicit none
      type(OptimizationBehavior_type), intent(in) :: this
      real(IRL_double), intent(in) :: a_parameter
      call F_OptimizationBehavior_setLambdaDecrease(this%c_object, a_parameter)
    end subroutine OptimizationBehavior_class_setLambdaDecrease

    subroutine OptimizationBehavior_class_setDelayJacobianAmt(this, a_parameter)
      implicit none
      type(OptimizationBehavior_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_parameter
      call F_OptimizationBehavior_setDelayJacobianAmt(this%c_object, a_parameter)
    end subroutine OptimizationBehavior_class_setDelayJacobianAmt

    subroutine OptimizationBehavior_class_setInitialAngle(this, a_parameter)
      implicit none
      type(OptimizationBehavior_type), intent(in) :: this
      real(IRL_double), intent(in) :: a_parameter
      call F_OptimizationBehavior_setInitialAngle(this%c_object, a_parameter)
    end subroutine OptimizationBehavior_class_setInitialAngle

    subroutine OptimizationBehavior_class_setInitialDistance(this, a_parameter)
      implicit none
      type(OptimizationBehavior_type), intent(in) :: this
      real(IRL_double), intent(in) :: a_parameter
      call F_OptimizationBehavior_setInitialDistance(this%c_object, a_parameter)
    end subroutine OptimizationBehavior_class_setInitialDistance

    subroutine OptimizationBehavior_class_setFinDiffAngle(this, a_parameter)
      implicit none
      type(OptimizationBehavior_type), intent(in) :: this
      real(IRL_double), intent(in) :: a_parameter
      call F_OptimizationBehavior_setFinDiffAngle(this%c_object, a_parameter)
    end subroutine OptimizationBehavior_class_setFinDiffAngle

    function OptimizationBehavior_class_getParameters(this) result(a_parameters)
      implicit none
      type(OptimizationBehavior_type), intent(in) :: this
      real(IRL_double), dimension(4) :: a_parameters
      call F_OptimizationBehavior_getParameters(this%c_object, a_parameters)
    end function OptimizationBehavior_class_getParameters


end module f_OptimizationBehavior_class
