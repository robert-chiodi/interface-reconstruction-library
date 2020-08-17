!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_constants.f90
!!
!! This file contains the Fortran interface to
!! IRL functions that deal with setting constants.

!> \brief This module contains mappings to the
!! IRL C interface that deal with setting global 
!! constants that are used in the IRL library.
module f_constants
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  interface
    subroutine F_setVFBounds(a_VF_low) &
      bind(C, name="c_setVFBounds")
      import
      implicit none
      real(C_DOUBLE), intent(in) :: a_VF_low
    end subroutine F_setVFBounds
  end interface

  interface
    subroutine F_setVFTolerance_IterativeDistanceFinding(a_tolerance) &
      bind(C, name="c_setVFTolerance_IterativeDistanceFinding")
      import
      implicit none
      real(C_DOUBLE), intent(in) :: a_tolerance
    end subroutine F_setVFTolerance_IterativeDistanceFinding
  end interface

  interface
    subroutine F_setMinimumVolToTrack(a_minimum_volume_to_track) &
      bind(C, name="c_setMinimumVolToTrack")
      import
      implicit none
      real(C_DOUBLE), intent(in) :: a_minimum_volume_to_track
    end subroutine F_setMinimumVolToTrack
  end interface

  interface
    subroutine F_setMinimumSAToTrack(a_minimum_surface_area_to_track) &
      bind(C, name="c_setMinimumSAToTrack")
      import
      implicit none
      real(C_DOUBLE), intent(in) :: a_minimum_surface_area_to_track
    end subroutine F_setMinimumSAToTrack
  end interface

contains

  subroutine setVFBounds(a_VF_low)
    implicit none
    real(IRL_double), intent(in) :: a_VF_low
    call F_setVFBounds(a_VF_low)
  end subroutine setVFBounds

  subroutine setVFTolerance_IterativeDistanceFinding(a_tolerance)
    implicit none
    real(IRL_double), intent(in) :: a_tolerance
    call F_setVFTolerance_IterativeDistanceFinding(a_tolerance)
  end subroutine setVFTolerance_IterativeDistanceFinding

  subroutine setMinimumVolToTrack(a_minimum_volume_to_track)
    implicit none
    real(IRL_double), intent(in) :: a_minimum_volume_to_track
    call F_setMinimumVolToTrack(a_minimum_volume_to_track)
  end subroutine setMinimumVolToTrack

  subroutine setMinimumSAToTrack(a_minimum_surface_area_to_track)
    implicit none
    real(IRL_double), intent(in) :: a_minimum_surface_area_to_track
    call F_setMinimumSAToTrack(a_minimum_surface_area_to_track)
  end subroutine setMinimumSAToTrack


end module f_constants
