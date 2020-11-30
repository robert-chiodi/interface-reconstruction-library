!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

! In this example, the ability to
! switch between cutting methods
! at runtime for the C/Fortran
! IRL interface is demonstrated.
! All three methods are capable of
! performing the same calculations,
! with various performance differences
! depending on the case configuration.
!
! Note, the ability to change cutting
! method at runtime is only available if the
! compile-time flag C_STATIC_CUTTING is NOT
! defined.

program main
  use irl_fortran_interface
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  real(DP), dimension(1:3,1:8) :: cube_pts
  type(SepVM_type) :: phase_moments_recursive_simplex
  type(SepVM_type) :: phase_moments_half_edge
  type(SepVM_type) :: phase_moments_iterative_simplex
  type(RectCub_type) :: cube
  type(PlanarSep_type) :: planar_separator

  real(DP), dimension(3) :: plane_normal
  real(DP) :: plane_distance

  ! Define unit-cubic cell
  cube_pts(1:3, 1) = (/ 0.5_DP, -0.5_DP, -0.5_DP/)
  cube_pts(1:3, 2) = (/ 0.5_DP,  0.5_DP, -0.5_DP/)
  cube_pts(1:3, 3) = (/ 0.5_DP,  0.5_DP,  0.5_DP/)
  cube_pts(1:3, 4) = (/ 0.5_DP, -0.5_DP,  0.5_DP/)
  cube_pts(1:3, 5) = (/-0.5_DP, -0.5_DP, -0.5_DP/)
  cube_pts(1:3, 6) = (/-0.5_DP,  0.5_DP, -0.5_DP/)
  cube_pts(1:3, 7) = (/-0.5_DP,  0.5_DP,  0.5_DP/)
  cube_pts(1:3, 8) = (/-0.5_DP, -0.5_DP,  0.5_DP/)

  ! Allocate cube storage and construct to be cell(:,:)
  call new(cube)
  call construct(cube,cube_pts)

  ! Construct the PlanarSeparator object in IRL.
  call new(planar_separator)

  ! Define interface reconstruction representing a x-z sheet
  ! centered -0.25 from cell center
  plane_normal = (/0.0_DP, 1.0_DP,0.0_DP/)
  plane_distance = -0.2_DP
  call addPlane(planar_separator,plane_normal,plane_distance)
  plane_normal = (/0.0_DP, -1.0_DP,0.0_DP/)
  plane_distance = 0.3_DP
  call addPlane(planar_separator,plane_normal,plane_distance)

  ! Perform cutting through IRL library
  ! After using getVolumeMoments_setMethod,
  ! subsequent calls to getMoments and
  ! getNormMoments will use
  ! that cutting method to calculate the
  ! volume moments.
  call new(phase_moments_recursive_simplex)
  call new(phase_moments_half_edge)
  call new(phase_moments_iterative_simplex)
  call getMoments_setMethod(0) ! Recursive Simplex Cutting
  call getNormMoments(cube, planar_separator, phase_moments_recursive_simplex)
  call getMoments_setMethod(1) ! Half Edge Cutting
  call getNormMoments(cube, planar_separator, phase_moments_half_edge)
  call getMoments_setMethod(2) ! Iterative Simplex Cutting
  call getNormMoments(cube, planar_separator, phase_moments_iterative_simplex)

  ! Print out the computed results.
  ! The phase moments are stored as internal to the PlanarSeparator (0)
  ! and external to the PlanarSeparator (1)
  write(*,'(A)')
  write(*,'(A)') 'Comparison between expected and computed results'
  write(*,'(A)') '================================================'
  write(*,'(A)') 'Volume between planes '
  write(*,'(A,F10.5)')  '   Expected                     : ', 0.1_DP
  write(*,'(A,F10.5)')  '   Computed (Recursive Simplex) : ',getVolume(phase_moments_recursive_simplex,0)
  write(*,'(A,F10.5)')  '   Computed (Half Edge)         : ',getVolume(phase_moments_half_edge,0)
  write(*,'(A,F10.5)')  '   Computed (Iterative Simplex) : ',getVolume(phase_moments_iterative_simplex,0)
  write(*,'(A)') 'Centroid for volume between planes '
  write(*,'(A,3F10.5)') '   Expected                     : ', 0.0_DP, -0.25_DP, 0.0_DP
  write(*,'(A,3F10.5)') '   Computed (Recursive Simplex) : ',getCentroid(phase_moments_recursive_simplex,0)
  write(*,'(A,3F10.5)') '   Computed (Half Edge)         : ',getCentroid(phase_moments_half_edge,0)
  write(*,'(A,3F10.5)') '   Computed (Iterative Simplex) : ',getCentroid(phase_moments_iterative_simplex,0)
  write(*,'(A)') 'Volume outside of planes '
  write(*,'(A,F10.5)')  '   Expected                     : ', 0.9_DP
  write(*,'(A,F10.5)')  '   Computed (Recursive Simplex) : ',getVolume(phase_moments_recursive_simplex,1)
  write(*,'(A,F10.5)')  '   Computed (Half Edge)         : ',getVolume(phase_moments_half_edge,1)
  write(*,'(A,F10.5)')  '   Computed (Iterative Simplex) : ',getVolume(phase_moments_iterative_simplex,1)
  write(*,'(A)') 'Centroid for volume outside of planes '
  write(*,'(A,3F10.5)') '   Expected                     : ', 0.0_DP, -(0.1_DP*(-0.25_DP))/0.9_DP, 0.0_DP
  write(*,'(A,3F10.5)') '   Computed (Recursive Simplex) : ',getCentroid(phase_moments_recursive_simplex,1)
  write(*,'(A,3F10.5)') '   Computed (Half Edge)         : ',getCentroid(phase_moments_half_edge,1)
  write(*,'(A,3F10.5)') '   Computed (Iterative Simplex) : ',getCentroid(phase_moments_iterative_simplex,1)
  write(*,'(A)')

end program main
