!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

  
! In this example a RectangularCuboid
! is cut by a two-plane PlanarSeparator that
! divides the RectangularCuboid into volumes internal
! and external to the separator.
! The computed results compared to the correct
! results are then printed to screen for comparison.
program main
  use irl_fortran_interface
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  real(DP), dimension(1:3,1:8) :: cube_pts
  type(SepVM_type) :: phase_moments
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

  ! Define unit-cubic cell using IRL RectangularCuboid class
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

  ! getNormalizedVolumeMoments is the main work horse of IRL.
  ! In general, it is capable of intersecting a
  ! polygon or polyhedron by a set of planes to return
  ! moments of the intersection between the polygon/polyhedron
  ! and the half-spaces of the planes.
  call new(phase_moments)
  call getNormMoments(cube, planar_separator, phase_moments)

  ! Print out the computed results.
  ! The phase moments are stored as internal to the PlanarSeparator (0)
  ! and external to the PlanarSeparator (1)
  write(*,'(A)')
  write(*,'(A)') 'Comparison between expected and computed results'
  write(*,'(A)') '================================================'
  write(*,'(A)') 'Volume between planes '
  write(*,'(A,F10.5)')  '   Expected: ', 0.1_DP
  write(*,'(A,F10.5)')  '   Computed: ',getVolume(phase_moments,0)
  write(*,'(A)') 'Centroid for volume between planes '
  write(*,'(A,3F10.5)') '   Expected: ', 0.0_DP, -0.25_DP, 0.0_DP
  write(*,'(A,3F10.5)') '   Computed: ',getCentroid(phase_moments,0)
  write(*,'(A)') 'Volume outside of planes '
  write(*,'(A,F10.5)')  '   Expected: ', 0.9_DP
  write(*,'(A,F10.5)')  '   Computed: ',getVolume(phase_moments,1)
  write(*,'(A)') 'Centroid for volume outside of planes '
  write(*,'(A,3F10.5)') '   Expected: ', 0.0_DP, -(0.1_DP*(-0.25_DP))/0.9_DP, 0.0_DP
  write(*,'(A,3F10.5)') '   Computed: ',getCentroid(phase_moments,1)
  write(*,'(A)')

end program main
