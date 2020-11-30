!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

! In this example, the moment-of-fluid
! reconstruction (MoF) method implemented inside
! IRL will be demonstrated. First, a PlanarSeparator
! will be directly specified and used to obtain
! SeperatedVolumeMoments for a cell. This will
! then be used to perform a MoF reconstruction.
! The plane from the MoF reconstruction
! is compared to the one initially given.
program main
  use irl_fortran_interface
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  real(DP), dimension(1:3,1:8) :: rectangular_cuboid_pts
  type(SepVM_type) :: phase_moments
  real(DP), dimension(1:3) :: normal
  type(RectCub_type) :: rectangular_cuboid
  type(PlanarSep_type) :: correct_planar_separator
  type(PlanarSep_type) :: found_planar_separator_2D, found_planar_separator_3D
  type(PlanarSep_type) :: found_planar_separator_2D_weighted, found_planar_separator_3D_weighted
  real(DP), dimension(4) :: plane

  ! Define unit-cubic cell
  rectangular_cuboid_pts(1:3, 1) = (/ 0.5_DP, -0.5_DP, -0.5_DP/)
  rectangular_cuboid_pts(1:3, 2) = (/ 0.5_DP,  0.5_DP, -0.5_DP/)
  rectangular_cuboid_pts(1:3, 3) = (/ 0.5_DP,  0.5_DP,  0.5_DP/)
  rectangular_cuboid_pts(1:3, 4) = (/ 0.5_DP, -0.5_DP,  0.5_DP/)
  rectangular_cuboid_pts(1:3, 5) = (/-0.5_DP, -0.5_DP, -0.5_DP/)
  rectangular_cuboid_pts(1:3, 6) = (/-0.5_DP,  0.5_DP, -0.5_DP/)
  rectangular_cuboid_pts(1:3, 7) = (/-0.5_DP,  0.5_DP,  0.5_DP/)
  rectangular_cuboid_pts(1:3, 8) = (/-0.5_DP, -0.5_DP,  0.5_DP/)

  ! Allocate RectangularCuboid storage and construct to be unit cell
  call new(rectangular_cuboid)
  call construct(rectangular_cuboid, rectangular_cuboid_pts)


  ! Allocate the three PlanarSeparators we'll use
  call new(correct_planar_separator)
  call new(found_planar_separator_2D)
  call new(found_planar_separator_3D)
  call new(found_planar_separator_2D_weighted)
  call new(found_planar_separator_3D_weighted)  

  ! Allocate SeparatedMoments<VolumeMoments> objects to store cutting results
  call new(phase_moments)


  ! Define interface reconstruction representing a slanted
  ! line across the cell
  normal = (/0.5_DP*sqrt(3.0_DP), 0.5_DP,0.0_DP/)
  call addPlane(correct_planar_separator,normal,-0.15_DP)

  ! Perform cutting through IRL library to obtain resulting SeparatedMoments<VolumeMoments>
  call getNormMoments(rectangular_cuboid, correct_planar_separator, phase_moments)

  ! Now perform a MoF reconstruction and store the result
  ! By default, equal weights will be applied in the optimization
  ! when trying to recover the phase centroids.
  call reconstructMOF2D(rectangular_cuboid, phase_moments, found_planar_separator_2D)
  call reconstructMOF3D(rectangular_cuboid, phase_moments, found_planar_separator_3D)
  
  ! We can also supply our own weights for MoF to use, allowing us
  ! to bias matching one centroid over the other. Here, let's just
  ! try to match the internal centroid.
  call reconstructMOF2D(rectangular_cuboid, phase_moments, 1.0_DP, 0.0_DP, &
       found_planar_separator_2D_weighted)
  call reconstructMOF3D(rectangular_cuboid, phase_moments, 1.0_DP, 0.0_DP, &
       found_planar_separator_3D_weighted)
  
  write(*, '(A)')
  write(*, '(A)') 'Comparison between given and computed results'
  write(*, '(A)') '================================================'
  write(*, '(A)') 'Default (equal) weights '
  plane = getPlane(correct_planar_separator, 0)
  write(*, '(A,4F10.5)') '   Given Plane   : ', plane(1:4)
  plane = getPlane(found_planar_separator_2D, 0)  
  write(*, '(A,4F10.5)') '   Computed (2D) : ', plane(1:4)
  plane = getPlane(found_planar_separator_3D, 0)  
  write(*, '(A,4F10.5)') '   Computed (3D) : ', plane(1:4)
  write(*, '(A)')
  write(*, '(A)') 'Set (unequal) weights '
  plane = getPlane(correct_planar_separator, 0)
  write(*, '(A,4F10.5)') '   Given Plane   : ', plane(1:4)
  plane = getPlane(found_planar_separator_2D_weighted, 0)  
  write(*, '(A,4F10.5)') '   Computed (2D) : ', plane(1:4)
  plane = getPlane(found_planar_separator_3D_weighted, 0)  
  write(*, '(A,4F10.5)') '   Computed (3D) : ', plane(1:4)
  write(*, '(A)')  


end program
