!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https:!mozilla.org/MPL/2.0/.

! In this example, the ELVIRA
! reconstruction method implemented inside
! IRL will be demonstrated. First, a PlanarSeparator
! will be directly specified and used to obtain
! SeperatedVolumeMoments for a 3x3 stencil. This will
! then be used to perform a ELVIRA reconstruction
! using the liquid volume fraction. The calculated interface
! is then shown to directly recover the planar interface
! that created the original volume fraction field.
! While this is shown for 2D, the exact same process
! can be done for 3D, with 27 members in the
! ELVIRANeighborhood due to the increase in
! stencil size.
program main
  use irl_fortran_interface
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  real(DP), dimension(1:3,1:8) :: rectangular_cuboid_pts
  real(DP), dimension(-1:1,-1:1) :: volume_fraction
  real(DP), dimension(1:3) :: normal
  real(DP), dimension(1:4) :: plane
  type(ELVIRANeigh_type) :: elvira_neighborhood
  type(RectCub_type), dimension(-1:1,-1:1) :: rectangular_cuboid
  type(PlanarSep_type) :: correct_planar_separator
  type(PlanarSep_type) :: found_planar_separator

  integer :: i, j

  ! Allocate storage for two PlanarSeparator and request objects from it
  call new(correct_planar_separator)
  call new(found_planar_separator)

  ! Define interface reconstruction representing a slanted
  ! line across the cell
  normal = (/0.5_DP*sqrt(3.0_DP), 0.5_DP,0.0_DP/)
  call addPlane(correct_planar_separator,normal,-0.15_DP)


  ! We are now going to setup the neighborhood of information
  ! that is needed to reconstruct with ELVIRA. Here, we will
  ! set the ELVIRANeighborhood to use 9 cells, and then
  ! construct the member variables and set them for the
  ! neighborhood. We will use the correct_planar_separator
  ! to generate the liquid volume fraction data for each cell,
  ! and then show that this can be directly recovered using
  ! ELVIRA.

  ! Allocate and obtain ELVIRANeighborhood object
  call new(elvira_neighborhood)

  ! Set its size 9 (since 3x3 stencil)
  call setSize(elvira_neighborhood,9)

  do j = -1,1
    do i = -1,1
      ! Set unit cell
      rectangular_cuboid_pts(1:3, 1) = (/ 0.5_DP, -0.5_DP, -0.5_DP/)
      rectangular_cuboid_pts(1:3, 2) = (/ 0.5_DP,  0.5_DP, -0.5_DP/)
      rectangular_cuboid_pts(1:3, 3) = (/ 0.5_DP,  0.5_DP,  0.5_DP/)
      rectangular_cuboid_pts(1:3, 4) = (/ 0.5_DP, -0.5_DP,  0.5_DP/)
      rectangular_cuboid_pts(1:3, 5) = (/-0.5_DP, -0.5_DP, -0.5_DP/)
      rectangular_cuboid_pts(1:3, 6) = (/-0.5_DP,  0.5_DP, -0.5_DP/)
      rectangular_cuboid_pts(1:3, 7) = (/-0.5_DP,  0.5_DP,  0.5_DP/)
      rectangular_cuboid_pts(1:3, 8) = (/-0.5_DP, -0.5_DP,  0.5_DP/)

      ! Shift the unit cell
      rectangular_cuboid_pts(1,:) = rectangular_cuboid_pts(1,:) + real(i,DP)
      rectangular_cuboid_pts(2,:) = rectangular_cuboid_pts(2,:) + real(j,DP)

      ! Allocate a rectangular cuboid and set its vertices
      call new(rectangular_cuboid(i,j))
      call construct(rectangular_cuboid(i,j),rectangular_cuboid_pts)

      ! Perform cutting through IRL library to obtain resulting volume
      ! internal to correct_planar_separator reconstruction.
      ! Since cell volumes are 1.0 in this case, the volume is also the volume fraction.
      call getNormMoments(rectangular_cuboid(i,j), correct_planar_separator, volume_fraction(i,j))

      ! Set the member in the ELVIRANeighborhood using the cell and corresponding volume fraction.
      ! For ELVIRA in 2D, the last variable in setMember should be -1.
      call setMember(elvira_neighborhood,rectangular_cuboid(i,j), volume_fraction(i,j), i, j, -1)

    end do
  end do

  ! With the whole neighborhood setup, we can now find the PlanarSeparator from ELVIRA.
  call reconstructELVIRA2D(elvira_neighborhood, found_planar_separator)

  ! Obtain the plane as a 4-element vector to compare to initial normal
  plane = getPlane(found_planar_separator,0)

  write(*,'(A)')
  write(*,'(A)') 'Comparison between given and computed normal'
  write(*,'(A)') '================================================'
  write(*,'(A)') 'Normal vectors '
  write(*,'(A,3F10.5)') '   Given    : ', normal(1:3)
  write(*,'(A,3F10.5)') '   Computed : ', plane(1:3)
  write(*,'(A)')
  write(*,'(A)') 'Plane distance '
  write(*,'(A,3F10.5)') '   Given    : ', -0.15_DP
  write(*,'(A, F10.5)') '   Computed : ', plane(4)
  write(*,'(A)')

end program main
