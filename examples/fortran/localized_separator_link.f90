!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

! This example covers how to create
! a LocalizedSeparatorLink network to
! cut a polyhedron by a plane (as a PlanarSeparator)
!  (localized inside PlanarLocalizers).
! This is done to directly obtain
! volumetric moments for the polyhedron
! internal/external to the separator (below/above the plane)
! for each localizer given.

! This will be shown using a Hexahedron
! that lays across a collection of
! localizers representing RectangularCuboid objects and
! mimicking a uniform Cartesian mesh.
! These PlanarLocalizer planes will be set directly, and in general
! are capable of representing any mesh made up of
! convex polyhedra.
program main
  use irl_fortran_interface
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer :: i,j,k

  real(DP), dimension(1:3,1:8) :: dodecahedron_volume_pts
  type(TagAccVM_SepVM_type) :: phase_moments_LocalizedSeparatorLink

  type(Dod_type) :: dodecahedron_volume
  type(PlanarLoc_type), dimension(-1:1,-1:1,-1:1) :: planar_localizer
  type(PlanarSep_type), dimension(-1:1,-1:1,-1:1) :: planar_separator
  type(LocSepLink_type), dimension(-1:1,-1:1,-1:1) :: localized_separator_link
  integer :: unique_id

  real(DP), dimension(3) :: plane_normal
  real(DP) :: plane_distance

  type(SepVM_type) :: separated_volume_moments

  ! Allocate the actual objects in IRL
  do k = -1,1
    do j = -1,1
      do i = -1,1
        ! These relate to implicit constructors,
        ! allocating the memory for the PlanarLocalizer
        ! and PlanarSeparator
        call new(planar_localizer(i,j,k))
        call new(planar_separator(i,j,k))

        ! Construct the LocalizedSeparatorLink objects,
        ! providing the PlanarLocalizer and PlanarSeparator
        ! it is comprised of.
        ! Internal to IRL, this is done using pointers, so changes to the
        ! objects in the planar_localizer and planar_separator
        ! arrays will be reflected in the localized_separator_link.
        call new(localized_separator_link(i,j,k),&
                 planar_localizer(i,j,k),&
                 planar_separator(i,j,k))

      end do
    end do
  end do

  ! Construct the PlanarLocalizer objects to
  ! represent a [-1.5 x 1.5]^3 uniform cubic grid.
  ! This is performed in the below helper function
  ! `contained` at the bottom of this program.
  call makeCubicPlanarLocalizer()

  ! Let's now setup the PlanarSeparators involved, which
  ! will separate the volume of the hexahedron 
  ! after it is localized in the region dictated by each PlanarLocalizer
  ! in the localized_separator and localized_separator_link.
  ! For the sake of simplicity,
  ! lets separate by a flat plane across the middle of
  ! the domain with a normal (0.0,1.0,0.0) and a distance of 0.0.
  ! In general, the PlanarSeparator can have any plane orientation,
  ! and require no continuity between the PlanarSeparators in
  ! the LocalizedSeparatorLink.
  plane_normal = (/0.0_DP,1.0_DP,0.0_DP/)
  plane_distance = 0.0_DP  
  do k = -1,1
     do j = -1,1
        do i = -1,1
           call addPlane(planar_separator(i,j,k),plane_normal,plane_distance)
        end do
     end do
  end do
  
  ! Specify the connectivity amongst the different
  ! LocalizedSeparatorLink objects. This provides information
  ! on which two LocalizedSeparatorLink objects each plane
  ! in the contained PlanarLocalizer separates. In a graph
  ! sense, this specifies a directed edge from one
  ! LocalizedSeparatorLink to another (which can be
  ! considered noes in the graph).
  call setupLinking(localized_separator_link)


  ! Specify the polyhedron that will be separated, in this case
  ! a Hexahedron that intersects all cells in our
  ! mimicked mesh.
  dodecahedron_volume_pts(1:3, 1) = (/ 0.25_DP, -1.00_DP, -1.00_DP/)
  dodecahedron_volume_pts(1:3, 2) = (/ 1.00_DP,  1.00_DP, -1.00_DP/)
  dodecahedron_volume_pts(1:3, 3) = (/ 1.00_DP,  1.00_DP,  1.00_DP/)
  dodecahedron_volume_pts(1:3, 4) = (/ 0.25_DP, -1.00_DP,  1.00_DP/)
  dodecahedron_volume_pts(1:3, 5) = (/-0.25_DP, -1.00_DP, -1.00_DP/)
  dodecahedron_volume_pts(1:3, 6) = (/-1.00_DP,  1.00_DP, -1.00_DP/)
  dodecahedron_volume_pts(1:3, 7) = (/-1.00_DP,  1.00_DP,  1.00_DP/)
  dodecahedron_volume_pts(1:3, 8) = (/-0.25_DP, -1.00_DP,  1.00_DP/)

  ! Allocate Dodecahedron storage and construct to use vertices dodecahedron_volume_pts
  call new(dodecahedron_volume)
  call construct(dodecahedron_volume, dodecahedron_volume_pts)

  ! Allocate AccumulatedVolumeMoments_SeparatedVolumeMoments which
  ! will store the moments for each localizer in the LocalizedSeparatorLink
  ! network.
  call new(phase_moments_LocalizedSeparatorLink)

  ! With linking setup, now perform all the calculation of
  ! localized SeparatedMoments<VolumeMoments> at once. To do this,
  ! only one link in the network needs to be passed. The network
  ! will then be traversed starting from that point.
  ! The moments for a specific PlanarLocalizer region can be obtained
  ! using its unique Id with the function getSepVMAtTag.
  call getNormMoments(dodecahedron_volume, localized_separator_link(-1,-1,-1), &
          phase_moments_LocalizedSeparatorLink)

  write(*,'(A)')
  write(*,'(A)') 'Comparison between Expected and Computed Moments'
  write(*,'(A)') '================================================================'
  write(*,'(A)') ' '
  write(*,'(A)')       'For Localizer(0,-1,-1) '
  write(*,'(A)')       'Internal volume          '
  write(*,'(A,F10.5)') '   Expected:   ', 11.0_DP/64.0_DP
  call getSepVMAtTag(phase_moments_LocalizedSeparatorLink,getId(localized_separator_link(0,-1,-1)),separated_volume_moments)
  write(*,'(A,F10.5)') '   Computed:   ', &
    getVolume(separated_volume_moments,0)
  write(*,'(A)')        'Internal volume centroid'
  write(*,'(A,3F10.5)') '   Expected:  ', (/0.0_DP, -96.0_DP/(11.0_DP*12.0_DP), -0.75_DP/)
  write(*,'(A,3F10.5)') '   Computed:  ', &
    getCentroid(separated_volume_moments,0)
  write(*,'(A)') ' '
  write(*,'(A)')        'For Localizer(0,0,0) '
  write(*,'(A)')        'Internal volume '
  write(*,'(A,F10.5)')  '   Expected:  ', 5.0_DP/32.0_DP + 1.0_DP/3.0_DP
  call getSepVMAtTag(phase_moments_LocalizedSeparatorLink,getId(localized_separator_link(0,0,0)),separated_volume_moments)
  write(*,'(A,F10.5)')  '   Computed:  ', &
    getVolume(separated_volume_moments,0)
  write(*,'(A)')        'Internal volume centroid'
  write(*,'(A,3F10.5)') '   Expected:  ', (/0.0_DP, -104.0_DP/423.0_DP, 0.0_DP/)
  write(*,'(A,3F10.5)') '   Computed:  ', &
    getCentroid(separated_volume_moments,0)

  write(*,'(A)')        'External volume         '
  write(*,'(A,F10.5)')  '   Expected:  ', 0.5_DP
  call getSepVMAtTag(phase_moments_LocalizedSeparatorLink,getId(localized_separator_link(0,0,0)),separated_volume_moments)
  write(*,'(A,F10.5)')  '   Computed:  ', &
    getVolume(separated_volume_moments,1)
  write(*,'(A)')        'External volume centroid'
  write(*,'(A,3F10.5)') '   Expected:  ', (/0.0_DP, 0.25_DP, 0.0_DP/)
  write(*,'(A,3F10.5)') '   Computed:  ', &
    getCentroid(separated_volume_moments,1)

contains

  ! Set the PlanarLocalizers to mimick the 3x3x3 uniform Cartesian mesh.
  ! Each PlanarLocalizer object will represent a single cell of
  ! dx = 1.0, with the 3x3x3 region centered at (0.0, 0.0, 0.0)
  subroutine makeCubicPlanarLocalizer()
    implicit none
    integer :: ii,jj,kk
    ! use of planar_localizer variables from main program

    real(DP), dimension(1:3) :: lower_pt ! Lower point of hex we're representing
    real(DP), dimension(1:3) :: upper_pt ! Upper point of hex we're representing

    do kk = -1, 1
       do jj = -1, 1
          do ii = -1, 1          
             ! Cell centers would be at i,j,k.
             lower_pt = real((/ii,jj,kk/),DP)
             upper_pt = real((/ii,jj,kk/),DP)
             
             ! Shift by 0.5*dx = 0.5 to get cell vertex positions.
             lower_pt = lower_pt - (/0.5_DP,0.5_DP,0.5_DP/)
             upper_pt = upper_pt + (/0.5_DP,0.5_DP,0.5_DP/)
             
             ! Add the six planes that represent the cube
             call addPlane(planar_localizer(ii,jj,kk),(/-1.0_DP, 0.0_DP, 0.0_DP/),-lower_pt(1)) ! x- face
             call addPlane(planar_localizer(ii,jj,kk),(/ 1.0_DP, 0.0_DP, 0.0_DP/), upper_pt(1)) ! x+ face
             
             call addPlane(planar_localizer(ii,jj,kk),(/ 0.0_DP,-1.0_DP, 0.0_DP/),-lower_pt(2)) ! y- face
             call addPlane(planar_localizer(ii,jj,kk),(/ 0.0_DP, 1.0_DP, 0.0_DP/), upper_pt(2)) ! y+ face
             
             call addPlane(planar_localizer(ii,jj,kk),(/ 0.0_DP, 0.0_DP,-1.0_DP/),-lower_pt(3)) ! z- face
             call addPlane(planar_localizer(ii,jj,kk),(/ 0.0_DP, 0.0_DP, 1.0_DP/), upper_pt(3)) ! z+ face
          end do
       end do
    end do
             
  end subroutine makeCubicPlanarLocalizer

  ! Setup the edges between the LocalizedSeparatorLink objects.
  ! Here, the plane order used in PlanarLocalizer needs to be respected,
  ! where the PlanarLocalizer objects were setup in the
  ! makeCubicPlanarLocalizer to be ordered as x-, x+, y-, y+, z-, z+.
  ! Neighbors therefore neeed to be specified in this way, so the
  ! 0th neighbor will be the x-1 neighbor. If this neighbor does not
  ! exist (e.g. for i = 0), a nullptr is supplied instead, which
  ! is used to indicate a lack of connectivity and termination of the graph.
  ! A unique id value must also be given to each LocalizedSeparatorLink,
  ! which is used to tag the moments returned during getNormalizedVolumeMoments
  ! with the respected region they came from.  
  subroutine setupLinking(localized_separator_link)
    implicit none
    type(LocSepLink_type), dimension(-1:1,-1:1,-1:1), intent(inout) :: localized_separator_link    
    integer :: ii,jj,kk

    do kk = -1,1
       do jj = -1,1
          do ii=-1,1
             unique_id = unique_id + 1
             call setId(localized_separator_link(ii,jj,kk),unique_id)    

             ! First face, x-
             if(ii-1 .lt. -1) then
                ! Boundary plane
                call setEdgeConnectivityNull(localized_separator_link(ii,jj,kk),0)
             else
                ! Plane has neighbor
                call setEdgeConnectivity(localized_separator_link(ii,jj,kk),0,localized_separator_link(ii-1,jj,kk))
             end if

             ! Second face, x+
             if(ii+1 .gt. 1) then
                ! Boundary plane
                call setEdgeConnectivityNull(localized_separator_link(ii,jj,kk),1)
             else
                ! Plane has neighbor
                call setEdgeConnectivity(localized_separator_link(ii,jj,kk),1,localized_separator_link(ii+1,jj,kk))
             end if

             ! Third face, y-
             if(jj-1 .lt. -1) then
                ! Boundary plane
                call setEdgeConnectivityNull(localized_separator_link(ii,jj,kk),2)
             else
                ! Plane has neighbor
                call setEdgeConnectivity(localized_separator_link(ii,jj,kk),2,localized_separator_link(ii,jj-1,kk))
             end if

             ! Fourth face, y+
             if(jj+1 .gt. 1) then
                ! Boundary plane
                call setEdgeConnectivityNull(localized_separator_link(ii,jj,kk),3)
             else
                ! Plane has neighbor
                call setEdgeConnectivity(localized_separator_link(ii,jj,kk),3,localized_separator_link(ii,jj+1,kk))
             end if

             ! Fifth face, z-
             if(kk-1 .lt. -1) then
                ! Boundary plane
                call setEdgeConnectivityNull(localized_separator_link(ii,jj,kk),4)
             else
                ! Plane has neighbor
                call setEdgeConnectivity(localized_separator_link(ii,jj,kk),4,localized_separator_link(ii,jj,kk-1))
             end if

             ! Sixth face, z+
             if(kk+1 .gt. 1) then
                ! Boundary plane
                call setEdgeConnectivityNull(localized_separator_link(ii,jj,kk),5)
             else
                ! Plane has neighbor
                call setEdgeConnectivity(localized_separator_link(ii,jj,kk),5,localized_separator_link(ii,jj,kk+1))
             end if
          end do
       end do
    end do

  end subroutine setupLinking

end program main
