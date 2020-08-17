!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

! This example demonstrates the ability
! to create a polygon with N-vertices and then
! return individual triangles from its surface.
! The integration of surface moments for localized
! regions is also demonstrated.
program main
  use irl_fortran_interface
  implicit none
  integer, parameter :: DP = kind(1.0d0)

  type(Poly_type) :: polygon
  type(Tri_type) :: triangle
  type(PlanarLoc_type) :: planar_localizer
  real(DP), dimension(3,5) :: polygon_pts

  integer :: n
  real(DP), dimension(4) :: existence_plane
  real(DP) :: full_area, localized_area, tmp_area
  real(DP) :: triangle_area, localized_triangle_area

  ! Initialize the different polygon objects
  call new(polygon)
  call new(triangle)


  ! Set the polygon to be a pentagon in the XZ Plane
  ! NOTE: There is an assumed connectivity between
  ! subsequent vertices and between the last vertex
  ! and the first one
  polygon_pts(:,1) = (/0.0_DP, 0.0_DP, 0.0_DP/)
  polygon_pts(:,2) = (/0.0_DP, 0.0_DP, 1.0_DP/)
  polygon_pts(:,3) = (/0.5_DP, 0.0_DP, 1.5_DP/)
  polygon_pts(:,4) = (/1.0_DP, 0.0_DP, 1.0_DP/)
  polygon_pts(:,5) = (/1.0_DP, 0.0_DP, 0.0_DP/)
  call construct(polygon, 5, polygon_pts)

  ! Calculate the plane the polygon exists on
  ! Following the right-hand-rule, this should be
  ! n = (0.0, 1.0, 0.0) and d = 0.0
  call calculateAndSetPlaneOfExistence(polygon)
  existence_plane =  getPlaneOfExistence(polygon)

  ! Surface moments for this polygon can be computed
  ! directly or the polygon can be used in
  ! getVolumeMoments to first perform
  ! truncation by a set of planes.
  full_area = calculateVolume(polygon)

  ! Localizer that removes all but the top part of the pentagon.
  call new(planar_localizer)
  call addPlane(planar_localizer, (/0.0_DP, 0.0_DP, -1.0_DP/), -1.0_DP)
  call getNormMoments(polygon, planar_localizer, localized_area)


  ! Individual triangles can also be obtained from
  ! polygon as well.
  triangle_area = 0.0_DP
  do n = 1, getNumberOfSimplicesInDecomposition(polygon)
    call getSimplexFromDecomposition(polygon, n-1, triangle)
    tmp_area = calculateVolume(triangle)
    triangle_area = triangle_area + tmp_area
  end do

  ! These triangles can also be localized using the PlanarLocalizer object
  localized_triangle_area = 0.0_DP
  do n = 1, getNumberOfSimplicesInDecomposition(polygon)
    call getSimplexFromDecomposition(polygon, n-1, triangle)
    call getNormMoments(triangle, planar_localizer, tmp_area)
    localized_triangle_area = localized_triangle_area + tmp_area
  end do

  write(*,'(A)')
  write(*,'(A)') '          Polygon usage demonstration           '
  write(*,'(A)') '================================================'
  write(*,'(A, 4F10.5)') ' Expected existence plane   : ', (/0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP/)
  write(*,'(A, 4F10.5)') ' Calculated existence plane : ', existence_plane(1:4)
  write(*,'(A)')
  write(*,'(A, F10.5)') ' Correct surface area    :  ', 1.0_DP + 2.0_DP / 8.0_DP
  write(*,'(A, F10.5)') ' Calculated surface area :  ', full_area
  write(*,'(A)')
  write(*,'(A, F10.5)') ' Correct surface area                   :  ', 1.0_DP + 2.0_DP / 8.0_DP
  write(*,'(A, F10.5)') ' Calculated surface area from triangles :  ', triangle_area
  write(*,'(A)')  
  write(*,'(A, F10.5)') ' Correct localized surface area              :  ', 2.0_DP / 8.0_DP
  write(*,'(A, F10.5)') ' Calculated localized surface area :  ', localized_area
  write(*,'(A)')
  write(*,'(A, F10.5)') ' Calculated localized surface area from triangles :  ', localized_triangle_area  

end program main

