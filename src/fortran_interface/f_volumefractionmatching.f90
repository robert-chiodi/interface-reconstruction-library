!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_volumefractionmatching.f90
!!
!! This file deals with setting the distances to
!! each plane in a planar reconstruction to match a given
!! volume fraction for the provided cell.

!> This module contains mappings to the
!! IRL C interface that deals with setting the distance
!! to each plane in a reconstruction to recreate
!! the volume fraction on the provided polyhedron.
module f_volumefractionmatching
  use f_DefinedTypes
  use f_RectCub_class
  use f_Hex_class
  use f_Tet_class
  use f_TriPrism_class
  use f_Pyrmd_class
  use f_PlanarSep_class
  use f_PlanarSepPathGroup_class
  implicit none

  interface matchVolumeFraction
    ! Set distances in PlanarSep to recreate volume fraction on a Rectangular Cuboid
    module procedure matchVolumeFraction_RectCub_PlanarSep
    ! Set distances in PlanarSep to recreate volume fraction on a Rectangular Cuboid
    ! using default tolerance value
    module procedure matchVolumeFraction_RectCub_PlanarSep_Default
    ! Set distances in PlanarSep to recreate volume fraction on a Hexahedron
    module procedure matchVolumeFraction_Hex_PlanarSep
    ! Set distances in PlanarSep to recreate volume fraction on a Hexahedron
    ! using default tolerance value
    module procedure matchVolumeFraction_Hex_PlanarSep_Default  
    ! Set distances in PlanarSep to recreate volume fraction on a Tet
    module procedure matchVolumeFraction_Tet_PlanarSep
    ! Set distances in PlanarSep to recreate volume fraction on a Tet
    ! using default tolerance value
    module procedure matchVolumeFraction_Tet_PlanarSep_Default    
    ! Set distances in PlanarSep to recreate volume fraction on a TriangularPrism
    module procedure matchVolumeFraction_TriPrism_PlanarSep
    ! Set distances in PlanarSep to recreate volume fraction on a TriangularPrism
    ! using default tolerance value
    module procedure matchVolumeFraction_TriPrism_PlanarSep_Default    
    ! Set distances in PlanarSep to recreate volume fraction on a Pyramid
    module procedure matchVolumeFraction_Pyrmd_PlanarSep
    ! Set distances in PlanarSep to recreate volume fraction on a Pyramid
    ! using default tolerance value
    module procedure matchVolumeFraction_Pyrmd_PlanarSep_Default   
  end interface matchVolumeFraction

  interface matchGroupVolumeFraction
    module procedure matchGroupVolumeFraction_Tet_PlanarSepPathGroup
    module procedure matchGroupVolumeFraction_Tet_PlanarSepPathGroup_Default    
    module procedure matchGroupVolumeFraction_Pyrmd_PlanarSepPathGroup
    module procedure matchGroupVolumeFraction_Pyrmd_PlanarSepPathGroup_Default
    module procedure matchGroupVolumeFraction_TriPrism_PlanarSepPathGroup
    module procedure matchGroupVolumeFraction_TriPrism_PlanarSepPathGroup_Default
    module procedure matchGroupVolumeFraction_Hex_PlanarSepPathGroup
    module procedure matchGroupVolumeFraction_Hex_PlanarSepPathGroup_Default        
  end interface

  interface
    subroutine F_matchVolumeFraction_RectCub_PlanarSep_DefTol(a_rectangular_cuboid, a_volume_fraction, &
         a_planar_separator) &
    bind(C, name="c_matchVolumeFraction_RectCub_PlanarSep_Default")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_RectCub) :: a_rectangular_cuboid ! Pointer to RectCub object
      real(C_DOUBLE) :: a_volume_fraction ! Volume fraction to recreate
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
    end subroutine F_matchVolumeFraction_RectCub_PlanarSep_DefTol
  end interface

  interface
    subroutine F_matchVolumeFraction_RectCub_PlanarSep(a_rectangular_cuboid, a_volume_fraction, a_planar_separator, &
         a_volume_fraction_tolerance) &
    bind(C, name="c_matchVolumeFraction_RectCub_PlanarSep")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_RectCub) :: a_rectangular_cuboid ! Pointer to RectCub object
      real(C_DOUBLE) :: a_volume_fraction ! Volume fraction to recreate
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
      real(C_DOUBLE) :: a_volume_fraction_tolerance ! Tolerance to recreate volume fraction within
    end subroutine F_matchVolumeFraction_RectCub_PlanarSep
  end interface

  interface
    subroutine F_matchVolumeFraction_Hex_PlanarSep_DefTol(a_hexahedron, a_volume_fraction, &
         a_planar_separator) &
    bind(C, name="c_matchVolumeFraction_Hex_PlanarSep_Default")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Hex) :: a_hexahedron ! Pointer to a_hexahedron object
      real(C_DOUBLE) :: a_volume_fraction ! Volume fraction to recreate
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
    end subroutine F_matchVolumeFraction_Hex_PlanarSep_DefTol
  end interface

  interface
    subroutine F_matchVolumeFraction_Tet_PlanarSep(a_tet, a_volume_fraction, a_planar_separator, &
         a_volume_fraction_tolerance) &
    bind(C, name="c_matchVolumeFraction_Tet_PlanarSep")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Tet) :: a_tet ! Pointer to a_hexahedron object
      real(C_DOUBLE) :: a_volume_fraction ! Volume fraction to recreate
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
      real(C_DOUBLE) :: a_volume_fraction_tolerance ! Tolerance to recreate volume fraction within
    end subroutine F_matchVolumeFraction_Tet_PlanarSep
  end interface
  
  interface
    subroutine F_matchVolumeFraction_Tet_PlanarSep_DefTol(a_tet, a_volume_fraction, &
         a_planar_separator) &
    bind(C, name="c_matchVolumeFraction_Tet_PlanarSep_Default")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Tet) :: a_tet ! Pointer to a_hexahedron object
      real(C_DOUBLE) :: a_volume_fraction ! Volume fraction to recreate
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
    end subroutine F_matchVolumeFraction_Tet_PlanarSep_DefTol
  end interface

  interface
    subroutine F_matchVolumeFraction_Hex_PlanarSep(a_hexahedron, a_volume_fraction, a_planar_separator, &
         a_volume_fraction_tolerance) &
    bind(C, name="c_matchVolumeFraction_Hex_PlanarSep")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Hex) :: a_hexahedron ! Pointer to a_hexahedron object
      real(C_DOUBLE) :: a_volume_fraction ! Volume fraction to recreate
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
      real(C_DOUBLE) :: a_volume_fraction_tolerance ! Tolerance to recreate volume fraction within
    end subroutine F_matchVolumeFraction_Hex_PlanarSep
  end interface
  
  interface
    subroutine F_matchVolumeFraction_TriPrism_PlanarSep_DefTol(a_triangular_prism, a_volume_fraction, &
         a_planar_separator) &
    bind(C, name="c_matchVolumeFraction_TriPrism_PlanarSep_Default")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_TriPrism) :: a_triangular_prism ! Pointer to a_triangular_prism object
      real(C_DOUBLE) :: a_volume_fraction ! Volume fraction to recreate
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
    end subroutine F_matchVolumeFraction_TriPrism_PlanarSep_DefTol
  end interface

  interface
    subroutine F_matchVolumeFraction_TriPrism_PlanarSep(a_triangular_prism, a_volume_fraction, &
         a_planar_separator, a_volume_fraction_tolerance) &
    bind(C, name="c_matchVolumeFraction_TriPrism_PlanarSep")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_TriPrism) :: a_triangular_prism ! Pointer to a_triangular_prism object
      real(C_DOUBLE) :: a_volume_fraction ! Volume fraction to recreate
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
      real(C_DOUBLE) :: a_volume_fraction_tolerance ! Tolerance to recreate volume fraction within
    end subroutine F_matchVolumeFraction_TriPrism_PlanarSep
  end interface
  
  interface
    subroutine F_matchVolumeFraction_Pyrmd_PlanarSep_DefTol(a_pyramid, a_volume_fraction, &
         a_planar_separator) &
    bind(C, name="c_matchVolumeFraction_Pyrmd_PlanarSep_Default")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Pyrmd) :: a_pyramid ! Pointer to a_pyramid object
      real(C_DOUBLE) :: a_volume_fraction ! Volume fraction to recreate
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
    end subroutine F_matchVolumeFraction_Pyrmd_PlanarSep_DefTol
  end interface

  interface
    subroutine F_matchVolumeFraction_Pyrmd_PlanarSep(a_pyramid, a_volume_fraction, &
         a_planar_separator, a_volume_fraction_tolerance) &
    bind(C, name="c_matchVolumeFraction_Pyrmd_PlanarSep")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Pyrmd) :: a_pyramid ! Pointer to a_pyramid object
      real(C_DOUBLE) :: a_volume_fraction ! Volume fraction to recreate
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
      real(C_DOUBLE) :: a_volume_fraction_tolerance ! Tolerance to recreate volume fraction within
    end subroutine F_matchVolumeFraction_Pyrmd_PlanarSep
  end interface

  interface
    subroutine F_matchGroupVolumeFraction_Tet_PlanarSepPathGroup_DefTol(a_cell, a_size, a_volume_fraction, &
                                                                        a_planar_separator) &
    bind(C, name="c_matchGroupVolumeFraction_Tet_PlanarSepPathGroup_Default")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Tet) :: a_cell
      integer(C_INT) :: a_size
      real(C_DOUBLE), dimension(*) :: a_volume_fraction ! dimension(1:a_size)
      type(c_PlanarSepPathGroup) :: a_planar_separator ! Pointer to PlanarSep object
    end subroutine F_matchGroupVolumeFraction_Tet_PlanarSepPathGroup_DefTol
  end interface

  interface
    subroutine F_matchGroupVolumeFraction_Tet_PlanarSepPathGroup(a_cell, a_size, a_volume_fraction, &
                                                                        a_planar_separator, a_volume_fraction_tolerance) &
    bind(C, name="c_matchGroupVolumeFraction_Tet_PlanarSepPathGroup")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Tet) :: a_cell
      integer(C_INT) :: a_size
      real(C_DOUBLE), dimension(*) :: a_volume_fraction ! dimension(1:a_size)
      type(c_PlanarSepPathGroup) :: a_planar_separator ! Pointer to PlanarSep object
      real(C_DOUBLE) :: a_volume_fraction_tolerance ! Tolerance to recreate volume fraction within      
    end subroutine F_matchGroupVolumeFraction_Tet_PlanarSepPathGroup
  end interface

  interface
    subroutine F_matchGroupVolumeFraction_Pyrmd_PlanarSepPathGroup_DefTol(a_cell, a_size, a_volume_fraction, &
                                                                        a_planar_separator) &
    bind(C, name="c_matchGroupVolumeFraction_Pyrmd_PlanarSepPathGroup_Default")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Pyrmd) :: a_cell
      integer(C_INT) :: a_size
      real(C_DOUBLE), dimension(*) :: a_volume_fraction ! dimension(1:a_size)
      type(c_PlanarSepPathGroup) :: a_planar_separator ! Pointer to PlanarSep object
    end subroutine F_matchGroupVolumeFraction_Pyrmd_PlanarSepPathGroup_DefTol
  end interface

  interface
    subroutine F_matchGroupVolumeFraction_Pyrmd_PlanarSepPathGroup(a_cell, a_size, a_volume_fraction, &
                                                                        a_planar_separator, a_volume_fraction_tolerance) &
    bind(C, name="c_matchGroupVolumeFraction_Pyrmd_PlanarSepPathGroup")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Pyrmd) :: a_cell
      integer(C_INT) :: a_size
      real(C_DOUBLE), dimension(*) :: a_volume_fraction ! dimension(1:a_size)
      type(c_PlanarSepPathGroup) :: a_planar_separator ! Pointer to PlanarSep object
      real(C_DOUBLE) :: a_volume_fraction_tolerance ! Tolerance to recreate volume fraction within      
    end subroutine F_matchGroupVolumeFraction_Pyrmd_PlanarSepPathGroup
  end interface  

  interface
    subroutine F_matchGroupVolumeFraction_TriPrism_PlanarSepPathGroup_DefTol(a_cell, a_size, a_volume_fraction, &
                                                                        a_planar_separator) &
    bind(C, name="c_matchGroupVolumeFraction_TriPrism_PlanarSepPathGroup_Default")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_TriPrism) :: a_cell
      integer(C_INT) :: a_size
      real(C_DOUBLE), dimension(*) :: a_volume_fraction ! dimension(1:a_size)
      type(c_PlanarSepPathGroup) :: a_planar_separator ! Pointer to PlanarSep object
    end subroutine F_matchGroupVolumeFraction_TriPrism_PlanarSepPathGroup_DefTol
  end interface

  interface
    subroutine F_matchGroupVolumeFraction_TriPrism_PlanarSepPathGroup(a_cell, a_size, a_volume_fraction, &
                                                                        a_planar_separator, a_volume_fraction_tolerance) &
    bind(C, name="c_matchGroupVolumeFraction_TriPrism_PlanarSepPathGroup")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_TriPrism) :: a_cell
      integer(C_INT) :: a_size
      real(C_DOUBLE), dimension(*) :: a_volume_fraction ! dimension(1:a_size)
      type(c_PlanarSepPathGroup) :: a_planar_separator ! Pointer to PlanarSep object
      real(C_DOUBLE) :: a_volume_fraction_tolerance ! Tolerance to recreate volume fraction within      
    end subroutine F_matchGroupVolumeFraction_TriPrism_PlanarSepPathGroup
  end interface  

  interface
    subroutine F_matchGroupVolumeFraction_Hex_PlanarSepPathGroup_DefTol(a_cell, a_size, a_volume_fraction, &
                                                                        a_planar_separator) &
    bind(C, name="c_matchGroupVolumeFraction_Hex_PlanarSepPathGroup_Default")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Hex) :: a_cell
      integer(C_INT) :: a_size
      real(C_DOUBLE), dimension(*) :: a_volume_fraction ! dimension(1:a_size)
      type(c_PlanarSepPathGroup) :: a_planar_separator ! Pointer to PlanarSep object
    end subroutine F_matchGroupVolumeFraction_Hex_PlanarSepPathGroup_DefTol
  end interface

  interface
    subroutine F_matchGroupVolumeFraction_Hex_PlanarSepPathGroup(a_cell, a_size, a_volume_fraction, &
                                                                        a_planar_separator, a_volume_fraction_tolerance) &
    bind(C, name="c_matchGroupVolumeFraction_Hex_PlanarSepPathGroup")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Hex) :: a_cell
      integer(C_INT) :: a_size
      real(C_DOUBLE), dimension(*) :: a_volume_fraction ! dimension(1:a_size)
      type(c_PlanarSepPathGroup) :: a_planar_separator ! Pointer to PlanarSep object
      real(C_DOUBLE) :: a_volume_fraction_tolerance ! Tolerance to recreate volume fraction within      
    end subroutine F_matchGroupVolumeFraction_Hex_PlanarSepPathGroup
  end interface  

contains

  subroutine matchVolumeFraction_RectCub_PlanarSep_Default(a_rectangular_cuboid, a_volume_fraction, &
       a_planar_separator)
    implicit none
    type(RectCub_type) :: a_rectangular_cuboid ! Pointer to RectCub object
    real(IRL_double) :: a_volume_fraction ! Volume fraction to recreate
    type(PlanarSep_type) :: a_planar_separator ! Pointer to PlanarSep object

    call F_matchVolumeFraction_RectCub_PlanarSep_DefTol(a_rectangular_cuboid%c_object, a_volume_fraction, &
          a_planar_separator%c_object)

  end subroutine matchVolumeFraction_RectCub_PlanarSep_Default

  subroutine matchVolumeFraction_RectCub_PlanarSep(a_rectangular_cuboid, a_volume_fraction, &
       a_planar_separator, a_volume_fraction_tolerance)
    implicit none
    type(RectCub_type) :: a_rectangular_cuboid ! Pointer to RectCub object
    real(IRL_double) :: a_volume_fraction ! Volume fraction to recreate
    type(PlanarSep_type) :: a_planar_separator ! Pointer to PlanarSep object
    real(IRL_double) :: a_volume_fraction_tolerance ! Tolerance to recreate volume fraction within

    call F_matchVolumeFraction_RectCub_PlanarSep(a_rectangular_cuboid%c_object, a_volume_fraction, &
          a_planar_separator%c_object, a_volume_fraction_tolerance)

  end subroutine matchVolumeFraction_RectCub_PlanarSep

  subroutine matchVolumeFraction_Hex_PlanarSep_Default(a_hexahedron, a_volume_fraction, a_planar_separator)
    implicit none
    type(Hex_type) :: a_hexahedron ! Pointer to Hex object
    real(IRL_double) :: a_volume_fraction ! Volume fraction to recreate
    type(PlanarSep_type) :: a_planar_separator ! Pointer to PlanarSep object

    call F_matchVolumeFraction_Hex_PlanarSep_DefTol(a_hexahedron%c_object, a_volume_fraction, &
          a_planar_separator%c_object)

  end subroutine matchVolumeFraction_Hex_PlanarSep_Default

  subroutine matchVolumeFraction_Hex_PlanarSep(a_hexahedron, a_volume_fraction, &
       a_planar_separator, a_volume_fraction_tolerance)
    implicit none
    type(Hex_type) :: a_hexahedron ! Pointer to Hex object
    real(IRL_double) :: a_volume_fraction ! Volume fraction to recreate
    type(PlanarSep_type) :: a_planar_separator ! Pointer to PlanarSep object
    real(IRL_double) :: a_volume_fraction_tolerance ! Tolerance to recreate volume fraction within

    call F_matchVolumeFraction_Hex_PlanarSep(a_hexahedron%c_object, a_volume_fraction, &
          a_planar_separator%c_object, a_volume_fraction_tolerance)

  end subroutine matchVolumeFraction_Hex_PlanarSep
  
  subroutine matchVolumeFraction_Tet_PlanarSep_Default(a_tet, a_volume_fraction, &
       a_planar_separator)
    implicit none
    type(Tet_type) :: a_tet ! Pointer to Tet object
    real(IRL_double) :: a_volume_fraction ! Volume fraction to recreate
    type(PlanarSep_type) :: a_planar_separator ! Pointer to PlanarSep object

    call F_matchVolumeFraction_Tet_PlanarSep_DefTol(a_tet%c_object, a_volume_fraction, &
          a_planar_separator%c_object)

  end subroutine matchVolumeFraction_Tet_PlanarSep_Default

  subroutine matchVolumeFraction_Tet_PlanarSep(a_tet, a_volume_fraction, &
       a_planar_separator, a_volume_fraction_tolerance)
    implicit none
    type(Tet_type) :: a_tet ! Pointer to Tet object
    real(IRL_double) :: a_volume_fraction ! Volume fraction to recreate
    type(PlanarSep_type) :: a_planar_separator ! Pointer to PlanarSep object
    real(IRL_double) :: a_volume_fraction_tolerance ! Tolerance to recreate volume fraction within

    call F_matchVolumeFraction_Tet_PlanarSep(a_tet%c_object, a_volume_fraction, &
          a_planar_separator%c_object, a_volume_fraction_tolerance)

  end subroutine matchVolumeFraction_Tet_PlanarSep
  
  subroutine matchVolumeFraction_TriPrism_PlanarSep_Default(a_triangular_prism, a_volume_fraction, &
       a_planar_separator)
    implicit none
    type(TriPrism_type) :: a_triangular_prism ! Pointer to TriPrism object
    real(IRL_double) :: a_volume_fraction ! Volume fraction to recreate
    type(PlanarSep_type) :: a_planar_separator ! Pointer to PlanarSep object

    call F_matchVolumeFraction_TriPrism_PlanarSep_DefTol(a_triangular_prism%c_object, a_volume_fraction, &
          a_planar_separator%c_object)

  end subroutine matchVolumeFraction_TriPrism_PlanarSep_Default

  subroutine matchVolumeFraction_TriPrism_PlanarSep(a_triangular_prism, a_volume_fraction, &
       a_planar_separator, a_volume_fraction_tolerance)
    implicit none
    type(TriPrism_type) :: a_triangular_prism ! Pointer to TriPrism object
    real(IRL_double) :: a_volume_fraction ! Volume fraction to recreate
    type(PlanarSep_type) :: a_planar_separator ! Pointer to PlanarSep object
    real(IRL_double) :: a_volume_fraction_tolerance ! Tolerance to recreate volume fraction within

    call F_matchVolumeFraction_TriPrism_PlanarSep(a_triangular_prism%c_object, a_volume_fraction, &
          a_planar_separator%c_object, a_volume_fraction_tolerance)

  end subroutine matchVolumeFraction_TriPrism_PlanarSep
  
  subroutine matchVolumeFraction_Pyrmd_PlanarSep_Default(a_pyramid, a_volume_fraction, &
       a_planar_separator)
    implicit none
    type(Pyrmd_type) :: a_pyramid ! Pointer to a_pyramid object
    real(IRL_double) :: a_volume_fraction ! Volume fraction to recreate
    type(PlanarSep_type) :: a_planar_separator ! Pointer to PlanarSep object

    call F_matchVolumeFraction_Pyrmd_PlanarSep_DefTol(a_pyramid%c_object, a_volume_fraction, &
          a_planar_separator%c_object)

  end subroutine matchVolumeFraction_Pyrmd_PlanarSep_Default

  subroutine matchVolumeFraction_Pyrmd_PlanarSep(a_pyramid, a_volume_fraction, &
       a_planar_separator, a_volume_fraction_tolerance)
    implicit none
    type(Pyrmd_type) :: a_pyramid ! Pointer to a_pyramid object
    real(IRL_double) :: a_volume_fraction ! Volume fraction to recreate
    type(PlanarSep_type) :: a_planar_separator ! Pointer to PlanarSep object
    real(IRL_double) :: a_volume_fraction_tolerance ! Tolerance to recreate volume fraction within

    call F_matchVolumeFraction_Pyrmd_PlanarSep(a_pyramid%c_object, a_volume_fraction, &
          a_planar_separator%c_object, a_volume_fraction_tolerance)

  end subroutine matchVolumeFraction_Pyrmd_PlanarSep

  subroutine matchGroupVolumeFraction_Tet_PlanarSepPathGroup_Default(a_cell, a_size, a_volume_fraction, &
       a_planar_separator)
    implicit none
    type(Tet_type) :: a_cell
    integer, intent(in) :: a_size
    real(IRL_double), dimension(:), intent(in) :: a_volume_fraction ! Volume fraction to recreate
    type(PlanarSepPathGroup_type) :: a_planar_separator ! Pointer to PlanarSepPathGroup object

    call F_matchGroupVolumeFraction_Tet_PlanarSepPathGroup_DefTol(a_cell%c_object, a_size, a_volume_fraction, &
          a_planar_separator%c_object)

  end subroutine matchGroupVolumeFraction_Tet_PlanarSepPathGroup_Default

  subroutine matchGroupVolumeFraction_Tet_PlanarSepPathGroup(a_cell, a_size, a_volume_fraction, a_planar_separator, &
       a_volume_fraction_tolerance)
    implicit none
    type(Tet_type) :: a_cell
    integer, intent(in) :: a_size
    real(IRL_double), dimension(:), intent(in) :: a_volume_fraction ! Volume fraction to recreate
    type(PlanarSepPathGroup_type) :: a_planar_separator ! Pointer to PlanarSepPathGroup object
    real(IRL_double) :: a_volume_fraction_tolerance ! Tolerance to recreate volume fraction within    

    call F_matchGroupVolumeFraction_Tet_PlanarSepPathGroup(a_cell%c_object, a_size, a_volume_fraction, &
          a_planar_separator%c_object, a_volume_fraction_tolerance)

  end subroutine matchGroupVolumeFraction_Tet_PlanarSepPathGroup

  subroutine matchGroupVolumeFraction_Pyrmd_PlanarSepPathGroup_Default(a_cell, a_size, a_volume_fraction, &
       a_planar_separator)
    implicit none
    type(Pyrmd_type) :: a_cell
    integer, intent(in) :: a_size
    real(IRL_double), dimension(:), intent(in) :: a_volume_fraction ! Volume fraction to recreate
    type(PlanarSepPathGroup_type) :: a_planar_separator ! Pointer to PlanarSepPathGroup object

    call F_matchGroupVolumeFraction_Pyrmd_PlanarSepPathGroup_DefTol(a_cell%c_object, a_size, a_volume_fraction, &
          a_planar_separator%c_object)

  end subroutine matchGroupVolumeFraction_Pyrmd_PlanarSepPathGroup_Default

  subroutine matchGroupVolumeFraction_Pyrmd_PlanarSepPathGroup(a_cell, a_size, a_volume_fraction, a_planar_separator, &
       a_volume_fraction_tolerance)
    implicit none
    type(Pyrmd_type) :: a_cell
    integer, intent(in) :: a_size
    real(IRL_double), dimension(:), intent(in) :: a_volume_fraction ! Volume fraction to recreate
    type(PlanarSepPathGroup_type) :: a_planar_separator ! Pointer to PlanarSepPathGroup object
    real(IRL_double) :: a_volume_fraction_tolerance ! Tolerance to recreate volume fraction within    

    call F_matchGroupVolumeFraction_Pyrmd_PlanarSepPathGroup(a_cell%c_object, a_size, a_volume_fraction, &
          a_planar_separator%c_object, a_volume_fraction_tolerance)

  end subroutine matchGroupVolumeFraction_Pyrmd_PlanarSepPathGroup

  subroutine matchGroupVolumeFraction_TriPrism_PlanarSepPathGroup_Default(a_cell, a_size, a_volume_fraction, &
       a_planar_separator)
    implicit none
    type(TriPrism_type) :: a_cell
    integer, intent(in) :: a_size
    real(IRL_double), dimension(:), intent(in) :: a_volume_fraction ! Volume fraction to recreate
    type(PlanarSepPathGroup_type) :: a_planar_separator ! Pointer to PlanarSepPathGroup object

    call F_matchGroupVolumeFraction_TriPrism_PlanarSepPathGroup_DefTol(a_cell%c_object, a_size, a_volume_fraction, &
          a_planar_separator%c_object)

  end subroutine matchGroupVolumeFraction_TriPrism_PlanarSepPathGroup_Default

  subroutine matchGroupVolumeFraction_TriPrism_PlanarSepPathGroup(a_cell, a_size, a_volume_fraction, a_planar_separator, &
       a_volume_fraction_tolerance)
    implicit none
    type(TriPrism_type) :: a_cell
    integer, intent(in) :: a_size
    real(IRL_double), dimension(:), intent(in) :: a_volume_fraction ! Volume fraction to recreate
    type(PlanarSepPathGroup_type) :: a_planar_separator ! Pointer to PlanarSepPathGroup object
    real(IRL_double) :: a_volume_fraction_tolerance ! Tolerance to recreate volume fraction within    

    call F_matchGroupVolumeFraction_TriPrism_PlanarSepPathGroup(a_cell%c_object, a_size, a_volume_fraction, &
          a_planar_separator%c_object, a_volume_fraction_tolerance)

  end subroutine matchGroupVolumeFraction_TriPrism_PlanarSepPathGroup

  subroutine matchGroupVolumeFraction_Hex_PlanarSepPathGroup_Default(a_cell, a_size, a_volume_fraction, &
       a_planar_separator)
    implicit none
    type(Hex_type) :: a_cell
    integer, intent(in) :: a_size
    real(IRL_double), dimension(:), intent(in) :: a_volume_fraction ! Volume fraction to recreate
    type(PlanarSepPathGroup_type) :: a_planar_separator ! Pointer to PlanarSepPathGroup object

    call F_matchGroupVolumeFraction_Hex_PlanarSepPathGroup_DefTol(a_cell%c_object, a_size, a_volume_fraction, &
          a_planar_separator%c_object)

  end subroutine matchGroupVolumeFraction_Hex_PlanarSepPathGroup_Default

  subroutine matchGroupVolumeFraction_Hex_PlanarSepPathGroup(a_cell, a_size, a_volume_fraction, a_planar_separator, &
       a_volume_fraction_tolerance)
    implicit none
    type(Hex_type) :: a_cell
    integer, intent(in) :: a_size
    real(IRL_double), dimension(:), intent(in) :: a_volume_fraction ! Volume fraction to recreate
    type(PlanarSepPathGroup_type) :: a_planar_separator ! Pointer to PlanarSepPathGroup object
    real(IRL_double) :: a_volume_fraction_tolerance ! Tolerance to recreate volume fraction within    

    call F_matchGroupVolumeFraction_Hex_PlanarSepPathGroup(a_cell%c_object, a_size, a_volume_fraction, &
          a_planar_separator%c_object, a_volume_fraction_tolerance)

  end subroutine matchGroupVolumeFraction_Hex_PlanarSepPathGroup
  
  
end module f_volumefractionmatching
