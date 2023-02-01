!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_Paraboloidarator_class.f90
!!
!! This file allows use of the IRL Paraboloidarator
!! class through a fortran interface.

!> \brief A fortran type class that allows the creation of
!! IRL's Paraboloidarator class along with enabling
!! some of its methods.
module f_Paraboloid_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  use f_ObjServer_Paraboloid_class
  use f_Poly_class
  use f_PlanarSep_class
  use f_RectCub_class
  use f_TriangulatedParaboloid_class
  implicit none

  type, public, bind(C) :: c_Paraboloid
    type(C_PTR), private :: object = C_NULL_PTR
    logical(C_BOOL), private :: is_owning  = .false.
  end type c_Paraboloid

  type, public :: Paraboloid_type
    type(c_Paraboloid) :: c_object
  contains
    final :: Paraboloid_class_delete
  end type Paraboloid_type

  interface new
    module procedure Paraboloid_class_new
    module procedure Paraboloid_class_newFromObjectAllocationServer
  end interface
  interface setDatum
    module procedure Paraboloid_class_setDatum
  end interface
  interface setReferenceFrame
    module procedure Paraboloid_class_setReferenceFrame
  end interface
  interface setAlignedParaboloid
    module procedure Paraboloid_class_setAlignedParaboloid
  end interface
  interface setParaboloidFromPolygon
    module procedure Paraboloid_class_setParaboloidFromPolygon
  end interface
  interface setParaboloidFromPlanarSep
    module procedure Paraboloid_class_setParaboloidFromPlanarSep
  end interface
  interface setParaboloidJibben
    module procedure Paraboloid_class_setParaboloidJibben
  end interface
  interface setParaboloidFull
    module procedure Paraboloid_class_setParaboloidFull
  end interface
  interface setParaboloidEmpty
    module procedure Paraboloid_class_setParaboloidEmpty
  end interface
  interface copy
    module procedure Paraboloid_class_copy
  end interface
  interface getDatum
    module procedure Paraboloid_class_getDatum
  end interface
  interface getReferenceFrame
    module procedure Paraboloid_class_getReferenceFrame
  end interface
  interface getAlignedParaboloid
    module procedure Paraboloid_class_getAlignedParaboloid
  end interface
  interface triangulateInsideCuboid
    module procedure Paraboloid_class_triangulateInsideCuboid
  end interface
  interface isFlipped
    module procedure Paraboloid_class_isFlipped
  end interface
  interface printToScreen
    module procedure Paraboloid_class_printToScreen
  end interface


  interface

    subroutine F_Paraboloid_new(this) &
      bind(C, name="c_Paraboloid_new")
      import
      implicit none
      type(c_Paraboloid) :: this
    end subroutine F_Paraboloid_new

    subroutine F_Paraboloid_newFromObjectAllocationServer(this, a_object_allocation_server) &
      bind(C, name="c_Paraboloid_newFromObjectAllocationServer")
      import
      implicit none
      type(c_Paraboloid) :: this
      type(c_ObjServer_Paraboloid) :: a_object_allocation_server ! Pointer to ObjectAllocationServer<Paraboloid>
    end subroutine F_Paraboloid_newFromObjectAllocationServer

    subroutine F_Paraboloid_delete(this) &
      bind(C, name="c_Paraboloid_delete")
      import
      implicit none
      type(c_Paraboloid) :: this
    end subroutine F_Paraboloid_delete

    subroutine F_Paraboloid_setDatum(this, a_datum) &
      bind(C, name="c_Paraboloid_setDatum")
      import
      implicit none
      type(c_Paraboloid) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_datum !  dimension(1:3)
    end subroutine F_Paraboloid_setDatum

    subroutine F_Paraboloid_setReferenceFrame(this, a_normal1, a_normal2, a_normal3) &
      bind(C, name="c_Paraboloid_setReferenceFrame")
      import
      implicit none
      type(c_Paraboloid) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_normal1 !  dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(in) :: a_normal2 !  dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(in) :: a_normal3 !  dimension(1:3)
    end subroutine F_Paraboloid_setReferenceFrame

    subroutine F_Paraboloid_setAlignedParaboloid(this, a_coeff_a, a_coeff_b) &
      bind(C, name="c_Paraboloid_setAlignedParaboloid")
      import
      implicit none
      type(c_Paraboloid) :: this
      real(C_DOUBLE), intent(in) :: a_coeff_a 
      real(C_DOUBLE), intent(in) :: a_coeff_b 
    end subroutine F_Paraboloid_setAlignedParaboloid

    subroutine F_Paraboloid_setParaboloidFromPolygon(this, a_poly) &
      bind(C, name="c_Paraboloid_setParaboloidFromPolygon")
      import
      implicit none
      type(c_Paraboloid) :: this
      type(c_Poly), intent(in) :: a_poly
    end subroutine F_Paraboloid_setParaboloidFromPolygon

    subroutine F_Paraboloid_setParaboloidFromPlanarSep(this, a_plane, a_pt_ref) &
      bind(C, name="c_Paraboloid_setParaboloidFromPlanarSep")
      import
      implicit none
      type(c_Paraboloid) :: this
      type(c_PlanarSep), intent(in) :: a_plane
      real(C_DOUBLE), dimension(*), intent(in) :: a_pt_ref
    end subroutine F_Paraboloid_setParaboloidFromPlanarSep

    subroutine F_Paraboloid_setParaboloidJibben(this, a_plane, a_cuboid, a_npoly, a_vfrac, a_nvert, a_vert_coords) &
      bind(C, name="c_Paraboloid_setParaboloidJibben")
      import
      implicit none
      type(c_Paraboloid) :: this
      type(c_PlanarSep), intent(in) :: a_plane
      type(c_RectCub), intent(in)  :: a_cuboid
      integer(C_INT), intent(in) :: a_npoly
      real(C_DOUBLE), dimension(*), intent(in) :: a_vfrac
      integer(C_INT), dimension(*), intent(in) :: a_nvert
      real(C_DOUBLE), dimension(*), intent(in) :: a_vert_coords
    end subroutine F_Paraboloid_setParaboloidJibben

    subroutine F_Paraboloid_setParaboloidFull(this) &
      bind(C, name="c_Paraboloid_setParaboloidFull")
      import
      implicit none
      type(c_Paraboloid) :: this
    end subroutine F_Paraboloid_setParaboloidFull

    subroutine F_Paraboloid_setParaboloidEmpty(this) &
      bind(C, name="c_Paraboloid_setParaboloidEmpty")
      import
      implicit none
      type(c_Paraboloid) :: this
    end subroutine F_Paraboloid_setParaboloidEmpty

    subroutine F_Paraboloid_copy(this, a_other_Paraboloid) &
      bind(C, name="c_Paraboloid_copy")
      import
      implicit none
      type(c_Paraboloid) :: this
      type(c_Paraboloid) :: a_other_Paraboloid
    end subroutine F_Paraboloid_copy

    subroutine F_Paraboloid_getDatum(this, a_datum) &
      bind(C, name="c_Paraboloid_getDatum")
      import
      implicit none
      type(c_Paraboloid) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_datum
    end subroutine F_Paraboloid_getDatum

    subroutine F_Paraboloid_getReferenceFrame(this, a_frame) &
      bind(C, name="c_Paraboloid_getReferenceFrame")
      import
      implicit none
      type(c_Paraboloid) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_frame
    end subroutine F_Paraboloid_getReferenceFrame

    subroutine F_Paraboloid_getAlignedParaboloid(this, a_aligned_paraboloid) &
      bind(C, name="c_Paraboloid_getAlignedParaboloid")
      import
      implicit none
      type(c_Paraboloid) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_aligned_paraboloid
    end subroutine F_Paraboloid_getAlignedParaboloid

    subroutine F_Paraboloid_triangulateInsideCuboid(this, a_cuboid, a_surface) &
      bind(C, name="c_Paraboloid_triangulateInsideCuboid")
      import
      implicit none
      type(c_Paraboloid) :: this
      type(c_RectCub), intent(in)  :: a_cuboid
      type(c_TriangulatedParaboloid), intent(out) :: a_surface
    end subroutine F_Paraboloid_triangulateInsideCuboid

    subroutine F_Paraboloid_printToScreen(this) &
      bind(C, name="c_Paraboloid_printToScreen")
      import
      implicit none
      type(c_Paraboloid) :: this
    end subroutine F_Paraboloid_printToScreen

  end interface


  contains

    subroutine Paraboloid_class_new(this)
      implicit none
      type(Paraboloid_type), intent(inout) :: this
      call F_Paraboloid_new(this%c_object)
    end subroutine Paraboloid_class_new

    subroutine Paraboloid_class_newFromObjectAllocationServer(this, a_object_allocation_server)
      implicit none
      type(Paraboloid_type), intent(inout) :: this
      type(ObjServer_Paraboloid_type), intent(in) :: a_object_allocation_server
      call F_Paraboloid_newFromObjectAllocationServer(this%c_object, a_object_allocation_server%c_object)
    end subroutine Paraboloid_class_newFromObjectAllocationServer

    impure elemental subroutine Paraboloid_class_delete(this)
      implicit none
      type(Paraboloid_type), intent(in) :: this
      call F_Paraboloid_delete(this%c_object)
    end subroutine Paraboloid_class_delete

    subroutine Paraboloid_class_setDatum(this, a_datum)
      implicit none
      type(Paraboloid_type), intent(in) :: this
      real(IRL_double), dimension(1:3), intent(in) :: a_datum
      call F_Paraboloid_setDatum(this%c_object, a_datum)
    end subroutine Paraboloid_class_setDatum

    subroutine Paraboloid_class_setReferenceFrame(this, a_normal1, a_normal2, a_normal3)
      implicit none
      type(Paraboloid_type), intent(in) :: this
      real(IRL_double), dimension(1:3), intent(in) :: a_normal1
      real(IRL_double), dimension(1:3), intent(in) :: a_normal2
      real(IRL_double), dimension(1:3), intent(in) :: a_normal3
      call F_Paraboloid_setReferenceFrame(this%c_object, a_normal1, a_normal2, a_normal3)
    end subroutine Paraboloid_class_setReferenceFrame

    subroutine Paraboloid_class_setAlignedParaboloid(this, a_coeff_a, a_coeff_b)
      implicit none
      type(Paraboloid_type), intent(in) :: this
      real(IRL_double), intent(in) :: a_coeff_a
      real(IRL_double), intent(in) :: a_coeff_b
      call F_Paraboloid_setAlignedParaboloid(this%c_object, a_coeff_a, a_coeff_b)
    end subroutine Paraboloid_class_setAlignedParaboloid

    subroutine Paraboloid_class_setParaboloidFromPolygon(this, a_poly)
      implicit none
      type(Paraboloid_type), intent(in) :: this
      type(Poly_type), intent(in) :: a_poly
      call F_Paraboloid_setParaboloidFromPolygon(this%c_object, a_poly%c_object)
    end subroutine Paraboloid_class_setParaboloidFromPolygon

    subroutine Paraboloid_class_setParaboloidFromPlanarSep(this, a_plane, a_pt_ref)
      implicit none
      type(Paraboloid_type), intent(in) :: this
      type(PlanarSep_type), intent(in) :: a_plane
      real(IRL_double), dimension(1:3), intent(in) :: a_pt_ref
      call F_Paraboloid_setParaboloidFromPlanarSep(this%c_object, a_plane%c_object, a_pt_ref)
    end subroutine Paraboloid_class_setParaboloidFromPlanarSep

    subroutine Paraboloid_class_setParaboloidJibben(this, a_plane, a_cuboid, a_npoly, a_vfrac, a_nvert, a_vert_coords)
      implicit none
      type(Paraboloid_type), intent(in) :: this
      type(PlanarSep_type), intent(in) :: a_plane
      type(RectCub_type), intent(in) :: a_cuboid
      integer(IRL_UnsignedIndex_t), intent(in) :: a_npoly
      real(IRL_double), dimension(:), intent(in) :: a_vfrac
      integer(IRL_UnsignedIndex_t), dimension(:), intent(in) :: a_nvert
      real(IRL_double), dimension(:), intent(in) :: a_vert_coords
      call F_Paraboloid_setParaboloidJibben(this%c_object, a_plane%c_object, a_cuboid%c_object, a_npoly, &
      a_vfrac, a_nvert, a_vert_coords)
    end subroutine Paraboloid_class_setParaboloidJibben

    subroutine Paraboloid_class_setParaboloidFull(this)
      implicit none
      type(Paraboloid_type), intent(in) :: this
      call F_Paraboloid_setParaboloidFull(this%c_object)
    end subroutine Paraboloid_class_setParaboloidFull

    subroutine Paraboloid_class_setParaboloidEmpty(this)
      implicit none
      type(Paraboloid_type), intent(in) :: this
      call F_Paraboloid_setParaboloidEmpty(this%c_object)
    end subroutine Paraboloid_class_setParaboloidEmpty

    subroutine Paraboloid_class_copy(this, a_other_Paraboloid)
      implicit none
      type(Paraboloid_type), intent(inout) :: this
      type(Paraboloid_type), intent(in) :: a_other_Paraboloid
      call F_Paraboloid_copy(this%c_object, a_other_Paraboloid%c_object)
    end subroutine Paraboloid_class_copy

    function Paraboloid_class_getDatum(this) result(a_datum)
      implicit none
      type(Paraboloid_type), intent(in) :: this
      real(IRL_double), dimension(1:3) :: a_datum
      call F_Paraboloid_getDatum(this%c_object, a_datum)
    end function Paraboloid_class_getDatum

    function Paraboloid_class_getReferenceFrame(this) result(a_frame)
      implicit none
      type(Paraboloid_type), intent(in) :: this
      real(IRL_double), dimension(1:9) :: a_frame
      call F_Paraboloid_getReferenceFrame(this%c_object, a_frame)
    end function Paraboloid_class_getReferenceFrame

    function Paraboloid_class_getAlignedParaboloid(this) result(a_aligned_paraboloid)
      implicit none
      type(Paraboloid_type), intent(in) :: this
      real(IRL_double), dimension(1:2) :: a_aligned_paraboloid
      call F_Paraboloid_getAlignedParaboloid(this%c_object, a_aligned_paraboloid)
    end function Paraboloid_class_getAlignedParaboloid

    subroutine Paraboloid_class_triangulateInsideCuboid(this, a_cuboid, a_surface)
      type(Paraboloid_type), intent(in) :: this
      type(RectCub_type), intent(in) :: a_cuboid
      type(TriangulatedParaboloid_type) :: a_surface
      call F_Paraboloid_triangulateInsideCuboid(this%c_object, a_cuboid%c_object, a_surface%c_object)
    end subroutine Paraboloid_class_triangulateInsideCuboid

    function Paraboloid_class_isFlipped(this) result(a_flipped)
      implicit none
      type(Paraboloid_type), intent(in) :: this
      logical(1) :: a_flipped
      a_flipped = .false.
      return
    end function Paraboloid_class_isFlipped

    subroutine Paraboloid_class_printToScreen(this)
      implicit none
      type(Paraboloid_type), intent(in) :: this
      call F_Paraboloid_printToScreen(this%c_object)
    end subroutine Paraboloid_class_printToScreen


end module f_Paraboloid_class
