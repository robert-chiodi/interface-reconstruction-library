!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_geometriccuttinghelpers.f90
!!
!! This file provides access to helper functions
!! often used during geometric cutting.

!> This module contains mappings to the
!! IRL C interface that provides access to functions
!! often used to geoemtric cutting operations. See
!! the C interface file
!! src/c_interface/c_geometric_cutting_helpers.h 
!! for more information.
module f_GeometricCuttingHelpers
  use f_DefinedTypes
  use f_PlanarSep_class
  use f_PlanarLoc_class
  use f_LocSepLink_class
  use f_LocSepGroupLink_class
  implicit none

  interface isPtInt
    ! Check if Pt is internal to a PlanarSeparator
    module procedure isPtInt_PlanarSep
    ! Check if Pt is internal to a PlanarLocalizer
    module procedure isPtInt_PlanarLoc
  end interface isPtInt

  interface locatePt
    module procedure locatePt_LocSepLink
    module procedure locatePt_LocSepGroupLink    
  end interface

  interface
    function F_isPtInt_PlanarSep(a_pt, a_separator) result (is_internal) &
    bind(C, name="c_isPtInt_PlanarSep")
      import
      implicit none
      real(C_DOUBLE), dimension(*), intent(in) :: a_pt ! dimension(1:3)
      type(c_PlanarSep) :: a_separator ! Pointer to PlanarSeparator
      logical(C_BOOL) :: is_internal ! Boolean, true if internal
    end function F_isPtInt_PlanarSep
  end interface

  interface
    function F_isPtInt_PlanarLoc(a_pt, a_localizer) result (is_internal) &
    bind(C, name="c_isPtInt_PlanarLoc")
      import
      implicit none
      real(C_DOUBLE), dimension(*), intent(in) :: a_pt ! dimension(1:3)
      type(c_PlanarLoc) :: a_localizer ! Pointer to PlanarLocalizer
      logical(C_BOOL) :: is_internal ! Boolean, true if internal
    end function F_isPtInt_PlanarLoc
  end interface

  interface
    function F_locatePt_LocSepLink(a_pt, a_locseplink) result(a_id) &
    bind(C, name="c_locatePt_LocSepLink")
      import
      implicit none
      real(C_DOUBLE), dimension(*), intent(in) :: a_pt ! dimension(1:3)
      type(c_LocSepLink) :: a_locseplink ! Pointer to LocalizedSeparatorLink
      integer(C_INT) :: a_id ! ID of reconstruction point is internal to
    end function F_locatePt_LocSepLink
  end interface  

  interface
    function F_locatePt_LocSepGroupLink(a_pt, a_locsepgrouplink) result(a_id) &
    bind(C, name="c_locatePt_LocSepGroupLink")
      import
      implicit none
      real(C_DOUBLE), dimension(*), intent(in) :: a_pt ! dimension(1:3)
      type(c_LocSepGroupLink) :: a_locsepgrouplink ! Pointer to LocalizedSeparatorGroupLink
      integer(C_INT) :: a_id ! ID of reconstruction point is internal to
    end function F_locatePt_LocSepGroupLink
  end interface  

  
contains

  function isPtInt_PlanarSep(a_pt, a_separator) result (is_internal)
    implicit none
      real(IRL_double), dimension(3), intent(in) :: a_pt
      type(PlanarSep_type), intent(in) :: a_separator
      logical(1) :: is_internal
      is_internal = F_isPtInt_PlanarSep(a_pt, a_separator%c_object)
  end function isPtInt_PlanarSep

  function isPtInt_PlanarLoc(a_pt, a_localizer) result (is_internal)
    implicit none
      real(IRL_double), dimension(3), intent(in) :: a_pt
      type(PlanarLoc_type), intent(in) :: a_localizer
      logical(1) :: is_internal
      is_internal = F_isPtInt_PlanarLoc(a_pt, a_localizer%c_object)
  end function isPtInt_PlanarLoc

  function locatePt_LocSepLink(a_pt, a_locseplink) result (a_id)
    implicit none
      real(IRL_double), dimension(3), intent(in) :: a_pt
      type(LocSepLink_type), intent(in) :: a_locseplink
      integer(IRL_UnsignedIndex_t) :: a_id
      a_id = F_locatePt_LocSepLink(a_pt, a_locseplink%c_object)
  end function locatePt_LocSepLink

  function locatePt_LocSepGroupLink(a_pt, a_locsepgrouplink) result (a_id)
    implicit none
      real(IRL_double), dimension(3), intent(in) :: a_pt
      type(LocSepGroupLink_type), intent(in) :: a_locsepgrouplink
      integer(IRL_UnsignedIndex_t) :: a_id
      a_id = F_locatePt_LocSepGroupLink(a_pt, a_locsepgrouplink%c_object)
  end function locatePt_LocSepGroupLink
  
end module f_GeometricCuttingHelpers
