!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file irl_fortran_interface.f90
!! 
!! This file serves to provide a single
!! include directive when using the
!! IRL fortran interface.

!> \brief This is just a master wrapper for
!! the entire IRL fortran interface. For information
!! about each module, view the documentation for the 
!! module itself.
module irl_fortran_interface
  use f_definedtypes
  use f_ObjServer_PlanarSep_class
  use f_ObjServer_PlanarLoc_class
  use f_ObjServer_LocSepLink_class
  use f_ObjServer_LocSep_class  
  use f_ObjServer_LocLink_class
  use f_SepVM_class
  use f_SepVM_d3_class
  use f_SepVol_class
  use f_VMAN_class
  use f_PlanarLoc_class
  use f_PlanarSep_class
  use f_LocLink_class
  use f_PlanarSepPath_class
  use f_PlanarSepPathGroup_class
  use f_LocSep_class  
  use f_LocSepLink_class
  use f_LocSepGroupLink_class
  use f_Dod_class
  use f_CapDod_class
  use f_CapDod_LLLL_class
  use f_CapDod_LLLT_class
  use f_CapDod_LTLT_class
  use f_CapDod_LLTT_class
  use f_CapDod_LTTT_class
  use f_CapDod_TTTT_class
  use f_CapOcta_LLL_class
  use f_CapOcta_LLT_class
  use f_CapOcta_LTT_class
  use f_CapOcta_TTT_class  
  use f_SymTet_class
  use f_SymPyrmd_class
  use f_SymTriPrism_class
  use f_SymHex_class
  use f_CapDod_d3_class
  use f_Poly24_class
  use f_Poly24_d3_class
  use f_RectCub_class
  use f_Hex_class
  use f_TriPrism_class
  use f_Octa_class
  use f_Pyrmd_class
  use f_Tet_class
  use f_Tri_class
  use f_Poly_class
  use f_DivPoly_class
  use f_listvm_vman_class
  use f_TagAccVM_SepVM_class
  use f_TagAccVM_VM_class
  use f_TagAccVM_SepVol_class
  use f_TagAccVM_Vol_class
  use f_TagAccVM2_Vol_class
  use f_TagAccListVM_VMAN_class
  use f_ELVIRANeigh_class
  use f_LVIRANeigh_Hex_class
  use f_LVIRANeigh_Tet_class
  use f_LVIRANeigh_RectCub_class
  use f_R2PNeigh_RectCub_class
  use f_Constants
  use f_GeometricCuttingHelpers
  use f_CutPolygon
  use f_volumefractionmatching
  use f_getMoments
  use f_ReconstructionInterface
  use f_Serializer
  use f_ByteBuffer_class
  use f_R2PWeighting_class
  use f_OptimizationBehavior_class

end module irl_fortran_interface
