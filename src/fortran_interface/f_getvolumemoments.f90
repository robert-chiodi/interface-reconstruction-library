!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_getvolumemoments.f90
!!
!! This file deals with subdivinding and integrating
!! volume moments for polyhedra.

!> This module contains mappings to the
!! IRL C interface that deal with intersecting polyhedron
!! volumes and integrating these volumes to obtain
!! volumetric moments.
module f_getMoments
  use f_DefinedTypes
  use f_Dod_class
  use f_CapDod_class
  use f_CapDod_LLLL_class
  use f_CapDod_LLLT_class
  use f_CapDod_LTLT_class
  use f_CapDod_LLTT_class
  use f_CapDod_LTTT_class
  use f_CapDod_TTTT_class
  use f_SymTet_class
  use f_SymPyrmd_class
  use f_SymTriPrism_class
  use f_SymHex_class  
  use f_CapDod_d3_class
  use f_Poly24_class
  use f_Poly24_d3_class
  use f_RectCub_class
  use f_Pyrmd_class
  use f_Hex_class
  use f_TriPrism_class
  use f_Octa_class
  use f_CapOcta_LLL_class
  use f_CapOcta_LLT_class
  use f_CapOcta_LTT_class
  use f_CapOcta_TTT_class    
  use f_Tet_class
  use f_Poly_class
  use f_Tri_class
  use f_SepVM_class
  use f_SepVM_d3_class
  use f_SepVol_class
  use f_TagAccVM_SepVM_class
  use f_TagAccVM_VM_class
  use f_TagAccVM_SepVol_class
  use f_TagAccVM_Vol_class
  use f_TagAccVM2_Vol_class  
  use f_PlanarSep_class
  use f_PlanarSepPathGroup_class
  use f_PlanarLoc_class
  use f_LocLink_class
  use f_TagAccListVM_VMAN_class
  use f_LocSepLink_class
  use f_LocSep_class  
  use f_LocSepGroupLink_class  
  implicit none

  interface getMoments_setMethod
    module procedure getMoments_setMethod
  end interface getMoments_setMethod

  ! Moments that have been normalized by volume
  interface getNormMoments
    ! Cut Dod by LocSepLink to get SeparatedMoments<VM>
    module procedure getNormMoments_Dod_LocSepLink_SepVM
    ! Cut Dod by LocSepLink to get Volume
    module procedure getNormMoments_Dod_LocSepLink_Vol
    ! Cut TriPrism by LocSepLink to get Volume
    module procedure getNormMoments_TriPrism_LocSepLink_Vol
    ! Cut Octa by LocSepLink to get Volume
    module procedure getNormMoments_Octa_LocSepLink_Vol
    ! Cut CapDod by LocSepLink to get SeparatedMoments<VM>
    module procedure getNormMoments_CapDod_LocSepLink_SepVM
    ! Cut CapDodWithDoubles3 by LocSepLink to get SeparatedMoments<VMAndDoubles<3>>
    module procedure getNormMoments_CapDod_d3_LocSepLink_SepVM_d3
    ! Cut Poly24_d3 by LocSepLink to get SeparatedMoments<VMAndDoubles<3>>
    module procedure getNormMoments_Poly24_LocSepLink_SepVM
    ! Cut Poly24_d3 by LocSepLink to get SeparatedMoments<VMAndDoubles<3>>
    module procedure getNormMoments_Poly24_d3_LocSepLink_SepVM_d3
    ! Cut Tet by LocSepLink to get SeparatedMoments<VM>
    module procedure getNormMoments_Tet_LocSepLink_SepVM
    ! Cut Dod by LocSepLink to get TagAccVM<SeparatedMoments<VM>>
    module procedure getNormMoments_Dod_LocSepLink_TagAccVM_SepVol
    ! Cut Tet by LocSepLink to get TagAccVM<SeparatedMoments<VM>>
    module procedure getNormMoments_Tet_LocSepLink_TagAccVM_SepVol        
    ! Cut Octa by LocSepLink to get TagAccVM<SeparatedMoments<VM>>
    module procedure getNormMoments_Octa_LocSepLink_TagAccVM_SepVol        
    ! Cut RectCub by PlanarSep to get Volume
    module procedure getNormMoments_RectCub_PlanarSep_Vol
    ! Cut Tet by PlanarSep to get Volume
    module procedure getNormMoments_Tet_PlanarSep_Vol
    ! Cut TriPrism  by PlanarSep to get Volume
    module procedure getNormMoments_TriPrism_PlanarSep_Vol
    ! Cut Pyrmd  by PlanarSep to get Volume
    module procedure getNormMoments_Pyrmd_PlanarSep_Vol            
    ! Cut Hex by PlanarSep to get Volume
    module procedure getNormMoments_Hex_PlanarSep_Vol
    ! Cut Hex by PlanarSep to get Volume
    module procedure getNormMoments_Hex_PlanarSep_SepVM
    ! Cut Dod by PlanarSep to get SeparatedMoments<VM>
    module procedure getNormMoments_Dod_PlanarSep_SepVM
    ! Cut CapDod by LocSepLink to get TagAccVM<SeparatedMoments<VM>>
    module procedure getNormMoments_CapDod_LocSepLink_TagAccVM_SepVM
    ! Cut Dod by LocSepLink to get TagAccVM<SeparatedMoments<VM>>
    module procedure getNormMoments_Dod_LocSepLink_TagAccVM_SepVM
    ! Cut RectCub by PlanarSep to get SeparatedMoments<VM>
    module procedure getNormMoments_RectCub_PlanarSep_SepVM
    ! Cut Tri by PlanarLoc to get Volume (Surface Area)
    module procedure getNormMoments_Tri_PlanarLoc_Vol
    ! Cut Poly by PlanarLoc to get Volume (Surface Area)
    module procedure getNormMoments_Poly_PlanarLoc_Vol
    ! Cut Tri by LocLink to get TagAccVM<VM>
    module procedure getNormMoments_Tri_LocLink_TagAccVM_VM
    ! Cut CapDod_LLLL by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapDod_LLLL_LocSepLink_TagAccVM_SepVol
    ! Cut CapDod_LLLT by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapDod_LLLT_LocSepLink_TagAccVM_SepVol
    ! Cut CapDod_LTLT by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapDod_LTLT_LocSepLink_TagAccVM_SepVol
    ! Cut CapDod_LLTT by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapDod_LLTT_LocSepLink_TagAccVM_SepVol
    ! Cut CapDod_LTTT by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapDod_LTTT_LocSepLink_TagAccVM_SepVol
    ! Cut CapDod_TTTT by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapDod_TTTT_LocSepLink_TagAccVM_SepVol
    ! Cut CapOcta_LLL by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapOcta_LLL_LocSepLink_TagAccVM_SepVol
    ! Cut CapOcta_LLT by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapOcta_LLT_LocSepLink_TagAccVM_SepVol
    ! Cut CapOcta_LTT by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapOcta_LTT_LocSepLink_TagAccVM_SepVol
    ! Cut CapOcta_TTT by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapOcta_TTT_LocSepLink_TagAccVM_SepVol 
    ! Cut SymTet by LocSepGroupLink to get TagAccVM2<Vol>
    module procedure getNormMoments_SymTet_LocSepGroupLink_TagAccVM2_Vol
    ! Cut SymPyrmd by LocSepGroupLink to get TagAccVM2<Vol>
    module procedure getNormMoments_SymPyrmd_LocSepGroupLink_TagAccVM2_Vol
    ! Cut SymTriPrism by LocSepGroupLink to get TagAccVM2<Vol>
    module procedure getNormMoments_SymTriPrism_LocSepGroupLink_TagAccVM2_Vol
    ! Cut SymHex by LocSepGroupLink to get TagAccVM2<Vol>
    module procedure getNormMoments_SymHex_LocSepGroupLink_TagAccVM2_Vol       
    ! Cut SymTet by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_SymTet_LocSepLink_TagAccVM_SepVol
    ! Cut SymPyrmd by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVol
    ! Cut SymTriPrism by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVol
    ! Cut SymHex by LocSepLink to get TagAccVM<SepVM>
    module procedure getNormMoments_SymHex_LocSepLink_TagAccVM_SepVol
    ! Cut SymTet by LocSepLink to get TagAccVM<SepVM>
    module procedure getNormMoments_SymTet_LocSepLink_TagAccVM_SepVM
    ! Cut SymPyrmd by LocSepLink to get TagAccVM<SepVM>
    module procedure getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVM
    ! Cut SymTriPrism by LocSepLink to get TagAccVM<SepVM>
    module procedure getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVM
    ! Cut SymHex by LocSepLink to get TagAccVM<SepVM>
    module procedure getNormMoments_SymHex_LocSepLink_TagAccVM_SepVM         
    ! Cut CapDod_LLLL by LocSepLink to get SepVol
    module procedure getNormMoments_CapDod_LLLL_LocSepLink_SepVol
    ! Cut CapDod_LLLT by LocSepLink to get SepVol
    module procedure getNormMoments_CapDod_LLLT_LocSepLink_SepVol
    ! Cut CapDod_LTLT by LocSepLink to get SepVol
    module procedure getNormMoments_CapDod_LTLT_LocSepLink_SepVol
    ! Cut CapDod_LLTT by LocSepLink to get SepVol
    module procedure getNormMoments_CapDod_LLTT_LocSepLink_SepVol
    ! Cut CapDod_LTTT by LocSepLink to get SepVol
    module procedure getNormMoments_CapDod_LTTT_LocSepLink_SepVol
    ! Cut CapDod_TTTT by LocSepLink to get SepVol
    module procedure getNormMoments_CapDod_TTTT_LocSepLink_SepVol
    ! Cut CapOcta_LLL by LocSepLink to get SepVol
    module procedure getNormMoments_CapOcta_LLL_LocSepLink_SepVol
    ! Cut CapOcta_LLT by LocSepLink to get SepVol
    module procedure getNormMoments_CapOcta_LLT_LocSepLink_SepVol
    ! Cut CapOcta_LTT by LocSepLink to get SepVol
    module procedure getNormMoments_CapOcta_LTT_LocSepLink_SepVol
    ! Cut CapOcta_TTT by LocSepLink to get SepVol
    module procedure getNormMoments_CapOcta_TTT_LocSepLink_SepVol        
    ! Cut Tet by LocSepLink to get SepVol
    module procedure getNormMoments_Tet_LocSepLink_SepVol
    ! Cut CapDod_LLLL by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapDod_LLLT by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapDod_LTLT by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapDod_LLTT by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapDod_LTTT by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapDod_TTTT by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapOcta_LLL by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapOcta_LLT by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapOcta_LTT by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapOcta_TTT by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapDod_LLLL by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM_Vol
    ! Cut CapDod_LLLT by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM_Vol
    ! Cut CapDod_LTLT by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM_Vol
    ! Cut CapDod_LLTT by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM_Vol
    ! Cut CapDod_LTTT by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM_Vol
    ! Cut CapDod_TTTT by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM_Vol
    ! Cut CapOcta_LLL by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM_Vol
    ! Cut CapOcta_LLT by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM_Vol
    ! Cut CapOcta_LTT by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM_Vol
    ! Cut CapOcta_TTT by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM_Vol
    ! Cut Tet by PlanarSepPathGroup to get TagAccVM_Vol
    module procedure getNormMoments_Tet_PlanarSepPathGroup_TagAccVM_Vol
    ! Cut Pyrmd by PlanarSepPathGroup to get TagAccVM_Vol
    module procedure getNormMoments_Pyrmd_PlanarSepPathGroup_TagAccVM_Vol
    ! Cut TriPrism by PlanarSepPathGroup to get TagAccVM_Vol
    module procedure getNormMoments_TriPrism_PlanarSepPathGroup_TagAccVM_Vol
    ! Cut Hex by PlanarSepPathGroup to get TagAccVM_Vol
    module procedure getNormMoments_Hex_PlanarSepPathGroup_TagAccVM_Vol
    ! Cut Polyhedra by LocalizedSeparator to get Separated Volumes
    module procedure getNormMoments_SymTet_LocSep_SepVol
    module procedure getNormMoments_SymPyrmd_LocSep_SepVol
    module procedure getNormMoments_SymTriPrism_LocSep_SepVol
    module procedure getNormMoments_SymHex_LocSep_SepVol
    module procedure getNormMoments_CapOcta_LLL_LocSep_SepVol
    module procedure getNormMoments_CapOcta_LLT_LocSep_SepVol
    module procedure getNormMoments_CapOcta_LTT_LocSep_SepVol
    module procedure getNormMoments_CapOcta_TTT_LocSep_SepVol
    module procedure getNormMoments_CapDod_LLLL_LocSep_SepVol
    module procedure getNormMoments_CapDod_LLLT_LocSep_SepVol
    module procedure getNormMoments_CapDod_LTLT_LocSep_SepVol
    module procedure getNormMoments_CapDod_LLTT_LocSep_SepVol
    module procedure getNormMoments_CapDod_LTTT_LocSep_SepVol
    module procedure getNormMoments_CapDod_TTTT_LocSep_SepVol      
  end interface getNormMoments

  ! Moments that are still weighted by volume
  interface getMoments
    ! Cut Dod by LocSepLink to get Volume
    module procedure getNormMoments_Dod_LocSepLink_Vol
    ! Cut TriPrism by LocSepLink to get Volume
    module procedure getNormMoments_TriPrism_LocSepLink_Vol
    ! Cut Octa by LocSepLink to get Volume
    module procedure getNormMoments_Octa_LocSepLink_Vol
    ! Cut CapDod by LocSepLink to get SeparatedMoments<VM>
    module procedure getMoments_CapDod_LocSepLink_SepVM
    ! Cut Dod by LocSepLink to get SeparatedMoments<VM>
    module procedure getMoments_Dod_LocSepLink_SepVM
    ! Cut Dod by LocSepLink to get TagAccVM<SeparatedMoments<VM>>
    module procedure getMoments_Dod_LocSepLink_TagAccVM_SepVM
    ! Cut Dod by LocSepLink to get TagAccVM<SeparatedMoments<VM>>
    module procedure getNormMoments_Dod_LocSepLink_TagAccVM_SepVol
    ! Cut Tet by LocSepLink to get TagAccVM<SeparatedMoments<VM>>
    module procedure getNormMoments_Tet_LocSepLink_TagAccVM_SepVol            
    ! Cut Octa by LocSepLink to get TagAccVM<SeparatedMoments<VM>>
    module procedure getMoments_Octa_LocSepLink_TagAccVM_SepVM
    ! Cut Octa by LocSepLink to get TagAccVM<SeparatedMoments<VM>>
    module procedure getNormMoments_Octa_LocSepLink_TagAccVM_SepVol    
    ! Cut CapDod by LocSepLink to get SeparatedMoments<VM>
    module procedure getMoments_Poly24_LocSepLink_SepVM
    ! Cut Tri by LocLink to get TagAccListVM<VMAN>
    module procedure getMoments_Tri_LocLink_TagAccListVM_VMAN
    ! Cut RectCub by PlanarSep to get Volume
    module procedure getNormMoments_RectCub_PlanarSep_Vol
    ! Cut Tet by PlanarSep to get Volume
    module procedure getNormMoments_Tet_PlanarSep_Vol
    ! Cut TriPrism  by PlanarSep to get Volume
    module procedure getNormMoments_TriPrism_PlanarSep_Vol
    ! Cut Pyrmd  by PlanarSep to get Volume
    module procedure getNormMoments_Pyrmd_PlanarSep_Vol                
    ! Cut Hex by PlanarSep to get Volume
    module procedure getNormMoments_Hex_PlanarSep_Vol
    ! Cut Tri by PlanarLoc to get Volume (Surface Area)
    module procedure getNormMoments_Tri_PlanarLoc_Vol
    ! Cut Poly by PlanarLoc to get Volume (Surface Area)
    module procedure getNormMoments_Poly_PlanarLoc_Vol
    ! Cut CapDod_LLLL by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapDod_LLLL_LocSepLink_TagAccVM_SepVol
    ! Cut CapDod_LLLT by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapDod_LLLT_LocSepLink_TagAccVM_SepVol
    ! Cut CapDod_LTLT by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapDod_LTLT_LocSepLink_TagAccVM_SepVol
    ! Cut CapDod_LLTT by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapDod_LLTT_LocSepLink_TagAccVM_SepVol
    ! Cut CapDod_LTTT by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapDod_LTTT_LocSepLink_TagAccVM_SepVol
    ! Cut CapDod_TTTT by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapDod_TTTT_LocSepLink_TagAccVM_SepVol
    ! Cut CapOcta_LLL by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapOcta_LLL_LocSepLink_TagAccVM_SepVol
    ! Cut CapOcta_LLT by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapOcta_LLT_LocSepLink_TagAccVM_SepVol
    ! Cut CapOcta_LTT by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapOcta_LTT_LocSepLink_TagAccVM_SepVol
    ! Cut CapOcta_TTT by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_CapOcta_TTT_LocSepLink_TagAccVM_SepVol
    ! Cut SymTet by LocSepGroupLink to get TagAccVM2<Vol>
    module procedure getNormMoments_SymTet_LocSepGroupLink_TagAccVM2_Vol
    ! Cut SymPyrmd by LocSepGroupLink to get TagAccVM2<Vol>
    module procedure getNormMoments_SymPyrmd_LocSepGroupLink_TagAccVM2_Vol
    ! Cut SymTriPrism by LocSepGroupLink to get TagAccVM2<Vol>
    module procedure getNormMoments_SymTriPrism_LocSepGroupLink_TagAccVM2_Vol
    ! Cut SymHex by LocSepGroupLink to get TagAccVM2<Vol>
    module procedure getNormMoments_SymHex_LocSepGroupLink_TagAccVM2_Vol           
    ! Cut SymTet by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_SymTet_LocSepLink_TagAccVM_SepVol
    ! Cut SymPyrmd by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVol
    ! Cut SymTriPrism by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVol
    ! Cut SymHex by LocSepLink to get TagAccVM<SepVol>
    module procedure getNormMoments_SymHex_LocSepLink_TagAccVM_SepVol
    ! Cut SymTet by LocSepLink to get TagAccVM<SepVM>
    module procedure getMoments_SymTet_LocSepLink_TagAccVM_SepVM
    ! Cut SymPyrmd by LocSepLink to get TagAccVM<SepVM>
    module procedure getMoments_SymPyrmd_LocSepLink_TagAccVM_SepVM
    ! Cut SymTriPrism by LocSepLink to get TagAccVM<SepVM>
    module procedure getMoments_SymTriPrism_LocSepLink_TagAccVM_SepVM
    ! Cut SymHex by LocSepLink to get TagAccVM<SepVM>
    module procedure getMoments_SymHex_LocSepLink_TagAccVM_SepVM             
    ! Cut CapDod_LLLL by LocSepLink to get SepVol
    module procedure getNormMoments_CapDod_LLLL_LocSepLink_SepVol
    ! Cut CapDod_LLLT by LocSepLink to get SepVol
    module procedure getNormMoments_CapDod_LLLT_LocSepLink_SepVol
    ! Cut CapDod_LTLT by LocSepLink to get SepVol
    module procedure getNormMoments_CapDod_LTLT_LocSepLink_SepVol
    ! Cut CapDod_LLTT by LocSepLink to get SepVol
    module procedure getNormMoments_CapDod_LLTT_LocSepLink_SepVol
    ! Cut CapDod_LTTT by LocSepLink to get SepVol
    module procedure getNormMoments_CapDod_LTTT_LocSepLink_SepVol
    ! Cut CapDod_TTTT by LocSepLink to get SepVol
    module procedure getNormMoments_CapDod_TTTT_LocSepLink_SepVol
    ! Cut CapOcta_LLL by LocSepLink to get SepVol
    module procedure getNormMoments_CapOcta_LLL_LocSepLink_SepVol
    ! Cut CapOcta_LLT by LocSepLink to get SepVol
    module procedure getNormMoments_CapOcta_LLT_LocSepLink_SepVol
    ! Cut CapOcta_LTT by LocSepLink to get SepVol
    module procedure getNormMoments_CapOcta_LTT_LocSepLink_SepVol
    ! Cut CapOcta_TTT by LocSepLink to get SepVol
    module procedure getNormMoments_CapOcta_TTT_LocSepLink_SepVol        
    ! Cut Tet by LocSepLink to get SepVol
    module procedure getNormMoments_Tet_LocSepLink_SepVol        
    ! Cut CapDod_LLLL by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapDod_LLLT by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapDod_LTLT by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapDod_LLTT by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapDod_LTTT by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapDod_TTTT by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapOcta_LLL by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapOcta_LLT by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapOcta_LTT by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapOcta_TTT by LocSepGroupLink to get TagAccVM<TagAccVM<Vol>>
    module procedure getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM2_Vol
    ! Cut CapDod_LLLL by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM_Vol
    ! Cut CapDod_LLLT by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM_Vol
    ! Cut CapDod_LTLT by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM_Vol
    ! Cut CapDod_LLTT by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM_Vol
    ! Cut CapDod_LTTT by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM_Vol
    ! Cut CapDod_TTTT by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM_Vol
    ! Cut CapOcta_LLL by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM_Vol
    ! Cut CapOcta_LLT by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM_Vol
    ! Cut CapOcta_LTT by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM_Vol
    ! Cut CapOcta_TTT by LocSepGroupLink to get TagAccVM_Vol
    module procedure getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM_Vol
    ! Cut Tet by PlanarSepPathGroup to get TagAccVM_Vol
    module procedure getNormMoments_Tet_PlanarSepPathGroup_TagAccVM_Vol
    ! Cut Pyrmd by PlanarSepPathGroup to get TagAccVM_Vol
    module procedure getNormMoments_Pyrmd_PlanarSepPathGroup_TagAccVM_Vol
    ! Cut TriPrism by PlanarSepPathGroup to get TagAccVM_Vol
    module procedure getNormMoments_TriPrism_PlanarSepPathGroup_TagAccVM_Vol
    ! Cut Hex by PlanarSepPathGroup to get TagAccVM_Vol
    module procedure getNormMoments_Hex_PlanarSepPathGroup_TagAccVM_Vol
    ! Cut Polyhedra by LocalizedSeparator to get Separated Volumes
    module procedure getNormMoments_SymTet_LocSep_SepVol
    module procedure getNormMoments_SymPyrmd_LocSep_SepVol
    module procedure getNormMoments_SymTriPrism_LocSep_SepVol
    module procedure getNormMoments_SymHex_LocSep_SepVol
    module procedure getNormMoments_CapOcta_LLL_LocSep_SepVol
    module procedure getNormMoments_CapOcta_LLT_LocSep_SepVol
    module procedure getNormMoments_CapOcta_LTT_LocSep_SepVol
    module procedure getNormMoments_CapOcta_TTT_LocSep_SepVol
    module procedure getNormMoments_CapDod_LLLL_LocSep_SepVol
    module procedure getNormMoments_CapDod_LLLT_LocSep_SepVol
    module procedure getNormMoments_CapDod_LTLT_LocSep_SepVol
    module procedure getNormMoments_CapDod_LLTT_LocSep_SepVol
    module procedure getNormMoments_CapDod_LTTT_LocSep_SepVol
    module procedure getNormMoments_CapDod_TTTT_LocSep_SepVol       
  end interface getMoments

  interface
    subroutine F_getMoments_setMethod(a_cutting_method) &
    bind(C, name="c_getMoments_setMethod")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      integer(C_INT) :: a_cutting_method
    end subroutine F_getMoments_setMethod
  end interface

  interface
    subroutine F_getNormMoments_Dod_LocSepLink_SepVM(a_Dod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Dod_LocSepLink_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Dod) :: a_Dod ! Pointer to Dod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVM) :: a_moments_to_return ! Where separated moments is returned to
    end subroutine F_getNormMoments_Dod_LocSepLink_SepVM
  end interface
  
  interface
    subroutine F_getNormMoments_Dod_LocSepLink_Vol(a_Dod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Dod_LocSepLink_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Dod) :: a_Dod ! Pointer to Dod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      real(C_DOUBLE) :: a_moments_to_return ! Where volume is returned to
    end subroutine F_getNormMoments_Dod_LocSepLink_Vol
  end interface
  
  interface
    subroutine F_getNormMoments_TriPrism_LocSepLink_Vol(a_triangular_prism, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_TriPrism_LocSepLink_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_TriPrism) :: a_triangular_prism ! Pointer to TriangularPrism object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      real(C_DOUBLE) :: a_moments_to_return ! Where volume is returned to
    end subroutine F_getNormMoments_TriPrism_LocSepLink_Vol
  end interface
  
  interface
    subroutine F_getNormMoments_Octa_LocSepLink_Vol(a_octahedron, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Octa_LocSepLink_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Octa) :: a_octahedron ! Pointer to Octahedron object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      real(C_DOUBLE) :: a_moments_to_return ! Where volume is returned to
    end subroutine F_getNormMoments_Octa_LocSepLink_Vol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LocSepLink_SepVM(a_Capped_Dod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LocSepLink_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod) :: a_Capped_Dod ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVM) :: a_moments_to_return ! Where separated moments is returned to
    end subroutine F_getNormMoments_CapDod_LocSepLink_SepVM
  end interface

  interface
    subroutine F_getNormMoments_CapDod_d3_LocSepLink_SepVM_d3(a_Capped_Dod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_d3_LocSepLink_SepVM_d3")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_d3) :: a_Capped_Dod ! Pointer to CapDodWithDoubles3 object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVM_d3) :: a_moments_to_return ! Where separated moments is returned to
    end subroutine F_getNormMoments_CapDod_d3_LocSepLink_SepVM_d3
  end interface

  interface
    subroutine F_getNormMoments_Poly24_LocSepLink_SepVM(a_polyhedron_24, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Poly24_LocSepLink_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Poly24) :: a_polyhedron_24 ! Pointer to Poly24 object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVM) :: a_moments_to_return ! Where separated moments is returned to
    end subroutine F_getNormMoments_Poly24_LocSepLink_SepVM
  end interface

  interface
    subroutine F_getNormMoments_Poly24_d3_LocSepLink_SepVM_d3(a_polyhedron_24, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Poly24_d3_LocSepLink_SepVM_d3")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Poly24_d3) :: a_polyhedron_24 ! Pointer to Poly24_d3 object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVM_d3) :: a_moments_to_return ! Where separated moments is returned to
    end subroutine F_getNormMoments_Poly24_d3_LocSepLink_SepVM_d3
  end interface

  interface
    subroutine F_getMoments_CapDod_LocSepLink_SepVM(a_Capped_Dod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getMoments_CapDod_LocSepLink_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod) :: a_Capped_Dod ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVM) :: a_moments_to_return ! Where separated moments is returned to
    end subroutine F_getMoments_CapDod_LocSepLink_SepVM
  end interface

  interface
    subroutine F_getMoments_Dod_LocSepLink_SepVM(a_Dod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getMoments_Dod_LocSepLink_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Dod) :: a_Dod ! Pointer to Dod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVM) :: a_moments_to_return ! Where separated moments is returned to
    end subroutine F_getMoments_Dod_LocSepLink_SepVM
  end interface

  interface
    subroutine F_getMoments_Poly24_LocSepLink_SepVM(a_polyhedron_24, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getMoments_Poly24_LocSepLink_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Poly24) :: a_polyhedron_24 ! Pointer to Poly24 object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVM) :: a_moments_to_return ! Where separated moments is returned to
    end subroutine F_getMoments_Poly24_LocSepLink_SepVM
  end interface

  interface
    subroutine F_getNormMoments_Tet_LocSepLink_SepVM(a_tet, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Tet_LocSepLink_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Tet) :: a_tet ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVM) :: a_moments_to_return ! Where separated moments is returned to
    end subroutine F_getNormMoments_Tet_LocSepLink_SepVM
  end interface

  interface
    subroutine F_getNormMoments_RectCub_PlanarSep_Vol(a_rectangular_cuboid, a_planar_separator, a_moments_to_return) &
    bind(C, name="c_getNormMoments_RectCub_PlanarSep_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_RectCub) :: a_rectangular_cuboid ! Pointer to RectangularCuboid object
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
      real(C_DOUBLE) :: a_moments_to_return ! Where volume is returned to
    end subroutine F_getNormMoments_RectCub_PlanarSep_Vol
  end interface

  interface
    subroutine F_getNormMoments_Tet_PlanarSep_Vol(a_tet, a_planar_separator, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Tet_PlanarSep_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Tet) :: a_tet ! Pointer to Tet object
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
      real(C_DOUBLE) :: a_moments_to_return ! Where volume is returned to
    end subroutine F_getNormMoments_Tet_PlanarSep_Vol
  end interface 

  interface
    subroutine F_getNormMoments_TriPrism_PlanarSep_Vol(a_tri_prism, a_planar_separator, a_moments_to_return) &
    bind(C, name="c_getNormMoments_TriPrism_PlanarSep_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_TriPrism) :: a_tri_prism ! Pointer to TriangularPrism object
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
      real(C_DOUBLE) :: a_moments_to_return ! Where volume is returned to
    end subroutine F_getNormMoments_TriPrism_PlanarSep_Vol
  end interface

  interface
    subroutine F_getNormMoments_Pyrmd_PlanarSep_Vol(a_pyramid, a_planar_separator, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Pyrmd_PlanarSep_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Pyrmd) :: a_pyramid ! Pointer to Pyramid object
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
      real(C_DOUBLE) :: a_moments_to_return ! Where volume is returned to
    end subroutine F_getNormMoments_Pyrmd_PlanarSep_Vol
  end interface
  
  interface
    subroutine F_getNormMoments_Hex_PlanarSep_Vol(a_hexahedron, a_planar_separator, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Hex_PlanarSep_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Hex) :: a_hexahedron ! Pointer to a_hexahedron object
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
      real(C_DOUBLE) :: a_moments_to_return ! Where volume is returned to
    end subroutine F_getNormMoments_Hex_PlanarSep_Vol
  end interface
  
  interface
    subroutine F_getNormMoments_Hex_PlanarSep_SepVM(a_hexahedron, a_planar_separator, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Hex_PlanarSep_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Hex) :: a_hexahedron ! Pointer to a_hexahedron object
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
      type(C_SepVM) :: a_moments_to_return ! Where volume is returned to
    end subroutine F_getNormMoments_Hex_PlanarSep_SepVM
  end interface    

  interface
    subroutine F_getNormMoments_Dod_PlanarSep_SepVM(a_Dod, a_planar_separator, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Dod_PlanarSep_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Dod) :: a_Dod ! Pointer to Dod object
      type(c_PlanarSep) :: a_planar_separator ! Pointer to PlanarSep object
      type(c_SepVM) :: a_moments_to_return ! Where separated moments is returned to
    end subroutine F_getNormMoments_Dod_PlanarSep_SepVM
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LocSepLink_TagAccVM_SepVM(a_Capped_Dod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LocSepLink_TagAccVM_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod) :: a_Capped_Dod ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVM) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_CapDod_LocSepLink_TagAccVM_SepVM
  end interface

  interface
    subroutine F_getNormMoments_Dod_LocSepLink_TagAccVM_SepVM(a_Dod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Dod_LocSepLink_TagAccVM_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Dod) :: a_Dod ! Pointer to Dod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVM) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_Dod_LocSepLink_TagAccVM_SepVM
  end interface

  interface
    subroutine F_getMoments_Dod_LocSepLink_TagAccVM_SepVM(a_Dod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getMoments_Dod_LocSepLink_TagAccVM_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Dod) :: a_Dod ! Pointer to Dod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVM) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getMoments_Dod_LocSepLink_TagAccVM_SepVM
  end interface

  interface
    subroutine F_getNormMoments_Dod_LocSepLink_TagAccVM_SepVol(a_Dod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Dod_LocSepLink_TagAccVM_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Dod) :: a_Dod ! Pointer to Dod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVol) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_Dod_LocSepLink_TagAccVM_SepVol
  end interface

  interface
    subroutine F_getNormMoments_Tet_LocSepLink_TagAccVM_SepVol(a_tet, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Tet_LocSepLink_TagAccVM_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Tet) :: a_tet ! Pointer to Dod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVol) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_Tet_LocSepLink_TagAccVM_SepVol
  end interface   
 
  interface
    subroutine F_getMoments_Octa_LocSepLink_TagAccVM_SepVM(a_octahedron, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getMoments_Octa_LocSepLink_TagAccVM_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Octa) :: a_octahedron ! Pointer to Dod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVM) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<Vol>> is stored
    end subroutine F_getMoments_Octa_LocSepLink_TagAccVM_SepVM
  end interface

  interface
    subroutine F_getNormMoments_Octa_LocSepLink_TagAccVM_SepVol(a_octahedron, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Octa_LocSepLink_TagAccVM_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Octa) :: a_octahedron ! Pointer to Dod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVol) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<Vol>> is stored
    end subroutine F_getNormMoments_Octa_LocSepLink_TagAccVM_SepVol
  end interface  
 
  interface
    subroutine F_getNormMoments_RectCub_PlanarSep_SepVM(a_rectangular_cuboid, a_planar_separator, a_moments_to_return) &
    bind(C, name="c_getNormMoments_RectCub_PlanarSep_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_RectCub) :: a_rectangular_cuboid ! Pointer to Dod object
      type(c_PlanarSep) :: a_planar_separator ! Pointer to LocSepLink object
      type(c_SepVM) :: a_moments_to_return ! Where separated moments is returned to
    end subroutine F_getNormMoments_RectCub_PlanarSep_SepVM
  end interface

  interface
    subroutine F_getNormMoments_Tri_LocLink_TagAccVM_VM(a_tri, a_localizer_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Tri_LocLink_TagAccVM_VM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Tri) :: a_tri ! Pointer to Tri object
      type(c_LocLink) :: a_localizer_link ! Pointer to LocLink object
      type(c_TagAccVM_VM) :: a_moments_to_return ! Pointer to TagAccVM<VM>
    end subroutine F_getNormMoments_Tri_LocLink_TagAccVM_VM
  end interface

  interface
    subroutine F_getNormMoments_Tri_PlanarLoc_Vol(a_tri, a_planar_localizer, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Tri_PlanarLoc_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Tri) :: a_tri ! Pointer to Tri object
      type(c_PlanarLoc) :: a_planar_localizer ! Pointer to PlanarLoc object
      real(C_DOUBLE) :: a_moments_to_return ! Pointer to double to to store volume
    end subroutine F_getNormMoments_Tri_PlanarLoc_Vol
  end interface

  interface
    subroutine F_getNormMoments_Poly_PlanarLoc_Vol(a_poly, a_planar_localizer, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Poly_PlanarLoc_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Poly) :: a_poly ! Pointer to Poly object
      type(c_PlanarLoc) :: a_planar_localizer ! Pointer to PlanarLoc object
      real(C_DOUBLE) :: a_moments_to_return ! Pointer to double to to store volume
    end subroutine F_getNormMoments_Poly_PlanarLoc_Vol
  end interface

  interface
    subroutine F_getMoments_Tri_LocLink_TagAccListVM_VMAN(a_tri, a_localizer_link, a_moments_to_return) &
    bind(C, name="c_getMoments_Tri_LocLink_TagAccListVM_VMAN")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Tri) :: a_tri ! Pointer to Tri object
      type(c_LocLink) :: a_localizer_link ! Pointer to LocLink object
      type(c_TagAccListVM_VMAN) :: a_moments_to_return ! Pointer to TagAccListVM<VMAN>
    end subroutine F_getMoments_Tri_LocLink_TagAccListVM_VMAN
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LLLL_LocSepLink_TagAccVM_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LLLL_LocSepLink_TagAccVM_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LLLL) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVol) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_CapDod_LLLL_LocSepLink_TagAccVM_SepVol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LLLT_LocSepLink_TagAccVM_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LLLT_LocSepLink_TagAccVM_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LLLT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVol) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_CapDod_LLLT_LocSepLink_TagAccVM_SepVol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LTLT_LocSepLink_TagAccVM_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LTLT_LocSepLink_TagAccVM_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LTLT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVol) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_CapDod_LTLT_LocSepLink_TagAccVM_SepVol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LLTT_LocSepLink_TagAccVM_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LLTT_LocSepLink_TagAccVM_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LLTT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVol) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_CapDod_LLTT_LocSepLink_TagAccVM_SepVol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LTTT_LocSepLink_TagAccVM_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LTTT_LocSepLink_TagAccVM_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LTTT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVol) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_CapDod_LTTT_LocSepLink_TagAccVM_SepVol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_TTTT_LocSepLink_TagAccVM_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_TTTT_LocSepLink_TagAccVM_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_TTTT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVol) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_CapDod_TTTT_LocSepLink_TagAccVM_SepVol
  end interface    

  interface
    subroutine F_getNormMoments_CapOcta_LLL_LocSepLink_TagAccVM_SepVol(a_CapOcta, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapOcta_LLL_LocSepLink_TagAccVM_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapOcta_LLL) :: a_CapOcta ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVol) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_CapOcta_LLL_LocSepLink_TagAccVM_SepVol
  end interface      

  interface
    subroutine F_getNormMoments_CapOcta_LLT_LocSepLink_TagAccVM_SepVol(a_CapOcta, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapOcta_LLT_LocSepLink_TagAccVM_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapOcta_LLT) :: a_CapOcta ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVol) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_CapOcta_LLT_LocSepLink_TagAccVM_SepVol
  end interface

  interface
     subroutine F_getNormMoments_CapOcta_LTT_LocSepLink_TagAccVM_SepVol(a_CapOcta, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapOcta_LTT_LocSepLink_TagAccVM_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapOcta_LTT) :: a_CapOcta ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVol) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_CapOcta_LTT_LocSepLink_TagAccVM_SepVol
  end interface

  interface
     subroutine F_getNormMoments_CapOcta_TTT_LocSepLink_TagAccVM_SepVol(a_CapOcta, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapOcta_TTT_LocSepLink_TagAccVM_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapOcta_TTT) :: a_CapOcta ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVol) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_CapOcta_TTT_LocSepLink_TagAccVM_SepVol
  end interface

  interface
     subroutine F_getNormMoments_SymTet_LocSepGroupLink_TagAccVM2_Vol(a_sym_tet, a_localized_separator_group_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_SymTet_LocSepGroupLink_TagAccVM2_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_SymTet) :: a_sym_tet ! Pointer to SymTet object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM2_Vol) :: a_moments_to_return ! Where TagAccVM2<Vol> is stored
    end subroutine F_getNormMoments_SymTet_LocSepGroupLink_TagAccVM2_Vol
  end interface

  interface
    subroutine F_getNormMoments_SymPyrmd_LocSepGroupLink_TagAccVM2_Vol(a_sym_pyramid, &
         a_localized_separator_group_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_SymPyrmd_LocSepGroupLink_TagAccVM2_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_SymPyrmd) :: a_sym_pyramid ! Pointer to SymPyrmd object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM2_Vol) :: a_moments_to_return ! Where TagAccVM2<Vol> is stored
    end subroutine F_getNormMoments_SymPyrmd_LocSepGroupLink_TagAccVM2_Vol
  end interface

  interface
    subroutine F_getNormMoments_SymTriPrism_LocSepGroupLink_TagAccVM2_Vol(a_sym_tri_prism, &
         a_localized_separator_group_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_SymTriPrism_LocSepGroupLink_TagAccVM2_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_SymTriPrism) :: a_sym_tri_prism ! Pointer to SymTriPrism object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM2_Vol) :: a_moments_to_return ! Where TagAccVM2<Vol> is stored
    end subroutine F_getNormMoments_SymTriPrism_LocSepGroupLink_TagAccVM2_Vol
  end interface

  interface
     subroutine F_getNormMoments_SymHex_LocSepGroupLink_TagAccVM2_Vol(a_sym_hex, a_localized_separator_group_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_SymHex_LocSepGroupLink_TagAccVM2_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_SymHex) :: a_sym_hex ! Pointer to SymHex object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM2_Vol) :: a_moments_to_return ! Where TagAccVM2<Vol> is stored
    end subroutine F_getNormMoments_SymHex_LocSepGroupLink_TagAccVM2_Vol
  end interface
  
  interface
     subroutine F_getNormMoments_SymTet_LocSepLink_TagAccVM_SepVol(a_sym_tet, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_SymTet_LocSepLink_TagAccVM_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_SymTet) :: a_sym_tet ! Pointer to SymTet object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVol) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_SymTet_LocSepLink_TagAccVM_SepVol
  end interface

  interface
    subroutine F_getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVol(a_sym_pyramid, a_localized_separator_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_SymPyrmd) :: a_sym_pyramid ! Pointer to SymPyrmd object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVol) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVol
  end interface

  interface
    subroutine F_getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVol(a_sym_tri_prism, a_localized_separator_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_SymTriPrism) :: a_sym_tri_prism ! Pointer to SymTriPrism object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVol) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVol
  end interface

  interface
     subroutine F_getNormMoments_SymHex_LocSepLink_TagAccVM_SepVol(a_sym_hex, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_SymHex_LocSepLink_TagAccVM_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_SymHex) :: a_sym_hex ! Pointer to SymHex object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVol) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_SymHex_LocSepLink_TagAccVM_SepVol
  end interface  

  interface
     subroutine F_getNormMoments_SymTet_LocSepLink_TagAccVM_SepVM(a_sym_tet, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_SymTet_LocSepLink_TagAccVM_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_SymTet) :: a_sym_tet ! Pointer to SymTet object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVM) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_SymTet_LocSepLink_TagAccVM_SepVM
  end interface

  interface
    subroutine F_getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVM(a_sym_pyramid, a_localized_separator_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_SymPyrmd) :: a_sym_pyramid ! Pointer to SymPyrmd object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVM) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVM
  end interface

  interface
    subroutine F_getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVM(a_sym_tri_prism, a_localized_separator_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_SymTriPrism) :: a_sym_tri_prism ! Pointer to SymTriPrism object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVM) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVM
  end interface

  interface
     subroutine F_getNormMoments_SymHex_LocSepLink_TagAccVM_SepVM(a_sym_hex, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_SymHex_LocSepLink_TagAccVM_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_SymHex) :: a_sym_hex ! Pointer to SymHex object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVM) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getNormMoments_SymHex_LocSepLink_TagAccVM_SepVM
  end interface
  
  interface
    subroutine F_getMoments_SymTet_LocSepLink_TagAccVM_SepVM(a_sym_tet, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getMoments_SymTet_LocSepLink_TagAccVM_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_SymTet) :: a_sym_tet ! Pointer to SymTet object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVM) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getMoments_SymTet_LocSepLink_TagAccVM_SepVM
  end interface

  interface
    subroutine F_getMoments_SymPyrmd_LocSepLink_TagAccVM_SepVM(a_sym_pyramid, a_localized_separator_link, &
         a_moments_to_return) &
    bind(C, name="c_getMoments_SymPyrmd_LocSepLink_TagAccVM_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_SymPyrmd) :: a_sym_pyramid ! Pointer to SymPyrmd object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVM) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getMoments_SymPyrmd_LocSepLink_TagAccVM_SepVM
  end interface

  interface
    subroutine F_getMoments_SymTriPrism_LocSepLink_TagAccVM_SepVM(a_sym_tri_prism, a_localized_separator_link, &
         a_moments_to_return) &
    bind(C, name="c_getMoments_SymTriPrism_LocSepLink_TagAccVM_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_SymTriPrism) :: a_sym_tri_prism ! Pointer to SymTriPrism object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVM) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getMoments_SymTriPrism_LocSepLink_TagAccVM_SepVM
  end interface

  interface
     subroutine F_getMoments_SymHex_LocSepLink_TagAccVM_SepVM(a_sym_hex, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getMoments_SymHex_LocSepLink_TagAccVM_SepVM")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_SymHex) :: a_sym_hex ! Pointer to SymHex object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_TagAccVM_SepVM) :: a_moments_to_return ! Where TagAccVM<SeparatedMoments<VM>> is stored
    end subroutine F_getMoments_SymHex_LocSepLink_TagAccVM_SepVM
  end interface
  
  interface
    subroutine F_getNormMoments_CapDod_LLLL_LocSepLink_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LLLL_LocSepLink_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LLLL) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVol) :: a_moments_to_return ! Where SepVol is stored
    end subroutine F_getNormMoments_CapDod_LLLL_LocSepLink_SepVol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LLLT_LocSepLink_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LLLT_LocSepLink_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LLLT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVol) :: a_moments_to_return ! Where SepVol is stored
    end subroutine F_getNormMoments_CapDod_LLLT_LocSepLink_SepVol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LTLT_LocSepLink_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LTLT_LocSepLink_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LTLT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVol) :: a_moments_to_return ! Where SepVol is stored
    end subroutine F_getNormMoments_CapDod_LTLT_LocSepLink_SepVol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LLTT_LocSepLink_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LLTT_LocSepLink_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LLTT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVol) :: a_moments_to_return ! Where SepVol is stored
    end subroutine F_getNormMoments_CapDod_LLTT_LocSepLink_SepVol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LTTT_LocSepLink_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LTTT_LocSepLink_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LTTT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVol) :: a_moments_to_return ! Where SepVol is stored
    end subroutine F_getNormMoments_CapDod_LTTT_LocSepLink_SepVol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_TTTT_LocSepLink_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_TTTT_LocSepLink_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_TTTT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVol) :: a_moments_to_return ! Where SepVol is stored
    end subroutine F_getNormMoments_CapDod_TTTT_LocSepLink_SepVol
  end interface    

  interface
    subroutine F_getNormMoments_CapOcta_LLL_LocSepLink_SepVol(a_CapOcta, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapOcta_LLL_LocSepLink_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapOcta_LLL) :: a_CapOcta ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVol) :: a_moments_to_return ! Where SepVol is stored
    end subroutine F_getNormMoments_CapOcta_LLL_LocSepLink_SepVol
  end interface      

  interface
    subroutine F_getNormMoments_CapOcta_LLT_LocSepLink_SepVol(a_CapOcta, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapOcta_LLT_LocSepLink_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapOcta_LLT) :: a_CapOcta ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVol) :: a_moments_to_return ! Where SepVol is stored
    end subroutine F_getNormMoments_CapOcta_LLT_LocSepLink_SepVol
  end interface

  interface
     subroutine F_getNormMoments_CapOcta_LTT_LocSepLink_SepVol(a_CapOcta, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapOcta_LTT_LocSepLink_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapOcta_LTT) :: a_CapOcta ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVol) :: a_moments_to_return ! Where SepVol is stored
    end subroutine F_getNormMoments_CapOcta_LTT_LocSepLink_SepVol
  end interface

  interface
     subroutine F_getNormMoments_CapOcta_TTT_LocSepLink_SepVol(a_CapOcta, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapOcta_TTT_LocSepLink_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapOcta_TTT) :: a_CapOcta ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVol) :: a_moments_to_return ! Where SepVol is stored
    end subroutine F_getNormMoments_CapOcta_TTT_LocSepLink_SepVol
  end interface    

   interface
     subroutine F_getNormMoments_Tet_LocSepLink_SepVol(a_tet, a_localized_separator_link, a_moments_to_return) &
    bind(C, name="c_getNormMoments_Tet_LocSepLink_SepVol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Tet) :: a_tet ! Pointer to CapDod object
      type(c_LocSepLink) :: a_localized_separator_link ! Pointer to LocSepLink object
      type(c_SepVol) :: a_moments_to_return ! Where SepVol is stored
    end subroutine F_getNormMoments_Tet_LocSepLink_SepVol
  end interface
  
  interface
    subroutine F_getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM2_Vol(a_CapDod, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM2_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LLLL) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM2_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM2_Vol
  end interface
 
  interface
    subroutine F_getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM2_Vol(a_CapDod, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM2_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LLLT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM2_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM2_Vol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM2_Vol(a_CapDod, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM2_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LTLT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM2_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM2_Vol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM2_Vol(a_CapDod, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM2_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LLTT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM2_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM2_Vol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM2_Vol(a_CapDod, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM2_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LTTT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM2_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM2_Vol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM2_Vol(a_CapDod, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM2_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_TTTT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM2_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM2_Vol
  end interface    

  interface
    subroutine F_getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM2_Vol(a_CapOcta, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM2_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapOcta_LLL) :: a_CapOcta ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM2_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM2_Vol
  end interface      

  interface
    subroutine F_getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM2_Vol(a_CapOcta, a_localized_separator_group_link, &
      a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM2_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapOcta_LLT) :: a_CapOcta ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM2_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM2_Vol
  end interface

  interface
    subroutine F_getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM2_Vol(a_CapOcta, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM2_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapOcta_LTT) :: a_CapOcta ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM2_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM2_Vol
  end interface

  interface
    subroutine F_getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM2_Vol(a_CapOcta, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM2_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapOcta_TTT) :: a_CapOcta ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM2_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM2_Vol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM_Vol(a_CapDod, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LLLL) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM_Vol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM_Vol(a_CapDod, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LLLT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM_Vol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM_Vol(a_CapDod, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LTLT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM_Vol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM_Vol(a_CapDod, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LLTT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM_Vol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM_Vol(a_CapDod, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_LTTT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM_Vol
  end interface

  interface
    subroutine F_getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM_Vol(a_CapDod, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapDod_TTTT) :: a_CapDod ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM_Vol
  end interface    

  interface
    subroutine F_getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM_Vol(a_CapOcta, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapOcta_LLL) :: a_CapOcta ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM_Vol
  end interface      

  interface
    subroutine F_getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM_Vol(a_CapOcta, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapOcta_LLT) :: a_CapOcta ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM_Vol
  end interface

  interface
    subroutine F_getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM_Vol(a_CapOcta, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapOcta_LTT) :: a_CapOcta ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM_Vol
  end interface

  interface
    subroutine F_getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM_Vol(a_CapOcta, a_localized_separator_group_link, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_CapOcta_TTT) :: a_CapOcta ! Pointer to CapDod object
      type(c_LocSepGroupLink) :: a_localized_separator_group_link ! Pointer to LocSepGroupLink object
      type(c_TagAccVM_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM_Vol
  end interface

  interface
    subroutine F_getNormMoments_Tet_PlanarSepPathGroup_TagAccVM_Vol(a_tet, a_planar_separator_path_group, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_Tet_PlanarSepPathGroup_TagAccVM_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Tet) :: a_tet ! Pointer to Tet object
      type(c_PlanarSepPathGroup) :: a_planar_separator_path_group ! Pointer to PlanarSepPathGroup object
      type(c_TagAccVM_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_Tet_PlanarSepPathGroup_TagAccVM_Vol
  end interface

  interface
    subroutine F_getNormMoments_Pyrmd_PlanarSepPathGroup_TagAccVM_Vol(a_pyramid, a_planar_separator_path_group, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_Pyrmd_PlanarSepPathGroup_TagAccVM_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Pyrmd) :: a_pyramid ! Pointer to Pyrmd object
      type(c_PlanarSepPathGroup) :: a_planar_separator_path_group ! Pointer to PlanarSepPathGroup object
      type(c_TagAccVM_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_Pyrmd_PlanarSepPathGroup_TagAccVM_Vol
  end interface

  interface
    subroutine F_getNormMoments_TriPrism_PlanarSepPathGroup_TagAccVM_Vol(a_tri_prism, a_planar_separator_path_group, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_TriPrism_PlanarSepPathGroup_TagAccVM_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_TriPrism) :: a_tri_prism ! Pointer to TriPrism object
      type(c_PlanarSepPathGroup) :: a_planar_separator_path_group ! Pointer to PlanarSepPathGroup object
      type(c_TagAccVM_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_TriPrism_PlanarSepPathGroup_TagAccVM_Vol
  end interface


  interface
    subroutine F_getNormMoments_Hex_PlanarSepPathGroup_TagAccVM_Vol(a_hexahedron, a_planar_separator_path_group, &
         a_moments_to_return) &
    bind(C, name="c_getNormMoments_Hex_PlanarSepPathGroup_TagAccVM_Vol")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_Hex) :: a_hexahedron ! Pointer to Hexahedron object
      type(c_PlanarSepPathGroup) :: a_planar_separator_path_group ! Pointer to PlanarSepPathGroup object
      type(c_TagAccVM_Vol) :: a_moments_to_return 
    end subroutine F_getNormMoments_Hex_PlanarSepPathGroup_TagAccVM_Vol
 end interface

  interface
     subroutine F_getNormMoments_SymTet_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return) &
          bind(C, name="c_getNormMoments_SymTet_LocSep_SepVol")
       use, intrinsic :: iso_c_binding
       import
       implicit none
       type(c_SymTet) :: a_poly
       type(c_LocSep) :: a_localized_separator
       type(C_SepVol) :: a_moments_to_return 
     end subroutine F_getNormMoments_SymTet_LocSep_SepVol
  end interface

  
  interface
     subroutine F_getNormMoments_SymPyrmd_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return) &
          bind(C, name="c_getNormMoments_SymPyrmd_LocSep_SepVol")
       use, intrinsic :: iso_c_binding
       import
       implicit none
       type(c_SymPyrmd) :: a_poly
       type(c_LocSep) :: a_localized_separator
       type(C_SepVol) :: a_moments_to_return 
     end subroutine F_getNormMoments_SymPyrmd_LocSep_SepVol
  end interface


  interface
     subroutine F_getNormMoments_SymTriPrism_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return) &
          bind(C, name="c_getNormMoments_SymTriPrism_LocSep_SepVol")
       use, intrinsic :: iso_c_binding
       import
       implicit none
       type(c_SymTriPrism) :: a_poly
       type(c_LocSep) :: a_localized_separator
       type(C_SepVol) :: a_moments_to_return 
     end subroutine F_getNormMoments_SymTriPrism_LocSep_SepVol
  end interface

  interface
     subroutine F_getNormMoments_SymHex_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return) &
          bind(C, name="c_getNormMoments_SymHex_LocSep_SepVol")
       use, intrinsic :: iso_c_binding
       import
       implicit none
       type(c_SymHex) :: a_poly
       type(c_LocSep) :: a_localized_separator
       type(C_SepVol) :: a_moments_to_return 
     end subroutine F_getNormMoments_SymHex_LocSep_SepVol
  end interface


  
  interface
     subroutine F_getNormMoments_CapOcta_LLL_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return) &
          bind(C, name="c_getNormMoments_CapOcta_LLL_LocSep_SepVol")
       use, intrinsic :: iso_c_binding
       import
       implicit none
       type(c_CapOcta_LLL) :: a_poly
       type(c_LocSep) :: a_localized_separator
       type(C_SepVol) :: a_moments_to_return 
     end subroutine F_getNormMoments_CapOcta_LLL_LocSep_SepVol
  end interface

  interface
     subroutine F_getNormMoments_CapOcta_LLT_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return) &
          bind(C, name="c_getNormMoments_CapOcta_LLT_LocSep_SepVol")
       use, intrinsic :: iso_c_binding
       import
       implicit none
       type(c_CapOcta_LLT) :: a_poly
       type(c_LocSep) :: a_localized_separator
       type(C_SepVol) :: a_moments_to_return 
     end subroutine F_getNormMoments_CapOcta_LLT_LocSep_SepVol
  end interface

  interface
     subroutine F_getNormMoments_CapOcta_LTT_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return) &
          bind(C, name="c_getNormMoments_CapOcta_LTT_LocSep_SepVol")
       use, intrinsic :: iso_c_binding
       import
       implicit none
       type(c_CapOcta_LTT) :: a_poly
       type(c_LocSep) :: a_localized_separator
       type(C_SepVol) :: a_moments_to_return 
     end subroutine F_getNormMoments_CapOcta_LTT_LocSep_SepVol
  end interface

  interface
     subroutine F_getNormMoments_CapOcta_TTT_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return) &
          bind(C, name="c_getNormMoments_CapOcta_TTT_LocSep_SepVol")
       use, intrinsic :: iso_c_binding
       import
       implicit none
       type(c_CapOcta_TTT) :: a_poly
       type(c_LocSep) :: a_localized_separator
       type(C_SepVol) :: a_moments_to_return 
     end subroutine F_getNormMoments_CapOcta_TTT_LocSep_SepVol
  end interface



    interface
     subroutine F_getNormMoments_CapDod_LLLL_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return) &
          bind(C, name="c_getNormMoments_CapDod_LLLL_LocSep_SepVol")
       use, intrinsic :: iso_c_binding
       import
       implicit none
       type(c_CapDod_LLLL) :: a_poly
       type(c_LocSep) :: a_localized_separator
       type(C_SepVol) :: a_moments_to_return 
     end subroutine F_getNormMoments_CapDod_LLLL_LocSep_SepVol
  end interface


    interface
     subroutine F_getNormMoments_CapDod_LLLT_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return) &
          bind(C, name="c_getNormMoments_CapDod_LLLT_LocSep_SepVol")
       use, intrinsic :: iso_c_binding
       import
       implicit none
       type(c_CapDod_LLLT) :: a_poly
       type(c_LocSep) :: a_localized_separator
       type(C_SepVol) :: a_moments_to_return 
     end subroutine F_getNormMoments_CapDod_LLLT_LocSep_SepVol
  end interface


    interface
     subroutine F_getNormMoments_CapDod_LTLT_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return) &
          bind(C, name="c_getNormMoments_CapDod_LTLT_LocSep_SepVol")
       use, intrinsic :: iso_c_binding
       import
       implicit none
       type(c_CapDod_LTLT) :: a_poly
       type(c_LocSep) :: a_localized_separator
       type(C_SepVol) :: a_moments_to_return 
     end subroutine F_getNormMoments_CapDod_LTLT_LocSep_SepVol
  end interface


    interface
     subroutine F_getNormMoments_CapDod_LLTT_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return) &
          bind(C, name="c_getNormMoments_CapDod_LLTT_LocSep_SepVol")
       use, intrinsic :: iso_c_binding
       import
       implicit none
       type(c_CapDod_LLTT) :: a_poly
       type(c_LocSep) :: a_localized_separator
       type(C_SepVol) :: a_moments_to_return 
     end subroutine F_getNormMoments_CapDod_LLTT_LocSep_SepVol
  end interface


    interface
     subroutine F_getNormMoments_CapDod_LTTT_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return) &
          bind(C, name="c_getNormMoments_CapDod_LTTT_LocSep_SepVol")
       use, intrinsic :: iso_c_binding
       import
       implicit none
       type(c_CapDod_LTTT) :: a_poly
       type(c_LocSep) :: a_localized_separator
       type(C_SepVol) :: a_moments_to_return 
     end subroutine F_getNormMoments_CapDod_LTTT_LocSep_SepVol
  end interface


    interface
     subroutine F_getNormMoments_CapDod_TTTT_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return) &
          bind(C, name="c_getNormMoments_CapDod_TTTT_LocSep_SepVol")
       use, intrinsic :: iso_c_binding
       import
       implicit none
       type(c_CapDod_TTTT) :: a_poly
       type(c_LocSep) :: a_localized_separator
       type(C_SepVol) :: a_moments_to_return 
     end subroutine F_getNormMoments_CapDod_TTTT_LocSep_SepVol
  end interface     

  
contains

  subroutine getMoments_setMethod(a_cutting_method)
    use, intrinsic :: iso_c_binding
    implicit none
    integer(IRL_UnsignedIndex_t), intent(in) :: a_cutting_method
    call F_getMoments_setMethod(a_cutting_method)
  end subroutine getMoments_setMethod

  subroutine getNormMoments_Dod_LocSepLink_SepVM(a_Dod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Dod_type), intent(in) :: a_Dod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVM_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Dod_LocSepLink_SepVM &
          (a_dod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_Dod_LocSepLink_SepVM
  
  subroutine getNormMoments_Dod_LocSepLink_Vol(a_Dod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Dod_type), intent(in) :: a_Dod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      real(IRL_double), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Dod_LocSepLink_Vol &
          (a_dod%c_object, a_localized_separator_link%c_object, a_moments_to_return)

  end subroutine getNormMoments_Dod_LocSepLink_Vol
  
  subroutine getNormMoments_TriPrism_LocSepLink_Vol(a_triangular_prism, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(TriPrism_type), intent(in) :: a_triangular_prism
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      real(IRL_double), intent(inout) :: a_moments_to_return

      call F_getNormMoments_TriPrism_LocSepLink_Vol &
          (a_triangular_prism%c_object, a_localized_separator_link%c_object, a_moments_to_return)

  end subroutine getNormMoments_TriPrism_LocSepLink_Vol
  
  subroutine getNormMoments_Octa_LocSepLink_Vol(a_octahedron, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Octa_type), intent(in) :: a_octahedron
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      real(IRL_double), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Octa_LocSepLink_Vol &
          (a_octahedron%c_object, a_localized_separator_link%c_object, a_moments_to_return)

  end subroutine getNormMoments_Octa_LocSepLink_Vol

  subroutine getNormMoments_CapDod_LocSepLink_SepVM(a_Capped_Dod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_type), intent(in) :: a_Capped_Dod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVM_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LocSepLink_SepVM &
          (a_capped_dod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LocSepLink_SepVM

  subroutine getNormMoments_CapDod_d3_LocSepLink_SepVM_d3(a_Capped_Dod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_d3_type), intent(in) :: a_Capped_Dod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVM_d3_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_d3_LocSepLink_SepVM_d3 &
          (a_capped_dod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_d3_LocSepLink_SepVM_d3

  subroutine getNormMoments_Poly24_LocSepLink_SepVM(a_polyhedron_24, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Poly24_type), intent(in) :: a_polyhedron_24
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVM_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Poly24_LocSepLink_SepVM &
          (a_polyhedron_24%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_Poly24_LocSepLink_SepVM

  subroutine getNormMoments_Poly24_d3_LocSepLink_SepVM_d3(a_polyhedron_24, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Poly24_d3_type), intent(in) :: a_polyhedron_24
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVM_d3_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Poly24_d3_LocSepLink_SepVM_d3 &
          (a_polyhedron_24%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_Poly24_d3_LocSepLink_SepVM_d3

  subroutine getMoments_CapDod_LocSepLink_SepVM(a_Capped_Dod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_type), intent(in) :: a_Capped_Dod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVM_type), intent(inout) :: a_moments_to_return

      call F_getMoments_CapDod_LocSepLink_SepVM &
          (a_capped_dod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getMoments_CapDod_LocSepLink_SepVM

  subroutine getMoments_Dod_LocSepLink_SepVM(a_Dod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Dod_type), intent(in) :: a_Dod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVM_type), intent(inout) :: a_moments_to_return

      call F_getMoments_Dod_LocSepLink_SepVM &
          (a_dod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getMoments_Dod_LocSepLink_SepVM

  subroutine getMoments_Poly24_LocSepLink_SepVM(a_polyhedron_24, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Poly24_type), intent(in) :: a_polyhedron_24
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVM_type), intent(inout) :: a_moments_to_return

      call F_getMoments_Poly24_LocSepLink_SepVM &
          (a_polyhedron_24%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getMoments_Poly24_LocSepLink_SepVM

  subroutine getNormMoments_Tet_LocSepLink_SepVM(a_tet, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Tet_type), intent(in) :: a_tet
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVM_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Tet_LocSepLink_SepVM &
          (a_tet%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_Tet_LocSepLink_SepVM

  subroutine getNormMoments_RectCub_PlanarSep_Vol(a_rectangular_cuboid, a_planar_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(RectCub_type), intent(in) :: a_rectangular_cuboid
      type(PlanarSep_type), intent(in) :: a_planar_separator
      real(IRL_double), intent(inout) :: a_moments_to_return

      call F_getNormMoments_RectCub_PlanarSep_Vol &
          (a_rectangular_cuboid%c_object, a_planar_separator%c_object, a_moments_to_return)

  end subroutine getNormMoments_RectCub_PlanarSep_Vol

  subroutine getNormMoments_Tet_PlanarSep_Vol(a_tet, a_planar_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Tet_type), intent(in) :: a_tet
      type(PlanarSep_type), intent(in) :: a_planar_separator
      real(IRL_double), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Tet_PlanarSep_Vol &
          (a_tet%c_object, a_planar_separator%c_object, a_moments_to_return)

  end subroutine getNormMoments_Tet_PlanarSep_Vol

  subroutine getNormMoments_TriPrism_PlanarSep_Vol(a_tri_prism, a_planar_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(TriPrism_type), intent(in) :: a_tri_prism
      type(PlanarSep_type), intent(in) :: a_planar_separator
      real(IRL_double), intent(inout) :: a_moments_to_return

      call F_getNormMoments_TriPrism_PlanarSep_Vol &
          (a_tri_prism%c_object, a_planar_separator%c_object, a_moments_to_return)

  end subroutine getNormMoments_TriPrism_PlanarSep_Vol

  subroutine getNormMoments_Pyrmd_PlanarSep_Vol(a_pyramid, a_planar_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Pyrmd_type), intent(in) :: a_pyramid
      type(PlanarSep_type), intent(in) :: a_planar_separator
      real(IRL_double), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Pyrmd_PlanarSep_Vol &
          (a_pyramid%c_object, a_planar_separator%c_object, a_moments_to_return)

  end subroutine getNormMoments_Pyrmd_PlanarSep_Vol
  
  subroutine getNormMoments_Hex_PlanarSep_Vol(a_hexahedron, a_planar_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Hex_type), intent(in) :: a_hexahedron
      type(PlanarSep_type), intent(in) :: a_planar_separator
      real(IRL_double), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Hex_PlanarSep_Vol &
          (a_hexahedron%c_object, a_planar_separator%c_object, a_moments_to_return)

  end subroutine getNormMoments_Hex_PlanarSep_Vol
  
  subroutine getNormMoments_Hex_PlanarSep_SepVM(a_hexahedron, a_planar_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Hex_type), intent(in) :: a_hexahedron
      type(PlanarSep_type), intent(in) :: a_planar_separator
      type(SepvM_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Hex_PlanarSep_SepvM &
          (a_hexahedron%c_object, a_planar_separator%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_Hex_PlanarSep_SepVM

  subroutine getNormMoments_Dod_PlanarSep_SepVM(a_Dod, a_planar_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Dod_type), intent(in) :: a_Dod
      type(PlanarSep_type), intent(in) :: a_planar_separator
      type(SepVM_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Dod_PlanarSep_SepVM &
          (a_dod%c_object, a_planar_separator%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_Dod_PlanarSep_SepVM

  subroutine getNormMoments_CapDod_LocSepLink_TagAccVM_SepVM(a_Capped_Dod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_type), intent(in) :: a_Capped_Dod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVM_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LocSepLink_TagAccVM_SepVM &
          (a_capped_dod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LocSepLink_TagAccVM_SepVM

  subroutine getNormMoments_Dod_LocSepLink_TagAccVM_SepVM(a_Dod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Dod_type), intent(in) :: a_Dod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVM_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Dod_LocSepLink_TagAccVM_SepVM &
          (a_dod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_Dod_LocSepLink_TagAccVM_SepVM

  subroutine getMoments_Dod_LocSepLink_TagAccVM_SepVM(a_Dod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Dod_type), intent(in) :: a_Dod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVM_type), intent(inout) :: a_moments_to_return

      call F_getMoments_Dod_LocSepLink_TagAccVM_SepVM &
          (a_dod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getMoments_Dod_LocSepLink_TagAccVM_SepVM

  subroutine getNormMoments_Dod_LocSepLink_TagAccVM_SepVol(a_Dod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Dod_type), intent(in) :: a_Dod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Dod_LocSepLink_TagAccVM_SepVol &
          (a_dod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_Dod_LocSepLink_TagAccVM_SepVol

  subroutine getNormMoments_Tet_LocSepLink_TagAccVM_SepVol(a_tet, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Tet_type), intent(in) :: a_tet
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Tet_LocSepLink_TagAccVM_SepVol &
          (a_tet%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_Tet_LocSepLink_TagAccVM_SepVol
  
  subroutine getMoments_Octa_LocSepLink_TagAccVM_SepVM(a_octahedron, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Octa_type), intent(in) :: a_octahedron
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVM_type), intent(inout) :: a_moments_to_return

      call F_getMoments_Octa_LocSepLink_TagAccVM_SepVM &
          (a_octahedron%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getMoments_Octa_LocSepLink_TagAccVM_SepVM

  subroutine getNormMoments_Octa_LocSepLink_TagAccVM_SepVol(a_octahedron, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Octa_type), intent(in) :: a_octahedron
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Octa_LocSepLink_TagAccVM_SepVol &
          (a_octahedron%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_Octa_LocSepLink_TagAccVM_SepVol

  
  subroutine getNormMoments_RectCub_PlanarSep_SepVM(a_rectangular_cuboid, a_planar_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(RectCub_type), intent(in) :: a_rectangular_cuboid
      type(PlanarSep_type), intent(in) :: a_planar_separator
      type(SepVM_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_RectCub_PlanarSep_SepVM &
          (a_rectangular_cuboid%c_object, a_planar_separator%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_RectCub_PlanarSep_SepVM

  subroutine getNormMoments_Tri_LocLink_TagAccVM_VM(a_tri, a_localizer_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Tri_type), intent(in) :: a_tri
      type(LocLink_type), intent(in) :: a_localizer_link
      type(TagAccVM_VM_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Tri_LocLink_TagAccVM_VM &
          (a_tri%c_object, a_localizer_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_Tri_LocLink_TagAccVM_VM

  subroutine getNormMoments_Tri_PlanarLoc_Vol(a_tri, a_planar_localizer, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Tri_type), intent(in) :: a_tri
      type(PlanarLoc_type), intent(in) :: a_planar_localizer
      real(IRL_double), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Tri_PlanarLoc_Vol &
          (a_tri%c_object, a_planar_localizer%c_object, a_moments_to_return)

  end subroutine getNormMoments_Tri_PlanarLoc_Vol

  subroutine getNormMoments_Poly_PlanarLoc_Vol(a_polygon, a_planar_localizer, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Poly_type), intent(in) :: a_polygon
      type(PlanarLoc_type), intent(in) :: a_planar_localizer
      real(IRL_double), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Poly_PlanarLoc_Vol &
          (a_polygon%c_object, a_planar_localizer%c_object, a_moments_to_return)

  end subroutine getNormMoments_Poly_PlanarLoc_Vol

  subroutine getMoments_TRI_LocLink_TagAccListVM_VMAN(a_tri, a_localizer_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Tri_type), intent(in) :: a_tri
      type(LocLink_type), intent(in) :: a_localizer_link
      type(TagAccListVM_VMAN_type), intent(inout) :: a_moments_to_return

      call F_getMoments_TRI_LocLink_TagAccListVM_VMAN &
          (a_tri%c_object, a_localizer_link%c_object, a_moments_to_return%c_object)

  end subroutine getMoments_TRI_LocLink_TagAccListVM_VMAN

  subroutine getNormMoments_CapDod_LLLL_LocSepLink_TagAccVM_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LLLL_type), intent(in) :: a_CapDod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LLLL_LocSepLink_TagAccVM_SepVol &
          (a_CapDod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LLLL_LocSepLink_TagAccVM_SepVol    

  subroutine getNormMoments_CapDod_LLLT_LocSepLink_TagAccVM_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LLLT_type), intent(in) :: a_CapDod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LLLT_LocSepLink_TagAccVM_SepVol &
          (a_CapDod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LLLT_LocSepLink_TagAccVM_SepVol

  subroutine getNormMoments_CapDod_LTLT_LocSepLink_TagAccVM_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LTLT_type), intent(in) :: a_CapDod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LTLT_LocSepLink_TagAccVM_SepVol &
          (a_CapDod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LTLT_LocSepLink_TagAccVM_SepVol

  subroutine getNormMoments_CapDod_LLTT_LocSepLink_TagAccVM_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LLTT_type), intent(in) :: a_CapDod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LLTT_LocSepLink_TagAccVM_SepVol &
          (a_CapDod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LLTT_LocSepLink_TagAccVM_SepVol

  subroutine getNormMoments_CapDod_LTTT_LocSepLink_TagAccVM_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LTTT_type), intent(in) :: a_CapDod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LTTT_LocSepLink_TagAccVM_SepVol &
          (a_CapDod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LTTT_LocSepLink_TagAccVM_SepVol

  subroutine getNormMoments_CapDod_TTTT_LocSepLink_TagAccVM_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_TTTT_type), intent(in) :: a_CapDod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_TTTT_LocSepLink_TagAccVM_SepVol &
          (a_CapDod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_TTTT_LocSepLink_TagAccVM_SepVol

  subroutine getNormMoments_CapOcta_LLL_LocSepLink_TagAccVM_SepVol(a_CapOcta, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_LLL_type), intent(in) :: a_CapOcta
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_LLL_LocSepLink_TagAccVM_SepVol &
          (a_CapOcta%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapOcta_LLL_LocSepLink_TagAccVM_SepVol    

  subroutine getNormMoments_CapOcta_LLT_LocSepLink_TagAccVM_SepVol(a_CapOcta, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_LLT_type), intent(in) :: a_CapOcta
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_LLT_LocSepLink_TagAccVM_SepVol &
          (a_CapOcta%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapOcta_LLT_LocSepLink_TagAccVM_SepVol

  subroutine getNormMoments_CapOcta_LTT_LocSepLink_TagAccVM_SepVol(a_CapOcta, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_LTT_type), intent(in) :: a_CapOcta
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_LTT_LocSepLink_TagAccVM_SepVol &
          (a_CapOcta%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapOcta_LTT_LocSepLink_TagAccVM_SepVol

  subroutine getNormMoments_CapOcta_TTT_LocSepLink_TagAccVM_SepVol(a_CapOcta, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_TTT_type), intent(in) :: a_CapOcta
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_TTT_LocSepLink_TagAccVM_SepVol &
          (a_CapOcta%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapOcta_TTT_LocSepLink_TagAccVM_SepVol        

  subroutine getNormMoments_SymTet_LocSepGroupLink_TagAccVM2_Vol(a_sym_tet, a_localized_separator_group_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymTet_type), intent(in) :: a_sym_tet
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM2_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_SymTet_LocSepGroupLink_TagAccVM2_Vol &
          (a_sym_tet%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_SymTet_LocSepGroupLink_TagAccVM2_Vol

  subroutine getNormMoments_SymPyrmd_LocSepGroupLink_TagAccVM2_Vol(a_sym_pyramid, a_localized_separator_group_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymPyrmd_type), intent(in) :: a_sym_pyramid
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM2_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_SymPyrmd_LocSepGroupLink_TagAccVM2_Vol &
          (a_sym_pyramid%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_SymPyrmd_LocSepGroupLink_TagAccVM2_Vol

  subroutine getNormMoments_SymTriPrism_LocSepGroupLink_TagAccVM2_Vol(a_sym_tri_prism, &
       a_localized_separator_group_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymTriPrism_type), intent(in) :: a_sym_tri_prism
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM2_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_SymTriPrism_LocSepGroupLink_TagAccVM2_Vol &
          (a_sym_tri_prism%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_SymTriPrism_LocSepGroupLink_TagAccVM2_Vol

  subroutine getNormMoments_SymHex_LocSepGroupLink_TagAccVM2_Vol(a_sym_hex, a_localized_separator_group_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymHex_type), intent(in) :: a_sym_hex
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM2_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_SymHex_LocSepGroupLink_TagAccVM2_Vol &
          (a_sym_hex%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_SymHex_LocSepGroupLink_TagAccVM2_Vol
  
  subroutine getNormMoments_SymTet_LocSepLink_TagAccVM_SepVol(a_sym_tet, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymTet_type), intent(in) :: a_sym_tet
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_SymTet_LocSepLink_TagAccVM_SepVol &
          (a_sym_tet%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_SymTet_LocSepLink_TagAccVM_SepVol

  subroutine getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVol(a_sym_pyramid, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymPyrmd_type), intent(in) :: a_sym_pyramid
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVol &
          (a_sym_pyramid%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVol

  subroutine getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVol(a_sym_tri_prism, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymTriPrism_type), intent(in) :: a_sym_tri_prism
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVol &
          (a_sym_tri_prism%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVol

  subroutine getNormMoments_SymHex_LocSepLink_TagAccVM_SepVol(a_sym_hex, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymHex_type), intent(in) :: a_sym_hex
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_SymHex_LocSepLink_TagAccVM_SepVol &
          (a_sym_hex%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_SymHex_LocSepLink_TagAccVM_SepVol

  subroutine getNormMoments_SymTet_LocSepLink_TagAccVM_SepVM(a_sym_tet, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymTet_type), intent(in) :: a_sym_tet
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVM_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_SymTet_LocSepLink_TagAccVM_SepVM &
          (a_sym_tet%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_SymTet_LocSepLink_TagAccVM_SepVM

  subroutine getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVM(a_sym_pyramid, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymPyrmd_type), intent(in) :: a_sym_pyramid
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVM_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVM &
          (a_sym_pyramid%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVM

  subroutine getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVM(a_sym_tri_prism, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymTriPrism_type), intent(in) :: a_sym_tri_prism
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVM_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVM &
          (a_sym_tri_prism%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVM

  subroutine getNormMoments_SymHex_LocSepLink_TagAccVM_SepVM(a_sym_hex, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymHex_type), intent(in) :: a_sym_hex
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVM_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_SymHex_LocSepLink_TagAccVM_SepVM &
          (a_sym_hex%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_SymHex_LocSepLink_TagAccVM_SepVM

  subroutine getMoments_SymTet_LocSepLink_TagAccVM_SepVM(a_sym_tet, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymTet_type), intent(in) :: a_sym_tet
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVM_type), intent(inout) :: a_moments_to_return

      call F_getMoments_SymTet_LocSepLink_TagAccVM_SepVM &
          (a_sym_tet%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getMoments_SymTet_LocSepLink_TagAccVM_SepVM

  subroutine getMoments_SymPyrmd_LocSepLink_TagAccVM_SepVM(a_sym_pyramid, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymPyrmd_type), intent(in) :: a_sym_pyramid
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVM_type), intent(inout) :: a_moments_to_return

      call F_getMoments_SymPyrmd_LocSepLink_TagAccVM_SepVM &
          (a_sym_pyramid%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getMoments_SymPyrmd_LocSepLink_TagAccVM_SepVM

  subroutine getMoments_SymTriPrism_LocSepLink_TagAccVM_SepVM(a_sym_tri_prism, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymTriPrism_type), intent(in) :: a_sym_tri_prism
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVM_type), intent(inout) :: a_moments_to_return

      call F_getMoments_SymTriPrism_LocSepLink_TagAccVM_SepVM &
          (a_sym_tri_prism%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getMoments_SymTriPrism_LocSepLink_TagAccVM_SepVM

  subroutine getMoments_SymHex_LocSepLink_TagAccVM_SepVM(a_sym_hex, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymHex_type), intent(in) :: a_sym_hex
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(TagAccVM_SepVM_type), intent(inout) :: a_moments_to_return

      call F_getMoments_SymHex_LocSepLink_TagAccVM_SepVM &
          (a_sym_hex%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getMoments_SymHex_LocSepLink_TagAccVM_SepVM
  
  subroutine getNormMoments_CapDod_LLLL_LocSepLink_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LLLL_type), intent(in) :: a_CapDod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LLLL_LocSepLink_SepVol &
          (a_CapDod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LLLL_LocSepLink_SepVol    

  subroutine getNormMoments_CapDod_LLLT_LocSepLink_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LLLT_type), intent(in) :: a_CapDod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LLLT_LocSepLink_SepVol &
          (a_CapDod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LLLT_LocSepLink_SepVol

  subroutine getNormMoments_CapDod_LTLT_LocSepLink_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LTLT_type), intent(in) :: a_CapDod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LTLT_LocSepLink_SepVol &
          (a_CapDod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LTLT_LocSepLink_SepVol

  subroutine getNormMoments_CapDod_LLTT_LocSepLink_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LLTT_type), intent(in) :: a_CapDod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LLTT_LocSepLink_SepVol &
          (a_CapDod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LLTT_LocSepLink_SepVol

  subroutine getNormMoments_CapDod_LTTT_LocSepLink_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LTTT_type), intent(in) :: a_CapDod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LTTT_LocSepLink_SepVol &
          (a_CapDod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LTTT_LocSepLink_SepVol

  subroutine getNormMoments_CapDod_TTTT_LocSepLink_SepVol(a_CapDod, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_TTTT_type), intent(in) :: a_CapDod
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_TTTT_LocSepLink_SepVol &
          (a_CapDod%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_TTTT_LocSepLink_SepVol

  subroutine getNormMoments_CapOcta_LLL_LocSepLink_SepVol(a_CapOcta, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_LLL_type), intent(in) :: a_CapOcta
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_LLL_LocSepLink_SepVol &
          (a_CapOcta%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapOcta_LLL_LocSepLink_SepVol    

  subroutine getNormMoments_CapOcta_LLT_LocSepLink_SepVol(a_CapOcta, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_LLT_type), intent(in) :: a_CapOcta
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_LLT_LocSepLink_SepVol &
          (a_CapOcta%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapOcta_LLT_LocSepLink_SepVol

  subroutine getNormMoments_CapOcta_LTT_LocSepLink_SepVol(a_CapOcta, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_LTT_type), intent(in) :: a_CapOcta
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_LTT_LocSepLink_SepVol &
          (a_CapOcta%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapOcta_LTT_LocSepLink_SepVol

  subroutine getNormMoments_CapOcta_TTT_LocSepLink_SepVol(a_CapOcta, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_TTT_type), intent(in) :: a_CapOcta
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_TTT_LocSepLink_SepVol &
          (a_CapOcta%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapOcta_TTT_LocSepLink_SepVol  

  subroutine getNormMoments_Tet_LocSepLink_SepVol(a_tet, a_localized_separator_link, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Tet_type), intent(in) :: a_tet
      type(LocSepLink_type), intent(in) :: a_localized_separator_link
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Tet_LocSepLink_SepVol &
          (a_tet%c_object, a_localized_separator_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_Tet_LocSepLink_SepVol  

  subroutine getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM2_Vol(a_CapDod, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LLLL_type), intent(in) :: a_CapDod
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM2_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM2_Vol &
          (a_CapDod%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM2_Vol    

  subroutine getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM2_Vol(a_CapDod, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LLLT_type), intent(in) :: a_CapDod
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM2_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM2_Vol &
          (a_CapDod%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM2_Vol

  subroutine getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM2_Vol(a_CapDod, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LTLT_type), intent(in) :: a_CapDod
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM2_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM2_Vol &
          (a_CapDod%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM2_Vol

  subroutine getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM2_Vol(a_CapDod, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LLTT_type), intent(in) :: a_CapDod
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM2_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM2_Vol &
          (a_CapDod%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM2_Vol

  subroutine getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM2_Vol(a_CapDod, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LTTT_type), intent(in) :: a_CapDod
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM2_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM2_Vol &
          (a_CapDod%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM2_Vol

  subroutine getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM2_Vol(a_CapDod, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_TTTT_type), intent(in) :: a_CapDod
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM2_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM2_Vol &
          (a_CapDod%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM2_Vol

  subroutine getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM2_Vol(a_CapOcta, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_LLL_type), intent(in) :: a_CapOcta
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM2_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM2_Vol &
          (a_CapOcta%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM2_Vol    

  subroutine getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM2_Vol(a_CapOcta, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_LLT_type), intent(in) :: a_CapOcta
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM2_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM2_Vol &
          (a_CapOcta%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM2_Vol

  subroutine getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM2_Vol(a_CapOcta, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_LTT_type), intent(in) :: a_CapOcta
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM2_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM2_Vol &
          (a_CapOcta%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM2_Vol

  subroutine getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM2_Vol(a_CapOcta, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_TTT_type), intent(in) :: a_CapOcta
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM2_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM2_Vol &
          (a_CapOcta%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM2_Vol        

  subroutine getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM_Vol(a_CapDod, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LLLL_type), intent(in) :: a_CapDod
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM_Vol &
          (a_CapDod%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM_Vol    

  subroutine getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM_Vol(a_CapDod, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LLLT_type), intent(in) :: a_CapDod
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM_Vol &
          (a_CapDod%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM_Vol

  subroutine getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM_Vol(a_CapDod, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LTLT_type), intent(in) :: a_CapDod
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM_Vol &
          (a_CapDod%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM_Vol

  subroutine getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM_Vol(a_CapDod, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LLTT_type), intent(in) :: a_CapDod
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM_Vol &
          (a_CapDod%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM_Vol

  subroutine getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM_Vol(a_CapDod, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LTTT_type), intent(in) :: a_CapDod
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM_Vol &
          (a_CapDod%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM_Vol

  subroutine getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM_Vol(a_CapDod, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_TTTT_type), intent(in) :: a_CapDod
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM_Vol &
          (a_CapDod%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM_Vol

  subroutine getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM_Vol(a_CapOcta, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_LLL_type), intent(in) :: a_CapOcta
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM_Vol &
          (a_CapOcta%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM_Vol    

  subroutine getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM_Vol(a_CapOcta, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_LLT_type), intent(in) :: a_CapOcta
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM_Vol &
          (a_CapOcta%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM_Vol

  subroutine getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM_Vol(a_CapOcta, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_LTT_type), intent(in) :: a_CapOcta
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM_Vol &
          (a_CapOcta%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM_Vol

  subroutine getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM_Vol(a_CapOcta, a_localized_separator_group_link, &
       a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_TTT_type), intent(in) :: a_CapOcta
      type(LocSepGroupLink_type), intent(in) :: a_localized_separator_group_link
      type(TagAccVM_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM_Vol &
          (a_CapOcta%c_object, a_localized_separator_group_link%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM_Vol  

  subroutine getNormMoments_Tet_PlanarSepPathGroup_TagAccVM_Vol(a_tet, a_planar_separator_path_group, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Tet_type), intent(in) :: a_tet
      type(PlanarSepPathGroup_type), intent(in) :: a_planar_separator_path_group
      type(TagAccVM_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Tet_PlanarSepPathGroup_TagAccVM_Vol &
          (a_tet%c_object, a_planar_separator_path_group%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_Tet_PlanarSepPathGroup_TagAccVM_Vol

  subroutine getNormMoments_Pyrmd_PlanarSepPathGroup_TagAccVM_Vol(a_pyramid, a_planar_separator_path_group, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Pyrmd_type), intent(in) :: a_pyramid
      type(PlanarSepPathGroup_type), intent(in) :: a_planar_separator_path_group
      type(TagAccVM_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Pyrmd_PlanarSepPathGroup_TagAccVM_Vol &
          (a_pyramid%c_object, a_planar_separator_path_group%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_Pyrmd_PlanarSepPathGroup_TagAccVM_Vol

  subroutine getNormMoments_TriPrism_PlanarSepPathGroup_TagAccVM_Vol(a_tri_prism, a_planar_separator_path_group, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(TriPrism_type), intent(in) :: a_tri_prism
      type(PlanarSepPathGroup_type), intent(in) :: a_planar_separator_path_group
      type(TagAccVM_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_TriPrism_PlanarSepPathGroup_TagAccVM_Vol &
          (a_tri_prism%c_object, a_planar_separator_path_group%c_object, a_moments_to_return%c_object)

  end subroutine getNormMoments_TriPrism_PlanarSepPathGroup_TagAccVM_Vol


  subroutine getNormMoments_Hex_PlanarSepPathGroup_TagAccVM_Vol(a_hexahedron, a_planar_separator_path_group, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(Hex_type), intent(in) :: a_hexahedron
      type(PlanarSepPathGroup_type), intent(in) :: a_planar_separator_path_group
      type(TagAccVM_Vol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_Hex_PlanarSepPathGroup_TagAccVM_Vol &
          (a_hexahedron%c_object, a_planar_separator_path_group%c_object, a_moments_to_return%c_object)

    end subroutine getNormMoments_Hex_PlanarSepPathGroup_TagAccVM_Vol


  subroutine getNormMoments_SymTet_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymTet_type), intent(in) :: a_poly
      type(LocSep_type), intent(in) :: a_localized_separator
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_SymTet_LocSep_SepVol &
          (a_poly%c_object, a_localized_separator%c_object, a_moments_to_return%c_object)

    end subroutine getNormMoments_SymTet_LocSep_SepVol

  subroutine getNormMoments_SymPyrmd_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymPyrmd_type), intent(in) :: a_poly
      type(LocSep_type), intent(in) :: a_localized_separator
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_SymPyrmd_LocSep_SepVol &
          (a_poly%c_object, a_localized_separator%c_object, a_moments_to_return%c_object)

    end subroutine getNormMoments_SymPyrmd_LocSep_SepVol

  subroutine getNormMoments_SymTriPrism_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymTriPrism_type), intent(in) :: a_poly
      type(LocSep_type), intent(in) :: a_localized_separator
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_SymTriPrism_LocSep_SepVol &
          (a_poly%c_object, a_localized_separator%c_object, a_moments_to_return%c_object)
      
    end subroutine getNormMoments_SymTriPrism_LocSep_SepVol

  subroutine getNormMoments_SymHex_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(SymHex_type), intent(in) :: a_poly
      type(LocSep_type), intent(in) :: a_localized_separator
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_SymHex_LocSep_SepVol &
          (a_poly%c_object, a_localized_separator%c_object, a_moments_to_return%c_object)

    end subroutine getNormMoments_SymHex_LocSep_SepVol

    
  subroutine getNormMoments_CapOcta_LLL_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_LLL_type), intent(in) :: a_poly
      type(LocSep_type), intent(in) :: a_localized_separator
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_LLL_LocSep_SepVol &
          (a_poly%c_object, a_localized_separator%c_object, a_moments_to_return%c_object)

    end subroutine getNormMoments_CapOcta_LLL_LocSep_SepVol

  subroutine getNormMoments_CapOcta_LLT_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_LLT_type), intent(in) :: a_poly
      type(LocSep_type), intent(in) :: a_localized_separator
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_LLT_LocSep_SepVol &
          (a_poly%c_object, a_localized_separator%c_object, a_moments_to_return%c_object)

    end subroutine getNormMoments_CapOcta_LLT_LocSep_SepVol
 
  subroutine getNormMoments_CapOcta_LTT_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_LTT_type), intent(in) :: a_poly
      type(LocSep_type), intent(in) :: a_localized_separator
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_LTT_LocSep_SepVol &
          (a_poly%c_object, a_localized_separator%c_object, a_moments_to_return%c_object)

    end subroutine getNormMoments_CapOcta_LTT_LocSep_SepVol

  subroutine getNormMoments_CapOcta_TTT_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapOcta_TTT_type), intent(in) :: a_poly
      type(LocSep_type), intent(in) :: a_localized_separator
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapOcta_TTT_LocSep_SepVol &
          (a_poly%c_object, a_localized_separator%c_object, a_moments_to_return%c_object)

    end subroutine getNormMoments_CapOcta_TTT_LocSep_SepVol


  subroutine getNormMoments_CapDod_LLLL_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LLLL_type), intent(in) :: a_poly
      type(LocSep_type), intent(in) :: a_localized_separator
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LLLL_LocSep_SepVol &
          (a_poly%c_object, a_localized_separator%c_object, a_moments_to_return%c_object)

    end subroutine getNormMoments_CapDod_LLLL_LocSep_SepVol

  subroutine getNormMoments_CapDod_LLLT_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LLLT_type), intent(in) :: a_poly
      type(LocSep_type), intent(in) :: a_localized_separator
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LLLT_LocSep_SepVol &
          (a_poly%c_object, a_localized_separator%c_object, a_moments_to_return%c_object)

    end subroutine getNormMoments_CapDod_LLLT_LocSep_SepVol

  subroutine getNormMoments_CapDod_LTLT_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LTLT_type), intent(in) :: a_poly
      type(LocSep_type), intent(in) :: a_localized_separator
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LTLT_LocSep_SepVol &
          (a_poly%c_object, a_localized_separator%c_object, a_moments_to_return%c_object)

    end subroutine getNormMoments_CapDod_LTLT_LocSep_SepVol

  subroutine getNormMoments_CapDod_LLTT_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LLTT_type), intent(in) :: a_poly
      type(LocSep_type), intent(in) :: a_localized_separator
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LLTT_LocSep_SepVol &
          (a_poly%c_object, a_localized_separator%c_object, a_moments_to_return%c_object)

    end subroutine getNormMoments_CapDod_LLTT_LocSep_SepVol

  subroutine getNormMoments_CapDod_LTTT_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_LTTT_type), intent(in) :: a_poly
      type(LocSep_type), intent(in) :: a_localized_separator
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_LTTT_LocSep_SepVol &
          (a_poly%c_object, a_localized_separator%c_object, a_moments_to_return%c_object)

    end subroutine getNormMoments_CapDod_LTTT_LocSep_SepVol

  subroutine getNormMoments_CapDod_TTTT_LocSep_SepVol(a_poly, a_localized_separator, a_moments_to_return)
    use, intrinsic :: iso_c_binding
    implicit none
      type(CapDod_TTTT_type), intent(in) :: a_poly
      type(LocSep_type), intent(in) :: a_localized_separator
      type(SepVol_type), intent(inout) :: a_moments_to_return

      call F_getNormMoments_CapDod_TTTT_LocSep_SepVol &
          (a_poly%c_object, a_localized_separator%c_object, a_moments_to_return%c_object)

    end subroutine getNormMoments_CapDod_TTTT_LocSep_SepVol    

  
end module f_getMoments
