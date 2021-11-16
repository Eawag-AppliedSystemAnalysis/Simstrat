!###############################################################################
!#                                                                             #
!# aed.h                                                                       #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2013 - 2021 -  The University of Western Australia               #
!#                                                                             #
!#   GLM is free software: you can redistribute it and/or modify               #
!#   it under the terms of the GNU General Public License as published by      #
!#   the Free Software Foundation, either version 3 of the License, or         #
!#   (at your option) any later version.                                       #
!#                                                                             #
!#   GLM is distributed in the hope that it will be useful,                    #
!#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
!#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
!#   GNU General Public License for more details.                              #
!#                                                                             #
!#   You should have received a copy of the GNU General Public License         #
!#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
!#                                                                             #
!###############################################################################
#ifndef _AED_H_
#define _AED_H_

#define AED_VERSION  "2.0.4a"

#define MAX_MODELS 40

!# aed_phytoplankton constants
#define MAX_PHYTO_TYPES 256
#define MAX_ZOOP_TYPES  256
#define MAX_ZOOP_PREY    10
#define MAX_PATHO_TYPES 256

!# for aed_geochemistry
#define MAX_GC_COMPONENTS 20
#define MAX_GC_MINERALS   20

!# for aed_vegetation
#define MAX_VEG_TYPES   256

!# for aed_macrophytes
#define MAX_ZONES       256

!# for aed_ass
#define MAX_ASS_PARAMS  20


!#define MISVAL -9999.
#define MISVAL misval_

#define INP_LINE_LEN 512
#define STR_LEN       32

!# for generic vertical settling/mobility approaches
#define _MOB_OFF_ 0
#define _MOB_CONST_ 1
#define _MOB_TEMP_  2
#define _MOB_STOKES_  3
#define _MOB_MOTILE_  4
#define _MOB_ATTACHED_  5

#ifdef SINGLE
#  define AED_REAL real(4)
#else
#  define AED_REAL real(8)
#endif

#define DOUBLETYPE double precision

#define _STATE_VAR_(id)   column(id)%cell(layer_idx)
#define _STATE_VAR_S_(id) column(id)%cell_sheet
#define _DIAG_VAR_(id)    column(id)%cell(layer_idx)
#define _DIAG_VAR_S_(id)  column(id)%cell_sheet

#define _FLUX_VAR_(id)    column(id)%flux_pel(layer_idx)
#define _FLUX_VAR_T_(id)  column(id)%flux_atm
#define _FLUX_VAR_B_(id)  column(id)%flux_ben

#endif
