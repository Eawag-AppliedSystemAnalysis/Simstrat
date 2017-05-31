!###############################################################################
!#                                                                             #
!# aed2.h                                                                      #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!###############################################################################
#ifndef _AED2_H_
#define _AED2_H_

#define AED2_VERSION  "1.1.0alpha8"

#define MAX_MODELS 40

!# aed2_phytoplankton constants
#define MAX_PHYTO_TYPES 256
#define MAX_ZOOP_TYPES  256
#define MAX_ZOOP_PREY    10
#define MAX_PATHO_TYPES 256

!# for aed2_geochemistry
#define MAX_GC_COMPONENTS 20
#define MAX_GC_MINERALS   20

!# for aed2_vegetation
#define MAX_VEG_TYPES   256

!# for aed2_macrophytes
#define MAX_ZONES       256

!#define MISVAL -9999.
#define MISVAL misval_

#define INP_LINE_LEN 512
#define STR_LEN       32
#define REACTION_START_CH  '['
#define REACTION_END_CH    ']'
#define PLUS               '+'
#define MINUS              '-'
#define EQUALS             '='
#define L_PAREN            '('
#define R_PAREN            ')'

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
