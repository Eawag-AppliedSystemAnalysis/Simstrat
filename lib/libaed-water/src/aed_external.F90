!###############################################################################
!#                                                                             #
!# aed_external.F90                                                            #
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
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created May 2020                                                            #
!#                                                                             #
!###############################################################################

#include "aed.h"

! when building in the library we want them all - but only use the declarations
! when building in the final program we only want those for libs we dont have
#ifdef LIBDEF
#define NO_BENTHIC
#define NO_RIPARIAN
#define NO_DEMO
#define NO_DEV
#endif


!+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#++#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#ifdef NO_BENTHIC
MODULE aed_benthic
!-------------------------------------------------------------------------------
   IMPLICIT NONE
   !#---------------------------------------------------------------------------
   PRIVATE   !# By default make everything private
   PUBLIC aed_new_ben_model, aed_print_version_ben
   !#---------------------------------------------------------------------------

CONTAINS
!===============================================================================

!###############################################################################
FUNCTION aed_new_ben_model(modelname) RESULT(model)
USE aed_core
!-------------------------------------------------------------------------------
   CHARACTER(*),INTENT(in) :: modelname
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
!BEGIN
   NULLIFY(model)
END FUNCTION aed_new_ben_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
SUBROUTINE aed_print_version_ben
!-------------------------------------------------------------------------------
!BEGIN
   print*,"    libaed-ben not included "
END SUBROUTINE aed_print_version_ben
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!===============================================================================
END MODULE aed_benthic
#endif

!+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#++#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#ifdef NO_RIPARIAN
MODULE aed_riparian
!-------------------------------------------------------------------------------
   IMPLICIT NONE
   !#---------------------------------------------------------------------------
   PRIVATE   !# By default make everything private
   PUBLIC aed_new_rip_model, aed_print_version_rip
   !#---------------------------------------------------------------------------

CONTAINS
!===============================================================================

!###############################################################################
FUNCTION aed_new_rip_model(modelname) RESULT(model)
USE aed_core
!-------------------------------------------------------------------------------
   CHARACTER(*),INTENT(in) :: modelname
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
!BEGIN
   NULLIFY(model)
END FUNCTION aed_new_rip_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
SUBROUTINE aed_print_version_rip
!-------------------------------------------------------------------------------
!BEGIN
   print*,"    libaed-riparian no included "
END SUBROUTINE aed_print_version_rip
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!===============================================================================
END MODULE aed_riparian
#endif

!+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#++#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#ifdef NO_DEMO
MODULE aed_demo
!-------------------------------------------------------------------------------
   IMPLICIT NONE
   !#---------------------------------------------------------------------------
   PRIVATE   !# By default make everything private
   PUBLIC aed_new_dmo_model, aed_print_version_dmo
   !#---------------------------------------------------------------------------

CONTAINS
!===============================================================================

!###############################################################################
FUNCTION aed_new_dmo_model(modelname) RESULT(model)
USE aed_core
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: modelname
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!
!-------------------------------------------------------------------------------
!BEGIN
   NULLIFY(model)
END FUNCTION aed_new_dmo_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
SUBROUTINE aed_print_version_dmo
!-------------------------------------------------------------------------------
!BEGIN
   print*,"    libaed-demo no included"
END SUBROUTINE aed_print_version_dmo
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!===============================================================================
END MODULE aed_demo
#endif

!+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#++#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#ifdef NO_DEV
!-------------------------------------------------------------------------------
MODULE aed_dev
   IMPLICIT NONE
   !#---------------------------------------------------------------------------
   PRIVATE   !# By default make everything private
   PUBLIC aed_new_dev_model, aed_print_version_dev
   !#---------------------------------------------------------------------------
CONTAINS
!===============================================================================

!###############################################################################
FUNCTION aed_new_dev_model(modelname) RESULT(model)
USE aed_core
!-------------------------------------------------------------------------------
   CHARACTER(*),INTENT(in) :: modelname
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
!BEGIN
   NULLIFY(model)
END FUNCTION aed_new_dev_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
SUBROUTINE aed_print_version_dev
!-------------------------------------------------------------------------------
!BEGIN
   print*,"    libaed-dev not included"
END SUBROUTINE aed_print_version_dev
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!===============================================================================
END MODULE aed_dev
#endif

!+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#++#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
