!###############################################################################
!#                                                                             #
!# aed_debug.h                                                                 #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2013 - 2020 -  The University of Western Australia               #
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
#ifndef _AED_DEBUG_H_
#define _AED_DEBUG_H_

!#define DO_DEBUG 1
#if DO_DEBUG
# define ASSERT(expr)  if (.not.(expr)) print *,"assert failed in ",__FILE__," at line ",__LINE__
# define _CheckAllocStatus(status) if (status /= 0) print *,"(de)allocation error in ",__FILE__," at line ",__LINE__
# define CheckAllocStatus(status,a,b) if (status /= 0) print *,"(de)allocation error in ",__FILE__," at line ",__LINE__
# define _CheckFileIOStatus(status) if (status /= 0) print *,"file io error in ",__FILE__," at line ",__LINE__
# define CheckFileIOStatus(a,b,status,c,d) if (status /= 0) print *,"file io error in ",__FILE__," at line ",__LINE__
#else
# define ASSERT(expr)
# define _CheckAllocStatus(status)
# define CheckAllocStatus(status,a,b)
# define _CheckFileIOStatus(status)
# define CheckFileIOStatus(a,b,status,c,d)
#endif

#endif
