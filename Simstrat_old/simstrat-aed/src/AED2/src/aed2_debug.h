!###############################################################################
!#                                                                             #
!# aed_debug.h                                                                 #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
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
