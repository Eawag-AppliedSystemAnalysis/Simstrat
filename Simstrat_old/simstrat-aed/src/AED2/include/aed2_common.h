!******************************************************************************
!*                                                                            *
!* aed2_common.h                                                              *
!*                                                                            *
!* Developed by :                                                             *
!*     AquaticEcoDynamics (AED) Group                                         *
!*     School of Earth & Environment                                          *
!* (C) The University of Western Australia                                    *
!*                                                                            *
!******************************************************************************
#ifndef _AED2_COMMON_H_
#define _AED2_COMMON_H_


MODULE aed2_common

 INTERFACE

   FUNCTION aed2_new_model(modelname) RESULT(model)
      CHARACTER(*),INTENT(in) :: modelname
   END FUNCTION aed2_new_model

   SUBROUTINE aed2_build_model(model, namlst, do_prefix)
      CLASS (aed2_model_data_t),POINTER :: model
      INTEGER,INTENT(in)      :: namlst
      LOGICAL,INTENT(in)      :: do_prefix
   END SUBROUTINE aed2_build_model

   SUBROUTINE aed2_define_model(modelname, namlst)
      CHARACTER(*),INTENT(in) :: modelname
      INTEGER,INTENT(in)      :: namlst
   END SUBROUTINE aed2_define_model

 END INTERFACE

END MODULE aed2_common

#endif
