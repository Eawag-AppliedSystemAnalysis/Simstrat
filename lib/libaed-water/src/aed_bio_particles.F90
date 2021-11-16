!###############################################################################
!#                                                                             #
!# aed_bio_particles.F90                                                       #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2018 - 2021 -  The University of Western Australia               #
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
!# Created Jul 2018                                                            #
!#                                                                             #
!###############################################################################

#include "aed.h"

!
MODULE aed_bio_particles
!-------------------------------------------------------------------------------
! aed_bio_particles --- test particle model
!
! The AED2 module test contains basic equations that have no dependencies
!-------------------------------------------------------------------------------
   USE aed_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_bio_particles_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_bio_particles_data_t
      !# Variable identifiers
      INTEGER :: id_ptm_01, id_ptm_02, id_ptm_03, id_ptm_04,   &
                 id_ptm_05, id_ptm_06, id_ptm_07, id_ptm_08,   &
                 id_ptm_09, id_ptm_10, id_ptm_11, id_ptm_12,   &
                 id_ptm_13, id_ptm_14, id_ptm_15, id_ptm_16,   &
                 id_ptm_17, id_ptm_18
      INTEGER :: id_ptm101, id_ptm102, id_ptm103, id_ptm104,   &
                 id_ptm105, id_ptm106, id_ptm107, id_ptm108,   &
                 id_ptm109, id_ptm110, id_ptm111, id_ptm112,   &
                 id_ptm113, id_ptm114, id_ptm115, id_ptm116,   &
                 id_ptm117, id_ptm118
      INTEGER :: id_ptm_00
      INTEGER :: id_d_oxy, id_d_dc, id_d_dn, id_d_dp
      INTEGER :: id_oxy,id_amm,id_nit,id_frp,id_doc,id_don,id_dop
      INTEGER :: id_lht, id_larea

      AED_REAL :: vvel_new, vvel_old, decay_rate_new, decay_rate_old
      AED_REAL :: X_dwww, X_cdw, X_nc, X_pc, mass_limit

      CONTAINS
         PROCEDURE :: define             => aed_define_bio_particles
!        PROCEDURE :: calculate          => aed_calculate_bio_particles
!        PROCEDURE :: calculate_benthic  => aed_calculate_benthic_bio_particles
!        PROCEDURE :: calculate_riparian => aed_calculate_riparian_bio_particles
!        PROCEDURE :: calculate_dry      => aed_calculate_dry_bio_particles
!        PROCEDURE :: equilibrate        => aed_equilibrate_bio_particles
!        PROCEDURE :: mobility           => aed_mobility_bio_particles
!        PROCEDURE :: light_extinction   => aed_light_extinction_bio_particles
!        PROCEDURE :: delete             => aed_delete_bio_particles
         PROCEDURE :: particle_bgc       => aed_particle_bgc_bio_particles
   END TYPE

   INTEGER, PARAMETER :: PTM_MASS   = 15
   INTEGER, PARAMETER :: PTM_VVEL   = 14
   INTEGER, PARAMETER :: PTM_BIRTH  = 17
   INTEGER, PARAMETER :: PTM_AGE    = 18
   INTEGER, PARAMETER :: PTM_STATUS = 19

   LOGICAL :: extra_diag = .false.
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_define_bio_particles(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and the variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_bio_particles_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER  :: status

!  %% NAMELIST    %%  /aed_bio_particles/
!  %% Last Checked 20/08/2021
   AED_REAL :: vvel_new = 0.
   AED_REAL :: vvel_old = 0.
   AED_REAL :: decay_rate_new = 0.
   AED_REAL :: decay_rate_old = 0.
   AED_REAL :: mass_limit = 10.
   AED_REAL :: X_cdw = 0.5
   AED_REAL :: X_nc = 0.1
   AED_REAL :: X_pc = 0.01
   AED_REAL :: X_dwww = 1.0

!  From Module Globals
!  LOGICAL  :: extra_diag = .false.      !## Obsolete Use diag_level = 10
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST    %%  /aed_bio_particles/

   NAMELIST /aed_bio_particles/ vvel_new, vvel_old, &
                           decay_rate_new, decay_rate_old, &
                           mass_limit, extra_diag, diag_level, &
                           X_dwww, X_cdw, X_nc, X_pc
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Initialise

   ! Read the namelist
   read(namlst,nml=aed_bio_particles,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_bio_particles'

   print *,"        aed_bio_particles initialization"

   IF ( extra_diag ) diag_level = 10

   ! Set module parameters
   data%vvel_new       = vvel_new/secs_per_day
   data%vvel_old       = vvel_old/secs_per_day
   data%decay_rate_new = decay_rate_new/secs_per_day
   data%decay_rate_old = decay_rate_old/secs_per_day
   data%mass_limit     = mass_limit
   data%X_dwww         = X_dwww
   data%X_cdw          = X_cdw
   data%X_nc           = X_nc
   data%X_pc           = X_pc

   ! Diagnostic outputs for particle properties
   data%id_ptm_00 = aed_define_diag_variable('total_count', '#', 'particle count')
   data%id_ptm_14 = aed_define_diag_variable('total_vvel', 'm/s', 'sum of particle vvel')
   data%id_ptm_15 = aed_define_diag_variable('total_mass', 'g', 'sum of particle mass')
   data%id_ptm_17 = aed_define_diag_variable('total_birth', 'day', 'sum of birth date')
   data%id_ptm_18 = aed_define_diag_variable('total_age', 'days', 'sum of particle age')
   data%id_ptm114 = aed_define_diag_variable('vvel', 'm/s', 'last particle vvel')
   data%id_ptm115 = aed_define_diag_variable('mass', 'g', 'last particle mass')
   data%id_ptm117 = aed_define_diag_variable('birth', 'day', 'last particle birth time')
   data%id_ptm118 = aed_define_diag_variable('age', 'days', 'last particle age')

   ! Junk properties
   IF( diag_level >= 10 ) THEN
    data%id_ptm_01 = aed_define_diag_variable('bioptm01', '', 'bio_particles 01')
    data%id_ptm_02 = aed_define_diag_variable('bioptm02', '', 'bio_particles 02')
    data%id_ptm_03 = aed_define_diag_variable('bioptm03', '', 'bio_particles 03')
    data%id_ptm_04 = aed_define_diag_variable('bioptm04', '', 'bio_particles 04')
    data%id_ptm_05 = aed_define_diag_variable('bioptm05', '', 'bio_particles 05')
    data%id_ptm_06 = aed_define_diag_variable('bioptm06', '', 'bio_particles 06')
    data%id_ptm_07 = aed_define_diag_variable('bioptm07', '', 'bio_particles 07')
    data%id_ptm_08 = aed_define_diag_variable('bioptm08', '', 'bio_particles 08')
    data%id_ptm_09 = aed_define_diag_variable('bioptm09', '', 'bio_particles 09')
    data%id_ptm_10 = aed_define_diag_variable('bioptm10', '', 'bio_particles 10')
    data%id_ptm_11 = aed_define_diag_variable('bioptm11', '', 'bio_particles 11')
    data%id_ptm_12 = aed_define_diag_variable('bioptm12', '', 'bio_particles 12')
    data%id_ptm_13 = aed_define_diag_variable('bioptm13', '', 'bio_particles 13')
    data%id_ptm_16 = aed_define_diag_variable('bioptm16', '', 'bio_particles 16')

    data%id_ptm101 = aed_define_diag_variable('bioptm101', '', 'bio_particles101')
    data%id_ptm102 = aed_define_diag_variable('bioptm102', '', 'bio_particles102')
    data%id_ptm103 = aed_define_diag_variable('bioptm103', '', 'bio_particles103')
    data%id_ptm104 = aed_define_diag_variable('bioptm104', '', 'bio_particles104')
    data%id_ptm105 = aed_define_diag_variable('bioptm105', '', 'bio_particles105')
    data%id_ptm106 = aed_define_diag_variable('bioptm106', '', 'bio_particles106')
    data%id_ptm107 = aed_define_diag_variable('bioptm107', '', 'bio_particles107')
    data%id_ptm108 = aed_define_diag_variable('bioptm108', '', 'bio_particles108')
    data%id_ptm109 = aed_define_diag_variable('bioptm109', '', 'bio_particles109')
    data%id_ptm110 = aed_define_diag_variable('bioptm110', '', 'bio_particles110')
    data%id_ptm111 = aed_define_diag_variable('bioptm111', '', 'bio_particles111')
    data%id_ptm112 = aed_define_diag_variable('bioptm112', '', 'bio_particles112')
    data%id_ptm113 = aed_define_diag_variable('bioptm113', '', 'bio_particles113')
    data%id_ptm116 = aed_define_diag_variable('bioptm116', '', 'bio_particles116')
   ENDIF

   ! Diagnostic outputs for fluxes into cell
   data%id_d_oxy = aed_define_diag_variable('oxy_flux', 'mmol O2/m3/day','oxygen consumption')
   data%id_d_dc  = aed_define_diag_variable('dc_flux', 'mmol DOC/m3/day','dissolved C release')
   data%id_d_dn  = aed_define_diag_variable('dn_flux', 'mmol N/m3/day','dissolved N release')
   data%id_d_dp  = aed_define_diag_variable('dp_flux', 'mmol P/m3/day','dissolved P release')

   ! Linked state variables
   data%id_oxy = aed_locate_variable('OXY_oxy')
   data%id_amm = aed_locate_variable('NIT_amm')
   data%id_nit = aed_locate_variable('NIT_nit')
   data%id_frp = aed_locate_variable('PHS_frp')
   data%id_doc = aed_locate_variable('OGM_doc')
   data%id_don = aed_locate_variable('OGM_don')
   data%id_dop = aed_locate_variable('OGM_dop')

   ! Environment variables
   data%id_larea = aed_locate_sheet_global('layer_area')
   data%id_lht = aed_locate_global('layer_ht')

END SUBROUTINE aed_define_bio_particles
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_particle_bgc_bio_particles(data,column,layer_idx,ppid,partcl)
!ARGUMENTS
   CLASS (aed_bio_particles_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   INTEGER,INTENT(inout) :: ppid
   AED_REAL,DIMENSION(:),INTENT(inout) :: partcl
!
!LOCALS
   INTEGER :: n
   AED_REAL :: oxy_flux
   AED_REAL :: decay, area, thickness

   AED_REAL, PARAMETER :: buoyancy_age = 86400.
   AED_REAL, PARAMETER :: DT           = 15.*60.
   AED_REAL, PARAMETER :: decay_rate   = 0.15
   AED_REAL, PARAMETER :: fres         = 0.7
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Check if we are in a new cell, to reset cumulative counters
   IF (ppid == 0) THEN
     _DIAG_VAR_(data%id_ptm_14) = zero_
     _DIAG_VAR_(data%id_ptm_15) = zero_
     _DIAG_VAR_(data%id_ptm_17) = zero_
     _DIAG_VAR_(data%id_ptm_18) = zero_
     _DIAG_VAR_(data%id_d_oxy) = zero_
     _DIAG_VAR_(data%id_d_dc)  = zero_
     _DIAG_VAR_(data%id_d_dn)  = zero_
     _DIAG_VAR_(data%id_d_dp)  = zero_

     IF( diag_level >= 10 ) THEN
      _DIAG_VAR_(data%id_ptm_01) = zero_
      _DIAG_VAR_(data%id_ptm_02) = zero_
      _DIAG_VAR_(data%id_ptm_03) = zero_
      _DIAG_VAR_(data%id_ptm_04) = zero_
      _DIAG_VAR_(data%id_ptm_05) = zero_
      _DIAG_VAR_(data%id_ptm_06) = zero_
      _DIAG_VAR_(data%id_ptm_07) = zero_
      _DIAG_VAR_(data%id_ptm_08) = zero_
      _DIAG_VAR_(data%id_ptm_09) = zero_
      _DIAG_VAR_(data%id_ptm_10) = zero_
      _DIAG_VAR_(data%id_ptm_11) = zero_
      _DIAG_VAR_(data%id_ptm_12) = zero_
      _DIAG_VAR_(data%id_ptm_13) = zero_
      _DIAG_VAR_(data%id_ptm_16) = zero_
     ENDIF
   ENDIF


   ! Increment the particle count for this cell and set to diagnostic
   ppid = ppid + 1
   _DIAG_VAR_(data%id_ptm_00) = ppid !,AED_REAL)   ! total number of particles within a cell

    ! Particle decay, changing with age
   thickness = _STATE_VAR_(data%id_lht)
   area      = _STATE_VAR_S_(data%id_larea)

   IF((partcl(PTM_AGE)-partcl(PTM_BIRTH))<buoyancy_age) THEN
     decay = partcl(PTM_MASS) * (DT*data%decay_rate_new)  ! g / timestep
   ELSE
     decay = partcl(PTM_MASS) * (DT*data%decay_rate_old)  ! g / timestep
   ENDIF
   partcl(PTM_MASS) = partcl(PTM_MASS) - decay
   IF ( partcl(PTM_MASS) <=  data%mass_limit ) partcl(PTM_STATUS) = -1

   oxy_flux = data%X_dwww * (1e3/12.) * (decay/DT) * data%X_cdw / (area*thickness)  ! mmol C / m3/ s

   ! Respiration of decaying particles
   _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) - fres * oxy_flux
   _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + fres * oxy_flux * data%X_nc
   _FLUX_VAR_(data%id_frp) = _FLUX_VAR_(data%id_frp) + fres * oxy_flux * data%X_pc
   _FLUX_VAR_(data%id_nit) = _FLUX_VAR_(data%id_nit) + zero_

   ! DOM leakage from particles during decay
   _FLUX_VAR_(data%id_doc) = _FLUX_VAR_(data%id_doc) + (1.-fres) * oxy_flux
   _FLUX_VAR_(data%id_don) = _FLUX_VAR_(data%id_don) + (1.-fres) * oxy_flux * data%X_nc
   _FLUX_VAR_(data%id_dop) = _FLUX_VAR_(data%id_dop) + (1.-fres) * oxy_flux * data%X_pc

   ! Cumulative oxygen and nutrient fluxes into a cell
   _DIAG_VAR_(data%id_d_oxy) = _DIAG_VAR_(data%id_d_oxy) - fres * oxy_flux * secs_per_day    ! O2
   _DIAG_VAR_(data%id_d_dc)  = _DIAG_VAR_(data%id_d_dc) - (1.-fres)* oxy_flux * secs_per_day ! DOC
   _DIAG_VAR_(data%id_d_dn)  = &
                            _DIAG_VAR_(data%id_d_dn) - oxy_flux * data%X_nc * secs_per_day    ! DON + NH4 + NO3
   _DIAG_VAR_(data%id_d_dp)  = _DIAG_VAR_(data%id_d_dp) - oxy_flux * data%X_pc * secs_per_day ! DOP + FRP


   ! Particle bouyancy, changing with age
   IF((partcl(PTM_AGE)-partcl(PTM_BIRTH))<buoyancy_age) THEN
     partcl(PTM_VVEL) = data%vvel_new
   ELSE
     partcl(PTM_VVEL) = data%vvel_old
   ENDIF

   ! Set diagnostics, summarising particles in this cell

   ! 1st, Cumulate particle properties (divide by particle number for average)
   _DIAG_VAR_(data%id_ptm_14) = _DIAG_VAR_(data%id_ptm_14) + partcl(PTM_VVEL)
   _DIAG_VAR_(data%id_ptm_15) = &
                  _DIAG_VAR_(data%id_ptm_15) + partcl(PTM_MASS) ! total particle mass within a cell
   _DIAG_VAR_(data%id_ptm_17) = _DIAG_VAR_(data%id_ptm_17) + partcl(PTM_BIRTH)
   _DIAG_VAR_(data%id_ptm_18) = &
                  _DIAG_VAR_(data%id_ptm_18) + (partcl(PTM_AGE)-partcl(PTM_BIRTH)) /secs_per_day
   IF( diag_level >= 10 ) THEN
    _DIAG_VAR_(data%id_ptm_01) = _DIAG_VAR_(data%id_ptm_01) + partcl(1)
    _DIAG_VAR_(data%id_ptm_02) = _DIAG_VAR_(data%id_ptm_02) + partcl(2)
    _DIAG_VAR_(data%id_ptm_03) = _DIAG_VAR_(data%id_ptm_03) + partcl(3)
    _DIAG_VAR_(data%id_ptm_04) = _DIAG_VAR_(data%id_ptm_04) + partcl(4)
    _DIAG_VAR_(data%id_ptm_05) = _DIAG_VAR_(data%id_ptm_05) + partcl(5)
    _DIAG_VAR_(data%id_ptm_06) = _DIAG_VAR_(data%id_ptm_06) + partcl(6)
    _DIAG_VAR_(data%id_ptm_07) = _DIAG_VAR_(data%id_ptm_07) + partcl(7)
    _DIAG_VAR_(data%id_ptm_08) = _DIAG_VAR_(data%id_ptm_08) + partcl(8)
    _DIAG_VAR_(data%id_ptm_09) = _DIAG_VAR_(data%id_ptm_09) + partcl(9)
    _DIAG_VAR_(data%id_ptm_10) = _DIAG_VAR_(data%id_ptm_10) + partcl(10)
    _DIAG_VAR_(data%id_ptm_11) = _DIAG_VAR_(data%id_ptm_11) + partcl(11)
    _DIAG_VAR_(data%id_ptm_12) = _DIAG_VAR_(data%id_ptm_12) + partcl(12)
    _DIAG_VAR_(data%id_ptm_13) = _DIAG_VAR_(data%id_ptm_13) + partcl(13)
    _DIAG_VAR_(data%id_ptm_16) = _DIAG_VAR_(data%id_ptm_16) + partcl(16)
   ENDIF

   ! 2nd, Set particle property (this will therefore remember last particle only)
   _DIAG_VAR_(data%id_ptm114) = partcl(PTM_VVEL)
   _DIAG_VAR_(data%id_ptm115) = partcl(PTM_MASS)
   _DIAG_VAR_(data%id_ptm117) = partcl(PTM_BIRTH)
   _DIAG_VAR_(data%id_ptm118) = (partcl(PTM_AGE)-partcl(PTM_BIRTH))/secs_per_day
   IF( diag_level >= 10 ) THEN
    _DIAG_VAR_(data%id_ptm101) = partcl(1)
    _DIAG_VAR_(data%id_ptm102) = partcl(2)
    _DIAG_VAR_(data%id_ptm103) = partcl(3)
    _DIAG_VAR_(data%id_ptm104) = partcl(4)
    _DIAG_VAR_(data%id_ptm105) = partcl(5)
    _DIAG_VAR_(data%id_ptm106) = partcl(6)
    _DIAG_VAR_(data%id_ptm107) = partcl(7)
    _DIAG_VAR_(data%id_ptm108) = partcl(8)
    _DIAG_VAR_(data%id_ptm109) = partcl(9)
    _DIAG_VAR_(data%id_ptm110) = partcl(10)
    _DIAG_VAR_(data%id_ptm111) = partcl(11)
    _DIAG_VAR_(data%id_ptm112) = partcl(12)
    _DIAG_VAR_(data%id_ptm113) = partcl(13)
    _DIAG_VAR_(data%id_ptm116) = partcl(16)
   ENDIF

   !print *,'ptm ',ppid,(partcl(PTM_AGE)-partcl(PTM_BIRTH))/secs_per_day,partcl(PTM_MASS), &
   !                                       decay,partcl(PTM_VVEL),_DIAG_VAR_(data%id_ptm_15)

END SUBROUTINE aed_particle_bgc_bio_particles
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



END MODULE aed_bio_particles
