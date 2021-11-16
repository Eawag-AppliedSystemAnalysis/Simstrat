!###############################################################################
!#                                                                             #
!# aed_util.F90                                                                #
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
!# Created June  2012                                                          #
!#                                                                             #
!###############################################################################

#include "aed.h"

!
MODULE aed_util
!-------------------------------------------------------------------------------
!
! aed_util --- shared utility functions for aed modules
!
!-------------------------------------------------------------------------------
   USE aed_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC find_free_lun, aed_real2str, aed_int2str, qsort, exp_integral
   PUBLIC water_viscosity
   PUBLIC aed_gas_piston_velocity, aed_oxygen_sat, aed_n2o_sat
   PUBLIC aed_bio_temp_function,fTemp_function, fSal_function
   PUBLIC PO4AdsorptionFraction, in_zone_set
   PUBLIC InitialTemp, SoilTemp
   PUBLIC make_dir_path, param_file_type, CSV_TYPE, NML_TYPE
!

INTEGER, PARAMETER :: CSV_TYPE = 1
INTEGER, PARAMETER :: NML_TYPE = 2

!===============================================================================
CONTAINS



!###############################################################################
INTEGER FUNCTION find_free_lun()
!-------------------------------------------------------------------------------
! find a free logical unit number
!-------------------------------------------------------------------------------
!LOCALS
    INTEGER :: lun
    LOGICAL :: opend
!
!-------------------------------------------------------------------------------
!BEGIN
   DO lun = 10,99
      inquire( unit = lun, opened = opend )
      IF ( .not. opend ) THEN
         find_free_lun = lun
         RETURN
      ENDIF
   ENDDO

   find_free_lun = -1
END FUNCTION find_free_lun
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_real2str(num, str)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL, INTENT(in) :: num
   CHARACTER(*), INTENT(out) :: str
!
!BEGIN
!-------------------------------------------------------------------------------
   WRITE(str, "(e12.4)") num
   str = TRIM(ADJUSTL(str))
END SUBROUTINE aed_real2str
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_int2str(num, str)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER, INTENT(in) :: num
   CHARACTER(*), INTENT(out) :: str
!
!BEGIN
!-------------------------------------------------------------------------------
   WRITE(str, *) num
   str = TRIM(ADJUSTL(str))
END SUBROUTINE aed_int2str
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION exp_integral(inp) RESULT(E_ib)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: inp
!
!LOCALS
   AED_REAL  :: E_ib !-- Outgoing
   INTEGER   :: j
   AED_REAL  :: ff
!
!-------------------------------------------------------------------------------
!BEGIN
   ff = -1e-9
   IF(ABS(inp-10.0) < 12.0) THEN
     IF(inp==0.0) THEN
       E_ib = inp
     ELSE
       j  = 10+2*IABS(INT(inp))
       ff = 1.0/(REAL(j+1)**2.0)
       DO WHILE(j/=0)
         ff = (ff*REAL(j)*inp+1.0)/REAL(j*j)
         j  = j-1
       ENDDO
       ff   = ff*inp+LOG(1.781072418*ABS(inp))
       E_ib = ff
     ENDIF
   ELSE
     j = 5 + 20 / IABS(INT(inp))
     ff = inp
     DO WHILE(j/=0)
       ff = (1.0/(1.0/ff-1.0/REAL(j)))+inp
       j = j-1
     ENDDO
     ff  = EXP(inp)/ff
     E_ib = ff
   ENDIF

END FUNCTION exp_integral
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
!#                                                                             #
!# A fortran implementation of the quicksort algorithm.                        #
!#                                                                             #
!###############################################################################
RECURSIVE SUBROUTINE qsort(RA,IA,start,end)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: RA(:)
   INTEGER,INTENT(inout) :: IA(:)
   INTEGER,INTENT(in) :: start,end
!
!LOCALS
  INTEGER :: p, l, r
!
!-------------------------------------------------------------------------------
!BEGIN
  IF ( start .LT. end ) THEN
     l=start+1
     r=end
     p = IA(start);

     DO WHILE(l<r)
        IF (cmp(IA(l), p) .LE. 0) THEN
           l=l+1;
        ELSEIF (cmp(IA(r), p) .GE. 0) THEN
           r=r-1
        ELSE
           CALL swap(IA(l), IA(r))
        ENDIF
     ENDDO
     IF (cmp(IA(l), p) .LT. 0 ) THEN
        CALL swap(IA(l), IA(start))
        l=l-1
     ELSE
        l=l-1
        CALL swap(IA(l), IA(start))
     ENDIF

     CALL qsort(RA,IA,start,l)
     CALL qsort(RA,IA,r,end)
  ENDIF

CONTAINS

   !############################################################################
   SUBROUTINE swap(a, b)
   !----------------------------------------------------------------------------
     INTEGER,intent(inout) :: a, b
     INTEGER t
   !----------------------------------------------------------------------------
   !BEGIN
     t = a
     a = b
     b = t
   END SUBROUTINE swap

   !############################################################################
   INTEGER FUNCTION cmp(l,r)
   !----------------------------------------------------------------------------
      INTEGER,INTENT(in)::l,r
   !----------------------------------------------------------------------------
   !BEGIN
      IF ( RA(l) .LT. RA(r) ) THEN
         cmp = -1
      ELSEIF ( RA(l) .EQ. RA(r) ) THEN
         cmp = 0
      ELSE
         cmp = 1
      ENDIF
   END FUNCTION cmp
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END SUBROUTINE qsort
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
LOGICAL FUNCTION in_zone_set(matz, active_zones)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: matz
   AED_REAL,INTENT(in) :: active_zones(:)
!
!LOCALS
   INTEGER :: i, l
   LOGICAL :: res
!BEGIN
!-------------------------------------------------------------------------------
   res = .FALSE.
   l = size(active_zones)
   DO i=1,l
      IF ( active_zones(i) == matz ) THEN
         res = .TRUE.
         EXIT
      ENDIF
   ENDDO

   in_zone_set = res
END FUNCTION in_zone_set
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
FUNCTION water_viscosity(temperature) RESULT(mu)
!-------------------------------------------------------------------------------
! Calculates the molecular viscosity of water for a given temperature
!
! From Table A.1b, FLUID_MECHANICS With Engineering Applications
! by Robert L. Daugherty and Joseph B. Franzini,
! however, note these values are common in most fluid mechanics texts.
! NOTE: N s / m^2  = kg / m / s
!
!  Temp (C)     Viscosity (N s / m^2) x 10^3
!  --------     ---------
!      0          1.781
!      5          1.518
!     10          1.307
!     15          1.139
!     20          1.002
!     25          0.890
!     30          0.798
!     40          0.653
!     50          0.547
!     60          0.466
!     70          0.404
!     80          0.354
!     90          0.315
!    100          0.282
!
!-------------------------------------------------------------------------------
!ARGUMENTS
  AED_REAL,INTENT(inout)  :: temperature
  AED_REAL :: mu
!
!LOCALS
!
!-------------------------------------------------------------------------------
!BEGIN
   !-- Check for non-sensical temperatures
   IF( temperature<zero_ ) temperature = 0.0
   IF( temperature>100.0 ) temperature = 100.0

   IF( temperature<=20.0 ) THEN
     ! 0C to 20C
     ! y = 0.0008 * x^2 - 0.0556 * x + 1.7789
     ! r^2 = 0.9999
     mu = 0.0008 * temperature**2. - 0.0556 * temperature + 1.7789

   ELSEIF(temperature <= 60) THEN
     ! 20C to 60C
     ! y = 0.0002 * x^2 - 0.0323 * x + 1.5471
     ! r^2 = 0.9997
     mu = 0.0002 * temperature**2. - 0.0323 * temperature + 1.5471
   ELSE
     ! 60C to 100C
     ! y = 0.00006 * x^2 - 0.0141 * x + 1.1026
     ! r^2 = 0.9995
     mu = 0.00006 * temperature**2. - 0.0141 * temperature + 1.1026
   ENDIF

   ! Now convert to units of: N s / m^2
   mu = mu / 1e3

END FUNCTION water_viscosity
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
PURE AED_REAL FUNCTION aed_gas_piston_velocity(wshgt,wind,tem,sal,vel,depth,  &
                                                LA,schmidt_model,piston_model)
!-------------------------------------------------------------------------------
! Atmospheric-surface water exchange piston velocity for O2, CO2, N2O, CH4 etc
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(IN)           :: wshgt,wind
   AED_REAL,INTENT(IN)           :: tem,sal
   AED_REAL,INTENT(IN),OPTIONAL  :: vel,depth
   AED_REAL,INTENT(IN),OPTIONAL  :: LA
   INTEGER, INTENT(IN),OPTIONAL  :: schmidt_model, piston_model
!
!LOCALS
   ! Temporary variables
   AED_REAL :: schmidt,k_wind,k_flow,temp,salt,hgtCorrx,a,x,windsp,vel_l
   INTEGER  :: schmidt_model_l,piston_model_l
   ! Parameters
   AED_REAL,PARAMETER :: roughlength = 0.000114  ! momn roughness length (m)
!
!-------------------------------------------------------------------------------
!BEGIN
   k_wind = 0.  ! default to zero

   !-----------------------------------------------
   ! Decide on models to apply
   schmidt_model_l = 2 ! default
   IF (PRESENT(schmidt_model)) schmidt_model_l = schmidt_model
   piston_model_l  = 1 ! default
   IF (PRESENT(piston_model)) piston_model_l = piston_model
   vel_l = zero_
   IF (PRESENT(vel)) vel_l = vel

   ! These parameterizations assume 10m windspeed, and must be scaled by hgtCorrx
   ! Adjust the windspeed if the sensor height is not 10m
   hgtCorrx =  LOG(10.00 / roughLength) / LOG(wshgt / roughLength)
   windsp = wind * hgtCorrx

   !-----------------------------------------------
   ! Compute k_wind
   IF (PRESENT(LA)) THEN

      ! New option for the calculation of k_wind. Note that this has a
      ! "lake area" (LA) variable included in it.

      ! Valchon & Prairie 2013: The ecosystem size and shape dependence of
      !           gas transfer velocity versus wind speed relationships in lakes
      ! k600 = 2.51 (±0.99) + 1.48 (±0.34) · U10 + 0.39 (±0.08) · U10 · log10 LA

      k_wind = 2.51 + 1.48*windsp  +  0.39*windsp*log10(LA)

   ELSE
      temp=tem
      salt=sal
      IF (temp < 0.0) temp = 0.0;   IF (temp > 38.0) temp = 38.0
      IF (salt < 0.0) salt = 0.0;   IF (salt > 75.0) salt = 75.0

      ! Schmidt, Sc
      ! control value : Sc = 590 at 20°C and 35 psu
      schmidt = 590.

      SELECT CASE (schmidt_model_l)
      CASE (1)
         schmidt = (0.9 + 0.1*salt/35.0)*(1953.4+temp*(-128.0+temp*(3.9918-temp*0.050091)))
      CASE (2)
         schmidt = (0.9 + salt/350.0)
         schmidt = schmidt * (2073.1 -125.62*temp +3.6276*temp*temp - 0.043219*temp*temp*temp)
      CASE (3)
         ! http://www.geo.uu.nl/Research/Geochemistry/kb/Knowledgebook/O2_transfer.pdf
         schmidt = (1.0 + 3.4e-3*salt)
         schmidt = schmidt * (1800.6 -120.1*temp +3.7818*temp*temp - 0.047608*temp*temp*temp)
      CASE (4)
         ! CH4 one from Arianto Santoso <abs11@students.waikato.ac.nz>
         schmidt = 2039.2 - (120.31*temp) + (3.4209*temp*temp) - (0.040437*temp*temp*temp)
         schmidt = schmidt / 600
      CASE (5)
         ! CH4 from Sturm et al. 2014 (ex Wanninkhof, 1992)
         schmidt = 1897.8 - (114.28*temp) + (3.2902*temp*temp) - (0.039061*temp*temp*temp)
      CASE (6)
         ! N2O from Sturm et al. 2014 (ex Wanninkhof, 1992)
         schmidt = 2055.6 - (137.11*temp) + (4.3173*temp*temp) - (0.054350*temp*temp*temp)
      END SELECT


      ! Gas transfer velocity (cm/hr)
      SELECT CASE (piston_model_l)
      CASE (1)
        ! k = a u^2 (Sc/660)^-x : Wanninkhof 1992
        a = 0.31
        x = 0.50
        IF( windsp <3.) x = 0.66
        k_wind = a * (windsp**2) * (schmidt/660.0)**(-x)
      CASE (2)
        ! k = a u^2 (Sc/660)^-x : Wanninkhof 2014
        a = 0.251
        x = 0.50
        k_wind = a * (windsp**2) * (schmidt/660.0)**(-x)
      CASE (3)
        ! k = a u^2 (Sc/600)^-x : Ho et al., 2011
        a = 0.26
        x = 0.50
        k_wind = a * (windsp**2) * (schmidt/600.0)**(-x)
      CASE (4)
        ! k = K (Sc/600)^-x : Ho et al., 2016
        a = 0.266
        x = 0.50
        k_wind = ((0.77*vel_l**x)*(depth**(-x)) + a*(windsp)**2) * (schmidt/600.0)**(-x)
      CASE (5)
        ! k = K (Sc/600)^-x : Raymond and Cole, 2001
        a = 1.91
        x = 0.50
        k_wind = 1.91*exp(0.35*windsp) * (schmidt/600.0)**(-x)
      CASE (6)
        ! k = K (Sc/600)^-x : Borge et al., 2004
        a = 2.58
        x = 0.50
        k_wind = (1.0 + (1.719*vel_l**x)*(depth**x) + a*windsp) * (schmidt/600.0)**(-x)
      CASE (7)
        ! k = K (Sc/600)^-x : Rosentreter et al., 2016 CO2
        x = 0.50
        k_wind = (-0.08 + 0.26*vel_l + 0.83*windsp +0.59*depth ) * (schmidt/600.0)**(-x)
      CASE (8)
        ! k = K (Sc/600)^-x : Rosentreter et al., 2016 CH4
        x = 0.50
        k_wind = (-1.07 + 0.36*vel_l + 0.99*windsp +0.87*depth ) * (schmidt/600.0)**(-x)
      CASE (9)
        ! k = K (Sc/600)^-x : Liss and Merlivat, 1986
        x = 0.50
        IF (windsp <3) THEN
          a = 0.17*windsp
        ELSE IF (windsp <13) THEN
          a = 2.85*windsp - 9.65
        ELSE
          a = 5.9*windsp - 49.3
        ENDIF
        k_wind = a * (schmidt/600.0)**(-x)
      END SELECT

   ENDIF

   ! convert to m/s
   k_wind = k_wind / 3.6e5

   !-----------------------------------------------
   ! Compute k_flow
   k_flow = zero_   ! The above options are including flow in estuary ones.

   !-----------------------------------------------
   ! piston velocity is the sum due to flow and wind
   aed_gas_piston_velocity = k_flow + k_wind

END FUNCTION aed_gas_piston_velocity
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
PURE AED_REAL FUNCTION aed_oxygen_sat(salt,temp)
!-------------------------------------------------------------------------------
!  Calculated saturated oxygen concentration at salinity and temperature
! Taken from Riley and Skirrow (1974)
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: salt,temp
!
!LOCALS
   AED_REAL :: Tabs
   AED_REAL :: buf1, buf2, buf3, sol_coeff
!
!-------------------------------------------------------------------------------
!BEGIN
   buf1 = zero_ ; buf2 = zero_ ; buf3 = zero_ ; sol_coeff = zero_

   Tabs = temp + 273.15
   buf1 = -173.4292 + 249.6339 * 100.0 / Tabs + 143.3483 * LOG(Tabs/100.0)
   buf2 = 21.8492 * Tabs / 100.0
   buf3 = salt * (-0.033096 + 0.014259 * Tabs / 100.0 - 0.0017 * (Tabs / 100.0)**2.0)
   sol_coeff = buf1 - buf2 + buf3

   aed_oxygen_sat = 1.42763 * exp(sol_coeff) !in g/m3

   !Convert to mmol/m3
   aed_oxygen_sat = (aed_oxygen_sat / 32.) * 1e3
END FUNCTION aed_oxygen_sat
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION aed_n2o_sat(salt,temp)
!-------------------------------------------------------------------------------
!  gsw_N2Osol_SP_pt                            solubility of N2O in seawater
!
!  USAGE:
!   N2Osol = gsw_N2Osol_SP_pt(SP,pt)
!
!  DESCRIPTION:
!   Calculates the nitrous oxide, N2O, concentration expected at equilibrium
!   with air at an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar)
!   including saturated water vapor  This function uses the solubility
!   coefficients as listed in Hamme and Emerson (2004).
!
!   Note that this algorithm has not been approved by IOC and is not work
!   from SCOR/IAPSO Working Group 127. It is included in the GSW
!   Oceanographic Toolbox as it seems to be oceanographic best practice.
!
!  INPUT:
!   salt  =  Practical Salinity  (PSS-78)                         [ unitless ]
!   temp  =  potential temperature (ITS-90) referenced               [ deg C ]
!          to one standard atmosphere (0 dbar).
!
!  OUTPUT:
!   N2Osol = solubility of N2O                                      [ mol/L ]
!
!  AUTHOR:  Rich Pawlowicz, Paul Barker and Trevor McDougall
!                                                       [ help@teos-10.org ]
!
!  VERSION NUMBER: 3.05 (27th January 2015)
!
!  REFERENCES:
!   IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
!    seawater - 2010: Calculation and use of thermodynamic properties.
!    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
!    UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
!
!   Weiss, R.F. and B.A. Price, 1980: Nitrous oxide solubility in water and
!    seawater. Mar. Chem., 8, 347-359.
!
!   The software is available from http://www.TEOS-10.org
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: salt,temp
!
!LOCALS
   AED_REAL :: x, y, y_100, pt68, ph2odP
   AED_REAL :: a0,a1,a2,a3,b1,b2,b3,m0,m1,m2,m3
!
!-------------------------------------------------------------------------------
!BEGIN

  x = salt     !  Note that salinity argument is Practical Salinity, this is
               !  beacuse the major ionic components of seawater related to Cl
               !  are what affect the solubility of non-electrolytes in seawater

  pt68 = temp*1.00024 ! pt68 is the potential temperature in degress C on
                      ! the 1968 International Practical Temperature Scale IPTS-68.
  y = pt68 + 273.15
  y_100 = y*1e-2

  !  The coefficents below are from Table 2 of Weiss and Price (1980)
  a0 = -165.8806
  a1 =  222.8743
  a2 =  92.0792
  a3 = -1.48425
  b1 = -0.056235
  b2 =  0.031619
  b3 = -0.0048472

  m0 = 24.4543
  m1 = 67.4509
  m2 = 4.8489
  m3 = 0.000544

  ph2odP = exp(m0 - m1*100.0/y - m2*log(y_100) - m3*x) !  Moist air correction at 1 atm.

  !aed_n2o_sat [mol/L] = (exp(a0 + a1*100.0/y + a2*log(y_100) + a3*y_100 + x*(b1 + y_100*(b2 + b3*y_100))))/(1.-ph2odP);
  aed_n2o_sat = (exp(a0 + a1*100.0/y + a2*log(y_100) + a3*y_100*y_100 + x*(b1 + y_100*(b2 + b3*y_100))))/(1.-ph2odP)

  !Convert to mmol/m3
  aed_n2o_sat = aed_n2o_sat * 1e6

END FUNCTION aed_n2o_sat
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_bio_temp_function(numg, theta, T_std, T_opt, T_max, aTn, bTn, kTn, name)
!-------------------------------------------------------------------------------
! Numerical solver for continuous temperature function based on CAEDYM method
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)       :: numg
   AED_REAL,INTENT(in)      :: theta(:)
   AED_REAL,INTENT(inout)   :: T_std(:), T_opt(:), T_max(:)
   AED_REAL,INTENT(out)     :: aTn(:), bTn(:), kTn(:)
   CHARACTER(64),INTENT(in) :: name(:)
!
!LOCALS
   AED_REAL :: Ts       ! Min. temperature where fT(Ts)=I (usually 1)
   AED_REAL :: To       ! Optimum temperature where d(fT(To))/dT=0
   AED_REAL :: Tm       ! Maximum temperature where fT(Tm)=0
   AED_REAL :: in       ! Constant for fT(Ts)=in
   AED_REAL :: v        ! Constant v
   AED_REAL :: k,a,b    ! Model constants
   AED_REAL :: G        ! Function fT()
   AED_REAL :: devG     ! Derivative of fT()
   AED_REAL :: a0,a1,a2 ! Dummies
   AED_REAL :: tol      ! Tolerance
   INTEGER :: group     ! Group counter
   INTEGER :: i         ! Counters
   AED_REAL,PARAMETER :: t20=20.0
   LOGICAL,PARAMETER :: curvef=.true. ! T : f(T)=v**(T-20) at T=Tsta
                                      ! F : f(T) = 1 at T=Tsta
   AED_REAL,ALLOCATABLE,DIMENSION(:,:) :: value
!
!-------------------------------------------------------------------------------
!BEGIN
    write(*,"(11X,'Solving temperature functions for phytoplankton - ')")
    write(*,"(11X,' using the form : f(T) = v^(T-20)-v^(k(T-a))+b')")

    tol   = 0.05

    DO group=1,numg

      ! Set the constants for the correct group
      v = theta(group)

      IF(v < 1.01) THEN
        print "(/,5X,'WARNING: theta_growth for group ',I2,' < 1.01',/)",group
      ENDIF

      Tm = T_max(group)
      Ts = T_std(group)
      To = T_opt(group)


      IF (Ts<0.0 .AND. To<0.0 .AND. Tm<0.0) THEN
        ! The user inputs the values of kTn, aTn and bTn directly
        kTn(group) = -Ts
        bTn(group) = -Tm
        aTn(group) = -To

        ALLOCATE(value(401,1))
        ! Calculate the temperature function using 0.1 deg C intervals

        DO i = 0,400
          b = REAL(i)/10.0
          value(i+1,1) = v**(b-20) - v**(kTn(group) * (b - aTn(group))) + bTn(group)
        ENDDO

        ! Find the values of Tsta, T_opt and T_max from the temp function
        a=0.0
        DO i=1,SIZE(value,1)
          b=REAL(i-1)/10.0
          IF(value(i,1)>0.0) THEN
            T_max(group) = b
          ENDIF
          IF(value(i,1)>a) THEN
            T_opt(group) = b
            a=value(i,1)
          ENDIF
          IF(value(i,1)>v**(b-20)-tol .and. value(i,1)<v**(b-20)+tol) THEN
            T_std(group) = b
          ENDIF
        ENDDO
        DEALLOCATE(value)


      ELSE
        in = 1.0
        a0 = v**(Ts-t20)
        a1 = v**(To-t20)
        a2 = v**(Tm-t20)

        ! Perform the iteration to find the constants.
        ! First approximation of k.
        k = 6.0
        i = 0
        G = tol + 1.0
        ! Do the iterations until -tol < G < tol
        DO WHILE((G <= -tol) .OR. (G >= tol))
          i=i+1
          IF(i==100) THEN  ! Increases the tolerance if more than 100
            i=0            ! iterations are performed.
            tol=tol+0.01
          ENDIF
          IF(curvef) THEN
            ! Use the condition f(T)=v**(T-20) at T=Tsta
            G = k * v**(k * To) * a2 - a1 * (v**(k * Tm) - v**(k * Ts))
            devG = v**(k * To) * a2 * (in + k * To * log(v)) - a1 * log(v) &
              * (Tm * v**(k * Tm) - Ts * v**(k * Ts))
          ELSE
            ! Use the condition f(T)=1 at T=Tsta
            G = k * v**(k * To) * (a0 - a2 - in) - a1 * (v**(k * Ts) &
              - v**(k * Tm))
            devG = (a0 - a2 - in) * v**(k * To) * (in + k * To * log(v)) - a1 &
              * log(v) * (Ts * v**(k * Ts) - Tm * v**(k * Tm))
          ENDIF
          ! Find the next iteration of k
          k = k - G / devG
        ENDDO

        ! Get the remaining model constants
        IF(k/=0.0) THEN
          a=-log(a1/(k*v**(k*To)))/(k*log(v))
          IF(curvef) THEN
            b=v**(k*(Ts-a))
          ELSE
            b=in+v**(k*(Ts-a))-a0
          ENDIF
        ELSE
          a=0.0
          b=0.0
        ENDIF

        ! Set the model constants to the calculated values
        kTn(group) = k
        aTn(group) = a
        bTn(group) = b
      ENDIF

      IF (kTn(group) < 0.1 .AND. bTn(group) > 100.0) THEN
         print *,'Cannot solve for fT for: ', name(group)
         STOP
      ENDIF

    ENDDO

END SUBROUTINE aed_bio_temp_function
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION fTemp_function(method,T_max,T_std,theta,aTn,bTn,kTn,temp) RESULT(fT)
!-------------------------------------------------------------------------------
! Generic temperature function for phytoplankton and zooplankton taking into
! account a decrease in production above T_opt.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)  :: method
   AED_REAL,INTENT(in) :: T_max, T_std,theta,aTn,bTn,kTn
   AED_REAL,INTENT(in) :: temp  ! Temperature
!
!LOCALS
   AED_REAL  :: fT        !-- Value of the temperature function
   AED_REAL,PARAMETER  :: tp = 20.0
!
!-------------------------------------------------------------------------------
!BEGIN
   fT = one_

   IF ( method /= 1 ) RETURN

   IF (temp > T_max) THEN
       fT = zero_
   ELSEIF ( temp < T_std ) THEN
       IF (ABS(temp-tp) > 1+MINEXPONENT(temp)/2) THEN
         fT = theta**(temp-tp)
       ENDIF
   ELSE
      IF (ABS(temp-tp) > 1 + MINEXPONENT(temp)/2 .AND. &
          ABS((kTn*(temp-aTn)) + bTn) > 1 + MINEXPONENT(temp)/2) THEN
        fT = theta**(temp-tp) - theta**(kTn*(temp - aTn)) + bTn
      ENDIF
   ENDIF
END FUNCTION fTemp_function
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE PO4AdsorptionFraction(PO4AdsorptionModel, &
                                 PO4tot_,            &
                                 ParticleConc_,      &
                                 Kpo4p,K,Qm,         &
                                 PO4dis,PO4par,      &
                                 thepH)
!-------------------------------------------------------------------------------
! Routine to compute fraction of PO4 adsorped to sediment/particulate concentration
!-------------------------------------------------------------------------------
    INTEGER,  INTENT(IN)  :: PO4AdsorptionModel
    AED_REAL, INTENT(IN)  :: PO4tot_, ParticleConc_
    AED_REAL, INTENT(IN)  :: Kpo4p,K,Qm
    AED_REAL, INTENT(OUT) :: PO4dis, PO4par
    AED_REAL, INTENT(IN), OPTIONAL :: thepH

!-------------------------------------------------------------------------------
!LOCALS
    AED_REAL :: buffer, f_pH, pH
    AED_REAL :: PO4tot, ParticleConc
    AED_REAL,PARAMETER :: one_e_neg_ten = 1e-10

!
!-------------------------------------------------------------------------------
!BEGIN
   PO4dis   = zero_
   PO4par   = zero_
   buffer   = zero_
   f_pH     = one_

   ! calculate the total possible PO4 for sorption, and solids
   PO4tot        = MAX(one_e_neg_ten, PO4tot_ )       ! Co in Chao (mg)
   ParticleConc  = MAX(one_e_neg_ten, ParticleConc_ ) ! s in Chao  (mg = mol/L * g/mol * mg/g)


   IF(PO4AdsorptionModel == 1) THEN
     !-----------------------------------------------------
     ! This is the model for PO4 sorption from Ji 2008:
     !
     ! Ji, Z-G. 2008. Hydrodynamics and Water Quality. Wiley Press.

     PO4par = (Kpo4p*ParticleConc) / (one_+Kpo4p*ParticleConc) * PO4tot
     PO4dis = one_ / (one_+Kpo4p*ParticleConc) * PO4tot


   ELSEIF(PO4AdsorptionModel == 2) THEN
     !-----------------------------------------------------
     ! This is the model for PO4 sorption from Chao et al. 2010:
     !
     ! Chao, X. et al. 2010. Three-dimensional numerical simulation of
     !   water quality and sediment associated processes with application
     !   to a Mississippi delta lake. J. Environ. Manage. 91 p1456-1466.
     IF(PRESENT(thepH)) THEN
       pH = MIN(MAX(2.0,thepH),12.0)

       ! -0.0094x2 + 0.0428x + 0.9574
       ! (ursula.salmon@uwa.edu.au: fPH for PO4 sorption to Fe in Mine Lakes)
       f_pH = -0.0094*pH*pH + 0.0428*pH + 0.9574
     ELSE
       f_pH = one_
     END IF

     ! calculate particulate fraction based on quadratic solution

     ! Chao Eq 16
     buffer = SQRT(((PO4tot+(1./K)-(ParticleConc*Qm*f_pH)))**2. + (4.*f_pH*ParticleConc*Qm/K))
     PO4par  = 0.5 * ((PO4tot+(1./K)+(ParticleConc*Qm*f_pH))  - buffer  )

     ! Check for stupid solutions
     IF(PO4par > PO4tot) PO4par = PO4tot
     IF(PO4par < zero_) PO4par = zero_

     ! Now set dissolved portion
     PO4dis = PO4tot - PO4par

   ELSE
     !-----------------------------------------------------
     ! No model is selected

     PO4dis = PO4tot
     PO4par = zero_

   END IF

END SUBROUTINE PO4AdsorptionFraction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
AED_REAL FUNCTION fSal_function(salinity, minS, Smin, Smax, maxS )
!-------------------------------------------------------------------------------
! Salinity tolerance of biotic variables
!
!-------------------------------------------------------------------------------
!ARGUMENTS
  AED_REAL :: salinity
  AED_REAL :: minS, Smin, Smax, maxS
!
!LOCALS
  AED_REAL :: tmp1,tmp2
!
!-------------------------------------------------------------------------------
!BEGIN
   fSal_function = zero_  !## CAB [-Wmaybe-uninitialized]

   !-- Check for non-sensical salinity values

   ! Salinity is non-limiting for sals > Smin and <Smax
   ! minS is salinity when organisms stop production (feeding/photosynthesis)
   ! Note the exception for FW species: Smin and minS are both zero.

   ! Salinity is within the tolerance; no limitation.
   ! If sal in bott cell is > min sal & < max sal set in WQcons, pf=1
   IF( salinity >= Smin .and. salinity <= Smax ) THEN
     fSal_function = one_
     RETURN
   ENDIF

   ! Salinity is greater than the upper bound
   ! maxS is set in caedym_globals at 45psu
   IF( salinity > Smax ) THEN
     fSal_function = (-salinity*salinity+2.0*Smax*salinity-  &
            2.0*Smax*maxS+maxS*maxS)/((Smax-maxS)*(Smax-maxS))
   ENDIF

   ! Salinity is less than the lower bound but greater than low cut
   ! If sal is < min set in WQcons but > clamLowCut set at clamCons.dat)!
   tmp1 = zero_
   tmp2 = zero_
   IF( salinity < Smin .AND. salinity > minS ) THEN
     tmp1 = salinity-minS
     tmp2 = Smin-minS
     fSal_function =  (2*tmp1/Smin-(tmp1*tmp1/(Smin*Smin)))/ &
            (2*tmp2/Smin-(tmp2*tmp2/(Smin*Smin)))
   ENDIF

   ! Salinity is less than the minS
   ! If sal < lowest sal (hardwired at start of fn), shells close
   IF( salinity <= minS ) fSal_function = zero_

   ! If lower bound and low cut are both zero then it is a freshwater species
   ! and we need to set f(S) to one
   IF( Smin==zero_ ) THEN
     IF( salinity <= Smin ) fSal_function =  one_
   ENDIF

   ! Ensure salinity function is not negative
   IF(fSal_function <= zero_) fSal_function = zero_

 END FUNCTION fSal_function
!-------------------------------------------------------------------------------



!===============================================================================
SUBROUTINE InitialTemp(m,depth,wv,topTemp,botTemp,nSPinUpDays,tNew)

   INTEGER,intent(in)   :: m
   AED_REAL,intent(in)  :: wv,depth(0:m+1)
   AED_REAL,intent(in)  :: topTemp,botTemp,nSPinUpDays
   AED_REAL,intent(out) :: tNew(0:m+1)

   INTEGER  :: i
   AED_REAL :: w(m+1),t(0:m+1),tn(0:m+1),k(0:m+1),cp(m),a(m+1),b(m),c(m),d(m),z(0:m+1)
   AED_REAL :: ti, dt, da, f, g, mc, c1, c2, c3, c4


   AED_REAL :: bd = 1.3

   k(0) = 20    ! boundary layer conductance in w/(m^2 k)

   ti = 0      ! ti is time of day
   dt = 3600   ! dt is time step (sec)
   da = 0      ! da is day number

   f = .6      !
   g = 1-f     !
   mc = .12;   ! clay fraction

   !z(0) = 0.00
   !z(1) = 0.001
   !z(i+1) = z(i) + 0.005 * 1.5**(i-1)

   z = depth

   t(:) = botTemp
   t(0:1) = topTemp
   tn(m+1) = t(m+1)

    c1 = .65-.78*bd+.6*bd*bd
    c2 = 1.06*bd
    c3 = 1+2.6/sqrt(mc)
    c4 = .3+.1*bd*bd

    DO i=1, m
      cp(i) = (2400000*bd/2.65+4180000*wv)*(z(i+1)-z(i-1))/(2*dt)
      k(i) = (c1+c2*wv-(c1-c4) * exp(-(c3*wv)**4))/(z(i+1)-z(i))
    ENDDO

    DO WHILE(da < nSpinUpDays)
        ti = ti+dt/3600
        IF ( ti > 24 ) THEN
           ti = ti-24
           da = da+1
        ENDIF
        tn(0) = 5.0 ! topTemp !ta+am*sin(.261799*(ti-6))
        DO i=1, m
            c(i) = -k(i)*f
            a(i+1) = c(i)
            b(i) = f*(k(i)+k(i-1))+cp(i)
            d(i) = g*k(i-1)*t(i-1)+(cp(i)-g*(k(i)+k(i-1)))*t(i)+g*k(i)*t(i+1)
        ENDDO
        d(1) = d(1)+k(0)*tn(0)*f
        d(m) = d(m)+k(m)*f*tn(m+1)
        DO i=1, m-1
            c(i) = c(i)/b(i)
            d(i) = d(i)/b(i)
            b(i+1) = b(i+1)-a(i+1)*c(i)
            d(i+1) = d(i+1)-a(i+1)*d(i)
        ENDDO
        tn(m) = d(m)/b(m)
        DO i=m-1, 1,-1
            tn(i) = d(i)-c(i)*tn(i+1)
        ENDDO
    !    print *,'depth ,temperature, k(i)'

        DO i=0, m+1
    !         print *,'s', z(i),tn(i),k(i)
             tNew(i) = tn(i)
        ENDDO
    ENDDO
    !print *,"Done InitialTemp"
END SUBROUTINE InitialTemp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!===============================================================================
SUBROUTINE SoilTemp(m,depth,wv,topTemp,temp,heatflux)
   INTEGER,intent(in) :: m
   AED_REAL,intent(in) :: wv(:), depth(:), topTemp
   AED_REAL,intent(inout) :: temp(m+1)
   AED_REAL,intent(out) :: heatflux
!
    AED_REAL :: w(m+1),t(0:m+1),tn(0:m+1),k(0:m+1),cp(m),a(m+1),b(m),c(m),d(m),z(0:m+1)
    INTEGER  :: i
    AED_REAL :: ti, dt, da, f, g, mc, c1, c2, c3, c4

    AED_REAL :: bd = 1.3


    k(0) = 20    ! boundary layer conductance in w/(m^2 k)

    ti = 0      ! ti is time of day
    dt = 3600   ! dt is time step (sec)
    da = 0      ! da is day number

    f = .6      !
    g = 1-f     !
    mc = .12;   ! clay fraction

    z = depth

    DO i=1, m
        t(i) = temp(i)
    ENDDO
    t(0:1) = topTemp
    t(m+1) = temp(m+1)
    tn(m+1) = temp(m+1)

    c1 = .65-.78*bd+.6*bd*bd
    c2 = 1.06*bd
    c3 = 1+2.6/sqrt(mc)
    c4 = .3+.1*bd*bd

    DO i=1, m
        cp(i) = (2400000*bd/2.65+4180000*wv(i))*(z(i+1)-z(i-1))/(2*dt)
        k(i) = (c1+c2*wv(i)-(c1-c4) * exp(-(c3*wv(i))**4))/(z(i+1)-z(i))
    ENDDO

    tn(0) = topTemp
    DO i=1, m
        c(i) = -k(i)*f
        a(i+1) = c(i)
        b(i) = f*(k(i)+k(i-1))+cp(i)
        d(i) = g*k(i-1)*t(i-1)+(cp(i)-g*(k(i)+k(i-1)))*t(i)+g*k(i)*t(i+1)
    ENDDO
    d(1) = d(1)+k(0)*tn(0)*f
    d(m) = d(m)+k(m)*f*tn(m+1)
    DO i=1, m-1
        c(i) = c(i)/b(i)
        d(i) = d(i)/b(i)
        b(i+1) = b(i+1)-a(i+1)*c(i)
        d(i+1) = d(i+1)-a(i+1)*d(i)
    ENDDO
    tn(m) = d(m)/b(m)
    DO i=m-1, 1,-1
        tn(i) = d(i)-c(i)*tn(i+1)
    ENDDO

    heatflux = k(0)*( g*(t(0)-t(1)) + f*(tn(0)-tn(1))) !"W/m2"
    !print *, 'heat flux =',heatflux
    !print *,'depth ,temperature, k(i)'

    DO i=1, m+1
    !   print *,'s', z(i),tn(i),k(i)
        temp(i) = tn(i)
    ENDDO

    !print *,"Done SoilTemp"
END SUBROUTINE SoilTemp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION make_dir_path(dir)
!-------------------------------------------------------------------------------
! Create the directory path as specified
!-------------------------------------------------------------------------------
#ifdef __INTEL_COMPILER
   USE ifport
#endif
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: dir
!LOCALS
   INTEGER :: len, i, sys
   CHARACTER(len=128) :: d
   LOGICAL :: res = .TRUE.
#  define DIRSEP "/"
!BEGIN
!-------------------------------------------------------------------------------
   len = LEN_TRIM(dir)
!print*,'making dir path at "',TRIM(dir),'"'
   d(1:128) = ' '
   DO i=1,len
      IF ( dir(i:i) == '/' ) THEN
         IF ( i > 1 ) THEN
          ! CALL execute_command_line("mkdir " // TRIM(d), exitstat=sys)
! print*,'making dir at "',TRIM(d),'"'
#ifdef __INTEL_COMPILER
             sys = system("mkdir " // TRIM(d))
#else
             CALL system("mkdir " // TRIM(d))
#endif
         ENDIF
         d(i:i) = DIRSEP
      ELSE
         d(i:i) = dir(i:i)
      ENDIF
   ENDDO
! MAKEDIRQQ is an intel fortran extension
! MAKEDIRQQ can create only one directory at a time. You cannot create a new
! directory and a subdirectory below it in a single command. MAKEDIRQQ does not
! translate path delimiters. You can use either slash (/) or backslash (\) as
! valid delimiters.
!  CALL MAKEDIRQQ(d)
!  if not intel ...
!  CALL SYSTEM("mkdir "//d)
!  but the f2008 standard introduces execute_command_line as a standard way
!  however it seems the ifort version we have been using doesnt support it?
   make_dir_path = res
END FUNCTION make_dir_path
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION param_file_type(fname)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: fname
!LOCALS
   INTEGER :: len, i, ic, j
   CHARACTER(len=4) :: ext
!BEGIN
!-------------------------------------------------------------------------------
   param_file_type = -1
   len = LEN_TRIM(fname)
   IF (fname(len-3:len-3) == '.' ) THEN
      ext = '   '
      DO i=1, 3
         ic = ichar(fname(len:len))
         j = 4 - i
         IF (ic >= 65 .AND. ic <= 90) THEN
            ext(j:j) = char(ic+32)
         ELSE
            ext(j:j) = char(ic)
         ENDIF
         len = len - 1
      ENDDO

      IF (ext == 'csv') THEN
         param_file_type = CSV_TYPE
      ELSEIF (ext == 'nml') THEN
         param_file_type = NML_TYPE
      ENDIF
   ENDIF
END FUNCTION param_file_type
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed_util
