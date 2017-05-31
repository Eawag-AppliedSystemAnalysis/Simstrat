!###############################################################################
!#                                                                             #
!# aed2_sedcandi.F90                                                           #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created Jun 2012                                                            #
!# Amended March 2014 (Dan)                                                    #
!###############################################################################

!------------------------------------------------------------------------------!
!               AED - Aquatic Ecodynamics Modelling Library                    !
!------------------------------------------------------------------------------!
!                                                                              !
!                MODULE: Sediment Diagenetic Model                             !
!                                                                              !
!  Module created based on code from the sediment diagenesis model CANDI into  !
!      AED.  Original CANDI description is covered below in the header         !
!                         CANDI was originally                                 !
!         written by BERNARD P. BOUDREAU (Dalhousie) in FORTRAN77.             !
!               Subsequent converison to F90 and some extensions               !
!                 by Roger Luff (C.CANDI). Extra additions for                 !
!                DOM, metals, and improved geochemistry by Matt Hipsey         !
!                                                                              !
!------------------------------------------------------------------------------!
!                                                                              !
!  This program is free software; you can redistribute it and/or               !
!  modify it under the terms of the GNU General Public License                 !
!  version 2, as published by the Free Software Foundation.                    !
!                                                                              !
!  This program is distributed in the hope that it will be useful,             !
!  but WITHOUT ANY WARRANTY; without even the implied warranty of              !
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               !
!  GNU General Public License for more details.                                !
!                                                                              !
!  You should have received a copy of the GNU General Public License           !
!  along with this program; if not, write to the Free Software                 !
!  Foundation, Inc., 51 Franklin Street, 5TH Floor, Boston, MA 02110-1301, USA !
!                                                                              !
!------------------------------------------------------------------------------!
!                                                                              !
!                                                                              !
!     CANDI:  CArbon and Nutrient DIagenesis.                                  !
!     ORIGINAL PROGRAM BY Dr. BERNARD P. BOUDREAU                              !
!                                                                              !
!     Solves a complete CARBON-OXIDANT-NUTRIENT/BY-PRODUCT Model               !
!     using finite differencing in the space variable, x, and using            !
!     the "stiff"-initial-value code VODE as a method-of-lines integrator.     !
!     VODE was programed by P.N. Brown, G.D. Byrne and A.C. Hindmarsh          !
!     (SIAM Journal on Scientific and Statistical Computing, 1989,             !
!     vol. 10, p. 1038-1051.)  This program also calculates pH and             !
!     dissolved protolytic species, i.e. CO2, HCO3-, CO3=, H2S, HS-, etc.      !
!                                                                              !
!                                                                              !
!     The basic model assumptions are:                                         !
!                                                                              !
!                                                                              !
!     - This program can model up to two types of reactive organic             !
!       matter;                                                                !
!     - Six possible oxidants available, i.e. O2, NO3, MnO2, Fe(III), SO4,     !
!       and organic matter itself (i.e. methanogenesis),                       !
!       These used in the traditional sequence;  Fe(III) is written as         !
!       Fe(OH)3 here simply for convenience, as any solid Fe(III)-             !
!       oxy-hydroxide is covered by this model;                                !
!     - Model considers the dissolved by-products: Mn2+,                       !
!       Fe2+, total NH4, total H2S, total PO4, and CH4;                        !
!     - It includes (re-)oxidation of Mn2+, Fe2+, total NH4, total H2S,        !
!       FeS, and CH4 by various reactions; this allows for internal cycling;   !
!     - It allows for the formation of FeS by precititation,                   !
!       precipitation of PO4, reaction of Fe(OH)3 with TH2S, and               !
!       its oxidation by O2, MnO2, and Fe(OH)3;                                !
!     - Any C:N:P ratio in the organic matter is allowed;                      !
!     - bioturbation can be set to any functionality the user desires,         !
!       including dependencies on depth, time, oxygen or organic matter        !
!       concentration.  For illustration purposes only, DB is currently        !
!       depth-dependent only and allowed to decrease in two ways:              !
!           a) like a Gaussian with depth (Christensen, 1982, JGR)             !
!           b) in an approximation to a classic two-layer model: i.e.          !
!                                                                              !
!     |                                                                        !
!     |                                                                        !
!     | ____________________________                                           !
!     |                             \                                          !
!   D |                              \                                         !
!   B |                               \                                        !
!     |                                \                                       !
!     |                                 \                                      !
!     |                                  \                                     !
!     |                                   \                                    !
!     |                                    \                                   !
!    0______________________________________\____                              !
!     0                             x1      x2                                 !
!                 DEPTH  ->                                                    !
!                                                                              !
!           i.e. as x1 -> x2, Guinasso and Schink's (1975) 2-layer model is    !
!           obtained exactly;                                                  !
!                                                                              !
!       The mass balance is significantly better if the Gaussian               !
!       functionality is employed;                                             !
!     - Irrigation can occur and is modelled as a nonlocal                     !
!       sink/source as per Boudreau (1984, JMR);                               !
!     - Transport processes for solutes are molecular diffusion,               !
!       advection (burial) and irrigation;                                     !
!     - Transport processes for solids are bioturbation and burial;            !
!     - Porosity can be depth and time variable, and any user                  !
!       supplied function can be incorporated;                                 !
!     - Ammonia and phospahate are subject to reversible linear                !
!       adsorption-desorption.                                                 !
!                                                                              !
!                                                                              !
!     DISCLAIMER:                                                              !
!                                                                              !
!           The creator of this program has used his best efforts in           !
!           preparing this code.  It is not absolutely guaranteed              !
!           to be error free.  The author/programmer makes no                  !
!           warrantees, expressed or implied, including without                !
!           limitation warrantees of merchantability or fitness                !
!           for any particular purpose.  No liability is accepted              !
!           in any event for any damages, including accidental or              !
!           consequential damages, lost of profits, costs of lost              !
!           candi or programming materials, or otherwise in                     !
!           connection with or arising out of the use of this program.         !
!                                                                              !
!                                                                              !
!           No program can solve all possible problems.  Numerical             !
!           methods have difficulties solving differential equations           !
!           when the reactions are fast and are, therefore, complete           !
!           over a depth interval that is small compared to the full           !
!           depth range being modelled, i.e. 0 < x < L.  This                  !
!           problem is called "stiffness".  The ODE solver VODE.f is           !
!           designed to work with this type of problem, but a user             !
!           can easily overwhelm even its remakable abilities.  If             !
!           carbon is not in mass balance according to the output of           !
!           this code, you may be asking the code to do too much.              !
!           Try using a shallower bottom depth, or dividing the region         !
!           into smaller regions, or solving the problem iteratively           !
!           (i.e. the long time scale reactions, then the fast ones,           !
!           then resolving), or try simplifying your  problem so that          !
!           the stiffness is removed.                                          !
!                                                                              !
!                                                                              !
!                                                                              !
!                                                                              !
!                                                                              !
!     POROSITY & TORTUOSITY:                                                   !
!                                                                              !
!     The porosity is described by a spatially dependent function in a         !
!     subroutine called SED(x).  Porosity is currently time independent.       !
!     While anything could be put there, this program                          !
!     currently assumes that:                                                  !
!                                                                              !
!                   p = (p0 - p00)*exp(-b*x) + p00                             !
!                                                                              !
!     where p0 = porosity at x = 0, p00 = porosity at great depth,             !
!     b = an attenuation constant for porosity.                                !
!                                                                              !
!     Tortuosity, T, is calculated in a subroutine function, and it can be     !
!     a function of porosity, depth and time.  Currently, the assumed form     !
!     is the general "universal" correlation proposed by Boudreau (submitted). !
!                                                                              !
!                          2                                                   !
!                         T  =  1 - 2 ln(p)                                    !
!                                                                              !
!                                                                              !
!                                                                              !
!     ADVECTION:                                                               !
!                                                                              !
!     The burial, w, and porewater, u, velocities are also calculated          !
!     in the subroutine SED(x), which uses the relationship:                   !
!                                                                              !
!                    w(x)  = [(1-p00) * w00]/(1-p(x))                          !
!                                                                              !
!     where "w00" = solids burial velocity at great depth (an input            !
!     parameter.  Therefore, one must know "w00".                              !
!                                                                              !
!                                                                              !
!                                                                              !
!                                                                              !
!                                                                              !
!                                                                              !
!     IRRIGATION:                                                              !
!                                                                              !
!     Irrigation is described by a nonlocal exchange (see subroutine           !
!     EXCH):                                                                   !
!                                                                              !
!                   irrigation =  alpha(x,t)*(Cw - C(x,t))                     !
!                                                                              !
!     where alpha(x,t) is the rate of irrigation at depth x at time t,         !
!     Cw is the concentration of the given solute in the overlying             !
!     water, and C(x,t) is its value in the porewater.                         !
!                                                                              !
!                                                                              !
!                                                                              !
!                                                                              !
!                                                                              !
!     REACTIONS:                                                               !
!                                                                              !
!     While CANDI comes with some (reasonably standard) diagenetic             !
!     reactions pre-programmed, these are in no way meant to be an             !
!     exclusive representation of either chemical diagenesis or the            !
!     reactions that can be solved by this code.  They are merely              !
!     examples, though quite relevant.  The code can be altered readily        !
!     to account for other reactions (and species).                            !
!                                                                              !
!     The first set of reactions included in the code are                      !
!     the standard organic matter decomposition reactions                      !
!     (excluding waters and hydrogens):                                        !
!                                                                              !
!     a) Oxygen Reduction:                                                     !
!                                                                              !
!     (CH2O) (NH3) (H3PO4) + (A+2B) O2 -> A CO2 + B HNO3 + C H3PO4             !
!           A     B       C                                                    !
!                                                                              !
!                                                                              !
!       b) Nitrate Reduction:                                                  !
!                                                                              !
!                        4         A       4A        2                         !
! (CH2O) (NH3) (H3PO4) + -A NO3 -> - CO2 + -- HCO3 + - A N2 + B NH3 + C H3PO4  !
!       A     B       C  5         5        5        5                         !
!                                                                              !
!                                                                              !
!       c) Manganese Oxide Reduction:                                          !
!                                                                              !
! (CH2O) (NH3) (H3PO4) + 2A MnO2 + 3A CO2                                      !
!       A     B       C                                                        !
!                                  -> 4A HCO3 + 2A Mn(II) + B NH3 + C H3PO4    !
!                                                                              !
!       d) Iron Oxy-Hydroxide Reduction:                                       !
!                                                                              !
! (CH2O) (NH3) (H3PO4) + 4A Fe(OH)3 + 7A CO2                                   !
!       A     B       C                                                        !
!                          -> 8A HCO3 + 4A Fe(II) + B NH3 + C H3PO4            !
!                                                                              !
!                                                                              !
!     e) Sulfate Reduction:                                                    !
!                                                                              !
!                         A                 A                                  !
! (CH2O) (NH3) (H3PO4)  + - SO4 -> A HCO3 + - H2S + B NH3 + C H3PO4            !
!       A     B       C   2                 2                                  !
!                                                                              !
!                                                                              !
!     f) Methanogenesis:                                                       !
!                                                                              !
!                           A       A                                          !
! (CH2O) (NH3) (H3PO4)   -> - CH4 + - CO2 + B NH3 + C H3PO4                    !
!       A     B       C     2       2                                          !
!                                                                              !
!                                                                              !
!                                                                              !
!     Secondary Reactions:                                                     !
!                                                                              !
!     The following secondary reactions are currently included in              !
!     the model; however, others can be added by simply editing the            !
!     code, and in particular, subroutine FEX.                                 !
!                                                                              !
!     The by-product Mn2+, Fe2+ and CH4 oxidation reactions are:               !
!                                                                              !
!           2+              -                                                  !
!       4 Fe   + O2 + 8 HCO3  + 2 H2O -> 4 Fe(OH)3 + 8 CO2                     !
!                                                                              !
!           2+              -                                                  !
!       2 Mn   + O2 + 4 HCO3   ->  2 MnO2  + 4 CO2 + 2 H2O                     !
!                                                                              !
!       CH4 + O2 ->   CO2 + 2H20                                               !
!                                                                              !
!                      =            -                                          !
!       CH4 + CO2 + SO4   ->  2 HCO3   +  H2S                                  !
!                                                                              !
!                                                                              !
!     There is also the posibility of coupled Fe-oxidation and                 !
!     manganese reduction:                                                     !
!                                                                              !
!                   2+                -                  2+                    !
!       MnO2  + 2 Fe   + 2 H2O  2 HCO3  -> 2 Fe(OH)3 + Mn   + 2 CO2            !
!                                                                              !
!                                                                              !
!       Dissolve sulfide and ammonia can be oxidized with O2:                  !
!                                                                              !
!                          -        =                                          !
!       H2S + 2 O2 + 2 HCO3  ->  SO4   + 2 CO2  + 2 H2O                        !
!                                                                              !
!          +                  -         -                                      !
!       NH4  +  2 O2  + 2 HCO3   ->  NO3  + 2 CO2 + 3 H2O                      !
!                                                                              !
!     or with manganese oxide:                                                 !
!                                                                              !
!                                   2+     o         -                         !
!       MnO2  +  H2S  + 2 CO2 ->  Mn    + S  + 2 HCO3                          !
!                                                                              !
!                                                                              !
!     or with iron oxy-hydroxide:                                              !
!                                                                              !
!                                         2+    o         -                    !
!       2 Fe(OH)3  + H2S  + 4 CO2 ->  2 Fe   + S  + 4 HCO3  + 2 H2O            !
!                                                                              !
!                                                                              !
!     In addition, solid FeS can be oxidized by any of the following           !
!     reactions :                                                              !
!                                     2+    o         -                        !
!       FeS + 2 Fe(OH)3  6 CO2 -> 3 Fe   + S  + 6 HCO3                         !
!                                                                              !
!                                            2+     2+      2-         -       !
!       FeS + 4 MnO2  + 8 CO2 + 4 H2O -> 4 Mn   + Fe   + SO4   + 8 HCO3        !
!                                                                              !
!                           2+      =                                          !
!       FeS  +  2 O2  ->  Fe   + SO4                                           !
!                                                                              !
!                                                                              !
!     Dissolved iron can be removed by FeS formation:                          !
!                                                                              !
!         2+     -        -                                                    !
!       Fe   + HS   + HCO3  ->  FeS + CO2 + H2O                                !
!                                                                              !
!                                                                              !
!                                                                              !
!                                                                              !
!     DIAGENETIC EQUATIONS solved are:                                         !
!                                                                              !
!                                                                              !
!     1) for each organic matter TYPE (Gi, solid):                             !
!                                                                              !
!                                                                              !
!       dGi     1    d             dGi       dpswGi                            !
!       --   =  -- [ -- { ps DB(x) --- }  -  ------]  -  ki Gi Rox     (1)     !
!       dt      ps   dx            dx          dx                              !
!                                                                              !
!                                                                              !
!     2) for oxygen  (O2, dissolved):                                          !
!                                                                              !
!       dO2   1      d   p  d2O2      dvpO2                                    !
!       --- = - [DO2 --( -  ---- ) -  -----] +  a(x){(O2)o - O2}               !
!       dt    p      dx  T   dx        dx                                      !
!                                                                              !
!            ps  (A+2B)                                                        !
!          - --  ------  {k1 G1 + k2 G2} RO2 - {2 kNHOx TNH4 + 2 kSOx TH2S     !
!            p     A                                                           !
!                                      ps                                      !
!       + kMnOx MnII + kFeOx FeII  + 2 -- kFeSOx FeS + KCH4Ox CH4} 02          !
!                                      p                               (2)     !
!                                                                              !
!                                                                              !
!                                                                              !
!      3) for nitrate (NO3, dissolved):                                        !
!                                                                              !
!       dNO3    1       d   p  d2NO3       dvpNO3                              !
!       ----  = - [DNO3 --( -  -----  ) -  ------] +  a(x){(NO3)o - NO3}       !
!        dt     p       dx  T   dx           dx                                !
!                                                                              !
!            ps                                                                !
!          - -- {k1 G1 + k2 G2} (4/5 RNO3 - B/A RO2) + kNHOx O2 TNH4           !
!            p                                                         (3)     !
!                                                                              !
!                                                                              !
!      4) MnO2 (solid):                                                        !
!                                                                              !
!       dMnO2      1    d             dMnO2       dpswMnO2                     !
!       -----   =  -- [ -- { ps DB(x) ----- }  -  --------]                    !
!        dt        ps   dx              dx           dx                        !
!                                                                              !
!                                    p                                         !
!        - 2 {k1 G1 + k2 G2} RMnO2 + -- 2 kMnOx O2 MnII - 4 kFeSMn4 FeS MnO2   !
!                                    ps                                        !
!                                                                              !
!        - kMnFe MnO2 FeII  -  kMnO2TS MnO2 TH2                        (4)     !
!                                                                              !
!                                                                              !
!                                                                              !
!     5) Fe(OH)3 (solid):                                                      !
!                                                                              !
!       dFeOH3      1    d            dFeOH3       dpswFeOH3                   !
!       ------  =  -- [ -- { ps DB(x) ------ }  -  ---------]                  !
!         dt       ps   dx              dx            dx                       !
!                                                                              !
!                                                                              !
!            -  4 {k1 G1 + k2 G2} RFeOH3                                       !
!                                                                              !
!              p                                                               !
!            + -- {4 kFeOx O2 FeII} - 2 kTSFe3 TH2S FeOH3 + kFeSOx FeS O2      !
!              ps                                                              !
!                                                                              !
!            - 2 kFeSFe3 FeS FeOH3  + 2 kMnFe MnO2 FeII                (5)     !
!                                                                              !
!                                                                              !
!     6) Sulfate (SO4, solute):                                                !
!                                                                              !
!       dSO4   1       d    p  d2SO4      dvpSO4                               !
!       ---- = - [DSO4 -- ( -  ----- ) -  ------]  +  a(x){(SO4)o - SO4}       !
!        dt    p       dx   T   dx          dx                                 !
!                                                                              !
!                ps                                                            !
!              - -- (1/2) {k1 G1 + k2 G2} RSO4  +  kSOx O2 TH2S                !
!                 p                                                            !
!                                                                              !
!                ps                                                            !
!              + -- kFeSMn4 FeS MnO2  - KCH4SO4 CH4 SO4                (6)     !
!                p                                                             !
!                                                                              !
!                                                                              !
!                                                                              !
!                                                                              !
!     7) Total Phosphate (TPO4, solute):                                       !
!                                                                              !
!       dTPO4      d    p DTPO4                 d2TPO4                         !
!       -----  = { --[( ------- + ps DB(x) KPA) ------ )                       !
!        dt        dx      T                      dx                           !
!                                                                              !
!              - (p v + ps w KPA) TPO4] - p a(x){(TPO4)o - TPO4}               !
!                                                                              !
!                   C                                                          !
!        + ps - (k1 G1 + k2 G2) Rox - p kPO4PPT (KSPPO4-TPO4) }/PP     (7)     !
!                   A                                                          !
!                                                                              !
!                                                                              !
!     8) Total Ammonia (TNH4, solute):                                         !
!                                                                              !
!       dTNH4      d    p DTNH4                 d2TNH4                         !
!       -----  = { --[( ------- + ps DB(x) KNA) ------ )                       !
!        dt        dx      T                      dx                           !
!                                                                              !
!              - (p v + ps w KNA) TNH4] - p a(x){(TNH4)o - TNH4}               !
!                                                                              !
!                    B                                                         !
!              + ps  - {k1 G1 + k2 G2} {RNO3 + RMnO2 + RFeOH + RSO4}           !
!                    A                                                         !
!                                                                              !
!                - p kNHOx O2 TNH4 }/PN                                (8)     !
!                                                                              !
!                                                                              !
!     9) Total Dissolved Sulfide (TH2S, solute):                               !
!                                                                              !
!       dTH2S    1        d   p  d2TH2S     dvpTH2S                            !
!       -----  = - [DTH2S --( -  ------ ) - -------] + a(x){(TH2S)o - TH2S}    !
!        dt      p        dx  T   dx          dx                               !
!                                                                              !
!                 ps 1                                                         !
!               + -- - {k1 G1 + k2 G2} RSO4  - kSOx O2 TH2S                    !
!                 p  2                                                         !
!                                                                              !
!                 ps                                                           !
!               - -- {kFeSPPT [FeII TH2S - KSPFeS] + kTSFe3 TH2S FeOH3         !
!                 p                                                            !
!                                                                              !
!               - kMnO2TS MnO2 TH2S} + KCH4SO4 CH4 SO4                 (9)     !
!                                                                              !
!                                                                              !
!     10) Dissolved Manganese (solute)                                         !
!                                                                              !
!       dMnII    1        d   p  d2MnII     dvpMnII                            !
!       -----  = - [DMnII --( -  ------ ) - -------] + a(x){(MnII)o - MnII}    !
!        dt      p        dx  T    dx         dx                               !
!                                                                              !
!                 ps                                                           !
!               + -- 2 {k1 G1 + k2 G2} RMnO2  -  2 kMnOx O2 MnII               !
!                 p                                                            !
!                                                                              !
!                 ps                                                           !
!               + --  {4 kFeSMn4 FeS MnO2 +  kMnFe MnO2 FeII                   !
!                 p                                                            !
!                                                                              !
!               +  kMnO2TS MnO2 TH2S}                                  (10)    !
!                                                                              !
!                                                                              !
!     11) Dissolved Iron (FeII, solute):                                       !
!                                                                              !
!       dFeII    1        d   p  d2FeII     dvpFeII                            !
!       -----  = - [DFeII --( -  ------ ) - -------] + a(x){(FeII)o - FeII}    !
!        dt      p        dx  T    dx         dx                               !
!                                                                              !
!                 ps                                                           !
!               + -- 4 {k1 G1 + k2 G2} RFeOH3  -  4 kFeOx O2 FeII              !
!                 p                                                            !
!                                                                              !
!                 ps                                                           !
!               + -- {- kFeSPPT [FeII TH2S - KSPFeS]                           !
!                 p                                                            !
!                                                                              !
!               + 3 kFeSFe3 FeS FeOH3 + kFeSMn4 FeS MnO2                       !
!                                                                              !
!               - kMnFe MnO2 FeII + 2 kTSFe3 FeOH3 TH2S}               (11)    !
!                                                                              !
!                                                                              !
!     12) FeS, Iron Monosulfide (solid):                                       !
!                                                                              !
!       dFeS      1    d             dFeS     dpswFeS                          !
!       ----   =  -- [ -- { ps DB(x) ---- } - ------- ]                        !
!        dt       ps   dx             dx        dx                             !
!                                                                              !
!               + kFeSPPT [FeII TH2S - KSPFeS]                                 !
!                                                                              !
!               - kFeSFe3 FeS FeOH3 - kFeSMn4 FeS MnO2                         !
!                                                                              !
!               - kFeSOx FeS O2  - kPYRIT FeS                          (12)    !
!                                                                              !
!                                                                              !
!     13) Total Dissolved CO2 (TCO2, solute):                                  !
!                                                                              !
!       dTCO2    1        d   p  d2TCO2     dvpTCO2                            !
!       -----  = - [DTCO2 --( -  ------ ) - -------] + a(x){(TCO2)o - TCO2}    !
!        dt      p        dx  T   dx          dx                               !
!                                                                              !
!                 ps                                                           !
!               + -- {k1 G1 + k2 G2} {RO2+ RNO3 + RMnO2 + RFeOH3               !
!                 p                                                            !
!                                                                              !
!               + RSO4 + 0.5 RCH4} + KCH4Ox CH4 O2  + KCH4SO4 CH4 SO4          !
!                                                                      (13)    !
!                                                                              !
!     14) Methane (CH4, solute):                                               !
!                                                                              !
!       dCH4    1        d   p  d2CH4     dvpCH4                               !
!       ----  = - [DCH4  --( -  ----- ) - ------] + a(x){(CH4)o - CH4}         !
!        dt     p        dx  T   dx         dx                                 !
!                                                                              !
!                ps                                                            !
!              + -- {k1 G1 + k2 G2} {0.5 RCH4}                                 !
!                p                                                             !
!                                                                              !
!              - KCH4Ox CH4 O2  -  KCH4SO4 CH4 SO4                     (14)    !
!                                                                              !
!                                                                              !
!     15) Dissolved Calcium (solute)                                           !
!                                                                              !
!       dCa      1      d   p  d2CaI     dvpCa                                 !
!       -----  = - [DCa --( -  ----- ) - -----] + a(x){(Ca)o - Ca}             !
!        dt      p      dx  T    dx        dx                                  !
!                                                                              !
!                                                                              !
!               +  kCalDis Ca (CO3)                                    (15)    !
!                                                                              !
!                                                                              !
!                                                                              !
!       where:                                                                     !
!                                                                              !
!               O2                                                             !
!     RO2 = --------                                                   (16)    !
!           KO2 + O2                                                           !
!                                                                              !
!                                                                              !
!                NO3       K'O2                                                !
!     RNO3 = ---------- ---------                                      (17)    !
!            KNO3 + NO3 K'O2 + O2                                              !
!                                                                              !
!                                                                              !
!              MnO2       K'O2      K'NO3                                      !
!     RMnO2 = -------  ---------  -----------                          (18)    !
!             (MnO2)o  K'O2 + O2  K'NO3 + NO3                                  !
!                                                                              !
!                FeOH3      K'O2      K'NO3           (MnO2)o                  !
!     RFeOH3 = --------  ---------  ----------- -------------------            !
!              (FeOH3)o  K'O2 + O2  K'NO3 + NO3 K'Mn(MnO2) + (MnO2)o           !
!                                                                              !
!                                                                      (19)    !
!                                                                              !
!                SO4        K'O2        K'NO3         (MnO2)o                  !
!     RSO4 = ----------  ---------  ----------- -------------------            !
!            KSO4 + SO4  K'O2 + O2  K'NO3 + NO3 K'Mn(MnO2) + (MnO2)o           !
!                                                                              !
!                                                      (FeOH3)o                !
!                                                 ---------------------        !
!                                                 K'Fe(FeOH3) + (FeOH3)o       !
!                                                                              !
!                                                                      (20)    !
!                                                                              !
!                 K'O2        K'NO3        (MnO2)o                             !
!     RCH4 =   ---------  ----------- -------------------                      !
!              K'O2 + O2  K'NO3 + NO3 10 (MnO2) + (MnO2)o                      !
!                                                                              !
!                                 (FeOH3)o           K'SO4                     !
!                          --------------------- -------------                 !
!                          10 (FeOH3) + (FeOH3)o  K'SO4 + SO4                  !
!                                                                              !
!                                                                      (21)    !
!                                                                              !
!                                                                              !
!       Rox = RO2 + RNO3 + RMnO2 + RFeOH + RSO4 + RCH4                 (22)    !
!                                                                              !
!                                                                              !
!     and where:                                                               !
!                                                                              !
!         DO2,DNO3,DMnII,DFeII,DSO4,DCH4  = molecular diffusion                !
!            coefficients for O2,NO3,Mn2+,Fe2+,SO4 AND CH4, respectively,      !
!            corrected for tortuosity, i.e. D = Do*p**N, where                 !
!            Do is the free-solution  diffusion coefficient, p is              !
!            porosity and N is an exponent;                                    !
!                                                                              !
!         DTNH4, DTH2S, DTPO4, DTCO2 = weigthed average molecular              !
!            diffusion coefficients to total NH4, H2S, PO4,                    !
!            and TCO2, respectively;                                           !
!                                                                              !
!         w = burial velocity;                                                 !
!                                                                              !
!         DB(x) = bioturbation coefficient;                                    !
!                                                                              !
!         k1,k2  = rate constants for organic matter decay;                    !
!                                                                              !
!         kNHOx, kSOx, kMnOx, kFeOx = oxidation rate constants                 !
!            for TNH4, TH2S, FeII and MnII;                                    !
!                                                                              !
!         kFeSOx, kFeSMn4, kFeSFe3 = oxidation rate constants for FeS;         !
!                                                                              !
!         kPYRIT =  rate constant for pyrite formation;                        !
!                                                                              !
!         kFeSPPT, kTSFe3  =  rate conastants for FeS formation;               !
!                                                                              !
!         a(x) = irrigation coefficient (function);                            !
!                                                                              !
!         KO2, KNO3, KSO4 = Monod constants for O2, NO3, SO4;                  !
!                                                                              !
!         K'O2, K'NO3, etc = inhibition constants due to O2, NO3, etc.;        !
!                                                                              !
!         p = porosity;  ps = 1 - p, i.e. solid volume fraction;               !
!                                                                              !
!         T  = tortuosity, a function of p(x);                                 !
!                                                                              !
!         A, B, C = stoichiometric coefficients for the                        !
!           organic matter being decomposed;                                   !
!                                                                              !
!         KPA, KNA = linear adsorption coefficients of TPO4 and TNH4;          !
!                                                                              !
!         PP = (p + ps KPA);  PN = (p + ps KNA);                               !
!                                                                              !
!         kCalDis = rate constant for Calcite dissolution                      !
!                                                                              !
!                                                                              !
!                                                                              !
!     BOUNDARY CONDITIONS:                                                     !
!                                                                              !
!     Two different types of boundary conditions at x=0 are allowed            !
!     for the organic matter:                                                  !
!                                                                              !
!     a) a known concentration:                                                !
!                                                                              !
!                  G  = G0    at  x = 0                                (24)    !
!                                                                              !
!                                                                              !
!     b) or a prescribed flux for each of the G:                               !
!                                                                              !
!                        dG                                                    !
!             - ps DB(0) --  + ps w G(0) = FGi   at  x = 0             (25)    !
!                        dx                                                    !
!                                                                              !
!                                                                              !
!     The FeS has b.c. (25) with FGi = 0.                                      !
!                                                                              !
!                                                                              !
!     Solutes have two possible boundary conditions at the sediment-water      !
!     interface:                                                               !
!                                                                              !
!     a) a known concentration:                                                !
!                                                                              !
!                       C = Co        at x = 0                         (26)    !
!                                                                              !
!     b)  a boundary layer regulated flux:                                     !
!                                                                              !
!                        dC      1                                             !
!             - p**(n+1) --  =  --- (Cw-Co)    at  x = 0               (27)    !
!                        dx     del                                            !
!                                                                              !
!     where "del" is the diffusive boundary layer thickness.                   !
!     You must choose one or the other for all species at once.                !
!                                                                              !
!                                                                              !
!     At the bottom of the model, i.e. x = L, two types of                     !
!     boundary conditions are available.  The first is that the                !
!     gradients in each species disappears, i.e.                               !
!                                                                              !
!                       dC                                                     !
!                       --  =  0   at  x = L                           (28)    !
!                       dx                                                     !
!                                                                              !
!     This approximates the condition expected at depth in a semi-             !
!     infinite sediment column.  However, a second posibility is               !
!     that the concentrations are exactly or approximately known               !
!     at some depth x = L.  This would be the case when modelling              !
!     candi sets taken over a limited depth-interval.  In this case,            !
!                                                                              !
!                      C = C        at   x = L                         (29)    !
!                           L                                                  !
!                                                                              !
!                                                                              !
!     UNITS:                                                                   !
!                                                                              !
!     While it would be very useful for this code to take any set              !
!     of consistent units, the practical reality is that one needs             !
!     to choose a standard set.  For this code these are:                      !
!                                                                              !
!           - LENGTHS are in CM                                                !
!           - TIME is in YEARS                                                 !
!           - CONCENTRATION of SOLUTES in Micro-MOLES/CC = mM of porewater     !
!           - CONCENTRATION of SOLID SPECIES in % OF SOLIDS                    !
!           - FLUX of Organic Matter in Micro-Mole/(cm^2)/yr                   !
!                                                                              !
!------------------------------------------------------------------------------!

#include "aed2.h"

!
MODULE aed2_sedcandi

   USE aed2_gctypes
   USE aed2_gcsolver,   ONLY: nComponents, allComponents
   USE aed2_vode,       ONLY: DVODE

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: ConfigureCANDI,                                                   &
             InitialiseCANDI,                                                  &
             doCANDI, dfluxes, aed2_sed_candi_t
!  PUBLIC :: Y,                                                                &
!            CANDIVarNames,                                                    &
!            iSTEADY, SALT, TEMP, PRES,                                        &
!            fg0, wvel, uvel ,                                                 &
!            ps, poros, kg0var,                                                &
!            rpar, hco3y, hsy, FEX,                                            &
!            w00, db0, p00, x2, FixedBottomConc, IBC2,                         &
!            o2y, po4ly, PO4sy,feohy, FeOHBy, PartFluxes, IAP, KIAP, QIAP,     &
!            OutputUnits, deltaT, isSolid, FINswitch,                          &
!            timeswitch, num_days, fluxon, fluxoff, substep, driverDT

!###############################################################################
!   TYPE(AEDConstType) :: consts
!    TYPE(AEDConstType) :: c
 CHARACTER (LEN=4), PARAMETER :: CO_I="CO_I"  ! Constant depth distrib.        !
 CHARACTER (LEN=4), PARAMETER :: LI_I="LI_I"  ! Linear depth dependence        !
 CHARACTER (LEN=4), PARAMETER :: EX_I="EX_I"  ! Exponential depth depend.      !
!CHARACTER (LEN=4), PARAMETER :: IN_I="IN_I"  ! Multiple depth linear int.     !
!CHARACTER (LEN=4), PARAMETER :: FI_I="FI_I"  ! File input                     !
!CHARACTER (LEN=4), PARAMETER :: D2_I="D2_I"  ! 2D candi input                  !
!CHARACTER (LEN=4), PARAMETER :: RI_I="RI_I"  ! Along-river interpolation      !

 LOGICAL :: disableWQUpdate=.false. ! Disables WQ update with WQ3F             !

!#CAB# End of bits I've added .....
!###############################################################################


!------------------------------------------------------------------------------!
! Module level declarations                                                    !
  REAL(SEDP), PARAMETER :: ZERO  = 0.0D+00
  REAL(SEDP), PARAMETER :: HALF  = 0.5D+00
  REAL(SEDP), PARAMETER :: ONE   = 1.0D+00
  REAL(SEDP), PARAMETER :: TWO   = 2.0D+00
  REAL(SEDP), PARAMETER :: THREE = 3.0D+00
  REAL(SEDP), PARAMETER :: FOUR  = 4.0D+00
  REAL(SEDP), PARAMETER :: FIVE  = 5.0D+00
  REAL(SEDP), PARAMETER :: EIGHT = 8.0D+00
  REAL(SEDP), PARAMETER :: NINE  = 9.0D+00
  REAL(SEDP), PARAMETER :: TEN   = 1.0D+01
  REAL(SEDP), PARAMETER :: HUN   = 1.0D+02
  REAL(SEDP), PARAMETER :: YEAR  = 3.1536D+07

  !-- Number of simulated variables
  INTEGER, PARAMETER :: Mainmod  = 8    ! Compulsory variables

!------------------------------------------------------------------------------!

 CONTAINS

!------------------------------------------------------------------------------!
! configureCANDI                                                               !
!                                                                              !
! Subroutine written to allow dynamic setup & configuration of CANDI module    !
! during setup of an AED run.                                                  !
!                                                                              !
! Sets grid size, simulated variables and allocates necessary space for        !
! local arrays                                                                 !
!------------------------------------------------------------------------------!
 FUNCTION ConfigureCANDI(dia,nSedLayers,nCANDIVarsToSim,Names,VarSolidStatus,nSedCols) RESULT(candi)
   !-- Incoming                            ! Where do these come in from? - Dan
   TYPE(AEDConstDiagenesisType),     INTENT(IN) :: dia
   INTEGER,                          INTENT(IN) :: nSedLayers
   INTEGER,                          INTENT(IN) :: nCANDIVarsToSim
   LOGICAL,            DIMENSION(:), INTENT(IN) :: VarSolidStatus
   CHARACTER(LEN=*),   DIMENSION(:), INTENT(IN) :: Names
   !-- Outgoing
   INTEGER,                          INTENT(OUT):: nSedCols
   !-- Local
   TYPE(aed2_sed_candi_t),POINTER :: candi
   INTEGER  :: i,j,nNonTracers
   REAL(SEDP) :: v

   !---------------------------------------------------------------------------!
   !-- Check nSedLayers is suitable
   !IF(nSedLayers < 3) THEN
   !  PRINT *,' AED cannot currently configure CANDI to have <3 vertical layers'
   !  STOP
   !END IF
   IF(nSedLayers > 999) THEN
     PRINT *,' AED cannot currently configure CANDI to have >999 vertical layers'
     STOP
   END IF

   ALLOCATE(candi)

   ALLOCATE(candi%bottom(0:nCANDIVarsToSim-1))
   ALLOCATE(candi%FixedBottomConc(0:nCANDIVarsToSim-1))
   ALLOCATE(candi%PartFluxes(0:nCANDIVarsToSim-1))

   !-- Set local size variables,
   !MAXNPTS  = nSedLayers
   candi%nSPECIES = nCANDIVarsToSim
   nSedCols = candi%nSPECIES
   print *,'**************************************************************'
   print *,'nCANDIVarsToSim',nCANDIVarsToSim
   print *,'nSpecies',candi%nSpecies
   print *,'nSedCols',nSedCols
   print *,'**************************************************************'
   !-- Set local constants
   CALL SetLocalConstants(candi,dia)

   !-- Set local simulated variable names
   ALLOCATE(candi%CANDIVarNames(0:candi%nSPECIES-1))
   candi%CANDIVarNames(0:candi%nSPECIES-1) = Names(1:nCANDIVarsToSim)
   print *,'nSPECIES', candi%nSPECIES

   !-- Set local "isSolid" array to AED version
   ALLOCATE(candi%isSolid(0:nCANDIVarsToSim-1))
   IF(SIZE(VarSolidStatus) == nCANDIVarsToSim) THEN
     candi%isSolid = VarSolidStatus
   ELSE
     STOP " Size Mismatch in configureCANDI - VarSolidStatus"
   END IF

   !---------------------------------------------------------------------------!
   !-- Set local array inidcies to their respective columns
   !---------------------------------------------------------------------------!

   !---------------------------------------------------------------------------!
   !-------- COMPULSORY COMPONENTS ------------
   !dp -- Mainmod = O2 NO3 SO4 PO4 NH4 CH4 HS HCO3
   !dp -- kgmod = OM pools of the three OM models
   !
   !
   !dp -- Mainmod
   candi%O2y    = GetAEDColNum(candi,"oxy")
   IF(candi%O2y < 0 .OR. candi%O2y > nCANDIVarsToSim-1) THEN
     WRITE(*,*) " Compulsory CANDI Variable not simulated: DO. STOPPING"
     STOP
   END IF
   candi%NO3y   = GetAEDColNum(candi,"nit")
   IF(candi%NO3y < 0 .OR. candi%NO3y > nCANDIVarsToSim-1) THEN
     WRITE(*,*) " Compulsory CANDI Variable not simulated: NO3. STOPPING"
     STOP
   END IF
   candi%SO4y   = GetAEDColNum(candi,"so4")
   IF(candi%SO4y < 0 .OR. candi%SO4y > nCANDIVarsToSim-1) THEN
     WRITE(*,*) " Compulsory CANDI Variable not simulated: SO4. STOPPING"
     STOP
   END IF
   candi%PO4ly  = GetAEDColNum(candi,"frp")
   IF(candi%PO4ly < 0 .OR. candi%PO4ly > nCANDIVarsToSim-1) THEN
     WRITE(*,*) " Compulsory CANDI Variable not simulated: PO4. STOPPING"
     STOP
   END IF
   candi%NH4y   = GetAEDColNum(candi,"amm")
   IF(candi%NH4y < 0 .OR. candi%NH4y > nCANDIVarsToSim-1) THEN
     WRITE(*,*) " Compulsory CANDI Variable not simulated: NH4. STOPPING"
     STOP
   END IF
   candi%CH4y   = GetAEDColNum(candi,"ch4")
   IF(candi%CH4y < 0 .OR. candi%CH4y > nCANDIVarsToSim-1) THEN
     WRITE(*,*) " Compulsory CANDI Variable not simulated: CH4. STOPPING"
     STOP
   END IF
   candi%HSy    = GetAEDColNum(candi,"h2s")
   IF(candi%HSy < 0 .OR. candi%HSy > nCANDIVarsToSim-1) THEN
     WRITE(*,*) " Compulsory CANDI Variable not simulated: H2S. STOPPING"
     STOP
   END IF
   candi%HCO3y  = GetAEDColNum(candi,"dic")
   IF(candi%HCO3y < 0 .OR. candi%HCO3y > nCANDIVarsToSim-1) THEN
     WRITE(*,*) " Compulsory CANDI Variable not simulated: DIC. STOPPING"
     STOP
   END IF

   !-------- kgmod ------------ OMModel1
   IF(candi%param%OMModel < 1 .OR. candi%param%OMModel > 3) THEN
     PRINT *, ('OMModel must be 1, 2 or 3. STOPPING')
     STOP
   END IF

   IF(candi%param%OMModel == 1) THEN
     !-- Old fashioned POM->CO2 style
     candi%POMLy        = GetAEDColNum(candi,"poml")
     candi%POMRy        = GetAEDColNum(candi,"pomr")
     candi%POMspecialy  = GetAEDColNum(candi,"pomspecial")
     candi%kgmod = 3

     ! Unused in this model
     candi%DOCLy      = -1
     candi%DONLy      = -1
     candi%DOPLy      = -1
     candi%POCLy      = -1
     candi%PONLy      = -1
     candi%POPLy      = -1
     candi%DOCRy      = -1
     candi%DONRy      = -1
     candi%DOPRy      = -1
     candi%POCRy      = -1
     candi%PONRy      = -1
     candi%POPRy      = -1

     candi%DOMRy      = -1
     candi%dhydy      = -1
     candi%OAcy       = -1
     candi%H2y        = -1
     candi%POM1y      = -1
     candi%POM2y      = -1
     candi%POM3y      = -1
     candi%POM4y      = -1

     candi%BMety      = -1
     candi%BSuly      = -1
     candi%BIroy      = -1
     candi%BMany      = -1
     candi%BDeny      = -1
     candi%BAery      = -1
     candi%BFery      = -1
     candi%Necromassy = -1

!-------- kgmod ------------ OMModel2
     ELSEIF(candi%param%OMmodel == 2) THEN
     !-- CAEDYM style Laurie era
     candi%DOCLy  = GetAEDColNum(candi,"docl")
     candi%DONLy  = GetAEDColNum(candi,"donl")
     candi%DOPLy  = GetAEDColNum(candi,"dopl")
     candi%POCLy  = GetAEDColNum(candi,"pocl")
     candi%PONLy  = GetAEDColNum(candi,"ponl")
     candi%POPLy  = GetAEDColNum(candi,"popl")
     candi%DOCRy  = GetAEDColNum(candi,"docr")
     IF(candi%DOCRy < 0 .OR. candi%DOCRy > nCANDIVarsToSim-1) THEN       ! Why is there a condition just for DOCR? - Dan
       candi%simRefOM = .FALSE.
       candi%DONRy  = -1        ! Set these to -1
       candi%DOPRy  = -1        ! as a warning that
       candi%POCRy  = -1        ! they are mistakenly
       candi%PONRy  = -1        ! activated.
       candi%POPRy  = -1
     ELSE
       candi%simRefOM = .TRUE.
       candi%DONRy  = GetAEDColNum(candi,"donr")
       candi%DOPRy  = GetAEDColNum(candi,"dopr")
       candi%POCRy  = GetAEDColNum(candi,"pocr")
       candi%PONRy  = GetAEDColNum(candi,"ponr")
       candi%POPRy  = GetAEDColNum(candi,"popr")
     END IF ! End if DOCR bad
     IF (candi%simRefOM) THEN
       candi%kgmod = 12
     ELSE
       candi%kgmod = 6
     END IF ! End if simRefOM
     candi%POMLy              = -1
     candi%POMRy              = -1
     candi%POMspecialy        = -1
     candi%DOMRy              = -1
     candi%dhydy              = -1
     candi%OAcy               = -1
     candi%H2y                = -1
     candi%POM1y              = -1
     candi%POM2y              = -1
     candi%POM3y              = -1
     candi%POM4y              = -1

     candi%BMety              = -1
     candi%BSuly              = -1
     candi%BIroy              = -1
     candi%BMany              = -1
     candi%BDeny              = -1
     candi%BAery              = -1
     candi%BFery              = -1
     candi%Necromassy         = -1

   !-------- kgmod ------------
   ELSEIF(candi%param%OMmodel == 3) THEN
     candi%POMLy     = GetAEDColNum(candi,"poml")              !Dan added
     candi%POMRy     = GetAEDColNum(candi,"pomr")              !Dan added
     candi%DOMRy     = GetAEDColNum(candi,"domr")              !Dan added
     candi%dhydy     = GetAEDColNum(candi,"dhyd")              !Dan added
     candi%OAcy      = GetAEDColNum(candi,"OAc" )              !Dan added
     candi%H2y       = GetAEDColNum(candi,"H2"  )              !Dan added
     candi%POM1y     = GetAEDColNum(candi,"POM1")              !Dan added
     candi%POM2y     = GetAEDColNum(candi,"POM2")              !Dan added
     candi%POM3y     = GetAEDColNum(candi,"POM3")              !Dan added
     candi%POM4y     = GetAEDColNum(candi,"POM4")              !Dan added

     candi%BMety     = GetAEDColNum(candi,"BMet")              !Dan added
     candi%BSuly     = GetAEDColNum(candi,"BSul")              !Dan added
     candi%BIroy     = GetAEDColNum(candi,"BIro")              !Dan added
     candi%BMany     = GetAEDColNum(candi,"BMan")              !Dan added
     candi%BDeny     = GetAEDColNum(candi,"BDen")              !Dan added
     candi%BAery     = GetAEDColNum(candi,"BAer")              !Dan added
     candi%BFery     = GetAEDColNum(candi,"BFer")              !Dan added
     candi%Necromassy= GetAEDColNum(candi,"Necromass")         !Dan added
     candi%POMspecialy =GetAEDColNum(candi,"pomspecial")       !Dan added
     candi%kgmod = 20

     ! Unused in this model
     candi%DOCLy      = -1
     candi%DONLy      = -1
     candi%DOPLy      = -1
     candi%POCLy      = -1
     candi%PONLy      = -1
     candi%POPLy      = -1
     candi%DOCRy      = -1
     candi%DONRy      = -1
     candi%DOPRy      = -1
     candi%POCRy      = -1
     candi%PONRy      = -1
     candi%POPRy      = -1

  !ELSE

   END IF ! End if OMModel = 1, 2 or 3
   nNonTracers = Mainmod + candi%kgmod
   !-------- (END COMPULSORY COMPONENTS) ------------



   !---------------------------------------------------------------------------!
   !-- OPTIONAL COMPONENTS

   !------------------------------------------
   !-- Mn and Fe
   candi%param%simMnFe = .FALSE.

   candi%MnO2y   = GetAEDColNum(candi,"mno2a")
   candi%MnO2By  = GetAEDColNum(candi,"mno2b")
   candi%FeOHy   = GetAEDColNum(candi,"feoh3a")
   candi%FeOHBy  = GetAEDColNum(candi,"feoh3b")
   candi%MnIIy   = GetAEDColNum(candi,"mnii")
   candi%MnIVy   = GetAEDColNum(candi,"mniv")
   candi%FeIIy   = GetAEDColNum(candi,"feii")
   candi%FeIIIy  = GetAEDColNum(candi,"feiii")

   candi%param%simMnFe = .TRUE.

   nNonTracers = nNonTracers + 7 ! Was 7, might have been wrong - Dan

   candi%mnfemod = 8 ! 8 Was 8, changed to 7 - Dan

   IF(.NOT. candi%param%simMnFe) THEN
     candi%MnO2y   =-1
     candi%MnO2By  =-1
     candi%FeOHy   =-1
     candi%FeOHBy  =-1
     candi%MnIIy   =-1
     candi%MnIVy   =-1                ! Dan added: might explain the discrepancy.
     candi%FeIIy   =-1
     candi%FeIIIy  =-1

     candi%mnfemod = 0
   END IF
 !Check whether user has requested "simMnFe" for mnfe kinetics ... do this later

   !------------------------------------------
   !-- FeS and FeS2 (Pyrite)
   candi%param%simFeS = .FALSE.

       ! FeS is simulated
       candi%param%simFeS = .TRUE.

       candi%FeSY    = GetAEDColNum(candi,"fes")
       candi%FeS2Y   = GetAEDColNum(candi,"fes2")

       nNonTracers = nNonTracers + 2

       candi%fesmod = 2

   IF(.NOT. candi%param%simFeS) THEN
     candi%FeSY    =-1
     candi%FeS2Y   =-1

     candi%fesmod  = 0
   END IF


   !------------------------------------------
   !-- Adsorption of PO4/NH4
   candi%PO4sY   = GetAEDColNum(candi,"pip")
   candi%NH4sY   = GetAEDColNum(candi,"pin")
   IF(candi%PO4sY < 0 .OR. candi%PO4sY > nCANDIVarsToSim-1) THEN
     candi%PO4sy   =-1
     candi%NH4sy   =-1
     candi%po4smod = 0
   ELSE
     ! Adsorbed PO4 & NH4 is simulated
     candi%simAdsp = .TRUE.
     nNonTracers = nNonTracers + 2
     candi%po4smod = 2
   END IF

   !------------------------------------------
   !-- Calcite
   candi%param%simCaCO3 = .FALSE.
       ! Calcite is simulated
       candi%param%simCaCO3 = .TRUE.
       candi%CaY     = GetAEDColNum(candi,"ca")
       candi%calY    = GetAEDColNum(candi,"caco3")
       candi%ARAY    =-1
       nNonTracers = nNonTracers + 2
       candi%CaCO3mod = 2
   IF(.NOT. candi%param%simCaCO3) THEN
     candi%CaY       =-1
     candi%ARAY      =-1
     candi%CALY      =-1
     candi%CaCO3mod  = 0
   END IF

   !------------------------------------------
   !-- Siderite
   candi%param%simFeCO3 = .FALSE.
       ! Siderite is simulated
       candi%param%simFeCO3 = .TRUE.
       candi%sidY    = GetAEDColNum(candi,"feco3")
       nNonTracers = nNonTracers + 1
       candi%FeCO3mod = 1
   IF(.NOT. candi%param%simFeCO3) THEN
     candi%sidY      =-1
     candi%FeCO3mod  = 0
   END IF

   !------------------------------------------
   !-- Rhodocrosite
   candi%param%simMnCO3 = .FALSE.
       ! Siderite is simulated
       candi%param%simMnCO3 = .TRUE.
       candi%rodY    = GetAEDColNum(candi,"mnco3")
       nNonTracers = nNonTracers + 1
       candi%MnCO3mod = 1
   IF(.NOT. candi%param%simFeCO3) THEN
     candi%rodY      =-1
     candi%MnCO3mod  = 0
   END IF

   !------------------------------------------
   !-- Metal (X) and its sulfide (XS)
   candi%param%simX  = .FALSE.
   candi%simXS = .FALSE.
   candi%simXO = .FALSE.

   IF(.NOT. candi%param%simX) THEN
     candi%XY      =-1
     candi%XSY     =-1
     candi%POXY    =-1
     candi%DOXY    =-1

     candi%xmod    = 0
   END IF

   !------------------------------------------
   !-- Legacy C.CANDI switches
   candi%simC12 = .FALSE.
   candi%simFeII = .FALSE.
   candi%simCaCO3C12 = .FALSE.
   candi%simCaXCO3 = .FALSE.
   ! Fix these pH models to OFF - to be removed
   candi%alkmod = 0
   candi%rmmod  = 0


   !------------------------------------------
   !-- Non-active (tracer) variable definition
   candi%TracerMod = nCANDIVarsToSim - nNonTracers
   ALLOCATE(candi%TracerY(candi%Tracermod))  ! records tracer col numbers in "Y"

   IF(candi%TracerMod>1) THEN
     candi%simTracer = .TRUE.
   ELSE
     candi%simTracer = .FALSE.
   END IF
   DO i = 1,candi%Tracermod
     candi%TracerY(i) = Mainmod+candi%kgmod+candi%mnfemod+candi%alkmod+candi%rmmod+candi%FEIImod &
                + candi%po4smod+candi%FESmod+candi%CaXCO3mod+candi%CaCO3mod+candi%FeCO3mod+candi%MnCO3mod+candi%xmod &
                + candi%C12mod+candi%CaCO3C12mod + (i-1)
   END DO


! H+
   candi%protony = (nCANDIVarsToSim-3)+1-1


   !---------------------------------------------------------------------------!
   !-- Now finalise configuration
   !---------------------------------------------------------------------------!

   candi%MAXNEQ   = (candi%nSPECIES*candi%MAXNPTS) + 100
   candi%MU       = candi%nSPECIES
   candi%ML       = candi%nSPECIES
   candi%LRW      = 22 + 11*candi%MAXNEQ+(3*candi%ML+2*candi%MU)*candi%MAXNEQ
   candi%LIW      = 30 + candi%MAXNEQ

   !---------------------------------------------------------------------------!
   !-- Now allocate space required for parameter arrays
   !---------------------------------------------------------------------------!
   ALLOCATE(candi%DIFFC(0:candi%nSPECIES-1,candi%MAXNPTS))
   ALLOCATE(candi%RGC(candi%maxnpts))
   ALLOCATE(candi%RGN(candi%maxnpts))
   ALLOCATE(candi%RGP(candi%maxnpts))
   ALLOCATE(candi%RGX(candi%maxnpts))
   ALLOCATE(candi%ROX(candi%maxnpts))

   ALLOCATE(candi%FO2(candi%maxnpts))
   ALLOCATE(candi%FNO3(candi%maxnpts))
   ALLOCATE(candi%FMnO2(candi%maxnpts))
   ALLOCATE(candi%FFeOH(candi%maxnpts))
   ALLOCATE(candi%FSO4(candi%maxnpts))
   ALLOCATE(candi%FMet(candi%maxnpts))
   ALLOCATE(candi%FBHyd(candi%maxnpts))
   ALLOCATE(candi%FDHyd(candi%maxnpts))
   ALLOCATE(candi%FOAc(candi%maxnpts))
   ALLOCATE(candi%FH2(candi%maxnpts))

   ALLOCATE(candi%RPOM1(candi%maxnpts))
   ALLOCATE(candi%RPOM2(candi%maxnpts))
   ALLOCATE(candi%RPOM3(candi%maxnpts))
   ALLOCATE(candi%RPOM4(candi%maxnpts))
   ALLOCATE(candi%RPOMspecial(candi%maxnpts))
   ALLOCATE(candi%RNecro(candi%maxnpts))
   ALLOCATE(candi%ROAc(candi%maxnpts))
   ALLOCATE(candi%RH2(candi%maxnpts))
   ALLOCATE(candi%RAerDHyd(candi%maxnpts))
   ALLOCATE(candi%RDenO2DHyd(candi%maxnpts))
   ALLOCATE(candi%RDenNO3DHyd(candi%maxnpts))
   ALLOCATE(candi%RFerDHyd(candi%maxnpts))
   ALLOCATE(candi%RDHyd(candi%maxnpts))

   ALLOCATE(candi%dGFerDHyd(candi%maxnpts))
   ALLOCATE(candi%FTFerDHyd(candi%maxnpts))
   ALLOCATE(candi%dGAerOAc(candi%maxnpts))
   ALLOCATE(candi%FTAerOAc(candi%maxnpts))
   ALLOCATE(candi%dGDenOAc(candi%maxnpts))
   ALLOCATE(candi%FTDenOAc(candi%maxnpts))
   ALLOCATE(candi%dGDenH2(candi%maxnpts))
   ALLOCATE(candi%FTDenH2(candi%maxnpts))
   ALLOCATE(candi%dGManOAc(candi%maxnpts))
   ALLOCATE(candi%FTManOAc(candi%maxnpts))
   ALLOCATE(candi%dGIroOAc(candi%maxnpts))
   ALLOCATE(candi%FTIroOAc(candi%maxnpts))
   ALLOCATE(candi%dGIroH2(candi%maxnpts))
   ALLOCATE(candi%FTIroH2(candi%maxnpts))
   ALLOCATE(candi%dGSulOAc(candi%maxnpts))
   ALLOCATE(candi%FTSulOAc(candi%maxnpts))
   ALLOCATE(candi%dGSulH2(candi%maxnpts))
   ALLOCATE(candi%FTSulH2(candi%maxnpts))
   ALLOCATE(candi%dGMetOAc(candi%maxnpts))
   ALLOCATE(candi%FTMetOAc(candi%maxnpts))
   ALLOCATE(candi%dGMetH2(candi%maxnpts))
   ALLOCATE(candi%FTMetH2(candi%maxnpts))

   ALLOCATE(candi%RdeathFer(candi%maxnpts))
   ALLOCATE(candi%RdeathAer(candi%maxnpts))
   ALLOCATE(candi%RdeathDen(candi%maxnpts))
   ALLOCATE(candi%RdeathMan(candi%maxnpts))
   ALLOCATE(candi%RdeathIro(candi%maxnpts))
   ALLOCATE(candi%RdeathSul(candi%maxnpts))
   ALLOCATE(candi%RdeathMet(candi%maxnpts))
   ALLOCATE(candi%RdeathTot(candi%maxnpts))

   ALLOCATE(candi%RAerOAc(candi%maxnpts))
   ALLOCATE(candi%RDenOAc(candi%maxnpts))
   ALLOCATE(candi%RManOAc(candi%maxnpts))
   ALLOCATE(candi%RIroOAc(candi%maxnpts))
   ALLOCATE(candi%RSulOAc(candi%maxnpts))
   ALLOCATE(candi%RMetOAc(candi%maxnpts))
   ALLOCATE(candi%RAerH2(candi%maxnpts))
   ALLOCATE(candi%RDenH2(candi%maxnpts))
   ALLOCATE(candi%RManH2(candi%maxnpts))
   ALLOCATE(candi%RIroH2(candi%maxnpts))
   ALLOCATE(candi%RSulH2(candi%maxnpts))
   ALLOCATE(candi%RMetH2(candi%maxnpts))
   ALLOCATE(candi%RNH4OX(candi%maxnpts))
   ALLOCATE(candi%RMnOX(candi%maxnpts))
   ALLOCATE(candi%RFeOX(candi%maxnpts))
   ALLOCATE(candi%RTSOX(candi%maxnpts))
   ALLOCATE(candi%RCH4OX(candi%maxnpts))
   ALLOCATE(candi%RFeSOX(candi%maxnpts))
   ALLOCATE(candi%RFeS2OX(candi%maxnpts))

   ALLOCATE(candi%RNH4NO2(candi%maxnpts))
   ALLOCATE(candi%RMnNO3(candi%maxnpts))
   ALLOCATE(candi%RFeNO3(candi%maxnpts))
   ALLOCATE(candi%RTSNO3(candi%maxnpts))
   ALLOCATE(candi%RFeMnA(candi%maxnpts))
   ALLOCATE(candi%RFeMnB(candi%maxnpts))
   ALLOCATE(candi%RTSMnA(candi%maxnpts))
   ALLOCATE(candi%RTSMnB(candi%maxnpts))
   ALLOCATE(candi%RFeSMnA(candi%maxnpts))
   ALLOCATE(candi%RFeSMnB(candi%maxnpts))
   ALLOCATE(candi%RTSFeA(candi%maxnpts))
   ALLOCATE(candi%RTSFeB(candi%maxnpts))
   ALLOCATE(candi%RFeSFeA(candi%maxnpts))
   ALLOCATE(candi%RFeSFeB(candi%maxnpts))
   ALLOCATE(candi%RCH4SO4(candi%maxnpts))

   ALLOCATE(candi%RMnAge(candi%maxnpts))
   ALLOCATE(candi%RFeOHAppt(candi%maxnpts))
   ALLOCATE(candi%RFeOHBppt(candi%maxnpts))
   ALLOCATE(candi%RFeAge(candi%maxnpts))
   ALLOCATE(candi%RFeSppt(candi%maxnpts))
   ALLOCATE(candi%RPyrite(candi%maxnpts))
   ALLOCATE(candi%RXSppt(candi%maxnpts))
   ALLOCATE(candi%RSidppt(candi%maxnpts))
   ALLOCATE(candi%RMnO2Appt(candi%maxnpts))
   ALLOCATE(candi%RMnO2Bppt(candi%maxnpts))
   ALLOCATE(candi%RRodppt(candi%maxnpts))
   ALLOCATE(candi%RCalppt(candi%maxnpts))
   ALLOCATE(candi%RPO4ads(candi%maxnpts))
   ALLOCATE(candi%RNH4ads(candi%maxnpts))

   ALLOCATE(candi%reac(0:candi%nSPECIES-1,candi%MAXNPTS))
   ALLOCATE(candi%Y(0:candi%nSPECIES-1,candi%MAXNPTS))
   ALLOCATE(candi%Ydot(0:candi%nSPECIES-1,candi%MAXNPTS))
   ALLOCATE(candi%Ytemp(0:candi%nSPECIES-1,candi%MAXNPTS))
   ALLOCATE(candi%IAP(0:candi%nSPECIES-1,candi%MAXNPTS))
   ALLOCATE(candi%KIAP(0:candi%nSPECIES-1,candi%MAXNPTS))        !Dan added
   ALLOCATE(candi%QIAP(0:candi%nSPECIES-1,candi%MAXNPTS))        !Dan added
   !ALLOCATE(candi%IAP(candi%MAXNPTS,1:candi%nSPECIES))
   ALLOCATE(candi%kg0var(candi%maxnpts))
   ALLOCATE(candi%kgpo4up(candi%maxnpts))
   ALLOCATE(candi%Btot(candi%maxnpts))
   ALLOCATE(candi%pps(candi%maxnpts)    )
   ALLOCATE(candi%psp(candi%maxnpts)    )
   ALLOCATE(candi%poros(candi%maxnpts)  )
   ALLOCATE(candi%t2(candi%maxnpts)     )
   ALLOCATE(candi%dpdx(candi%maxnpts)   )
   ALLOCATE(candi%dt2dx(candi%maxnpts)  )
   ALLOCATE(candi%uvel(candi%maxnpts)   )
   ALLOCATE(candi%wvel(candi%maxnpts)   )
   ALLOCATE(candi%bioturb(candi%maxnpts))
   ALLOCATE(candi%dbdx(candi%maxnpts)   )
   ALLOCATE(candi%ps(candi%maxnpts)     )
   ALLOCATE(candi%cirrig(candi%maxnpts) )
   ALLOCATE(candi%dff(candi%maxnpts)    )
   ALLOCATE(candi%dudx(candi%maxnpts)   )
   ALLOCATE(candi%poros_bg(candi%maxnpts) )
   ALLOCATE(candi%poros_dot(candi%maxnpts))

   ALLOCATE(candi%top_bound(0:candi%nSPECIES-1))
   ALLOCATE(candi%dh(candi%maxnpts))
   ALLOCATE(candi%dh2(candi%maxnpts))
   ALLOCATE(candi%dhr(candi%MAXNPTS))
   ALLOCATE(candi%dhf(candi%MAXNPTS))
   ALLOCATE(candi%RPAR(candi%maxnpts))
   ALLOCATE(candi%rtol(candi%lrw))
   ALLOCATE(candi%atol(candi%liw))
           WRITE(*,"(6X,'Monkeys6')")
   !---------------------------------------------------------------------------!
   !-- Set the conservation of mass check varilable to zero
   candi%lost   = 0.0
   candi%oldday = 1

   !-- Set the tolerances for VODE
   candi%rtol  = 1.0D-6
   candi%atol  = 1.0D-9

   !-- Set grid dimensions
   CALL SetGrid(candi,candi%param%xl)

   !---------------------------------------------------------------------------!
   !-- Report Config
   !---------------------------------------------------------------------------!
   WRITE(*,"(8X,'AED Sediment Diagenesis Model Configuration: ')")
   WRITE(*,"(9X,'--- CANDI-AED v1.0 --- ')")
   WRITE(*,"(10X,'Refractory OM,             simRefOM   : ',L1)")candi%simRefOM
   WRITE(*,"(10X,'Mn and Fe,                 simMnFe    : ',L1)")candi%param%simMnFe
   WRITE(*,"(10X,'FeS and FeS2,              simFeS     : ',L1)")candi%param%simFeS
   WRITE(*,"(10X,'Calcite,                   simCaCO3   : ',L1)")candi%param%simCaCO3
   WRITE(*,"(10X,'Siderite,                  simFeCO3   : ',L1)")candi%param%simFeCO3
   WRITE(*,"(10X,'Rhodchrosite,              simMnCO3   : ',L1)")candi%param%simMnCO3
   WRITE(*,"(10X,'Trace metal species,       simX       : ',L1)")candi%param%simX
   WRITE(*,"(10X,'Metal sulfide,             simXS      : ',L1)")candi%simXS
   WRITE(*,"(10X,'Biologically active metal, simXO      : ',L1)")candi%simXO
   WRITE(*,"(10X,'PO4 and NH4 asdsorption,   simAdsp    : ',L1)")candi%simAdsp
   WRITE(*,"(10X,'C12 fraction,              simC12     : ',L1)")candi%simC12
   WRITE(*,"(10X,'Mn/Fe in CaCO3,            simCaXCO3  : ',L1)")candi%simCaXCO3
   WRITE(*,"(10X,'C12 calcite fraction,      simCaCO3C12: ',L1)")candi%simCaCO3C12
   WRITE(*,"(10X,'Tracer/Inactive component, simTracer  : ',L1,',',I3)")candi%simTracer, candi%TracerMod
   WRITE(*,"(10X,'Kinetic Reaction Mode,     rxn_mode   : ',I3,/)")candi%param%rxn_mode

   WRITE(*,"(10X,'OM Model                   OMModel    : ',I1)")candi%param%OMModel
   WRITE(*,"(10X,'FBIO                       FBIOswitch : ',I1)")candi%param%FBIOswitch
   WRITE(*,"(10X,'FT                         FTswitch   : ',I1)")candi%param%FTswitch
   WRITE(*,"(10X,'FTem                       FTemswitch : ',I1)")candi%param%FTemswitch
   WRITE(*,"(10X,'FIN                        FINswitch  : ',I1)")candi%param%FINswitch
   WRITE(*,"(10X,'FIn O2 only                FInO2OnlySwitch:',I1)")candi%param%FInO2OnlySwitch
   WRITE(*,"(10X,'FOM                        FOMswitch  : ',I1)")candi%param%FOMswitch
   WRITE(*,"(10X,'OM Approach                OMapproach : ',I1,/)")candi%param%OMapproach
   WRITE(*,"(10X,'Bacteria solid?            Bsolid     : ',I1,/)")candi%param%Bsolidswitch

   WRITE(*,"(10X,'Compulsory variables       mainmod    : ',I1)")mainmod
   WRITE(*,"(10X,'OM phases involved         kgmod      : ',I2)")candi%kgmod
   WRITE(*,"(10X,'Mn and Fe                  mnfemod    : ',I2)")candi%mnfemod
   WRITE(*,"(10X,'FeSmod                     fesmod     : ',I2)")candi%fesmod
   WRITE(*,"(10X,'Phosphate                  po4smod    : ',I2)")candi%po4smod
   WRITE(*,"(10X,'Calcite                    caco3      : ',I2)")candi%caco3mod
   WRITE(*,"(10X,'Siderite                   feco3mod   : ',I2)")candi%feco3mod
   WRITE(*,"(10X,'Rhodchrosite               mnco3      : ',I2)")candi%mnco3mod
   WRITE(*,"(10X,'Non-reactive species       tracermod  : ',I2)")candi%tracermod
   WRITE(*,"(10X,'Reactive species          nNonTracers : ',I2)")nNonTracers
   WRITE(*,"(10X,'Total species              nSpecies   : ',I2)")candi%nSpecies
   WRITE(*,"(10X,'Does tracermod + nNonTracers = nSpecies?',2/)")


   WRITE(*,"(9X,'Configuration Successful ')")

 END FUNCTION ConfigureCANDI
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!rl   ****************************************************************
!rl   Early diagenese model C. CANDI from R. Luff (1996-2003), based on
!rl   CANDI from  B.P. Boudreau.
!rl   ****************************************************************
!rl   AED VERSION
!rl   ****************************************************************
!rl   FILE: INIT.f
!rl   AIM : Definition of boundary conditions, initial conditions for
!rl         concentration profiles and model internal parameters.
!rl   ****************************************************************
!------------------------------------------------------------------------------!
 SUBROUTINE InitialiseCANDI(candi,bottomConcs,sedimentConcs)
 !SUBROUTINE InitialiseCANDI(bottomwaterConcs,sedimentConcs)
    !-- Incoming
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   REAL(SEDP), DIMENSION(0:candi%nSPECIES-1), INTENT(IN) :: bottomConcs
    !REAL(SEDP), DIMENSION(0:candi%nSPECIES-1), INTENT(IN) :: bottomwaterConcs
   REAL(SEDP), DIMENSION(0:candi%nSPECIES-1), INTENT(IN) :: sedimentConcs
   !-- Local
   INTEGER  :: i
   REAL(SEDP) :: tief
   REAL(SEDP) :: help1, help2
   REAL(SEDP) :: alpha1
   REAL(SEDP) :: fracfacorg

   candi%restartCANDI = .FALSE.
   !-- Calculate the porosity (P) and the Velocity (U,W) at as function of depth
   CALL SedProperties(candi)

   !#ifdef kgvar
   !kg0var(:)=kg0*exp(-rpar(:)**2.0/(2.0*x2)) &
   !          +0.16*(x2/w00+rpar(:)/wvel(:))**(-0.95)
   !#endif
   !#ifdef po4solid
   !kgpo4up(:)= kpo4up*exp(-rpar(:)**2.0/(2.0*x2))
   !#endif

   IF(candi%BC%IBC2 == 1) THEN
     ! Fixed surface OM concentrations - no flux rate
     candi%Fg0 = 0.0
     candi%Fg1 = 0.0
     candi%Fg2 = 0.0
   END IF

   !-- Set initial conditions
   IF (.NOT. candi%restartCANDI) THEN
     ! for Corg &1 and &2 and FeS
     fracfacorg = 1.0
     alpha1 = SQRT(candi%param%docl2dic/candi%param%db0)

     candi%Y(:,:)    = ZERO
     candi%Ydot(:,:) = ZERO

     DO i = 0,candi%nSPECIES-1
       ! All layers
       candi%Y(i,:) = sedimentConcs(i)

       ! Top layer reserved for bottom water condition
       candi%Y(i,1) = bottomConcs(i)
       !Y(i,1) = bottomwaterConcs(i)
       END DO

     !##IF(OMImitMethod == "LI_I") THEN
     !##  !-- Exponential decrease
     !##  CALL SetSedInitialCondition("LI_I",Y(DOCLy,:),sedimentConcs(DOCLy), OM_min*sedimentConcs(DOCLy),InitMinDepth)
     !##  CALL SetSedInitialCondition("LI_I",Y(POCLy,:),sedimentConcs(POCLy), OM_min*sedimentConcs(POCLy),InitMinDepth)
     !##  CALL SetSedInitialCondition("LI_I",Y(DOCRy,:),sedimentConcs(DOCRy), OM_min*sedimentConcs(DOCRy),InitMinDepth)
     !##  CALL SetSedInitialCondition("LI_I",Y(POCRy,:),sedimentConcs(POCRy), OM_min*sedimentConcs(POCRy),InitMinDepth)
     !##
     !##  CALL SetSedInitialCondition("LI_I",Y(DONLy,:),sedimentConcs(DONLy), OM_min*sedimentConcs(DONLy),InitMinDepth)
     !##  CALL SetSedInitialCondition("LI_I",Y(PONLy,:),sedimentConcs(PONLy), OM_min*sedimentConcs(PONLy),InitMinDepth)
     !##  CALL SetSedInitialCondition("LI_I",Y(DONRy,:),sedimentConcs(DONRy), OM_min*sedimentConcs(DONRy),InitMinDepth)
     !##  CALL SetSedInitialCondition("LI_I",Y(PONRy,:),sedimentConcs(PONRy), OM_min*sedimentConcs(PONRy),InitMinDepth)
     !##
     !##  CALL SetSedInitialCondition("LI_I",Y(DOPLy,:),sedimentConcs(DOPLy), OM_min*sedimentConcs(DOPLy),InitMinDepth)
     !##  CALL SetSedInitialCondition("LI_I",Y(POPLy,:),sedimentConcs(POPLy), OM_min*sedimentConcs(POPLy),InitMinDepth)
     !##  CALL SetSedInitialCondition("LI_I",Y(DOPRy,:),sedimentConcs(DOPRy), OM_min*sedimentConcs(DOPRy),InitMinDepth)
     !##  CALL SetSedInitialCondition("LI_I",Y(POPRy,:),sedimentConcs(POPRy), OM_min*sedimentConcs(POPRy),InitMinDepth)
     !##
     !##ELSE IF(OMImitMethod == "CO_I" .OR. OMImitMethod == "EX_I" ) THEN
     !##  !-- Constant
     !##  CALL SetSedInitialCondition(OMImitMethod,Y(DOCLy,:),sedimentConcs(DOCLy))
     !##  CALL SetSedInitialCondition(OMImitMethod,Y(DONLy,:),sedimentConcs(DONLy))
     !##  CALL SetSedInitialCondition(OMImitMethod,Y(DOPLy,:),sedimentConcs(DOPLy))
     !##
     !##  CALL SetSedInitialCondition(OMImitMethod,Y(POCLy,:),sedimentConcs(POCLy))
     !##  CALL SetSedInitialCondition(OMImitMethod,Y(PONLy,:),sedimentConcs(PONLy))
     !##  CALL SetSedInitialCondition(OMImitMethod,Y(POPLy,:),sedimentConcs(POPLy))
     !##
     !##  CALL SetSedInitialCondition(OMImitMethod,Y(DOCRy,:),sedimentConcs(DOCRy))
     !##  CALL SetSedInitialCondition(OMImitMethod,Y(DONRy,:),sedimentConcs(DONRy))
     !##  CALL SetSedInitialCondition(OMImitMethod,Y(DOPRy,:),sedimentConcs(DOPRy))
     !##
     !##  CALL SetSedInitialCondition(OMImitMethod,Y(POCRy,:),sedimentConcs(POCRy))
     !##  CALL SetSedInitialCondition(OMImitMethod,Y(PONRy,:),sedimentConcs(PONRy))
     !##  CALL SetSedInitialCondition(OMImitMethod,Y(POPRy,:),sedimentConcs(POPRy))
     !##
     !##END IF

     ! Hardcoded mixed-profile types for CdA sed candi
     !OM_min = 0.80
     !InitMinDepth = 100.0 !cm

     IF(candi%param%OMModel ==1) THEN
       CALL SetSedInitialCondition(candi,"LI_I",candi%Y(candi%POMRy,:),sedimentConcs(candi%POMRy), candi%BC%OM_minR*sedimentConcs(candi%POMRy),candi%BC%InitMinDepthR)
       CALL SetSedInitialCondition(candi,"EX_I",candi%Y(candi%POMLy,:),sedimentConcs(candi%POMLy))
       CALL SetSedInitialCondition(candi,"EX_I",candi%Y(candi%POMspecialy,:),sedimentConcs(candi%POMspecialy))

     ELSEIF(candi%param%OMModel ==2) THEN
       CALL SetSedInitialCondition(candi,"LI_I",candi%Y(candi%DOCRy,:),sedimentConcs(candi%DOCRy), candi%BC%OM_minR*sedimentConcs(candi%DOCRy),candi%BC%InitMinDepthR)
       CALL SetSedInitialCondition(candi,"LI_I",candi%Y(candi%POCRy,:),sedimentConcs(candi%POCRy), candi%BC%OM_minR*sedimentConcs(candi%POCRy),candi%BC%InitMinDepthR)

       CALL SetSedInitialCondition(candi,"LI_I",candi%Y(candi%DONRy,:),sedimentConcs(candi%DONRy), candi%BC%OM_minR*sedimentConcs(candi%DONRy),candi%BC%InitMinDepthR)
       CALL SetSedInitialCondition(candi,"LI_I",candi%Y(candi%PONRy,:),sedimentConcs(candi%PONRy), candi%BC%OM_minR*sedimentConcs(candi%PONRy),candi%BC%InitMinDepthR)

       CALL SetSedInitialCondition(candi,"LI_I",candi%Y(candi%DOPRy,:),sedimentConcs(candi%DOPRy), candi%BC%OM_minR*sedimentConcs(candi%DOPRy),candi%BC%InitMinDepthR)
       CALL SetSedInitialCondition(candi,"LI_I",candi%Y(candi%POPRy,:),sedimentConcs(candi%POPRy), candi%BC%OM_minR*sedimentConcs(candi%POPRy),candi%BC%InitMinDepthR)

       CALL SetSedInitialCondition(candi,"EX_I",candi%Y(candi%DOCLy,:),sedimentConcs(candi%DOCLy))
       CALL SetSedInitialCondition(candi,"EX_I",candi%Y(candi%DONLy,:),sedimentConcs(candi%DONLy))
       CALL SetSedInitialCondition(candi,"EX_I",candi%Y(candi%DOPLy,:),sedimentConcs(candi%DOPLy))

       CALL SetSedInitialCondition(candi,"EX_I",candi%Y(candi%POCLy,:),sedimentConcs(candi%POCLy))
       CALL SetSedInitialCondition(candi,"EX_I",candi%Y(candi%PONLy,:),sedimentConcs(candi%PONLy))
       CALL SetSedInitialCondition(candi,"EX_I",candi%Y(candi%POPLy,:),sedimentConcs(candi%POPLy))

     ELSEIF(candi%param%OMModel ==3) THEN
       CALL SetSedInitialCondition(candi,"LI_I",candi%Y(candi%POMRy,:),sedimentConcs(candi%POMRy), candi%BC%OM_minR*sedimentConcs(candi%POMRy),candi%BC%InitMinDepthR)
       !CALL SetSedInitialCondition(candi,"LI_I",Y(DOMRy,:),sedimentConcs(DOMRy), OM_minR*sedimentConcs(DOMRy),InitMinDepthR)

       CALL SetSedInitialCondition(candi,"EX_I",candi%Y(candi%DHydy,:),sedimentConcs(candi%dhydy))
       CALL SetSedInitialCondition(candi,"EX_I",candi%Y(candi%POMLy,:),sedimentConcs(candi%POMLy))

       !CALL SetSedInitialCondition(candi,"EX_I",candi%Y(candi%OAcy,:),sedimentConcs(candi%OAcy))
       !CALL SetSedInitialCondition(candi,"EX_I",candi%Y(candi%H2y,:),sedimentConcs(candi%H2y))
       !CALL SetSedInitialCondition(candi,"EX_I",candi%Y(candi%POM1y,:),sedimentConcs(candi%POM1y))
       !CALL SetSedInitialCondition(candi,"EX_I",candi%Y(candi%POM2y,:),sedimentConcs(candi%POM2y))
       !CALL SetSedInitialCondition(candi,"EX_I",candi%Y(candi%POM3y,:),sedimentConcs(candi%POM3y))
       !CALL SetSedInitialCondition(candi,"EX_I",candi%Y(candi%POM4y,:),sedimentConcs(candi%POM4y))

     END IF


     candi%FixedBottomConc(:) = candi%Y(:,1)

!     DO i = 1,maxnpts
!       IF(rpar(i) > InitSedDepth) THEN
!         Y(:,i) = ZERO
!       END IF
!     END DO
   ELSE
     !-- Restart
     ! CALL readauf(11)
   END IF
 END SUBROUTINE InitialiseCANDI
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! SetLocalConstants                                                            !
!                                                                              !
! Subroutine written to allow dynamic setup & configuration of CANDI module    !
! during setup of a AED run.                                                   !
!                                                                              !
!------------------------------------------------------------------------------!
 SUBROUTINE SetLocalConstants(candi,dia)
    IMPLICIT NONE
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi

    TYPE(AEDConstDiagenesisType), INTENT(IN) :: dia

    INTEGER  :: openStatus,i

    candi%time = dia%time

    candi%stcoef = dia%stcoef
!   candi%stcoef%fracCPL       = dia%stcoef%fracCPL               !Dan added
!   candi%stcoef%fracCPR       = dia%stcoef%fracCPR               !Dan added
!   candi%stcoef%fracCPspecial = dia%stcoef%fracCPspecial         !Dan added
!   candi%stcoef%fracNPL       = dia%stcoef%fracNPL               !Dan added
!   candi%stcoef%fracNPR       = dia%stcoef%fracNPR               !Dan added
!   candi%stcoef%fracNPspecial = dia%stcoef%fracNPspecial         !Dan added
!   candi%stcoef%fracPPL       = dia%stcoef%fraCPPL               !Dan added
!   candi%stcoef%fracPPR       = dia%stcoef%fracPPR               !Dan added
!   candi%stcoef%fracPPspecial = dia%stcoef%fracPPspecial         !Dan added

    candi%param = dia%param
    !# now changed those that were different
    candi%kNH4NO2       = ZERO                !dia%param%kNH4NO2
    candi%param%kPO4ads       = dia%param%kapo4
    candi%param%kNH4ads       = dia%param%kanh4

!# This lot can be removed when stuff works
#if 0
    candi%param%db0           = dia%param%db0
    candi%param%imix          = dia%param%imix
    candi%param%xs            = dia%param%xs
    candi%param%x1            = dia%param%x1
    candi%param%x2            = dia%param%x2
    candi%param%irrg          = dia%param%irrg
    candi%param%alpha0        = dia%param%alpha0
    candi%param%xirrig        = dia%param%xirrig
    candi%param%ventflow      = dia%param%ventflow
    candi%param%w00           = dia%param%w00
    candi%param%p0            = dia%param%p0
    candi%param%p00           = dia%param%p00
    candi%param%bp            = dia%param%bp
    candi%param%torteq        = dia%param%torteq
    candi%param%an            = dia%param%an
    candi%param%aa            = dia%param%aa
    candi%param%ab            = dia%param%ab
    candi%param%xl            = dia%param%xl
    candi%param%maxnpts       = dia%param%maxnpts  ! Deepest layer of sediment found in fabm.nml, for example, 50
    candi%param%OMModel       = dia%param%OMModel         ! Dan added
    candi%param%OMapproach    = dia%param%OMapproach      ! Dan added
    candi%param%FTemswitch    = dia%param%FTemswitch      ! Dan added
    candi%param%Bsolidswitch  = dia%param%Bsolidswitch    ! Dan added
    candi%param%FTswitch      = dia%param%FTswitch        ! Dan added
    candi%param%FBIOswitch    = dia%param%FBIOswitch      ! Dan added
    candi%param%FINswitch     = dia%param%FINswitch       ! Dan added
    candi%param%FInO2OnlySwitch = dia%param%FInO2OnlySwitch       ! Dan added
    candi%param%FOMswitch     = dia%param%FOMswitch       ! Dan added

    candi%param%pomr2dic       = dia%param%pomr2dic             ! Dan added
    candi%param%poml2dic       = dia%param%poml2dic             ! Dan added
    candi%param%pomspecial2dic = dia%param%pomspecial2dic       ! Dan added

    candi%param%docl2dic      = dia%param%docl2dic
    candi%param%donl2din      = dia%param%donl2din
    candi%param%dopl2dip      = dia%param%dopl2dip
    candi%param%pocl2docl     = dia%param%pocl2docl
    candi%param%ponl2donl     = dia%param%ponl2donl
    candi%param%popl2dopl     = dia%param%popl2dopl
    candi%param%docr2docl     = dia%param%docr2docl
    candi%param%donr2donl     = dia%param%donr2donl
    candi%param%dopr2dopl     = dia%param%dopr2dopl
    candi%param%pocr2docr     = dia%param%pocr2docr
    candi%param%ponr2donr     = dia%param%ponr2donr
    candi%param%popr2dopr     = dia%param%popr2dopr
    candi%param%pocvr2docr    = dia%param%pocvr2docr
    candi%param%ponvr2donr    = dia%param%ponvr2donr
    candi%param%popvr2dopr    = dia%param%popvr2dopr
      !OM Model 3 rate constants
    candi%param%domr2pomr     = dia%param%domr2pomr             !Dan added
    candi%param%poml2doml     = dia%param%poml2doml             !Dan added
   !candi%param%pomr2dic      = dia%param%pomr2dic              !Dan added
    candi%param%domr2dic      = dia%param%domr2dic              !Dan added

    candi%param%fracOAc      = dia%param%fracOAc
    candi%param%fracH2       = dia%param%fracH2
    candi%param%e            = dia%param%e
    candi%param%F            = dia%param%F
    candi%param%n            = dia%param%n
    candi%param%dPsi         = dia%param%dPsi
    candi%param%fuse         = dia%param%fuse
    candi%param%CellWeight   = dia%param%CellWeight
    candi%param%Tiny         = dia%param%Tiny
    candi%param%Temporary_proton = dia%param%Temporary_proton
    candi%param%KDHyd        = dia%param%KDHyd
    candi%param%KOAc         = dia%param%KOAc
    candi%param%KH2          = dia%param%KH2
    candi%param%kgrowthFer   = dia%param%kgrowthFer
    candi%param%kgrowthAer   = dia%param%kgrowthAer
    candi%param%kgrowthDen   = dia%param%kgrowthDen
    candi%param%kgrowthMan   = dia%param%kgrowthMan
    candi%param%kgrowthIro   = dia%param%kgrowthIro
    candi%param%kgrowthSul   = dia%param%kgrowthSul
    candi%param%kgrowthMet   = dia%param%kgrowthMet

    candi%param%kdeathFer    = dia%param%kdeathFer
    candi%param%kdeathAer    = dia%param%kdeathAer
    candi%param%kdeathDen    = dia%param%kdeathDen
    candi%param%kdeathMan    = dia%param%kdeathMan
    candi%param%kdeathIro    = dia%param%kdeathIro
    candi%param%kdeathSul    = dia%param%kdeathSul
    candi%param%kdeathMet    = dia%param%kdeathMet

    candi%param%kHyd1        = dia%param%kHyd1
    candi%param%kHyd2        = dia%param%kHyd2
    candi%param%kHyd3        = dia%param%kHyd3
    candi%param%kHyd4        = dia%param%kHyd4
    candi%param%kHydN        = dia%param%kHydN

     ! Stoichiometric coefficients
      ! OM model 1
    ! OM model 2
      !Nothing.
      ! OM model 3
    candi%param%FTR           = dia%param%FTR                   !Dan added
    candi%param%FTT           = dia%param%FTT                   !Dan added
    candi%param%deltaGATP     = dia%param%deltaGATP             !Dan added
    candi%param%BMax          = dia%param%BMax

    candi%param%YDHyAer      = dia%param%YDHyAer
    candi%param%YDHyFer      = dia%param%YDHyFer
    candi%param%YDenDHy      = dia%param%YDenDHy
    candi%param%YAerOAc      = dia%param%YAerOAc
    candi%param%YDenOAc      = dia%param%YDenOAc
    candi%param%YDenH2       = dia%param%YDenH2
    candi%param%YManOAc      = dia%param%YManOAc
    candi%param%YIroOAc      = dia%param%YIroOAc
    candi%param%YIroH2       = dia%param%YIroH2
    candi%param%YSulOAc      = dia%param%YSulOAc
    candi%param%YSulH2       = dia%param%YSulH2
    candi%param%YMetOAc      = dia%param%YMetOAc
    candi%param%YMetH2       = dia%param%YMetH2

    candi%param%dG0FerDHyd   = dia%param%dG0FerDHyd
    candi%param%dG0AerDHy    = dia%param%dG0AerDHy
    candi%param%dG0AerOAc    = dia%param%dG0AerOAc
    candi%param%dG0DenDHy    = dia%param%dG0DenDHy
    candi%param%dG0DenOAc    = dia%param%dG0DenOAc
    candi%param%dG0DenH2     = dia%param%dG0DenH2
    candi%param%dG0ManOAc    = dia%param%dG0ManOAc
    candi%param%dG0IroOAc    = dia%param%dG0IroOAc
    candi%param%dG0IroH2     = dia%param%dG0IroH2
    candi%param%dG0SulOAc    = dia%param%dG0SulOAc
    candi%param%dG0SulH2     = dia%param%dG0SulH2
    candi%param%dG0MetOAc    = dia%param%dG0MetOAc
    candi%param%dG0MetH2     = dia%param%dG0MetH2

    candi%stcoef%fracCDHyd    = dia%stcoef%fracCDHyd
    candi%stcoef%fracNDHyd    = dia%stcoef%fracNDHyd
    candi%stcoef%fracPDHyd    = dia%stcoef%fracPDHyd
    candi%stcoef%fracCOAc     = dia%stcoef%fracCOAc
    candi%stcoef%fracNOAc     = dia%stcoef%fracNOAc
    candi%stcoef%fracPOAc     = dia%stcoef%fracPOAc
    candi%stcoef%fracCH2      = dia%stcoef%fracCH2
    candi%stcoef%fracNH2      = dia%stcoef%fracNH2
    candi%stcoef%fracPH2      = dia%stcoef%fracPH2

    candi%param%kO2           = dia%param%ko2
    candi%param%kpO2          = dia%param%kpo2
    candi%param%lO2           = dia%param%lo2
    candi%param%lpO2          = dia%param%lpo2
    candi%param%kNO3          = dia%param%kno3
    candi%param%kpNO3         = dia%param%kpno3
    candi%param%lNO3          = dia%param%lno3
    candi%param%lpNO3         = dia%param%lpno3
    candi%param%kMnO2         = dia%param%kmn
    candi%param%kpMnO2        = dia%param%kpmn
    candi%param%lMnO2         = dia%param%lmn
    candi%param%lpMnO2        = dia%param%lpmn
    candi%param%kFeOH         = dia%param%kfe
    candi%param%kpFeOH        = dia%param%kpfe
    candi%param%lFeOH         = dia%param%lfe
    candi%param%lpFeOH        = dia%param%lpfe
    candi%param%kSO4          = dia%param%kso4
    candi%param%kpSO4         = dia%param%kpso4
    candi%param%lSO4          = dia%param%lso4
    candi%param%lpSO4         = dia%param%lpso4

    candi%param%kNH4OX        = dia%param%knh4ox
    candi%param%kMnOX         = dia%param%kmnox
    candi%param%kFeOX         = dia%param%kfeox
    candi%param%kTSOx         = dia%param%ktsox
    candi%param%kCH4OX        = dia%param%kch4ox
    candi%param%kFeSOX        = dia%param%kfesox
    candi%param%kFeS2OX       = dia%param%kfeS2ox

    candi%kNH4NO2       = ZERO                !dia%param%kNH4NO2
    candi%param%kMnNO3        = dia%param%kmnno3
    candi%param%kFeNO3        = dia%param%kfeno3  !dia%param%kFeNO3
    candi%param%ktsno3        = dia%param%ktsno3
    candi%param%kmnfe         = dia%param%kmnfe
    candi%param%ktsmn         = dia%param%ktsmn
    candi%param%kfesmn        = dia%param%kfesmn
    candi%param%ktsfe         = dia%param%ktsfe
    candi%param%kfesfe        = dia%param%kfesfe
    candi%param%kch4so4       = dia%param%kch4so4
    candi%param%kMnAge        = dia%param%kMnAge
!   candi%param%kMnAge        = dia%param%kfe2ox        !dia%param%kMnAge
    candi%param%kFeOHAppt     = dia%param%kFeOHAppt     !dia%param%kFeOHppt
    candi%param%kFeOHBppt     = dia%param%kFeOHBppt     !dia%param%kFeOHppt
    candi%param%kFeAge        = dia%param%kFeAge        !dia%param%kFeAge
!   candi%param%kFeAge        = dia%param%kfe1no3       !dia%param%kFeAge
    candi%param%kfesppt       = dia%param%kfesppt
    candi%param%kpyrite       = dia%param%kpyrite
    candi%param%kXSppt        = dia%param%kXSppt
    candi%param%kSidppt       = dia%param%kSidppt
    candi%param%kCalppt       = dia%param%kCalppt
    candi%param%kRodppt       = dia%param%kRodppt
    candi%param%kMnO2Appt     = dia%param%kMnO2Appt
    candi%param%kMnO2Bppt     = dia%param%kMnO2Bppt

    candi%param%kPO4ads       = dia%param%kapo4
    candi%param%kNH4ads       = dia%param%kanh4

    candi%param%rxn_mode   = dia%param%rxn_mode
#endif

    candi%BC = dia%BC

!   candi%BC%ibc2          = dia%BC%ibc2
!   candi%BC%ibbc          = dia%BC%ibbc
!   candi%BC%startSteady   = dia%BC%startSteady
!   candi%BC%flux_scale    = dia%BC%flux_scale
!   candi%BC%POMVR         = dia%BC%POMVR
!   candi%BC%OMInitMethodL = dia%BC%OMInitMethodL
!   candi%BC%OM_topL       = dia%BC%OM_topL
!   candi%BC%OM_minL       = dia%BC%OM_minL
!   candi%BC%OM_cfL        = dia%BC%OM_cfL
!   candi%BC%InitMinDepthL = dia%BC%InitMinDepthL
!   candi%BC%OMInitMethodR = dia%BC%OMInitMethodR
!   candi%BC%OM_topR       = dia%BC%OM_topR
!   candi%BC%OM_minR       = dia%BC%OM_minR
!   candi%BC%OM_cfR        = dia%BC%OM_cfR
!   candi%BC%InitMinDepthR = dia%BC%InitMinDepthR
!   candi%BC%InitSedDepth  = dia%BC%InitSedDepth
!   candi%BC%OutputUnits   = dia%BC%OutputUnits

    candi%mk = dia%Xmk   ! Metal X stoichiometry in MnO2
    candi%fl = dia%Xfl   ! Metal X stoichiometry in FeOH
    candi%fm = dia%param%Xfm   ! Metal X stoichiometry in FeS

    !IF(TRIM(candi%OMInitMethod) == 'E') THEN
      candi%BC%OMInitMethodL = "EX_I"
      candi%BC%OMInitMethodR = "LI_I"
    !END IF

    IF(TRIM(candi%BC%OMInitMethodL) == "C") candi%BC%OMInitMethodL = "CO_I"
    IF(TRIM(candi%BC%OMInitMethodR) == "L") candi%BC%OMInitMethodR = "LI_I"

    candi%FilePar     = "results/candi_aed/candi.log"
    candi%FileRat     = "results/candi_aed/rates.sed"
    candi%FileRat2    = "results/candi_aed/rates2.sed"
    candi%FileRatOAc  = "results/candi_aed/OAcrates.sed"
    candi%FileDanPar  = "results/candi_aed/Dan.sed"

    candi%iprint   = 0
    candi%job      = 1
    candi%xldouble = 5.0
    candi%hofmu    = 0
    candi%writestep= 365
    candi%iso4     = 1

    candi%sc = 106
    candi%sn = 16
    candi%sp = 1

    candi%corginp = 0

    IF (candi%param%imix <= 1) THEN
      candi%nld = zero
      candi%lld = one
    ELSE
      candi%nld = one
      candi%lld = zero
    END IF

    IF(candi%param%irrg == 0) THEN
      candi%param%alpha0 = zero
      candi%param%xirrig = zero
    END IF

    candi%npt = candi%maxnpts                         !Deepest sediment layer

    !--   Use explicit=0, use VODE-solver=1 (usevode)
    candi%usevode = 1
    candi%tout    = 0.5

 !   deltat = 1./365.25  !(assumed daily)
 !   deltat = 1./(365.25*24.)  !(assumed hourly)

    candi%deltat = (1./(365.25 * (86400./(candi%time%driverDT *candi%time%substep))))
    print *, 'deltat', candi%deltat

    !--   1=read  setupfile, 0=start with zero initialisation (restartCANDI)
    candi%restartCANDI = .FALSE.

    !--   1=write setupfile after the end of calculation (aufflagout)
    candi%aufflagout = 0
    !-- STEADY STATE CALCULATION? 1=YES 0=NO
    candi%iSTEADY = 0

    candi%param%kpO2    = candi%param%kO2
    candi%param%kpNO3   = candi%param%kNO3
    candi%param%kpMnO2  = candi%param%kMnO2
    candi%param%kpFeOH  = candi%param%kFeOH
    candi%param%kpSO4   = candi%param%kSO4

    ! TYPE OF BOTTOM BOUNDARY CONDITION: 1=NO GRADIENTS 2=KNOWN CONCENTRATION
    candi%BC%ibbc          = dia%BC%ibbc

    !     BOTTOM BOUNDARY CONDITION
    !     NOTE:  IF X = L IS BELOW THE DEPTH OF BIOTURBATION, THE
    !            SOLID CONCENTRATIONS ARE NOT USED AND CAN BE ENTERED
    !            AS ZEROS.
    candi%bottom(:) = one


    candi%BC%POMVR = candi%BC%POMVR * 2650 * 1000 / 1000 * 1000


    candi%param%Xname =  dia%param%Xname


    !-- If user defined deposition rates, then load in
    IF (candi%BC%IBC2==2 .OR. candi%BC%IBC2==3) THEN
!    IF (IBC2==2 .OR. IBC2==3 .OR. IBC2==10) THEN
      PRINT *,'       Prescribed particulate fluxes set from default_vals ... '

      IF (candi%BC%IBC2==3) THEN
        disableWQUpdate = .TRUE.
        candi%BC%IBC2=2
      END IF

      candi%PartFluxes = 0.00

    END IF
    IF (candi%BC%IBC2==10) THEN
      PRINT *,'       Prescribed bottom water condition / flux varying in time '
      PRINT *,'        ...> read from aed2_sediment_swibc.dat '
    END IF

 END SUBROUTINE SetLocalConstants
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Routine to set the initial vertical profiles for OM                          !
!------------------------------------------------------------------------------!
 SUBROUTINE SetSedInitialCondition(candi,mode,Ain,t_th,b_th,minDepth)
   !-- Incoming
    TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
    CHARACTER (LEN=4), INTENT(IN)       :: mode  !-- Type of profile          ?!
                                                 !  CO_I : constant profile    !
                                                 !  LI_I : linear profile      !
                                                 !  EX_I : exponential profile !
    REAL (SEDP), INTENT(IN)             :: t_th  ! Scalar value at the surface !
    REAL (SEDP), INTENT(IN), OPTIONAL   :: b_th  ! Scalar value at the bottom  !

    REAL (SEDP), INTENT(IN), OPTIONAL   :: minDepth
    !-- Outgoing
    REAL (SEDP), DIMENSION(:)           :: Ain   ! Array to initialise         !
    !-- Local
    REAL (SEDP) :: d1
    INTEGER  :: i, min_i
    !-------------------------------------------------------------------!
    !-- Linear interpolation between points
    IF(mode == LI_I) THEN
      min_i = MINLOC( ABS(candi%rpar(:)-minDepth),1 )
      IF(min_i>candi%npt-1)min_i=candi%npt-1
      Ain = t_th
      DO i = 2,min_i
        d1     = candi%rpar(i)/candi%rpar(min_i)
        Ain(i) = t_th - (t_th - b_th)*d1
      END DO
      DO i = min_i+1,candi%npt
        Ain(i) = b_th
      END DO
    !-- Exponential decrease in initial concentration
    ELSE IF(mode == EX_I) THEN
      Ain = t_th
      DO i = 2,candi%npt
        d1     = 0.5 * candi%rpar(i)+candi%rpar(i-1)
        Ain(i) = t_th*EXP(-candi%BC%OM_cfL*candi%rpar(i))
      END DO
    !-- Constant initial value
    ELSE IF(mode == CO_I) THEN
      Ain = t_th
    END IF
 END SUBROUTINE SetSedInitialCondition
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! ****************************************************************
! Early diagenese model C. CANDI from R. Luff (1996-2003), based on
! CANDI from  B.P. Boudreau.
! ****************************************************************
! VERSION 3.0 Roger Luff (rluff@gmx.de) GEOMAR, Kiel, Germany
! ****************************************************************
! The model solves a complete CARBON-OXIDANT-NUTRIENT/BY-PRODUCT Model
! using finite differencing in the space variable, x, and using
! the "stiff"-initial-value code VODE as a method-of-lines integrator.
! This program also calculates pH and dissolved protolytic species,
! i.e. CO2, HCO3-, CO3=, H2S, HS-, etc. using 3 different approaches.
! Program has been tested on different plattforms with different
! compilers, e.g. SUN/SOLARIS, LINUX, DOS (gnu g77), AIX, CRAY.
!
! See CANDI.TXT for a general overview of the model and the
! related publications.
! ****************************************************************
!------------------------------------------------------------------------------!
 SUBROUTINE doCANDI(candi,steadiness)
   !-- Incoming
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   REAL(SEDP), INTENT(INOUT)         :: steadiness
   !-- Local
   REAL(SEDP)   :: startT
   REAL(SEDP)   :: RWORK(candi%LRW)
   INTEGER  :: IWORK(candi%LIW)
   INTEGER  :: ITASK,ISTATE,IOPT,ITOL,IPAR(1)
   INTEGER  :: i, info, j, mf
   INTEGER  :: neq
   REAL(SEDP)   :: endday
   REAL(SEDP)   :: ssum, RPARfex(candi%maxneq)
   REAL(SEDP)   :: told, tcrit,RSW
   REAL(SEDP)   :: v !, dfluxes
   REAL(SEDP)   :: Yfex(candi%maxneq), Ydotfex(candi%maxneq)
   REAL(SEDP)   :: cflux(7)

   !-- Reset
!  Ytemp(:,:) = ZERO
   candi%Ytemp      = ZERO
   startT     = ZERO

   !-- Reset top boundary values
   candi%top_bound(:) = candi%Y(:,1)

   !-- Get molecular diffusion coefficients for solutes
   CALL difcoef2(V,candi%DF,candi%SALT,candi%TEMP,candi%PRES)
   CALL setdifcoef(candi,candi%DF)
   OPEN(15,FILE=TRIM(candi%Filepar),STATUS="UNKNOWN")

   !-- Reformat the 2D Y array into 1D form for the VODE solver
   neq = candi%npt*candi%nSPECIES
   CALL Copy2Dinto1D(neq,yfex,ydotfex,candi%nSPECIES,candi%y,candi%ydot)

   !-- Pre-condition the solver, DVODE
   ITOL = 1
   endday = candi%tout*365.0
   IF (candi%iSTEADY == 1) THEN
     TCRIT = candi%TOUT
     ITASK = 5
   ELSE
     ! endday = 0.0
     TCRIT = 0.0
     ITASK = 4    ! Seems ITASK =4 is critical to get this to do anything !MATT
   END IF

   IWORK = 0
   RWORK = ZERO
   ISTATE   = 1
   IWORK(1) = candi%ML
   IWORK(2) = candi%MU
   IWORK(6) = 3000
   RWORK(1) = TCRIT
   IOPT     = 1
   MF       = 25
   steadiness   = 9999.0
   !-- main solver loop:
   1000   CONTINUE
   IF (candi%iSTEADY == 0) THEN
     TCRIT    = startT + candi%deltaT
     RWORK(1) = TCRIT
   END IF
   RWORK(7) = 0.0
   told = startT
   CALL DVODE(candi, FEX,NEQ,Yfex,startT,TCRIT,ITOL,candi%RTOL,candi%ATOL,ITASK,ISTATE,  &
              IOPT,RWORK,candi%LRW,IWORK,candi%LIW,JEX,MF,RPARfex,IPAR)

   IF (iSTATE /= 2) CALL myerror (candi,istate,0,"          ")
   !#ifdef porosvar
   !!calculate the changes in porosity due to caco3 diss or prec
   !    CALL porosfunc(initflag)
   !#endif
   !-- To prevent the accumulation of small negatives, set small
   !-- negative values to zero
   DO i = 1,neq
     IF (ABS(Yfex(i)) <  1e-10 ) THEN
       candi%lost = candi%lost + Yfex(i)
       Yfex(i) = ZERO
     END IF
   END DO

   candi%day = startT*365.0
   CALL CheckSteadyStatus(candi,steadiness)
   !CALL mass(steadiness,6)
   8080     FORMAT(9(f20.10))
   goto 4000
   3000     continue
   !  #ifdef caco3
   !    #ifdef porosvar
   !      !rl finnish calculation when poros is small (0.02)
   !      if (porosmin  <   0.02) day = endday
   !    #endif
   !  #endif

   IF (candi%day <  endday .AND. steadiness >= 0.02) goto 1000
   4000   continue
   IF(candi%BC%IBC2 == 2 .or. candi%BC%IBC2 == 10) THEN

           ! If FLUX: First layer(1) no calculation here
 !    IF(OMModel ==1) THEN
 !      y(POMLy,1)      = y(POMLy,2)
 !      y(POMRy,1)      = y(POMRy,2)
 !      y(POMspecialy,1) = y(POMspecialy,2)
 !    ELSEIF(OMModel ==2) THEN
 !      y(DOCLy,1) = y(DOCLy,2)
 !      y(POCLy,1) = y(POCLy,2)
 !      IF(simRefOM) THEN
 !        y(DOCRy,1) = y(DOCRy,2)
 !        y(POCRy,1) = y(POCRy,2)
 !      ENDIF
 !    ELSEIF(OMModel ==3) THEN
 !      y(DOMRy,1) = y(DOMRy,2)
 !      y(POMLy,1) = y(POMLy,2)
 !      y(POMRy,1) = y(POMRy,2)
 !    ENDIF End if OMModel = 1, 2 or 3
   END IF ! IBC2 = 2 or 10 or End if endday, steadiness
   IF(candi%iSTEADY /= 0) THEN
     CALL mass(candi,steadiness,6)
   END IF ! End if iSteady
   CALL WriteRates(candi,15)
   CLOSE(15)
   CALL WriteRates2(candi,15)
   CLOSE (15)
   CALL WriteOAcRates(candi,15)
   CLOSE (15)
   RETURN
 END SUBROUTINE doCANDI
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!rl   ****************************************************************
!rl   Early diagenese model C. CANDI from R. Luff (1996-2003), based on
!rl   CANDI from  B.P. Boudreau.
!rl   ****************************************************************
!rl   VERSION 3.0 Roger Luff (rluff@gmx.de) GEOMAR, Kiel, Germany
!rl   ****************************************************************
!rl   FILE: FEX.f
!rl   AIM : Definition of the diagenetic equations to be solved at each
!rl         depth. Spatial derivatives have been replaced by finite
!rl         diffu, and the resulting ODEs are solved by the method-
!rl         of-lines using the integrator VODE.
!rl   ****************************************************************
!------------------------------------------------------------------------------!
 SUBROUTINE FEX(candi_, neq, t, yfex, ydotfex, rparfex, ipar)
   !-- Incoming
   TYPE(aed2_sed_candi_t),TARGET,INTENT(inout) :: candi_
   INTEGER,    INTENT(IN)                  :: neq
   INTEGER,    INTENT(IN)                  :: ipar
   REAL(SEDP), INTENT(IN)                  :: rparfex(*)
   REAL(SEDP), INTENT(IN)                  :: t
   REAL(SEDP), INTENT(INOUT)               :: yfex(neq)
   REAL(SEDP), INTENT(INOUT)               :: ydotfex(neq)
   !-- Local
   REAL(SEDP) :: maxdot
   REAL(SEDP) :: fracfacflux
   REAL(SEDP) :: diffu, adv, irrig
   INTEGER  :: i, j
   INTEGER  :: maxspec,maxlayer
   TYPE(aed2_sed_candi_t),POINTER :: candi

   candi => candi_

   IF(candi%usevode == 1) THEN
     !-- Re-populate 2D Y array from 1D array so we can manipulate
     CALL Copy1Dinto2D(neq,yfex,candi%nSPECIES,candi%y)
   END IF

   !   If candi is in connection with a WC-model, calculate the first
   !   sediment layer for the solutes with the water-concentrations of
   !   the WC-model.

   !#ifdef porosvar
   !psp(:) = ps(:)/poros(:)
   !pps(:) = poros(:)/ps(:)
   !#endif
   !-- Calculate RATE constants
   CALL RATES(candi)

   !-- Calculate REACTION terms
   CALL REACTION(candi)

   !-- Now perform advection, diffusion, irrigation etc.
   DO i = 0,candi%nSPECIES-1

     IF(candi%isSolid(i)) THEN
       CALL solmotion(candi,i,candi%PartFluxes(i))
     ELSE
       CALL liqmotion(candi,i)
     END IF
   END DO
   ! Add the reaction term
   candi%ydot(:,:) = candi%ydot(:,:) + candi%reac(:,:)
   !set top condition
   !ydot(:,1) = 0.0
   !check bottom condition
   candi%ydot(:,candi%npt)= candi%ydot(:,candi%npt)*candi%bottom(:)

   IF (candi%BC%ibbc == 2) THEN
     ! IBBC =2 implies fixed bottom concentration.
     DO i = 0,candi%nSPECIES-1
       !IF(isSolid(i)) THEN
         candi%ydot(i,candi%npt) = 0.0 !FixedBottomConc(i)*0.33
       !END IF !End if isSolid
     END DO
   END IF ! End if ibbc=2
   !#ifdef porosvar
     !!rl   calculate dporos/dt for all species
     !ydot(o2y,:)=ydot(o2y,:)-poros_dot(:)*y(o2y,:)/poros(:)
     !ydot(no3y,:)=ydot(no3y,:)-poros_dot(:)*y(no3y,:)/poros(:)
     !ydot(so4y,:)=ydot(so4y,:)-poros_dot(:)*y(so4y,:)/poros(:)
     !ydot(po4ly,:)=ydot(po4ly,:)-poros_dot(:)*y(po4ly,:)/poros(:)
     !ydot(nh4y,:)=ydot(nh4y,:)-poros_dot(:)*y(nh4y,:)/poros(:)
     !ydot(hsy,:)=ydot(hsy,:)-poros_dot(:)*y(hsy,:)/poros(:)
     !ydot(hco3y,:)=ydot(hco3y,:)-poros_dot(:)*y(hco3y,:)/poros(:)
     !ydot(ch4y,:)=ydot(ch4y,:)-poros_dot(:)*y(ch4y,:)/poros(:)
     !ydot(cor0y,:)=ydot(cor0y,:)+poros_dot(:)*y(cor0y,:)/ps(:)
     !ydot(tbohy,:)=ydot(tbohy,:)-poros_dot(:)*y(tbohy,:)/poros(:)
     !#ifdef mnfe
     !ydot(mno2y,:)=ydot(mno2y,:)+poros_dot(:)*y(mno2y,:)/ps(:)
     !ydot(feohy,:)=ydot(feohy,:)+poros_dot(:)*y(feohy,:)/ps(:)
     !ydot(mniiy,:)=ydot(mniiy,:)-poros_dot(:)*y(mniiy,:)/poros(:)
     !ydot(feiiy,:)=ydot(feiiy,:)-poros_dot(:)*y(feiiy,:)/poros(:)
     !#endif
     !#ifndef kgvar
     !ydot(cor1y,:)=ydot(cor1y,:)+poros_dot(:)*y(cor1y,:)/ps(:)
     !ydot(cor2y,:)=ydot(cor2y,:)+poros_dot(:)*y(cor2y,:)/ps(:)
     !#endif
     !#ifndef eqbb
     !ydot(co2y,:)=ydot(co2y,:)-poros_dot(:)*y(co2y,:)/poros(:)
     !#endif
     !#ifdef eqrm
     !ydot(h2sy,:)=ydot(h2sy,:)-poros_dot(:)*y(h2sy,:)/poros(:)
     !ydot(co3y,:)=ydot(co3y,:)-poros_dot(:)*y(co3y,:)/poros(:)
     !#endif
     !ydot(cay,:)=ydot(cay,:)-poros_dot(:)*y(cay,:)/poros(:)
     !ydot(caly,:)=ydot(caly,:)+poros_dot(:)*y(caly,:)/ps(:)
     !ydot(aray,:)=ydot(aray,:)+poros_dot(:)*y(aray,:)/ps(:)
     !#ifdef po4solid
     !ydot(po4sy,:)=ydot(po4sy,:)+poros_dot(:)*y(po4sy,:)/ps(:)
     !#endif
     !#ifdef fes
     !ydot(fesy,:)=ydot(fesy,:)+poros_dot(:)*y(fesy,:)/ps(:)
     !#endif
     !#ifdef feii
     !ydot(feii1y,:)=ydot(feii1y,:)+poros_dot(:)*y(feii1y,:)/ps(:)
     !ydot(feii2y,:)=ydot(feii2y,:)+poros_dot(:)*y(feii2y,:)/ps(:)
     !#endif
     !#ifdef tracer
     !ydot(pby,:)=ydot(pby,:)+poros_dot(:)*y(pby,:)/ps(:)
     !ydot(ashy,:)=ydot(ashy,:)+poros_dot(:)*y(ashy,:)/ps(:)
     !ydot(cly,:)=ydot(cly,:)-poros_dot(:)*y(cly,:)/poros(:)
     !#endif

     !#ifdef caxco3
     !ydot(caxmny,:)=ydot(caxmny,:)+poros_dot(:)*y(caxmny,:)/ps(:)
     !ydot(caxfey,:)=ydot(caxfey,:)+poros_dot(:)*y(caxfey,:)/ps(:)
     !#endif

     !#ifdef c12
     !ydot(co2c12y,:)=ydot(co2c12y,:)-poros_dot(:) * y(co2c12y,:)/poros(:)
     !ydot(co3c12y,:)=ydot(co3c12y,:)-poros_dot(:) * y(co3c12y,:)/poros(:)
     !ydot(hco3c12y,:)=ydot(hco3c12y,:)-poros_dot(:) * y(hco3c12y,:)/poros(:)
     !ydot(ch4c12y,:)=ydot(ch4c12y,:)-poros_dot(:) * y(ch4c12y,:)/poros(:)
     !ydot(org0c12y,:)=ydot(org0c12y,:)+poros_dot(:) * y(org0c12y,:)/ps(:)
     !ydot(org1c12y,:)=ydot(org1c12y,:)+poros_dot(:) * y(org1c12y,:)/ps(:)
     !ydot(org2c12y,:)=ydot(org2c12y,:)+poros_dot(:) * y(org2c12y,:)/ps(:)
     !#endif

     !#ifdef caco3c12
     !ydot(arac12y,:)=ydot(arac12y,:)+poros_dot(:)*y(arac12y,:)/ps(:)
     !ydot(calc12y,:)=ydot(calc12y,:)+poros_dot(:)*y(calc12y,:)/ps(:)
     !#endif

   !#endif

   !  Find out which species has the largest changes between two timesteps
   maxdot = 0.0
   DO i=1,candi%npt
     DO j=0,candi%nSPECIES-1
       IF (ABS(candi%ytemp(j,i)-candi%ydot(j,i)) > maxdot) THEN
         maxdot   = ABS(candi%ytemp(j,i)-candi%ydot(j,i))
         maxspec  = j
         maxlayer = i
       END IF
       candi%ytemp(j,i) = candi%ydot(j,i)
     END DO
   END DO
   !PRINT *, 'ytemp-ydot:', maxdot, maxspec, maxlayer

   IF(candi%usevode == 1) THEN
     !-- Fill 1D YFEX arrays from 2D arrays for the VODE solver
     CALL Copy2Dinto1D(neq,yfex,ydotfex,candi%nSPECIES,candi%y,candi%ydot)
   END IF

   RETURN
 END SUBROUTINE FEX
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!rl   ****************************************************************
!rl   Early diagenese model C. CANDI from R. Luff (1996-2003), based on
!rl   CANDI from  B.P. Boudreau.
!rl   ****************************************************************
!rl   VERSION 3.0 Roger Luff (rluff@gmx.de) GEOMAR, Kiel, Germany
!rl   ****************************************************************
!rl   FILE: JEX.f
!rl   AIM : A dummy routine needed in calling integrator VODE.
!rl         etc.) from file STEUER.DAT.
!rl   ****************************************************************
!------------------------------------------------------------------------------!
 SUBROUTINE JEX (candi_, neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
   TYPE(aed2_sed_candi_t),INTENT(in) :: candi_
   INTEGER,    INTENT(INOUT)         :: neq
   REAL(SEDP), INTENT(INOUT)         :: t
   REAL(SEDP), INTENT(INOUT)         :: y(neq)
   INTEGER,    INTENT(INOUT)         :: ml
   INTEGER,    INTENT(INOUT)         :: mu
   INTEGER,    INTENT(INOUT)         :: nrpd
   REAL(SEDP), INTENT(INOUT)         :: pd(nrpd,neq)
   REAL(SEDP), INTENT(INOUT)         :: rpar
   INTEGER,    INTENT(INOUT)         :: ipar

   RETURN
 END SUBROUTINE JEX
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!rl   ****************************************************************
!rl   Early diagenese model C. CANDI from R. Luff (1996-2003), based on
!rl   CANDI from  B.P. Boudreau.
!rl   ****************************************************************
!rl   VERSION 3.0 Roger Luff (rluff@gmx.de) GEOMAR, Kiel, Germany
!rl   ****************************************************************
!rl   FILE: MOTION.f
!rl   AIM : A set of routines to calculate the transport within the
!rl         sediment (Advection irrigation and Diffusion) for solute
!rl         and solid species.
!
!rl   calculates the diffusion and advection for solutes at the upper
!rl   boundary on even and uneven grids
!rl   see Boudreau 1997 p. 325 ff
!rl   ****************************************************************
!------------------------------------------------------------------------------!
 SUBROUTINE LiqMotion(candi,speci)
   !-- Incoming
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   INTEGER, INTENT(IN) :: speci
   !-- Local
   INTEGER  :: i
   REAL(SEDP)          :: diff,adv,irrig
   REAL(SEDP)          :: ds,ddsdx
   ds = candi%DIFFC(speci,2)/candi%t2(2)
   ddsdx = -candi%dh(3)*(candi%DIFFC(speci,2)/candi%t2(1))/(candi%dh(2)*(candi%dh(2)+candi%dh(3)))  &
           +(candi%dh(3)-candi%dh(2))*(candi%DIFFC(speci,2)/candi%t2(2))/(candi%dh(2)*candi%dh(3))  &
           +candi%dh(2)*(candi%DIFFC(speci,2)/candi%t2(3))/(candi%dh(3)*(candi%dh(2)+candi%dh(3)))
   diff = ds*(2.0*candi%top_bound(speci)/(candi%dh(2)*(candi%dh(2)+candi%dh(3)))                    &
            -2.0*candi%y(speci,2)/(candi%dh(2)*candi%dh(3))                                         &
            +2.0*candi%y(speci,3)/(candi%dh(3)*(candi%dh(2)+candi%dh(3))))
   adv = (-candi%dh(3)*candi%top_bound(speci)/(candi%dh(2)*(candi%dh(2)+candi%dh(3)))               &
         +(candi%dh(3)-candi%dh(2))*candi%y(speci,2)/(candi%dh(2)*candi%dh(3))                      &
         +candi%dh(2)*candi%y(speci,3)/(candi%dh(3)*(candi%dh(2)+candi%dh(3))))                     &
         *(ddsdx+1.0/candi%poros(2)*ds*candi%dpdx(2)-candi%uvel(2)) -candi%dudx(2)*candi%y(speci,2) &
         -1.0/candi%poros(2)*candi%dpdx(2)*candi%uvel(2)*candi%y(speci,2)
   irrig = candi%cirrig(2)*(candi%top_bound(speci)-candi%y(speci,2))
   candi%ydot(speci,2) = diff + adv + irrig
   DO i = 3,candi%npt-1
     ds = candi%DIFFC(speci,i)/candi%t2(i)
     ddsdx = -candi%dh(i+1)*(candi%DIFFC(speci,i)/candi%t2(i-1))/(candi%dh(i)*(candi%dh(i)+candi%dh(i+1))) &
           +(candi%dh(i+1)-candi%dh(i))*(candi%DIFFC(speci,i)/candi%t2(i))/(candi%dh(i)*candi%dh(i+1))     &
           +candi%dh(i)*(candi%DIFFC(speci,i)/candi%t2(i+1))/(candi%dh(i+1)*(candi%dh(i)+candi%dh(i+1)))
     diff = ds*(2.0*candi%y(speci,i-1)/(candi%dh(i)*(candi%dh(i)+candi%dh(i+1)))                          &
              -2.0*candi%y(speci,i)/(candi%dh(i)*candi%dh(i+1))                                           &
              +2.0*candi%y(speci,i+1)/(candi%dh(i+1)*(candi%dh(i)+candi%dh(i+1))))
     adv = (-candi%dh(i+1)*candi%y(speci,i-1)/(candi%dh(i)*(candi%dh(i)+candi%dh(i+1)))                   &
           +(candi%dh(i+1)-candi%dh(i))*candi%y(speci,i)/(candi%dh(i)*candi%dh(i+1))                      &
           +candi%dh(i)*candi%y(speci,i+1)/(candi%dh(i+1)*(candi%dh(i)+candi%dh(i+1))))                   &
           *(ddsdx+1.0/candi%poros(i)*ds*candi%dpdx(i)-candi%uvel(i)) -candi%dudx(i)*candi%y(speci,i)     &
           -1.0/candi%poros(i)*candi%dpdx(i)*candi%uvel(i)*candi%y(speci,i)
     irrig = candi%cirrig(i)*(candi%top_bound(speci)-candi%y(speci,i))
     candi%ydot(speci,i) = diff + adv + irrig
   ENDDO
   candi%ydot(speci,candi%npt) = candi%DIFFC(speci,candi%npt)/candi%t2(candi%npt) * TWO             &
                     * (-candi%y(speci,candi%npt)+candi%y(speci,candi%npt-1))/candi%dh2(candi%npt)
   RETURN
 END SUBROUTINE LiqMotion
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!rl   calculates the diffusion and advection for solids on even and uneven grids
!rl   see Boudreau 1997 p. 325 ff
!------------------------------------------------------------------------------!
 SUBROUTINE SolMotion(candi,speci,flux)
   !-- Incoming
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   INTEGER,    INTENT(IN)   :: speci
   REAL(SEDP), INTENT(IN)   :: flux
   !-- Local
   INTEGER  :: i
   REAL(SEDP)   :: dcdx,psdiff,dbdiff,diff,adv,dumm,psadv


   IF (flux > zero) THEN
     ! Eq 103 in Boudreau 1996
     dumm = candi%y(speci,3) +candi%tdh*(flux-candi%ps(2)*candi%wvel(2)*candi%y(speci,2))/candi%bioturb(2)/candi%ps(2)
   ELSE
     dumm = candi%y(speci,1)
   END IF


 !IF(speci==9) print *,'Dumm',speci,ps(2)*wvel(2)*y(speci,2)


   dcdx = (-candi%dh(3)*dumm/(candi%dh(2)*(candi%dh(2)+candi%dh(3)))                                   &
         +(candi%dh(3)-candi%dh(2))*candi%y(speci,2)/(candi%dh(2)*candi%dh(3))                               &
         +candi%dh(2)*candi%y(speci,3)/(candi%dh(3)*(candi%dh(2)+candi%dh(3))))

   !diff = bioturb(2)*(y(speci,3)-two*y(speci,2)+dumm)/dh2(2)
   diff = candi%bioturb(2)*(2.0*dumm/(candi%dh(2)*(candi%dh(2)+candi%dh(3)))                           &
            -2.0*candi%y(speci,2)/(candi%dh(2)*candi%dh(3))                                      &
            +2.0*candi%y(speci,3)/(candi%dh(3)*(candi%dh(2)+candi%dh(3))))

   adv  = -candi%wvel(2)*dcdx
   psadv= candi%dpdx(2)*candi%wvel(2)*candi%y(speci,2)/candi%ps(2)
   psdiff  = candi%bioturb(2)/candi%ps(2)*dcdx*(-candi%dpdx(2))
   dbdiff  = dcdx*candi%dbdx(2)
   candi%ydot(speci,2) = diff+adv+psadv+psdiff+dbdiff

   DO i = 3,candi%npt-1
     IF(candi%bioturb(i) > 1.0E-05) THEN
       dcdx = (-candi%dh(i+1)*candi%y(speci,i-1)/(candi%dh(i)*(candi%dh(i)+candi%dh(i+1)))                   &
             +(candi%dh(i+1)-candi%dh(i))*candi%y(speci,i)/(candi%dh(i)*candi%dh(i+1))                       &
             +candi%dh(i)*candi%y(speci,i+1)/(candi%dh(i+1)*(candi%dh(i)+candi%dh(i+1))))
     ELSE
       dcdx = (candi%y(speci,i)-candi%y(speci,i-1))/candi%dh(i-1)
     END IF

     psdiff  = candi%bioturb(i)/candi%ps(i)*dcdx*(-candi%dpdx(i))
     dbdiff  = dcdx*candi%dbdx(i)
     diff = candi%bioturb(i)*(2.0*candi%y(speci,i-1)/(candi%dh(i)*(candi%dh(i)+candi%dh(i+1)))               &
                       -2.0*candi%y(speci,i)/(candi%dh(i)*candi%dh(i+1))                         &
                       +2.0*candi%y(speci,i+1)/(candi%dh(i+1)*(candi%dh(i)+candi%dh(i+1))))
     adv  = - candi%wvel(i)*dcdx
     psadv= candi%dpdx(i)*candi%wvel(i)*candi%y(speci,i)/candi%ps(i)
     !rl   adv  =  Wvel(i)*(y(speci,i-1)-y(speci,i))/dhr(i)
     candi%ydot(speci,i) = diff+adv+psadv+psdiff+dbdiff
   ENDDO

   diff = (-candi%y(speci,candi%npt)+candi%y(speci,candi%npt-1)) *two*candi%bioturb(candi%npt)/candi%dh2(candi%npt)
   adv  = candi%wvel(candi%npt)*(candi%y(speci,candi%npt-1)-candi%y(speci,candi%npt))/candi%dhr(candi%npt)

   candi%ydot(speci,candi%npt) = diff+adv

   RETURN
 END SUBROUTINE SolMotion
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!rl   calculates the irrigation, for flux calculations
!------------------------------------------------------------------------------!
 REAL(SEDP) FUNCTION FIrrig(candi,speci,layer)
   !-- Incoming
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   INTEGER, INTENT(IN)                  :: speci
   INTEGER, INTENT(IN)                  :: layer

   firrig = candi%cirrig(layer)*(candi%top_bound(speci) - candi%y(speci,layer))

   RETURN
 END FUNCTION FIrrig
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!rl   ****************************************************************
!rl   Early diagenese model C. CANDI from R. Luff (1996-2003), based on
!rl   CANDI from  B.P. Boudreau.
!rl   ****************************************************************
!rl   VERSION 3.0 Roger Luff (rluff@gmx.de) GEOMAR, Kiel, Germany
!rl   ****************************************************************
!rl   FILE: REACTION.f
!rl   AIM : Definition/calculation of the reactions for all species in
!rl   the model.
!dp   (Subroutine RATES is below)
!dp   CONTENTS:
!dp   * C:N:P stoichiometry function
!dp   * Three organic matter model options:
!dp     - POM --> CO2
!dp     - POM --> DOM --> CO2
!dp     - POM --> Hydrolysis Products --> Fermentation Products
!dp                               \      /        |
!dp                                V    V         V
!dp                                  CO2    <--  DOMR
!dp   * Organic X
!     * REACTION TERMS:
!     To aid the user in both understanding and modifying this
!     code, here is a list of the variables used to represent
!     the reaction terms included in the model and the equations
!     used below: The numbers represent the formula numbers from Boudreau (1996)

!  (01)   RO2     = oxygen dependence of oxic organic matter decay rate
!  (02)   RNO3    = nitrate dependence of nitrate reduction
!  (03)   RMNO2   = manganese oxide dependence of MnO2 reduction
!  (04)?  RFEOH   = iron oxide dependence of Fe(OH)3 reduction
!  (05)   RSO4    = sulfate dependence of sulfate reduction
!  (06)   RFESPPT = rate of FeS precipitation
!  (07)   RPO4PPT = rate of PO4 precipitation
!  (08)   RFEOX   = rate of Fe(II) oxidation by O2
!  (09)   RMNOX   = rate of Mn(II) oxidation by O2
!  (10)   RSOX    = rate of TH2S oxidation by O2
!  (11)   RNHOX   = rate of TNH4 oxidation by O2
!  (12)   RCH4OX  = rate of CH4 oxidation by reaction with O2
!  (13)   RCH4SO4 = rate of CH4 oxidation by reaction with SO4
!  (14)   RMNFE   = rate of Fe(II) oxidation by reaction with MnO2
!  (15)   RMNO2TS = rate of TH2S oxidation by reaction with MnO2
!  (16)   RTSFE3  = rate of TH2S oxidation by reaction with Fe(OH)3
!  (17)   RFESFE3 = rate of FeS oxidation by reaction with Fe(OH)3
!  (18)   RFESMN4 = rate of FeS oxidation by reaction with MnO2
!  (19)   RFESOX  = rate of FeS oxidation by reaction with O2
!         RCH4    = sulfate inhibition on methanogenesis
!         RNH     = RNO3 + RMNO2 + RFEOH + RSO4 + RCH4
!         ROX     = RO2 + RNH
!         RG      = rate of organic matter decay without inhibition
!         RCALDIS = rate of calcite dissolution
!         RCALPPT = rate of calcite precipitation
!         RARADIS = rate of aragonite dissolution
!         RARAPPT = rate of aragonite precipitation
!         KPB     = rate constant of PB-210 decay
!mh  ???  RNO3TS  = rate of TH2S oxidation by reaction with NO3
!mh       RFe1OX  = rate of Fe(II)-oxidation in clays by O2, reactive fraction
!mh       RFe2OX  = rate of Fe(II)-oxidation in clays by O2, non-reactive fraction
!mh       RFe1NO3 = rate of Fe(II)-oxidation in clays by NO3, reactive fraction
!mh       RFe2NO3 = rate of Fe(II)-oxidation in clays by NO3, non-reactive fraction
!mh       RMnNO3  = rate of Mn(II) oxidation by NO3
!rl   ****************************************************************
!------------------------------------------------------------------------------!
 SUBROUTINE REACTION(candi)
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   !-- Local
   INTEGER  :: i

   candi%reac(:,:) = ZERO
  ! reac(O2y,2:npt) = -0.9
  !      RETURN
          ! SUBROUTINE REACTION
   DO i = 2,candi%npt ! Search "Bunyip"

     !--------------------------------------
     ! Estimate stoichiometry based on available OM.
     IF(candi%Y(candi%DOPLy,i) >1e-8)THEN
       candi%sc = candi%Y(candi%DOCLy,i)/candi%Y(candi%DOPLy,i)
       candi%sn = candi%Y(candi%DONLy,i)/candi%Y(candi%DOPLy,i)
       candi%sp = ONE
     ELSE
       candi%sc = 106.0
       candi%sn = 16.0
       candi%sp = ONE
     END IF
                                                        ! SUBROUTINE REACTION
     !Dan starts here
     !--------------------------------------
     IF ( candi%param%OMModel == 1 ) THEN !Simple OM --> CO2              !OM Model 1
       candi%reac(candi%POMLy,i) = - candi%param%poml2dic*candi%y(candi%POMLy,i) !*rox(i)
       candi%reac(candi%POMRy,i) = - candi%param%pomr2dic*candi%y(candi%POMRy,i) !*rox(i)
       candi%reac(candi%POMspecialy,i) = - candi%param%pomspecial2dic*candi%y(candi%POMspecialy,i) !*rox(i)
! 20141127 I am removing rox from the reaction just to test it.
        !Result: yes, rox should be turned off here. Boudreau 1996 p 481 reaction 29 is wrong. - Dan
     candi%reac(candi%POCLy,i) = 0.0
     candi%reac(candi%DOCLy,i) = 0.0
     candi%reac(candi%PONLy,i) = 0.0
     candi%reac(candi%DONLy,i) = 0.0
     candi%reac(candi%POPLy,i) = 0.0
     candi%reac(candi%DOPLy,i) = 0.0
     candi%reac(candi%POCRy,i) = 0.0
     candi%reac(candi%DOCRy,i) = 0.0
     candi%reac(candi%PONRy,i) = 0.0
     candi%reac(candi%DONRy,i) = 0.0
     candi%reac(candi%POPRy,i) = 0.0
     candi%reac(candi%DOPRy,i) = 0.0
     candi%reac(candi%dhydy,i) = 0.0

     candi%reac(candi%POMLy,i) = 0.0
     candi%reac(candi%POMRy,i) = 0.0
     candi%reac(candi%POMspecialy,i) = 0.0
     candi%reac(candi%POMRy,i) = 0.0
     candi%reac(candi%DOMRy,i) = 0.0
     candi%reac(candi%POM1y,i) =  0.0
     candi%reac(candi%POM2y,i) =  0.0
     candi%reac(candi%POM3y,i) =  0.0
     candi%reac(candi%POM4y,i) =  0.0
     candi%reac(candi%OAcy,i)  =  0.0
     candi%reac(candi%H2y,i)   =  0.0
     candi%reac(candi%BFery,i) =  0.0
     candi%reac(candi%BAery,i) =  0.0
     candi%reac(candi%BDeny,i) =  0.0
     candi%reac(candi%BMany,i) =  0.0
     candi%reac(candi%BIroy,i) =  0.0
     candi%reac(candi%BSuly,i) =  0.0
     candi%reac(candi%BMety,i) =  0.0
     candi%reac(candi%Necromassy,i) =  0.0
!       ELSEIF ( OMModel == 2 ) THEN !Laurie Era four pools                   ! OM Model 2
     ELSEIF ( candi%param%OMModel == 2 ) THEN !Laurie Era four pools
     ! Organic Group #1
     candi%reac(candi%POCLy,i) = - candi%param%pocl2docl*candi%y(candi%POCLy,i)*candi%rox(i)
     candi%reac(candi%DOCLy,i) = - candi%param%docl2dic*candi%y(candi%DOCLy,i)*candi%rox(i)             &
                               + candi%param%pocl2docl*candi%y(candi%POCLy,i)*candi%rox(i)            &
                               + candi%param%docr2docl*candi%y(candi%DOCRy,i)*candi%rox(i)
     candi%reac(candi%PONLy,i) = - candi%param%ponl2donl*candi%y(candi%PONLy,i)*candi%rox(i)
     candi%reac(candi%DONLy,i) = - candi%param%donl2din*candi%y(candi%DONLy,i)*candi%rox(i)             &
                               + candi%param%ponl2donl*candi%y(candi%PONLy,i)*candi%rox(i)            &
                               + candi%param%donr2donl*candi%y(candi%DONRy,i)*candi%rox(i)
     candi%reac(candi%POPLy,i) = - candi%param%popl2dopl*candi%y(candi%POPLy,i)*candi%rox(i)
     candi%reac(candi%DOPLy,i) = - candi%param%dopl2dip*candi%y(candi%DOPLy,i)*candi%rox(i)             &
                               + candi%param%popl2dopl*candi%y(candi%POPLy,i)*candi%rox(i)            &
                               + candi%param%dopr2dopl*candi%y(candi%DOPRy,i)*candi%rox(i)
     ! Organic Group #2
     candi%reac(candi%POCRy,i) = - candi%param%pocr2docr*candi%y(candi%POCRy,i)*candi%rox(i)
     candi%reac(candi%DOCRy,i) = - candi%param%docr2docl*candi%y(candi%DOCRy,i)*candi%rox(i)            &
                               + candi%param%pocr2docr*candi%y(candi%POCRy,i)*candi%rox(i)            &
                               + candi%param%pocvr2docr*(candi%BC%POMVR*0.5/12.0)*candi%rox(i)
     candi%reac(candi%PONRy,i) = - candi%param%ponr2donr*candi%y(candi%PONRy,i)*candi%rox(i)
     candi%reac(candi%DONRy,i) = - candi%param%donr2donl*candi%y(candi%DONRy,i)*candi%rox(i)            &
                               + candi%param%ponr2donr*candi%y(candi%PONRy,i)*candi%rox(i)            &
                               + candi%param%pocvr2docr*(candi%BC%POMVR*0.075/14.0)*candi%rox(i)
     candi%reac(candi%POPRy,i) = - candi%param%popr2dopr*candi%y(candi%POPRy,i)*candi%rox(i)
     candi%reac(candi%DOPRy,i) = - candi%param%dopr2dopl*candi%y(candi%DOPRy,i)*candi%rox(i)            &
                               + candi%param%popr2dopr*candi%y(candi%POPRy,i)*candi%rox(i)            &
                               + candi%param%pocvr2docr*(candi%BC%POMVR*0.005/30.9)*candi%rox(i)

     candi%reac(candi%POMLy,i) = 0.0
     candi%reac(candi%POMRy,i) = 0.0
     candi%reac(candi%POMspecialy,i) = 0.0
     candi%reac(candi%dhydy,i) = 0.0
     candi%reac(candi%POMRy,i) = 0.0
     candi%reac(candi%DOMRy,i) = 0.0
     candi%reac(candi%POM1y,i) =  0.0
     candi%reac(candi%POM2y,i) =  0.0
     candi%reac(candi%POM3y,i) =  0.0
     candi%reac(candi%POM4y,i) =  0.0
     candi%reac(candi%OAcy,i)  =  0.0
     candi%reac(candi%H2y,i)   =  0.0
     candi%reac(candi%BFery,i) =  0.0
     candi%reac(candi%BAery,i) =  0.0
     candi%reac(candi%BDeny,i) =  0.0
     candi%reac(candi%BMany,i) =  0.0
     candi%reac(candi%BIroy,i) =  0.0
     candi%reac(candi%BSuly,i) =  0.0
     candi%reac(candi%BMety,i) =  0.0
     candi%reac(candi%Necromassy,i) =  0.0

     ELSEIF ( candi%param%OMModel == 3 ) THEN                            ! OM Model 3
     candi%reac(candi%POM1y,i)       = - candi%RPOM1(i)
     candi%reac(candi%POM2y,i)       = - candi%RPOM2(i)
     candi%reac(candi%POM3y,i)       = - candi%RPOM3(i)
     candi%reac(candi%POM4y,i)       = - candi%RPOM4(i)
     candi%reac(candi%POMspecialy,i) = - candi%RPOMspecial(i)
     candi%reac(candi%dhydy,i) = +candi%RPOM1(i)+candi%RPOM2(i)+candi%RPOM3(i)+candi%RPOM4(i)+candi%RNecro(i)&  ! *rox(i)
                     -candi%RAerDHyd(i)-candi%RDenO2DHyd(i)-candi%RDenNO3DHyd(i)-candi%RFerDHyd(i)
     candi%reac(candi%OAcy,i)  = + candi%RFerDHyd(i)*candi%stcoef%fracOAc - candi%ROAc(i)
     candi%reac(candi%H2y,i)   = + candi%RFerDHyd(i)*candi%stcoef%fracH2  - candi%RH2(i)
     candi%reac(candi%BFery,i)     =candi%RFerDHyd(i)*candi%param%YDHyFer/candi%param%CellWeight -candi%RdeathFer(i)
     candi%reac(candi%BAery,i)     =candi%RAerDHyd(i)*candi%param%YDHyAer/candi%param%CellWeight+candi%RAerOAc(i)*candi%param%YAerOAc/candi%param%CellWeight-candi%RdeathAer(i)

     !print *,"RAerDHyd",RAerDHyd(:)
     !print *,"YDHyAer",YDHyAer
     !print *,"CellWeight",CellWeight
     !print *,"RAerOAc",RAerOAc(:)
     !print *,"YAerOAc",YAerOAc
     !print *,"RdeathAer",RdeathAer
     candi%reac(candi%BDeny,i)     = candi%RDenO2DHyd(i)*candi%param%YDenDHy/candi%param%CellWeight+candi%RDenNO3DHyd(i)*candi%param%YDenDHy/candi%param%CellWeight  &
                   + candi%RDenOAc(i)*candi%param%YDenOAc/candi%param%CellWeight+candi%RDenH2(i)*candi%param%YDenH2/candi%param%CellWeight       -candi%RdeathDen(i)
     candi%reac(candi%BMany,i)     = candi%RManOAc(i)*candi%param%YManOAc/candi%param%CellWeight                               -candi%RdeathMan(i)
     candi%reac(candi%BIroy,i)     = candi%RIroOAc(i)*candi%param%YIroOAc/candi%param%CellWeight+candi%RIroH2(i)*candi%param%YIroH2/candi%param%CellWeight   -candi%RdeathIro(i)
     candi%reac(candi%BSuly,i)     = candi%RSulOAc(i)*candi%param%YSulOAc/candi%param%CellWeight+candi%RSulH2(i)*candi%param%YSulH2/candi%param%CellWeight   -candi%RdeathSul(i)
     candi%reac(candi%BMety,i)     = candi%RMetOAc(i)*candi%param%YMetOAc/candi%param%CellWeight+candi%RMetH2(i)*candi%param%YMetH2/candi%param%CellWeight   -candi%RdeathMet(i)
     candi%reac(candi%Necromassy,i)= candi%RdeathTot(i)*candi%param%fuse-candi%RNecro(i)
     candi%reac(candi%protony,i)   = 0. !1. !y(protony,i)
                                        ! SUBROUTINE REACTION
     candi%reac(candi%POCLy,i) = 0.0
     candi%reac(candi%DOCLy,i) = 0.0
     candi%reac(candi%PONLy,i) = 0.0
     candi%reac(candi%DONLy,i) = 0.0
     candi%reac(candi%POPLy,i) = 0.0
     candi%reac(candi%DOPLy,i) = 0.0
     candi%reac(candi%POCRy,i) = 0.0
     candi%reac(candi%DOCRy,i) = 0.0
     candi%reac(candi%PONRy,i) = 0.0
     candi%reac(candi%DONRy,i) = 0.0
     candi%reac(candi%POPRy,i) = 0.0
     candi%reac(candi%DOPRy,i) = 0.0
     !candi%reac(candi%POMspecialy,i) = 0.0

END IF !End if OMmodel = 3
!     !--------------------------------------

!End Dan
     ! Organic X
     IF(candi%param%simX) THEN
      candi%reac(candi%POXy,i)    = - candi%param%pocl2docl*candi%y(candi%POXy,i)*candi%rox(i)
      candi%reac(candi%DOXy,i)    = - candi%param%docl2dic*candi%y(candi%DOXy,i)*candi%rox(i) + candi%param%pocl2docl*candi%y(candi%POXy,i)*candi%rox(i)
      candi%reac(candi%Xy,i)      =   candi%param%docl2dic*candi%y(candi%DOXy,i)*candi%rox(i)
     END IF
     !--------------------------------------
     ! O2
     candi%reac(candi%O2y,i) = - candi%rgC(i)*candi%FO2(i)        &
                             - TWO*candi%RNH4OX(i)             &
                             - HALF*candi%RMnOX(i)             &
                             - 0.25*candi%RFeOX(i)             &
                             - TWO*candi%RTSOX(i)              &
                             - candi%RCH4OX(i)                 &
                             - TWO*candi%psp(i)*candi%RFeSOX(i) &
                             - 3.5*candi%psp(i)*candi%RFeS2OX(i)
     !--------------------------------------
     ! NO3
     candi%reac(candi%NO3y,i) = - 0.8*candi%rgC(i)*candi%FNO3(i) & !/sc*(-(REAL(0.8,SEDP)*sc+REAL(0.6,SEDP)*sn)*rno3(i)+sn*ro2(i)) &
                              + candi%rnh4ox(i)           &
                              - candi%rnh4no2(i)          &
                              - TWO*candi%rmnno3(i)       &
                              - candi%rfeno3(i)           &
                              - candi%rtsno3(i)
     !--------------------------------------
     !N2
     candi%reac(candi%N2y,i) = + 0.4*candi%rgC(i)*candi%FNO3(i)

     !--------------------------------------
     ! SO4
     candi%reac(candi%SO4y,i) = - HALF*candi%rgC(i)*candi%FSO4(i)     &
                              + candi%rtsox(i)                     &
                              + candi%psp(i)*candi%rfesox(i)        &
                              + candi%psp(i)*TWO*candi%rfes2ox(i)   &
                              + candi%psp(i)*candi%rtsmnA(i)        &
                              + candi%psp(i)*candi%rtsmnB(i)        &
                              + candi%psp(i)*candi%rtsfeA(i)        &
                              + candi%psp(i)*candi%rtsfeB(i)        &
                              - candi%rch4so4(i)                   &
                              + 2.5*candi%rtsno3(i)                &
                            ! + candi%rfesno3(i)                   &
                              + candi%rfesmnA(i) * candi%psp(i)     &
                              + candi%rfesfeA(i) * candi%psp(i)     &
                              + candi%rfesmnB(i) * candi%psp(i)     &
                              + candi%rfesfeB(i) * candi%psp(i)
                                                          ! SUBROUTINE REACTION

     !--------------------------------------
     ! TPO4
     candi%reac(candi%PO4ly,i) = + candi%rgP(i) & !*candi%ROX(i) & I'm turning these off because I think Boudreau is wrong - Dan 20141127
                     !- candi%psp(i)*candi%RPO4ads(i)     &
                     - candi%psp(i)*candi%RPO4ads(i)

     candi%reac(candi%PO4sy,i) = candi%RPO4ads(i)  !       &
                     !- candi%psp(i)*candi%po4sdis(i)

     !--------------------------------------
     ! TNH4
     candi%reac(candi%NH4y,i) = + candi%rgN(i) & !*rox(i)  &I'm turning these off because I think Boudreau is wrong - Dan 20141127
                    - candi%rnh4ox(i)                    &
                    - candi%rnh4no2(i)                   &
                    !- candi%psp(i)*candi%nh4sppt(i)      &
                    - candi%psp(i)*candi%RNH4ads(i)

     candi%reac(candi%NH4sy,i)= candi%RNH4ads(i)           !&
                    !- candi%psp(i)*candi%nh4sdis(i)

     !--------------------------------------
     ! CH4
     candi%reac(candi%CH4y,i) = + candi%rgC(i)*HALF*candi%FMet(i)  &
                    - candi%RCH4OX(i)                           &
                    - candi%RCH4SO4(i)

     !--------------------------------------
     ! TH2S
     candi%reac(candi%HSy,i)  = + candi%rgC(i)*candi%FSO4(i)*HALF  &
                              - candi%RTSOX(i)                  &
                              - candi%psp(i)*candi%RTSMnA(i)     &
                              - candi%psp(i)*candi%RTSMnB(i)     &
                              - candi%psp(i)*candi%RTSFeA(i)     &
                              - candi%psp(i)*candi%RTSFeB(i)     &
                              + candi%RCH4SO4(i)                &
                              - candi%RFeSppt(i)                &
                              - candi%RPyrite(i)                &
                              - 2.5*candi%RTSNO3(i)             &
                              - candi%RXSppt(i)

                    ! SUBROUTINE REACTION
     !--------------------------------------
     ! TCO2
     candi%reac(candi%HCO3y,i) = + candi%rgC(i)*(candi%FO2(i)+candi%FNO3(i)+candi%FMnO2(i)+candi%FFeOH(i)+candi%FSO4(i)) &
                     + candi%rgC(i)*HALF*candi%FMet(i)    &
                     + candi%RCH4OX(i)                   &
                     + candi%RCH4SO4(i)                  &
                     - candi%psp(i)*candi%RSidppt(i)      &
                     - candi%psp(i)*candi%RCalppt(i)     !&
                     !- candi%psp(i)*candi%rarappt(i)     &
                     !+ candi%psp(i)*candi%rfeco3dis(i)   &
                     !+ candi%psp(i)*candi%rcaldis(i)     &
                     !+ candi%psp(i)*candi%raradis(i)
                             ! SUBROUTINE REACTION
     !--------------------------------------
     ! TX
     IF(candi%param%simX) THEN
      candi%reac(candi%Xy,i) = + candi%rgX(i)*candi%ROX(i)                             &
                             + candi%psp(i)*TWO*candi%rgC(i)*candi%mk*candi%FMnO2(i)   &
                             + candi%psp(i)*FOUR*candi%rgC(i)*candi%fl*candi%FFeOH(i)  &
                             - candi%mk*candi%RMnOX(i)                               &
                      !      - candi%fl*candi%RFeOHppt(i)                            &
                             + candi%fm*candi%RFeSOX(i)                              &
                             + candi%fm*candi%RFeS2OX(i)                             &
                             - candi%psp(i)*(TWO*candi%fl-candi%mk)*candi%RFeMnA(i)    &
                             - candi%psp(i)*(TWO*candi%fl-candi%mk)*candi%RFeMnB(i)    &
                             + candi%mk*candi%RMnNO3(i)                              &
                             + FOUR*candi%mk*candi%RTSMnA(i)                         &
                             + FOUR*candi%mk*candi%RTSMnB(i)                         &
                             + candi%psp(i)*(candi%fm+FOUR*candi%mk)*candi%RFeSMnA(i)  &
                             + candi%psp(i)*(candi%fm+FOUR*candi%mk)*candi%RFeSMnB(i)  &
                             + candi%psp(i)*EIGHT*candi%fl*candi%RTSFeA(i)            &
                             + candi%psp(i)*EIGHT*candi%fl*candi%RTSFeB(i)            &
                             + candi%psp(i)*(candi%fm+EIGHT*candi%fl)*candi%RFeSFeA(i) &
                             + candi%psp(i)*(candi%fm+EIGHT*candi%fl)*candi%RFeSFeB(i) &
                             - candi%fm*candi%RFeSppt(i)                             &
                             - candi%RXSppt(i)
     END IF
                                                          ! SUBROUTINE REACTION
     !--------------------------------------
     ! Ca
     candi%reac(candi%Cay,i)   = - candi%RCalppt(i)


     !--------------------------------------
     ! MnO2, FeOH3, MnII & FeII
     IF(candi%param%simMnFe) THEN

       !--------------------------------------
       ! MnO2(A)
       candi%reac(candi%MnO2y,i) = - TWO*candi%rgC(i)*candi%FMnO2(i)       &
                                 + candi%pps(i)*candi%RMnOX(i)           &
                                 + candi%pps(i)*FIVE*candi%RMnNO3(i)     &
                                 - candi%RFeMnA(i)                      &
                                 - FOUR*candi%RTSMnA(i)                 &
                                 - FOUR*candi%RFeSMnA(i)                &
                                 - candi%RMnAge(i)                      &
                                 + candi%psp(i)*candi%RMnO2Appt(i)

       !--------------------------------------
       ! MnO2(B)
       candi%reac(candi%MnO2By,i)= - candi%RFeMnB(i)                 &
                                 - FOUR*candi%RTsMnB(i)            &
                                 - FOUR*candi%RFeSMnB(i)           &
                                 + candi%RMnAge(i)                 &
                                 + candi%psp(i)*candi%RMnO2Bppt(i)
       !--------------------------------------
       ! FeOH3(A)
       candi%reac(candi%FeOHy,i) = - FOUR*candi%rgC(i)*candi%FFeOH(i) &
                              !  + candi%pps(i)*FOUR*candi%RFeOX(i) & ! Moved to Fe3+
                                 + candi%pps(i)*candi%RFeOHAppt(i)  &
                                 + candi%RFeMnA(i)                 &
                                 + candi%RFeMnB(i)                 &
                                 - EIGHT*candi%RTSFeA(i)           &
                                 - EIGHT*candi%RFeSFeA(i)          &
                                 - candi%RFeAge(i)

       !--------------------------------------
       ! FeOH3(B)
       candi%reac(candi%FeOHBy,i)= - EIGHT*candi%RTSFeB(i)           &
                       - EIGHT*candi%RFeSFeB(i)          &
                   !    + candi%pps(i)*candi%RFeOHAppt(i)       &      ! Check that it should be RFeOHA, not B
                       + candi%RFeAge(i)
                                       ! SUBROUTINE REACTION
       !--------------------------------------
       ! MnII
       candi%reac(candi%MnIIy,i) = + TWO*candi%rgC(i)*candi%FMnO2(i)       &
                                 - candi%RMnOX(i)                       &
                                 + candi%psp(i)*candi%RFeMnA(i)          &
                                 + candi%psp(i)*candi%RFeMnB(i)          &
                                 + candi%psp(i)*FOUR*candi%RTSMnA(i)     &
                                 + candi%psp(i)*FOUR*candi%RTSMnB(i)     &
                                 + candi%psp(i)*FOUR*candi%RFeSMnA(i)    &
                                 + candi%psp(i)*FOUR*candi%RFeSMnB(i)    &
                                 - FIVE*candi%RMnNO3(i)                 &
                                 - candi%psp(i)*candi%RMnO2Appt(i)       &
                                 - candi%psp(i)*candi%RMnO2Bppt(i)
       !--------------------------------------
       ! FeII
       candi%reac(candi%FeIIy,i) = + FOUR *candi%rgC(i)*candi%FFeOH(i)      &
                       - FOUR*candi%RFeOX(i)                  &
                       + candi%psp(i)*candi%RFeSOX(i)          &
                       + candi%psp(i)*candi%RFeS2OX(i)         &
                       - candi%psp(i)*TWO*candi%RFeMnA(i)      &
                       - candi%psp(i)*TWO*candi%RFeMnB(i)      &
                       + candi%psp(i)*EIGHT*candi%RTSFeA(i)    &
                       + candi%psp(i)*EIGHT*candi%RTSFeB(i)    &
                       - candi%psp(i)*candi%RFeSppt(i)         &
                       - candi%psp(i)*candi%RSidppt(i)         &
                       + FIVE*candi%RFeNO3(i)                 &
                       + candi%psp(i)*candi%RFeSMnA(i)         &
                       + candi%psp(i)*NINE*candi%RFeSFeA(i)    &
                       + candi%psp(i)*candi%RFeSMnB(i)         &
                       + candi%psp(i)*NINE*candi%RFeSFeB(i)

       !--------------------------------------
       ! FeIII
       candi%reac(candi%FeIIIy,i) = + FOUR*candi%RFeOX(i)
                                                                 ! SUBROUTINE REACTION
       !--------------------------------------
       ! FeS & FeS2
       IF(candi%param%simFeS) THEN
         candi%reac(candi%FeSy,i) = + candi%pps(i)*candi%RFeSppt(i)   &
                                  - candi%RFeSFeA(i)               &
                                  - candi%RFeSMnA(i)               &
                                  - candi%RFeSFeB(i)               &
                                  - candi%RFeSMnB(i)               &
                                  - candi%RFeSOx(i)                &
                                  - candi%RPyrite(i)

         candi%reac(candi%FeS2y,i)= + candi%RPyrite(i)               &
                        - candi%RFeS2OX(i)

         IF(candi%simXS) candi%reac(candi%XSy,i)  = + candi%pps(i)*candi%RXSppt(i)

       END IF


       !--------------------------------------
       ! FeCO3 (Siderite)
       IF(candi%param%simFeCO3) THEN
         candi%reac(candi%Sidy,i) = + candi%pps(i)*candi%RSidppt(i)

       END IF

       !--------------------------------------
       ! MnCO3 (Rhodcrosite)
       IF(candi%param%simMnCO3) THEN
         candi%reac(candi%rody,i) = + candi%pps(i)*candi%RRodppt(i)

       END IF

     END IF
                                                          ! SUBROUTINE REACTION
     !--------------------------------------
     ! CaCO3 (Calcite)
     IF(candi%param%simCaCO3) THEN
       candi%reac(candi%caly,i) =  candi%pps(i)*candi%RCalppt(i) !      &
                      !+ candi%psp(i)*candi%rcalppt(i)
     END IF

   END DO ! Search "Bunyip"
!print*, 'OK so far'
!            pause

   RETURN
 END SUBROUTINE REACTION
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!rl   ****************************************************************
!rl   Early diagenese model C. CANDI from R. Luff (1996-2003), based on
!rl   CANDI from  B.P. Boudreau.
!rl   ****************************************************************
!rl   VERSION 3.0 Roger Luff (rluff@gmx.de) GEOMAR, Kiel, Germany
!rl   ****************************************************************
!rl   FILE: RATES.f
!rl   AIM : Calculation of the factors and rates for the diagenetic equations of
!rl         all reactions.
!rl   ****************************************************************
!dp   (Subroutine REACTION is above)
!dp   CONTENTS:
!dp   * FBIO
!dp   * FTem
!dp   * FTEA - Approach 1 / Approach 2
!dp   * FIN - Approach 1 / Approach 2
!dp   * Rates for each of the reactions in the Reactions section above
!dp   * Secondary rates, including options to turn on or off Mn, Fe, FeS processes
!dp   * rg C:N:P, rox
!------------------------------------------------------------------------------!
 SUBROUTINE RATES(candi)
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   !-- Local
   INTEGER        :: i, j
   !INTEGER (Sedp) :: FTemswitch, FBIOswitch, FTSwitch
   REAL (Sedp)    :: FBIO_O2, FBIO_NO3, FBIO_MnO2, FBIO_FeOH, FBIO_SO4, FBIO_CH4
   REAL (Sedp)    :: FTem_O2, FTem_NO3, FTem_MnO2, FTem_FeOH, FTem_SO4, FTem_Met
   REAL (Sedp)    :: FTEA_O2, FTEA_NO3, FTEA_MnO2, FTEA_FeOH, FTEA_SO4, FTEA_CH4                         ! Dan added
   REAL (Sedp)    :: FIN_O2, FIN_NO3, FIN_MnO2, FIN_FeOH, FIN_SO4, FIN_CH4                               ! Dan added
   REAL (Sedp)    :: FTterm
   !REAL (Sedp)    :: FTMetOAc, FTMetH2, FTSulOAc, FTSulH2, FTIroOAc, FTIroH2, FTManOAc, FTManH2
   !REAL (Sedp)    :: FTDenOAc, FTDenH2, FTAerOAc, FTAerH2, FTFerOAc, FTFerH2
   !REAL (Sedp)    :: FDHyd, FH2, FOAc, FBHyd
   !-----------------------------
   !--- PRIMARY REDOX REACTIONS - SUBROUTINE RATES
!Dan starts here
!Firstly, make the limiting Factors
!FOM (Substrate limitation)

!pause
DO i=1, candi%npt !
      IF (candi%param%FOMswitch == 1) THEN
          candi%FDHyd = 1.00
          candi%FOAc  = 1.00
          candi%FH2   = 1.00
      ELSE IF (candi%param%FOMswitch == 2) THEN
          candi%FDHyd = candi%y(candi%dhydy,i)/( candi%y(candi%dhydy,i) + candi%param%KDHyd  )
          candi%FOAc  = candi%y(candi%OAcy, i)/( candi%y(candi%OAcy, i) + candi%param%KOAc   )
          candi%FH2   = candi%y(candi%H2y,  i)/( candi%y(candi%H2y,  i) + candi%param%KH2    )
      END IF ! End if FDOMswitch
END DO

!FTem (Temperature factor) - SUBROUTINE RATES
      IF (candi%param%FTemswitch == 1) THEN
          FTem_O2   = 1.000000000000000
          FTem_NO3  = 1.000000000000000
          FTem_MnO2 = 1.000000000000000
          FTem_FeOH = 1.000000000000000
          FTem_SO4  = 1.000000000000000
          FTem_Met  = 1.000000000000000
 !         print *, 'Tem switch = 1'
      ELSE IF (candi%param%FTemswitch == 2) THEN
          FTem_O2  = 1.000000000000000
          FTem_NO3 = 1.000000000000000
          FTem_MnO2= 1.000000000000000
          FTem_FeOH= 1.000000000000000
          FTem_SO4 = 1.000000000000000
          FTem_Met = 1.000000000000000
!          print *, 'Tem switch = 2'
      END IF !End if Temperature switch = 1 or 2
!---------------------------------------------------------------
!----------        Thermodynamic limit on Rox
!          If FTswitch = 1, FT switch is effectively off         - SUBROUTINE RATES
      IF (candi%param%FTswitch == 1) THEN
              candi%FTMetOAc  = 1
              candi%FTMetH2   = 1
              candi%FTSulOAc  = 1
              candi%FTSulH2   = 1
              candi%FTIroOAc  = 1
              candi%FTIroH2   = 1
              candi%FTManOAc  = 1
              candi%FTDenOAc  = 1
              candi%FTDenH2   = 1
              candi%FTAerOAc  = 1
              candi%FTFerDHyd = 1
!          If FTswitch = 1, FT switch is effectively off         - SUBROUTINE RATES
        ELSEIF (candi%param%FTswitch == 2 .AND. candi%param%OMModel < 3) THEN
              print*,'Error: Switch OMModel to 3 or FTSwitch to 1, please'
              STOP
      ELSEIF (candi%param%FTswitch == 2 .AND. candi%param%OMModel == 3) THEN
         candi%dGmp         =  candi%param%n*candi%param%F*candi%param%dPsi
         !---------- FT_Fer start ---------------------
         DO i= 1, candi%npt ! Search "Shark"
            candi%param%Temporary_proton = 10**(-(candi%y(candi%protony,i)))
            IF (candi%y(candi%DHydy,i)<candi%param%Tiny) THEN
                  candi%FTFerDHyd = 0.
            ELSE
               IF (candi%y(candi%hco3y,i)<candi%param%Tiny .OR. candi%y(candi%OAcy,i)<candi%param%Tiny .OR. candi%y(candi%H2y,i)<candi%param%Tiny) THEN
                  candi%FTFerDHyd = 1.
               ELSE
                  candi%dGFerDHyd = candi%param%dG0FerDHyd + candi%param%FTR*candi%param%FTT*log( (candi%y(candi%hco3y,i)**2)*(candi%param%Temporary_proton**4)*(candi%y(candi%OAcy,i)**2)*(candi%y(candi%H2y,i)**4) / ((candi%y(candi%DHydy,i)**1) ) )
                  candi%FTFerDHyd = 1/( candi%param%e** ((candi%dGAerOAc+candi%dGmp)/(candi%param%FTR*candi%param%FTT) ) + 1 )

               END IF
            END IF
         END DO ! Shark
!---------- FT_Fer end ---------------------

!---------- FT_O2 start ---------------------
DO i= 1, candi%npt ! Search "Albatross"
candi%param%Temporary_proton = 10**(-(candi%y(candi%protony,i)))
 IF (candi%y(candi%OAcy,i)<candi%param%Tiny .OR. candi%y(candi%O2y,i)<candi%param%Tiny) THEN
 candi%FTAerOAc = 0.
 ELSE
 IF (candi%y(candi%hco3y,i)<candi%param%Tiny)THEN
 candi%FTAerOAc = 1.
 ELSE
 candi%dGAerOAc = candi%param%dG0AerOAc + candi%param%FTR*candi%param%FTT*log(  (candi%y(candi%hco3y,i)**2) * (candi%param%Temporary_proton**1) / ( (candi%y(candi%OAcy,i)**1) * (candi%y(candi%O2y,i)**2) ) )
 !print*,"dGAerOAc",dGAerOAc
 candi%FTAerOAc = 1/( candi%param%e** ((candi%dGAerOAc+candi%dGmp)/(candi%param%FTR*candi%param%FTT) ) + 1 )
 END IF
 END IF
END DO ! Albatross
!---------- FT_O2 end ---------------------

!---------- FT_NO3 start ---------------------
!print*,"N2",y(N2y,:)
DO i = 1, candi%npt ! Fruit bat
candi%param%Temporary_proton = 10**(-(candi%y(candi%protony,i)))
 IF (candi%y(candi%OAcy,i)<candi%param%Tiny .OR. candi%y(candi%NO3y,i)<candi%param%Tiny)THEN
         candi%FTDenOAc = 0.
 ELSE
 IF (candi%y(candi%N2y,i)<candi%param%Tiny)THEN
         candi%FTDenOAc = 1.
 ELSE
 candi%dGDenOAc = candi%param%dG0DenOAc + candi%param%FTR*candi%param%FTT*log(  (candi%y(candi%N2y,i)**0.8) / ( (candi%y(candi%OAcy,i)**1) * (candi%y(candi%NO3y,i)**1.6) * (candi%param%Temporary_proton**0.6) ) )
 candi%FTDenOAc = 1/( candi%param%e** ((candi%dGDenOAc+candi%dGmp)/(candi%param%FTR*candi%param%FTT) ) + 1 )
 END IF
 END IF
 IF (candi%y(candi%H2y,i)<candi%param%Tiny .OR. candi%y(candi%NO3y,i)<candi%param%Tiny)THEN
         candi%FTDenH2 = 0.
         ELSE
 IF (candi%y(candi%N2y,i)<candi%param%Tiny)THEN
        candi%FTDenH2 = 1.
        ELSE
 candi%dGDenH2 = candi%param%dG0DenH2 + candi%param%FTR*candi%param%FTT*log(  (candi%y(candi%N2y,i)**0.2) / ( (candi%y(candi%H2y,i)**1) * (candi%y(candi%NO3y,i)**0.25) * (candi%param%Temporary_proton**0.4) ) )
 candi%FTDenH2 = 1/( candi%param%e** ((candi%dGDenH2+candi%dGmp)/(candi%param%FTR*candi%param%FTT) ) + 1 )
 END IF
 END IF
END DO ! Search "Fruit bat"
!---------- FT_NO3 end ---------------------

!---------- FT_MnO2 start ---------------------
DO i = 1, candi%npt ! Search "Yabbies"
candi%param%Temporary_proton = 10**(-(candi%y(candi%protony,i)))
 IF (candi%y(candi%OAcy,i)<candi%param%Tiny.OR.candi%y(candi%MnO2y,i)<candi%param%Tiny) THEN
 candi%FTManOAc = 0.
 ELSE
 IF (candi%y(candi%mniiy,i)<candi%param%Tiny.OR.candi%y(candi%hco3y,i)<candi%param%Tiny) THEN
 candi%FTManOAc = 1.
 ELSE
 candi%dGManOAc = candi%param%dG0ManOAc + candi%param%FTR*candi%param%FTT*log(  (candi%y(candi%mniiy,i)**4) * (candi%y(candi%hco3y,i)**2) / ( (candi%y(candi%OAcy,i)**1) * (candi%y(candi%MnO2y,i)**4) * (candi%param%Temporary_proton**7) ) )
 candi%FTManOAc = 1/( candi%param%e** ((candi%dGManOAc+candi%dGmp)/(candi%param%FTR*candi%param%FTT) ) + 1 )
 END IF
 END IF
END DO ! Search "Yabbies"
!---------- FT_MnO2 end ---------------------

!---------- FT_FeOH start ---------------------
DO i = 1, candi%npt ! Search "Emu"
candi%param%Temporary_proton = 10**(-(candi%y(candi%protony,i)))
 IF (candi%y(candi%OAcy,i)<candi%param%Tiny.OR. candi%y(candi%feohy,i)<candi%param%Tiny) THEN
 candi%FTIroOAc = 0.
 ELSE
 IF (candi%y(candi%feiiy,i)<candi%param%Tiny) THEN
 candi%FTIroOAc = 1.
 ELSE
 candi%dGIroOAc = candi%param%dG0IroOAc + candi%param%FTR*candi%param%FTT*log(  (candi%y(candi%feiiy,i)**8) / ( (candi%y(candi%OAcy,i)**1) * (candi%y(candi%feohy,i)**2) * (candi%param%Temporary_proton**4) ) )
 candi%FTIroOAc = 1/( candi%param%e** ((candi%dGIroOAc+candi%dGmp)/(candi%param%FTR*candi%param%FTT) ) + 1 )
 END IF
 END IF
 IF (candi%y(candi%H2y,i)<candi%param%Tiny .OR. candi%y(candi%feohy,i)<candi%param%Tiny) THEN
 candi%FTIroH2 = 0.
 IF (candi%y(candi%feiiy,i)<candi%param%Tiny) THEN
 candi%FTIroH2 = 1.
 ELSE
 candi%dGIroH2 = candi%param%dG0IroH2 + candi%param%FTR*candi%param%FTT*log(  (candi%y(candi%feiiy,i)**2) / ( (candi%y(candi%H2y,i)**1) * (candi%y(candi%feohy,i)**2) * (candi%param%Temporary_proton)**4 ) )
 candi%FTIroH2 = 1/( candi%param%e** ((candi%dGIroH2+candi%dGmp)/(candi%param%FTR*candi%param%FTT) ) + 1 )
 END IF
 END IF
 END DO ! Search "Emu"
!---------- FT_FeOH end ---------------------

!---------- FT_SO4 start ---------------------
DO i = 1, candi%npt ! Search "Galahs"
candi%param%Temporary_proton = 10**(-(candi%y(candi%protony,i)))
 IF (candi%y(candi%OAcy,i)<candi%param%Tiny.OR.candi%y(candi%so4y,i)<candi%param%Tiny) THEN
 candi%FTSulOAc = 0.
 ELSE
 IF (candi%y(candi%hsy,i)<candi%param%Tiny) THEN
 candi%FTSulOAc = 1.
 ELSE
 candi%dGSulOAc = candi%param%dG0SulOAc + candi%param%FTR*candi%param%FTT*log(  (candi%y(candi%hsy,i)**1) * (candi%y(candi%hco3y,i)**2) / ( (candi%y(candi%OAcy,i)**1) * (candi%y(candi%so4y,i)**1)  ) )
 candi%FTSulOAc = 1/( candi%param%e** ((candi%dGSulOAc+candi%dGmp)/(candi%param%FTR*candi%param%FTT) ) + 1 )
 END IF
 END IF
 IF (candi%y(candi%H2y,i)<candi%param%Tiny.OR.candi%y(candi%so4y,i)<candi%param%Tiny) THEN
 candi%FTSulH2 = 0.
 ELSE
 IF (candi%y(candi%hsy,i)<candi%param%Tiny) THEN
 candi%FTSulH2 = 1.
 ELSE
 candi%dGSulH2 = candi%param%dG0SulH2 + candi%param%FTR*candi%param%FTT*log(  (candi%y(candi%hsy,i)**0.25) / ( (candi%y(candi%H2y,i)**1) * (candi%y(candi%so4y,i)**0.25) * (candi%param%Temporary_proton)**0.25 ) )
 candi%FTSulH2 = 1/( candi%param%e** ((candi%dGSulH2+candi%dGmp)/(candi%param%FTR*candi%param%FTT) ) + 1 )
 END IF
 END IF
END DO ! Search "Galahs"
!---------- FT_SO4 end ---------------------

!---------- FT_Met start ---------------------
DO i = 1, candi%npt ! Search "Donkey"
candi%param%Temporary_proton = 10**(-(candi%y(candi%protony,i)))
 IF (candi%y(candi%OAcy,i)<candi%param%Tiny) THEN
 candi%FTMetOAc = 0.
 ELSE
 IF (candi%y(candi%CH4y,i)<candi%param%Tiny.OR.candi%y(candi%hco3y,i)<candi%param%Tiny) THEN
 candi%FTMetOAc = 1.
 ELSE
 candi%dGMetOAc = candi%param%dG0MetOAc + candi%param%FTR*candi%param%FTT*log(  (candi%y(candi%CH4y,i)**1) * (candi%y(candi%hco3y,i)**1) / ( (candi%y(candi%OAcy,i)**1)  ) )
 candi%FTMetOAc = 1/( candi%param%e** ((candi%dGMetOAc+candi%dGmp)/(candi%param%FTR*candi%param%FTT) ) + 1 )
 END IF
 END IF
 IF (candi%y(candi%H2y,i)<candi%param%Tiny.OR.candi%y(candi%hco3y,i)<candi%param%Tiny) THEN
 candi%FTMetH2 = 0.
 ELSE
 IF (candi%y(candi%CH4y,i)<candi%param%Tiny) THEN
 candi%FTMetH2 = 1.
 ELSE
 candi%dGMetH2 = candi%param%dG0MetH2 + candi%param%FTR*candi%param%FTT*log(  (candi%y(candi%CH4y,i)**0.25) / ( (candi%y(candi%H2y,i)**1) * (candi%y(candi%hco3y,i)**0.25) * (candi%param%Temporary_proton)**0.25 ) )
 candi%FTMetH2 = 1/( candi%param%e** ((candi%dGMetH2+candi%dGmp)/(candi%param%FTR*candi%param%FTT) ) + 1 )
 END IF
 END IF
END DO ! Search "Donkey"
!---------- FT_Met end ---------------------

        END IF ! End if Thermodynamic switch = 1 or 2
!print *, ('FT_SO4sum'),FT_SO4sum
  !-- Define kinetic factors FTEA and FIN - SUBROUTINE RATES
  !-----------------A1-----------------------A1A1A1A1A1A1A1A1A1A1A1A1 - SUBROUTINE RATES
  !FTEA (Monod limitation at low concentration of oxidants)
IF(candi%param%OMapproach < 1 .OR. candi%param%OMapproach > 2) THEN
     PRINT *, ('OMapproach must be 1 or 2. STOPPING')
     STOP
END IF  ! End if OM approach \= 1 or 2

! Start Define FIN (Monod inhibition of an oxidation process by another oxidant)
DO i=1, candi%npt ! Search "Wombat"
        IF(candi%param%FINswitch ==1) THEN
          FIN_O2   = 1.0000000000000000
          FIN_NO3  = 1.0000000000000000
          FIN_MnO2 = 1.0000000000000000
          FIN_FeOH = 1.0000000000000000
          FIN_SO4  = 1.0000000000000000
        !PRINT *,('FIN 1 Approach 1 or 2 : Wombats')
          ELSE
        IF(candi%param%FINswitch == 2 .AND. candi%param%OMapproach == 1) THEN
           !FIN_O2   = ( candi%param%kpo2  /(candi%param%kpo2+candi%y(candi%o2y,i))     )
           FIN_NO3  = ( candi%param%kpno3 /(candi%param%kpno3+candi%y(candi%no3y,i))   )
           FIN_MnO2 = ( candi%param%kpmno2/(candi%param%kpmno2+candi%y(candi%mno2y,i)) )
           FIN_FeOH = ( candi%param%kpfeoh/(candi%param%kpfeoh+candi%y(candi%feohy,i)) )
           FIN_SO4  = ( candi%param%kpso4 /(candi%param%kpso4+candi%y(candi%so4y,i))   )
        !PRINT *,('FIN 2 Approach 1 : Badgers')
        ELSE
        IF(candi%param%FINswitch ==2 .AND. candi%param%OMapproach == 2) THEN
          !FIN_O2   = (1.0 - (candi%y(candi%o2y,i)/candi%param%lpo2)    )
          FIN_NO3  = (1.0 - (candi%y(candi%NO3y,i)/candi%param%lpNO3)  )
          FIN_MnO2 = (1.0 - (candi%y(candi%MnO2y,i)/candi%param%lpMnO2))
          FIN_FeOH = (1.0 - (candi%y(candi%FeOHy,i)/candi%param%lpFeOH))
          FIN_SO4  = (1.0 - (candi%y(candi%SO4y,i)/candi%param%lpSO4)  )
        !PRINT *,('FIN 2 Approach 2 : Cows')
          ELSE
         STOP

        END IF ! End if FINswitch == 2 and OMapproach == 2
        END IF ! End if FINswitch == 2 and OMapproach == 1
        END IF ! End if FINswitch == 1
!END DO
       IF (candi%param%FInO2OnlySwitch==1)THEN
               FIN_O2 = 1.000
       ELSE
       !DO i=1, npt
       IF (candi%param%FInO2OnlySwitch==2)THEN
               IF(candi%param%OMapproach == 1) THEN
               FIN_O2 = ( candi%param%kpo2/(candi%param%kpo2+candi%y(candi%o2y,i)) )
               ELSEIF(candi%param%OMapproach == 1) THEN
               FIN_O2 = (1.0 - (candi%y(candi%o2y,i)/candi%param%lpo2))
               END IF ! End if OMapproach
       ELSE
       IF (candi%param%FInO2OnlySwitch>2 .OR. candi%param%FInO2OnlySwitch<1)THEN
       print*,"Set FInOnlySwitch to 1 or 2, please."
              print*,"FInO2OnlySwitch",candi%param%FInO2OnlySwitch
       STOP
       END IF ! End if FInOnlySwitch
       END IF
       END IF
END DO  ! Search "Wombat"

! End Define FIN (Monod inhibition of an oxidation process by another oxidant)

! Start FTEA and FIN calculations for ROX equations
DO i=1, candi%npt  !Search "Elephant"
   IF ( candi%param%OMapproach == 1 ) THEN
                   FTEA_O2   = ( candi%y(candi%o2y,i)  /(candi%param%kO2+candi%y(candi%o2y,i))      )
                   FTEA_NO3  = ( candi%y(candi%NO3y,i) /(candi%param%kNO3+candi%y(candi%NO3y,i))   )
                   FTEA_MnO2 = ( candi%y(candi%MnO2y,i)/(candi%param%kMnO2+candi%y(candi%MnO2y,i)))
                   FTEA_FeOH = ( candi%y(candi%FeOHy,i)/(candi%param%kFeOH+candi%y(candi%FeOHy,i)))
                   FTEA_SO4  = ( candi%y(candi%SO4y,i) /(candi%param%kSO4+candi%y(candi%SO4y,i))   )
  !Rate limiting factors for Approach 1
  !You multiply these by rgC to get the fraction of the organic matter oxidation that is consuming each TEA
          !-- Aerobic respiration
          candi%FO2(i)   = FTEA_O2  *FTem_O2  * candi%FTAerOAc(i)
          !-- Denitrification
          candi%FNO3(i)  = FTEA_NO3 *FTem_NO3 *FIN_O2 * candi%FTDenOAc(i) * candi%FTDenH2(i)
          !-- Manganese reduction
          candi%FMnO2(i) = FTEA_MnO2*FTem_MnO2*FIN_O2*FIN_NO3 * candi%FTManOAc(i)
          !-- Iron reduction
          candi%FFeOH(i) = FTEA_FeOH*FTem_FeOH*FIN_O2*FIN_NO3*FIN_MnO2 * candi%FTIroOAc(i) * candi%FTIroH2(i)
          !-- Sulfate reduction
          candi%FSO4(i)  = FTEA_SO4 *FTem_SO4 *FIN_O2*FIN_NO3*FIN_MnO2*FIN_FeOH * candi%FTSulOAc(i) * candi%FTSulH2(i)
          !-- Methanogenesis
          candi%FMet(i)  =           FTem_Met *FIN_O2*FIN_NO3*FIN_MnO2*FIN_FeOH*FIN_SO4 * candi%FTMetOAc(i) * candi%FTMetH2(i)
END IF !End if Approach 1
END DO ! Search "Elephant"

!-----------------A2-----------------------A2A2A2A2A2A2A2A2A2A2A2 - SUBROUTINE RATES
DO i=1,candi%npt ! Search "Tiger"
IF ( candi%param%OMapproach == 2 ) THEN
   !DO i=1,npt
      !FTEA (Monod limitation at low concentration of oxidant)
          FTEA_O2   = candi%y(candi%o2y,i)  /candi%param%lo2
          FTEA_NO3  = candi%y(candi%no3y,i) /candi%param%lno3
          FTEA_MnO2 = candi%y(candi%mno2y,i)/candi%param%lmno2
          FTEA_FeOH = candi%y(candi%feohy,i)/candi%param%lfeoh
          FTEA_SO4  = candi%y(candi%so4y,i) /candi%param%lso4
!You multiply these by rgC to get the fraction of the organic matter oxidation that is consuming each TEA
!print *,'i',i
!PRINT *,'A2 TEST ',lo2,lpo2,lno3,lpno3,lmno2,lpmno2,lfeoh,lpfeoh,lso4,lpso4
!PRINT *,'FTi     ',FTEA_O2,FIN_O2,FTEA_NO3,FIN_NO3,FTEA_MnO2,FIN_MnO2, FTEA_FeOH,FIN_FeOH,FTEA_SO4,FIN_SO4
      !-- Aerobic respiration
          IF ( candi%y(candi%o2y,i)>candi%param%lO2 ) THEN ! If O2 concentration is high ...
           candi%FO2(i)  = 1.0000000000000000*FTem_O2 ! ... then aerobic respiration goes at 1*kOM ...
           candi%FNO3(i) = 0.0000000000000000 !... and everything else is inhibited by O2.
           candi%FMnO2(i)= 0.0000000000000000
           candi%FFeOH(i)= 0.0000000000000000
           candi%FSO4(i) = 0.0000000000000000
           candi%FMet(i) = 0.0000000000000000
           ELSE
      !-- Denitrification
          IF ( candi%y(candi%NO3y,i)>candi%param%lNO3 ) THEN ! If NO3 concentration is high and if O2 concentration is low ...
          candi%FO2(i)  = FTEA_O2*FTem_O2 ! ... 1: then aerobic respiration is O2 limited ...
          candi%FNO3(i) = 1.0000000000000000*FTem_NO3*FIN_O2 ! ... denitrification goes at 1*kOM ...
          candi%FMnO2(i)= 0.0000000000000000 !... and everything else is inhibited by NO3.
          candi%FFeOH(i)= 0.0000000000000000
          candi%FSO4(i) = 0.0000000000000000
          candi%FMet(i) = 0.0000000000000000
          ELSE
      !-- Manganese reduction
            IF ( candi%y(candi%MnO2y,i)>candi%param%lMnO2  ) THEN       !If Mn concentration is high and if NO3 is low ...
             candi%FO2(i)  = FTEA_O2*FTem_O2 ! ... monkeys 2: then aerobic respiration is O2 limited ...
             candi%FNO3(i) = FTem_NO3*FTEA_NO3*FIN_O2 ! ... then denitrification is NO3 limited,
             candi%FMnO2(i)= 1.0000000000000000*FTem_MnO2*FIN_NO3*FIN_O2 ! ... manganese reduction goes at 1*kOM ...
             candi%FFeOH(i)= 0.0000000000000000 !... and everything else is inhibited by MnO2.
             candi%FSO4(i) = 0.0000000000000000
             candi%FMet(i) = 0.0000000000000000
            ELSE
      !-- Iron reduction
               IF ( candi%y(candi%FeOHy,i)>candi%param%lFeOH ) THEN     !If iron concentration is high and MnO2 concentration is low ...
               candi%FO2(i)  = FTEA_O2*FTem_O2 ! ... monkeys 3: then aerobic respiration is O2 limited ...
               candi%FNO3(i) = FTem_NO3*FTEA_NO3*FIN_O2 ! ... then denitrification is NO3 limited,
               candi%FMnO2(i)= FTem_MnO2*FTEA_MnO2*FIN_NO3*FIN_O2 !... then manganese reduction is manganese limited ...
               candi%FFeOH(i)= 1.0000000000000000*FTem_FeOH*FIN_MnO2*FIN_NO3*FIN_O2! ... iron reduction goes at 1*kOM ...
               candi%FSO4(i) = 0.0000000000000000 !... and the others are inhibited by FeOH.
               candi%FMet(i) = 0.0000000000000000
               ELSE
     !-- Sulfate reduction
                 IF ( candi%y(candi%SO4y,i)>candi%param%lSO4 ) THEN     !If SO4 concentration is high and iron concentration is low ...
                 candi%FO2(i)  = FTEA_O2*FTem_O2 ! ... monkeys 4: then aerobic respiration is O2 limited ...
                 candi%FNO3(i) = FTem_NO3*FTEA_NO3*FIN_O2 ! ... then denitrification is NO3 limited,
                 candi%FMnO2(i)= FTem_MnO2*FTEA_MnO2*FIN_NO3*FIN_O2 !... then manganese reduction is manganese limited ...
                 candi%FFeOH(i)= FTem_FeOH*FTEA_FeOH*FIN_MnO2*FIN_NO3*FIN_O2 !... then the iron reduction rate is iron limited ...
                 candi%FSO4(i) = 1.0000000000000000*FTem_SO4*FIN_FeOH*FIN_MnO2*FIN_NO3*FIN_O2 ! ... SO4 reduction goes at 1*kOM ...
                 candi%FMet(i) = 0.0000000000000000 !... and methanogenesis is inhibited by SO4.
                 ELSE
      !-- Methanogenesis
                 ! If SO4 concentration is low ...
                 candi%FO2(i)  = FTEA_O2*FTem_O2 ! ... monkeys 5: then aerobic respiration is O2 limited ...
                 candi%FNO3(i) = FTem_NO3*FTEA_NO3*FIN_O2 ! ... then denitrification is NO3 limited,
                 candi%FMnO2(i)= FTem_MnO2*FTEA_MnO2*FIN_NO3*FIN_O2 !... then manganese reduction is manganese limited ...
                 candi%FFeOH(i)= FTem_FeOH*FTEA_FeOH*FIN_MnO2*FIN_NO3*FIN_O2 !... then the iron reduction rate is iron limited ...
                 candi%FSO4(i) = FTem_SO4*FTEA_SO4*FIN_FeOH*FIN_MnO2*FIN_NO3*FIN_O2 !... then sulfate reduction is SO4 limited.
                 candi%FMet(i) = FTem_Met*FIN_SO4*FIN_FeOH*FIN_MnO2*FIN_NO3*FIN_O2 ! If methanogenesis happens at all, it is inhibited by all the other oxidants. If the others are low, then it reacts at close to 1*kOM.
                 END IF !-- Sulfate is high
              END IF !-- Iron is high
            END IF ! -- Manganese is high
          END IF ! -- Nitrate is high
        END IF ! -- End if oxygen is high
!PRINT *,'R ',RO2(i),RNO3(i),RMnO2(i),RFeOH(i),RSO4(i),RCH4(i)
!PAUSE
!END DO
END IF !End if Approach 2
END DO ! Search "Tiger"
! Make sure the FTEA don't go negative
DO i=1,candi%npt ! Search "Penguin"
        IF (FTEA_O2<0.0000000000000000000000000000000000) THEN
        FTEA_O2 = 1.00000E-20
        END IF ! End if FTEA_02 <0
!END DO !
!DO i=1,candi%npt
        IF (FTEA_NO3<0.0000000000000000000000000000000000) THEN
        FTEA_NO3 = 1.00000E-20
        END IF ! End if FTEA_02 <0
!END DO !
!DO i=1,npt
        IF (FTEA_MnO2<0.0000000000000000000000000000000000) THEN
        FTEA_MnO2 = 1.00000E-20
        END IF ! End if FTEA_02 <0
!END DO !
!DO i=1,npt
        IF (FTEA_FeOH<0.0000000000000000000000000000000000) THEN
        FTEA_FeOH = 1.00000E-20
        END IF ! End if FTEA_02 <0
!END DO !
!DO i=1,npt
        IF (FTEA_SO4<0.0000000000000000000000000000000000) THEN
        FTEA_SO4 = 1.00000E-20
        END IF ! End if FTEA_02 <0
END DO ! Search "Penguin"

! Make sure the rates don't go negative
DO i=1,candi%npt ! Search "Ferret"
        IF (candi%FO2(i)<0.00000000000000000000000000000000000000000) THEN
        candi%FO2(i) = 1.00000E-20
        END IF !
!END DO !
!DO i=1,npt
        IF (candi%FNO3(i)<0.00000000000000000000000000000000000000000) THEN
        candi%FNO3(i) = 1.00000E-20
        END IF !
!END DO !
!DO i=1,npt
        IF (candi%FMnO2(i)<0.00000000000000000000000000000000000000000) THEN
        candi%FMnO2(i) = 1.00000E-20
        END IF !
!END DO !
!DO i=1,npt
        IF (candi%FFeOH(i)<0.00000000000000000000000000000000000000000) THEN
        candi%FFeOH(i) = 1.00000E-20
        END IF !
!END DO !
!DO i=1,npt
        IF (candi%FSO4(i)<0.00000000000000000000000000000000000000000) THEN
        candi%FSO4(i) = 1.00000E-20
        END IF !
END DO ! Search "Ferret"

! Start FTEA and FIN calculations for ROX equations
   !-----------------------------
   !--- SECONDARY REDOX REACTIONS - SUBROUTINE RATES

   !-- NH4 oxidation by O2
   candi%rnh4ox(:)  = candi%param%knh4ox  * candi%y(candi%o2y,:)*candi%y(candi%nh4y,:)
     DO i=1,candi%npt ! Search keyword "Eagle"
     IF (candi%rnh4ox(i) < 0.000000000000000000000000000000000000000) THEN
     candi%rnh4ox(i) = 1.000000E-20
     END IF ! End if RNH4Ox<0
     END DO ! Search keyword "Eagle"
   !-- H2S oxidation by O2
   candi%rtsox(:)   = candi%param%ktsox   * candi%y(candi%o2y,:)*candi%y(candi%hsy,:)
     DO i=1,candi%npt
     IF (candi%rtsox(i) < 0.000000000000000000000000000000000000000) THEN
     candi%rtsox(i) = 1.000000E-20
     END IF
     END DO
   !-- CH4 oxidation by O2
   candi%rch4ox(:)  = candi%param%kch4ox  * candi%y(candi%o2y,:)*candi%y(candi%ch4y,:)
     !DO i=1,npt !Causes problems
     !IF (rch4ox(i) < 0.000000) THEN  !Causes problems
     !rch4ox(i) = 1.000000E-40  !Causes problems
     !END IF !Causes problems
     !END DO !Causes problems
   !-- H2S oxidation by NO3
   candi%RTSNO3(:)  = candi%param%kTSNO3  * candi%y(candi%no3y,:)*candi%y(candi%hsy,:)
     DO i=1,candi%npt !
     IF (candi%rtsNO3(i) < 0.000000000000000000000000000000000000000) THEN
     candi%rtsNO3(i) = 1.000000E-20
     END IF ! End if rtsNO3<0
     END DO
   !-- NH4 oxidation by NO2
   candi%rnh4no2(:) = candi%knh4no2 * candi%y(candi%o2y,:)*candi%y(candi%nh4y,:)
      DO i=1, candi%npt
     IF (candi%rnh4no2(i) < 0.000000000000000000000000000000000000000) THEN
     candi%rnh4no2(i) = 1.000000E-20
     END IF ! End if rnh4no2<0
     END DO
   !-- CH4 oxidation by SO4
   candi%rch4so4(:) = candi%param%kch4so4 * candi%y(candi%so4y,:)*candi%y(candi%ch4y,:)
     DO i=1, candi%npt
     IF (candi%rch4so4(i) < 0.000000000000000000000000000000000000000) THEN
     candi%rch4so4(i) = 1.000000E-20
     END IF ! End if rnh4no2<0                                         !  - SUBROUTINE RATES
     END DO
!DO i=1,npt ! Search keyword "Taipan"
   IF(candi%param%simMnFe) THEN

     !-- Mn2+ oxidation by O2
     candi%rmnox(:)  = candi%param%kmnox  * candi%y(candi%o2y,:)*candi%y(candi%mniiy,:)
     DO i=1, candi%npt
     IF (candi%rmnox(i) < 0.0000000000000000000000000000000000000) THEN
     candi%rmnox(i) = 1.000000E-20
     END IF
     END DO !
     !-- Fe2+ oxidation by O2
     candi%rfeox(:)  = candi%param%kfeox  * candi%y(candi%o2y,:)*candi%y(candi%feiiy,:)
     DO i=1, candi%npt
     IF (candi%rfeox(i) < 0.0000000000000000000000000000000000000) THEN
     candi%rfeox(i) = 1.000000E-20
     END IF !
     END DO
     !-- Fe2+ oxidation by MnO2A & MnO2B
     candi%rfemnA(:)  = candi%param%kmnfe  * candi%y(candi%mno2y,:)*candi%y(candi%feiiy,:)
     DO i=1, candi%npt
     IF (candi%rfemnA(i) < 0.0000000000000000000000000000000000000) THEN
     candi%rfemnA(i) = 1.000000E-20
     END IF !
     END DO
     candi%rfemnB(:)  = candi%param%kmnfe  * candi%y(candi%mno2By,:)*candi%y(candi%feiiy,:)
!     DO i=1, npt ! Causes problems
!     IF (rfemnB(i) < 0.0000) THEN ! Causes problems
!     rfemnB(i) = 1.000000E-40 ! Causes problems
!     END IF !  ! Causes problems
!     END DO ! Causes problems
     !-- Fe2+ oxidation by NO3
     candi%rfeno3(:) = candi%param%kfeno3 * candi%y(candi%no3y,:)*candi%y(candi%feiiy,:)
     DO i=1,candi%npt ! Search keyword "Water rat"
     IF (candi%rfeno3(i) < 0.00000000000000000000000000000000000000) THEN
     candi%rfeno3(i) = 1.000000E-20
     END IF ! End if RfeNO3<0
     END DO ! Search keyword "Water rat"
     !-- H2S oxidation by MnO2A + MnO2B
     candi%rtsmnA(:)  = candi%param%ktsmn  * candi%y(candi%mno2y,:)*candi%y(candi%hsy,:)
     DO i=1, candi%npt
     IF (candi%rtsmnA(i) < 0.0000000000000000000000000000000000) THEN
     candi%rtsmnA(i) = 1.000000E-20
     END IF !
     END DO
     candi%rtsmnB(:)  = candi%param%ktsmn  * candi%y(candi%mno2By,:)*candi%y(candi%hsy,:)
     DO i=1, candi%npt
     IF (candi%rtsmnB(i) < 0.0000000000000000000000000000000000) THEN
     candi%rtsmnB(i) = 1.000000E-20
     END IF !
     END DO
     !-- Mn2+ oxidation by NO3
     candi%rmnno3(:) = candi%param%kmnno3 * candi%y(candi%no3y,:)*candi%y(candi%mniiy,:)
     DO i=1,candi%npt ! Search keyword "Osprey"
     IF (candi%rmnno3(i) < 0.0000000000000000000000000000000000000) THEN
     candi%rmnno3(i) = 1.000000E-20
     END IF ! End if RMnNO3<0
     END DO !"Osprey"
     !-- H2S oxidation by Fe(OH)3A & Fe(OH)3B
     candi%rtsfeA(:)  = candi%param%ktsfe  * candi%y(candi%hsy,:)*candi%y(candi%feohy,:)
     DO i=1, candi%npt
     IF (candi%rtsfeA(i) < 0.0000000000000000000000000000000000000) THEN
     candi%rtsfeA(i) = 1.000000E-20
     END IF !
     END DO
     candi%rtsfeB(:)  = candi%param%ktsfe  * candi%y(candi%hsy,:)*candi%y(candi%feohBy,:)
     DO i=1, candi%npt
     IF (candi%rtsfeB(i) < 0.0000000000000000000000000000000000000) THEN
     candi%rtsfeB(i) = 1.000000E-20
     END IF !
     END DO
     !-- MnO2A ageing
     candi%rmnage(:) = candi%param%kmnage * candi%y(candi%mno2y,:)
     DO i=1, candi%npt
     IF (candi%rmnage(i) < 0.00000000000000000000000000000000000000000) THEN
     candi%rmnage(i) = 1.000000E-20
     END IF !
     END DO
                                          !  - SUBROUTINE RATES
                                          !  - If sim Mn Fe
     !-- Fe(OH)3A ageing
     candi%rfeage(:) = candi%param%kfeage * candi%y(candi%feohy,:)
     DO i=1, candi%npt
     IF (candi%rfeage(i) < 0.00000000000000000000000000000000000000000) THEN
     candi%rfeage(i) = 1.000000E-20
     END IF !
     END DO
    !-- Fe(OH)3A precipitation - based on Tufano 2009
     IF(candi%param%rxn_mode==0)THEN
             candi%RFeOHAppt(:) = 0.00
     ELSEIF(candi%param%rxn_mode==1)THEN
             candi%RFeOHAppt(:) = 0.00
     ELSEIF(candi%param%rxn_mode==2) THEN
       ! kfeohppt [Fe3+]
       candi%RFeOHAppt(:) = candi%param%kFeOHAppt * candi%y(candi%feiiiy,:)
       !RFeOHAppt(:) = kFeOHAppt * y(feohy,:)
     ELSEIF(candi%param%rxn_mode==3) THEN
       WHERE (candi%IAP(candi%FeOHy,:)== 0.00)
         candi%RFeOHAppt(:) = candi%param%kFeOHAppt * candi%y(candi%feiiiy,:)
       ELSEWHERE ( (ONE / candi%IAP(candi%FeOHy,:)) < ONE )
         ! Precipitation
         candi%RFeOHAppt(:) = candi%param%kfeohAppt * ( ONE - (ONE / (candi%IAP(candi%FeOHy,:))) )
       ELSEWHERE
         ! Dissolution
         candi%RFeOHAppt(:) = candi%param%kfeohAppt * (-ONE) * ( ONE - (candi%IAP(candi%FeOHy,:)) )
       END WHERE
       !  print *, ('FeOHAppt'),RFeOHAppt
       !  print *, ('IAP FeOHA'),IAP(FeOHy,:)
     END IF ! End if rxn mode = 0, 1, 2 or 3

    !-- Fe(OH)3B precipitation
     IF(candi%param%rxn_mode==0)THEN
             candi%RFeOHBppt(:) = 0.00
     ELSEIF(candi%param%rxn_mode==1)THEN
             candi%RFeOHBppt(:) = 0.00
     ELSEIF(candi%param%rxn_mode==2) THEN
       ! kfeohbppt [Fe3+]
       candi%RFeOHBppt(:) = candi%param%kFeOHBppt * candi%y(candi%feiiiy,:)
       !RFeOHBppt(:) = kFeOHBppt * y(feohy,:)
       !RFeOHBppt(:) = 1.
     ELSEIF(candi%param%rxn_mode==3) THEN
           WHERE (candi%IAP(candi%FeOHBy,:)== 0.00)
           candi%RFeOHBppt(:) = candi%param%kFeOHBppt * candi%y(candi%feiiiy,:)
       ELSEWHERE ( (ONE / candi%IAP(candi%FeOHBy,:)) < ONE )
         ! Precipitation
         candi%RFeOHBppt(:) = 0.
       ELSEWHERE
         ! Dissolution
         candi%RFeOHBppt(:) = candi%param%kfeohbppt * (-ONE) * ( ONE - (candi%IAP(candi%FeOHy,:)) )
       END WHERE
          END IF ! End if rxn mode = 0, 1, 2 or 3
     !-- MnO2A precipitation
     IF(candi%param%rxn_mode==0)THEN
             candi%RMnO2Appt(:) = 0.00
     ELSEIF(candi%param%rxn_mode==1)THEN
             candi%RMnO2Appt(:) = 0.00
     ELSEIF(candi%param%rxn_mode==2) THEN
             candi%RMnO2Appt(:) = candi%param%kMnO2Appt * candi%y(candi%mniiy,:)
     ELSEIF(candi%param%rxn_mode==3) THEN
      !kMnO2Appt [Mn2+]?
       WHERE (candi%IAP(candi%MnO2y,:) == 0.00 )
             candi%RMnO2Appt(:) = candi%param%kMnO2Appt * candi%y(candi%mniiy,:)
       ELSEWHERE ( (ONE / candi%IAP(candi%MnO2y,:)) < ONE )
       ! Precipitation
         candi%RMnO2Appt(:) = candi%param%kMnO2Appt * ( ONE - (ONE / (candi%IAP(candi%MnO2y,:))) )
       ELSEWHERE
         ! Dissolution
         candi%RMnO2Appt(:) = candi%param%kMnO2Appt * (-ONE) * ( ONE - (candi%IAP(candi%MnO2y,:)) )
       END WHERE
        ! print *, ('MnO2Appt'),RMnO2Appt
        ! print *, ('IAP MnO2A'),IAP(MnO2y,:)
     END IF ! End if rxn mode = 0, 1, 2 or 3

          !-- MnO2B precipitation
     IF(candi%param%rxn_mode==0)THEN
             candi%RMnO2Bppt(:) = 0.00
     ELSEIF(candi%param%rxn_mode==1)THEN
             candi%RMnO2Bppt(:) = 0.00
     ELSEIF(candi%param%rxn_mode==2) THEN
       candi%RMnO2Bppt(:) = candi%param%kMnO2Bppt * candi%y(candi%mniiy,:)
     ELSEIF(candi%param%rxn_mode==3) THEN
      WHERE (candi%IAP(candi%MnO2By,:)==0.00)
       candi%RMnO2Bppt(:) = candi%param%kMnO2Bppt * candi%y(candi%mniiy,:)
      ELSEWHERE ( (ONE / candi%IAP(candi%MnO2By,:)) < ONE )
         ! Precipitation
         candi%RMnO2Bppt(:) = candi%param%kMnO2Bppt * ( ONE - (ONE / (candi%IAP(candi%MnO2By,:))) )
       ELSEWHERE
         ! Dissolution
         candi%RMnO2Bppt(:) = candi%param%kMnO2Bppt * (-ONE) * ( ONE - (candi%IAP(candi%MnO2By,:)) )
       END WHERE

     END IF ! End if rxn mode = 0, 1, 2 or 3
                                          !  - SUBROUTINE RATES
                                          !  - If sim Mn Fe
     !-- FeCO3 precipitation
     IF(candi%param%rxn_mode==0)THEN
             candi%RSidppt(:) = 0.00
     ELSEIF(candi%param%rxn_mode==1)THEN
             candi%RSidppt(:) = 0.00
     ELSEIF(candi%param%rxn_mode==2) THEN
             ! kSidppt [Fe2+][CO3]
       candi%RSidppt(:) = candi%param%kSidppt * candi%y(candi%feiiy,:)*candi%y(candi%hco3y,:)
        !       RSidppt(:) = 1.
     ELSEIF(candi%param%rxn_mode==3) THEN
       WHERE (candi%IAP(candi%Sidy,:)==0.00)
             candi%RSidppt(:) = candi%param%kSidppt * candi%y(candi%feiiy,:)*candi%y(candi%hco3y,:)
       ELSEWHERE ( (ONE / candi%IAP(candi%Sidy,:)) < ONE )
         ! Precipitation
         candi%RSidppt(:) = candi%param%kSidppt * ( ONE - (ONE / (candi%IAP(candi%Sidy,:))) )
       ELSEWHERE
         ! Dissolution
         candi%RSidppt(:) = candi%param%kSidppt * (-ONE) * ( ONE - (candi%IAP(candi%Sidy,:)) )
       END WHERE

       END IF ! End if rxn mode = 0, 1, 2 or 3
                                          !  - SUBROUTINE RATES
                                          !  - If sim Mn Fe
    !-- MnCO3 precipitation
     IF(candi%param%rxn_mode==0)THEN
             candi%RRodppt(:) = 0.00
     ELSEIF(candi%param%rxn_mode==1)THEN
             candi%RRodppt(:) = 0.00
     ELSEIF(candi%param%rxn_mode==2) THEN
           candi%RRodppt(:) = candi%param%kRodppt * candi%y(candi%mniiy,:)*candi%y(candi%hco3y,:)
       !               RRodppt(:) = 1.
     ELSEIF(candi%param%rxn_mode==3) THEN
       WHERE (candi%IAP(candi%Rody,:)==0.00)
             candi%RRodppt(:) = candi%param%kRodppt * candi%y(candi%mniiy,:)*candi%y(candi%hco3y,:)
       ELSEWHERE ( (ONE / candi%IAP(candi%Rody,:)) < ONE )
         ! Precipitation
         candi%RRodppt(:) = candi%param%kRodppt * ( ONE - (ONE / (candi%IAP(candi%Rody,:))) )
       ELSEWHERE
         ! Dissolution
         candi%RRodppt(:) = candi%param%kRodppt * (-ONE) * ( ONE - (candi%IAP(candi%Rody,:)) )
       END WHERE
     END IF ! End if rxn mode = 0, 1, 2 or 3
   ELSE                 !  - Else if not Sim Mn Fe
                        !  - SUBROUTINE RATES
     candi%rmnox(:)     = 0.0
     candi%rfeox(:)     = 0.0
     candi%rfemnA(:)    = 0.0
     candi%rfemnB(:)    = 0.0
     candi%rfeno3(:)    = 0.0
     candi%rtsmnA(:)    = 0.0
     candi%rtsmnB(:)    = 0.0
     candi%rmnno3(:)    = 0.0
     candi%rtsfeA(:)    = 0.0
     candi%rtsfeB(:)    = 0.0
     candi%rmnage(:)    = 0.0
     candi%rfeage(:)    = 0.0
     candi%rfeohappt(:) = 0.0
     candi%rfeohbppt(:) = 0.0
     candi%RSidppt(:)   = 0.0
     candi%RRodppt(:)   = 0.0
     candi%RMnO2Appt(:) = 0.0
     candi%RMnO2Bppt(:) = 0.0

   END IF               !  - End if Sim Mn Fe
                        !  - SUBROUTINE RATES
   IF(candi%param%simFeS .AND. candi%param%simMnFe) THEN

     !-- FeS oxidation by O2
     candi%rfesox(:)  = candi%param%kfesox  * candi%y(candi%o2y,:)*candi%y(candi%fesy,:)
     DO i=1, candi%npt
     IF (candi%rfesox(i) < 0.0) THEN
     candi%rfesox(i) = 0.00 !1.0000000E-14
     END IF !
     END DO
     !-- FeS2 oxidation by O2
     candi%rfes2ox(:) = candi%param%kfes2ox * candi%y(candi%o2y,:)*candi%y(candi%fes2y,:)
     DO i=1, candi%npt
     IF (candi%rfes2ox(i) < 0.0) THEN
     candi%rfes2ox(i) = 0.00 ! 1.0000000E-14
     END IF !
     END DO
     !-- FeS oxidation by Fe(OH)3
     candi%rfesfeA(:)  = candi%param%kfesfe * candi%y(candi%fesy,:)*candi%y(candi%feohy,:)
     DO i=1, candi%npt
     IF (candi%rfesfeA(i) < 0.0) THEN
     candi%rfesfeA(i) = 0.00 ! 1.0000000E-14
     END IF !
     END DO
     candi%rfesfeB(:)  = candi%param%kfesfe * candi%y(candi%fesy,:)*candi%y(candi%feohBy,:)
     DO i=1, candi%npt
     IF (candi%rfesfeB(i) < 0.0) THEN
     candi%rfesfeB(i) =  0.00 !1.0000000E-14
     END IF !
     END DO
     !-- FeS oxidation by MnO2
     candi%rfesmnA(:)  = candi%param%kfesmn * candi%y(candi%fesy,:)*candi%y(candi%mno2y,:)
     DO i=1, candi%npt
     IF (candi%rfesmnA(i) < 0.0) THEN
     candi%rfesmnA(i) = 0.00 ! 1.0000000E-14
     END IF !
     END DO
     candi%rfesmnB(:)  = candi%param%kfesmn * candi%y(candi%fesy,:)*candi%y(candi%mno2By,:)
     DO i=1, candi%npt
     IF (candi%rfesmnB(i) < 0.0) THEN
     candi%rfesmnB(i) =  0.00 !1.0000000E-14
     END IF !
     END DO
 !                                      IF(simFeS .AND. simMnFe) THEN
     !-- FeS precipitation
     IF(candi%param%rxn_mode==0) THEN
       ! kFeSppt [Fe2+][H2S]
       !RFeSppt(:) = kFeSppt * y(feiiy,:)*y(hsy,:)
               candi%RFeSppt(:) = 0.00
     ELSEIF(candi%param%rxn_mode==1) THEN
               candi%RFeSppt(:) = 0.00
     ELSEIF(candi%param%rxn_mode==2) THEN
            candi%RFeSppt(:) = candi%param%kFeSppt * candi%y(candi%feiiy,:)*candi%y(candi%hsy,:)
     ELSEIF(candi%param%rxn_mode==3) THEN
       WHERE (candi%IAP(candi%FeSy,:)==0.00)
           candi%RFeSppt(:) = candi%param%kFeSppt * candi%y(candi%feiiy,:)*candi%y(candi%hsy,:)
       ELSEWHERE ( (ONE / candi%IAP(candi%FeSy,:)) < ONE )
         ! Precipitation
                 ! kFeSppt (IAP/Ksp -1), where IAP = [Fe2+][H2S]/[H]
         candi%RFeSppt(:) = candi%param%kFeSppt * ( ONE - (ONE / (candi%IAP(candi%FeSy,:))) )
       ELSEWHERE
         ! Dissolution
         candi%RFeSppt(:) = candi%param%kFeSppt * (-ONE) * ( ONE - (candi%IAP(candi%FeSy,:)) )

       END WHERE

     END IF ! If rxn mode = 0, 1, 2 or 3

     !-- FeS transformation to FeS2
     candi%RPyrite(:)   = candi%param%kPyrite * candi%y(candi%fesy,:) *candi%y(candi%hsy,:) !??

     !!-- S0 dispropoertionation
     !rsdispro(:)   = rksdispro * y(s0y,:) * (1. - (y(hsy,:)/H2Sstop))
!                                    !  - simFeS .AND. simMnFe
                                     !  - SUBROUTINE RATES
     IF (candi%param%simX) THEN

!       !-- XS precipitation
!       IF(param%rxn_mode==1) THEN
         ! kxsppt [X][H2S]
         candi%RXSppt(:)   = candi%param%kxsppt * candi%y(candi%xy,:) *candi%y(candi%hsy,:)
!       ELSEIF(param%rxn_mode==2)THEN
!         ! kxsppt (IAP/Ksp -1), where IAP = [X][H2S]/[H]^2
!         WHERE (IAP(:,xsy+1) > ONE)
!           RXSppt(:) = kxsppt * ( IAP(:,xsy+1)-ONE )
!         ELSEWHERE
!           RXSppt(:) = ZERO
!         END WHERE
!       ELSEIF(param%rxn_mode==0) THEN
!               RXSppt(:) = 1.00000E-20
!      END IF !End if rxn mode = 0, 1 or 2

     END IF ! End if SimX
                        !  - SUBROUTINE RATES
   ELSE

     candi%RFeSOX(:)  = 0.0
     candi%RFeS2Ox(:) = 0.0
     candi%RFesFeA(:) = 0.0
     candi%RFesMnA(:) = 0.0
     candi%RFesFeB(:) = 0.0
     candi%RFesMnB(:) = 0.0
     candi%RFeSppt(:) = 0.0
     candi%RPyrite(:) = 0.0
     candi%RXSppt(:)  = 0.0
     !candi%rsdispro(:)= 0.0

   END IF ! End if simFeS .AND. simMnFe
                        !  - SUBROUTINE RATES
   ! CaCO3 (Calcite)
   IF(candi%param%simCaCO3) THEN
    !-- CaCO3 precipitation
     IF(candi%param%rxn_mode==0)THEN
     candi%RCalppt(:) = 0.00
     ELSEIF(candi%param%rxn_mode==1)THEN
     candi%RCalppt(:) = 0.00
     ELSEIF(candi%param%rxn_mode==2) THEN
       ! kCalppt [Ca2+][CO3]
       candi%RCalppt(:) = candi%param%kcalppt * candi%y(candi%cay,:)*candi%y(candi%hco3y,:)
        !       RCalppt(:) = 1.
     ELSEIF(candi%param%rxn_mode==3) THEN
        WHERE (candi%IAP(candi%Caly,:)==0.00)
        candi%RCalppt(:) = candi%param%kcalppt * candi%y(candi%cay,:)*candi%y(candi%hco3y,:)
        ELSEWHERE ( (ONE / candi%IAP(candi%Caly,:)) < ONE )
         ! Precipitation
         candi%RCalppt(:) = candi%param%kCalppt * ( ONE - (ONE / (candi%IAP(candi%Caly,:))) )
       ELSEWHERE
         ! Dissolution
         candi%RCalppt(:) = candi%param%kCalppt * (-ONE) * ( ONE - (candi%IAP(candi%Caly,:)) )
       END WHERE
     END IF ! End if rxn mode = 0, 1 2 or 3
   END IF ! End if simulate CaCO3
                                                   !  - SUBROUTINE RATES
!END DO ! Search keyword "Taipan"
   !-- C. CANDI ---
   !IF(simC12) THEN
   !  rch4ox(:) = kch4ox*y(o2y,:)*(y(ch4y,:)+y(ch4c12y,:))
   !ELSE
   !  rch4ox(:) = kch4ox*y(o2y,:)*y(ch4y,:)
   !END IF
   !IF(simC12) THEN
   !  rch4so4(:) = kch4so4*y(so4y,:)*(y(ch4y,:)+y(ch4c12y,:))
   !ELSE
   !  rch4so4(:) = kch4so4*y(so4y,:)*y(ch4y,:)
   !END IF
   !
   !IF(simFeII) THEN
   !  rfe1ox(:)  = kfe1ox*y(o2y,:)*y(feii1y,:)
   !  rfe2ox(:)  = kfe2ox*y(o2y,:)*y(feii2y,:)
   !  rfe1no3(:) = kfe1no3*y(no3y,:)*y(feii1y,:)
   !  rfe2no3(:) = kfe2no3*y(no3y,:)*y(feii2y,:)
   !END IF
   !
   !rg(:) = kg0*y(cor0y,:) + kg1*y(cor1y,:) + kg2*y(cor2y,:)
   !------------
!Dan adding another thing here                        !  - SUBROUTINE RATES
!Set the rates of organic matter oxidation
IF (candi%param%OMModel == 1) THEN
   candi%rgC(:) = candi%param%poml2dic*candi%y(candi%POMLy,:)*candi%stcoef%fracCPL &
   + candi%param%pomr2dic*candi%y(candi%POMRy,:)*candi%stcoef%fracCPR + candi%param%pomspecial2dic*candi%y(candi%POMspecialy,:)*candi%stcoef%fracCPspecial
   candi%rgN(:) = candi%param%poml2dic*candi%y(candi%POMLy,:)*candi%stcoef%fracNPL &
   + candi%param%pomr2dic*candi%y(candi%POMRy,:)*candi%stcoef%fracNPR + candi%param%pomspecial2dic*candi%y(candi%POMspecialy,:)*candi%stcoef%fracNPspecial
   candi%rgP(:) = candi%param%poml2dic*candi%y(candi%POMLy,:)*candi%stcoef%fracPPL &
   + candi%param%pomr2dic*candi%y(candi%POMRy,:)*candi%stcoef%fracPPR + candi%param%pomspecial2dic*candi%y(candi%POMspecialy,:)*candi%stcoef%fracPPspecial
   !print *, 'rgC', rgC(1:4)
  ! print  *, 'fracCPL', fracCPL
  ! print  *, 'fracCPR', fracCPR
  ! print  *, 'fracCPspecial', fracCPspecial
   !print  *, 'fracCP', fracCPL
   !print  *, 'fracCPL', fracCPL
   !print  *, 'fracCPL', fracCPL

  ELSEIF (candi%param%OMModel == 2) THEN
     candi%rgC(:) = candi%param%docl2dic*candi%y(candi%DOCLy,:)
     candi%rgN(:) = candi%param%donl2din*candi%y(candi%DONLy,:)
     candi%rgP(:) = candi%param%dopl2dip*candi%y(candi%DOPLy,:)
    ELSEIF (candi%param%OMModel == 3) THEN
DO i = 1, candi%npt
     candi%Btot(i)     = candi%y(candi%BAery,i)+candi%y(candi%BDeny,i)+candi%y(candi%BMany,i)+candi%y(candi%BIroy,i)+candi%y(candi%BSuly,i)+candi%y(candi%BMety,i)+candi%y(candi%BFery,i)
     candi%FBHyd       = candi%Btot / ( candi%Btot + candi%param%BMax )
     candi%RPOM1       = candi%param%kHyd1*candi%y(candi%POM1y,i)*candi%FBHyd
     candi%RPOM2       = candi%param%kHyd2*candi%y(candi%POM2y,i)*candi%FBHyd
     candi%RPOM3       = candi%param%kHyd3*candi%y(candi%POM3y,i)*candi%FBHyd
     candi%RPOM4       = candi%param%kHyd4*candi%y(candi%POM4y,i)*candi%FBHyd
     candi%RNecro      = candi%param%kHydN*candi%y(candi%Necromassy,i)
     candi%RPOMspecial = candi%param%pomspecial2dic*candi%y(candi%POMspecialy,i)*candi%FBHyd
     candi%RAerDHyd    = candi%param%kgrowthAer*candi%y(candi%BAery,i)*candi%FDHyd(i)*FTEA_O2
     candi%RDenO2DHyd  = candi%param%kgrowthDen*candi%y(candi%BDeny,i)*candi%FDHyd(i)*FTEA_O2
     candi%RDenNO3DHyd = candi%param%kgrowthDen*candi%y(candi%BDeny,i)*candi%FDHyd(i)*FTEA_NO3          *FIN_O2
     candi%RFerDHyd    = candi%param%kgrowthFer*candi%y(candi%BFery,i)*candi%FDHyd(i)                   *FIN_O2*candi%FTFerDHyd
     candi%RAerOAc     = candi%param%kgrowthAer*candi%y(candi%BAery,i)*candi%FOAc    *FTEA_O2 * candi%FTAerOAc
     candi%RDenOAc     = candi%param%kgrowthDen*candi%y(candi%BDeny,i)*candi%FOAc    *FTEA_NO3* candi%FTDenOAc*FIN_O2
     candi%RDenH2      = candi%param%kgrowthDen*candi%y(candi%BDeny,i)*candi%FH2     *FTEA_NO3* candi%FTDenH2 *FIN_O2
     candi%RManOAc     = candi%param%kgrowthMan*candi%y(candi%BMany,i)*candi%FOAc    *FTEA_MnO2*candi%FTManOAc*FIN_O2*FIN_NO3
     candi%RIroOAc     = candi%param%kgrowthIro*candi%y(candi%BIroy,i)*candi%FOAc    *FTEA_FeOH*candi%FTIroOAc*FIN_O2*FIN_NO3*FIN_MnO2
     candi%RIroH2      = candi%param%kgrowthIro*candi%y(candi%BIroy,i)*candi%FH2     *FTEA_FeOH*candi%FTIroH2 *FIN_O2*FIN_NO3*FIN_MnO2
     candi%RSulOAc     = candi%param%kgrowthSul*candi%y(candi%BSuly,i)*candi%FOAc    *FTEA_SO4 *candi%FTSulOAc*FIN_O2*FIN_NO3*FIN_MnO2*FIN_FeOH
     candi%RSulH2      = candi%param%kgrowthSul*candi%y(candi%BSuly,i)*candi%FH2     *FTEA_SO4 *candi%FTSulH2 *FIN_O2*FIN_NO3*FIN_MnO2*FIN_FeOH
     candi%RMetOAc     = candi%param%kgrowthMet*candi%y(candi%BMety,i)*candi%FOAc              *candi%FTMetOAc*FIN_O2*FIN_NO3*FIN_MnO2*FIN_FeOH*FIN_SO4
     candi%RMetH2      = candi%param%kgrowthMet*candi%y(candi%BMety,i)*candi%FH2               *candi%FTMetH2 *FIN_O2*FIN_NO3*FIN_MnO2*FIN_FeOH*FIN_SO4
     candi%RDHyd       = candi%RAerDHyd(i)+candi%RDenO2DHyd(i)+candi%RDenNO3DHyd(i)+candi%RFerDHyd(i)
     candi%ROAc        = candi%RAerOAc(i) +candi%RDenOAc(i)+candi%RManOAc(i)+candi%RIroOAc(i)+candi%RSulOAc(i)+candi%RMetOAc(i)
     candi%RH2         =            +candi%RDenH2(i)            +candi%RIroH2(i) +candi%RSulH2(i) +candi%RMetH2(i)
     candi%RdeathFer(i)   = candi%param%kdeathFer * candi%y(candi%BFery,i)
     candi%RdeathAer(i)   = candi%param%kdeathAer * candi%y(candi%BAery,i)
     candi%RdeathDen(i)   = candi%param%kdeathDen * candi%y(candi%BDeny,i)
     candi%RdeathMan(i)   = candi%param%kdeathMan * candi%y(candi%BMany,i)
     candi%RdeathIro(i)   = candi%param%kdeathIro * candi%y(candi%BIroy,i)
     candi%RdeathSul(i)   = candi%param%kdeathSul * candi%y(candi%BSuly,i)
     candi%RdeathMet(i)   = candi%param%kdeathMet * candi%y(candi%BMety,i)
     candi%RdeathTot(i)   = candi%RdeathAer(i)+candi%RdeathDen(i)+candi%RdeathMan(i)+candi%RdeathIro(i)+candi%RdeathSul(i)+candi%RdeathMet(i)+candi%RdeathFer(i)
END DO !
     candi%rgC(:) = candi%RDHyd(:)*candi%stcoef%fracCDHyd + candi%ROAc(:)*candi%stcoef%fracCOAc + candi%RH2(:)*candi%stcoef%fracCH2
     candi%rgN(:) = candi%RDHyd(:)*candi%stcoef%fracNDHyd + candi%ROAc(:)*candi%stcoef%fracNOAc + candi%RH2(:)*candi%stcoef%fracNH2
     candi%rgP(:) = candi%RDHyd(:)*candi%stcoef%fracPDHyd + candi%ROAc(:)*candi%stcoef%fracPOAc + candi%RH2(:)*candi%stcoef%fracPH2
     !END IF ! End if OMmodel = 3
  !END IF ! End if OMmodel = 2
END IF ! OMmodel = 1, 2 or 3

DO i=1,candi%npt
IF (candi%rgC(i)<0.00000000000000000000000000000000000000000000) THEN
        candi%rgC(i)=1.0000E-14
END IF
IF (candi%rgN(i)<0.00000000000000000000000000000000000000000000) THEN
        candi%rgN(i)=1.0000E-14
END IF
IF (candi%rgP(i)<0.00000000000000000000000000000000000000000000) THEN
        candi%rgP(i)=1.0000E-14
END IF
END DO

   IF(candi%param%simMnFe) THEN
     candi%rox(:) = candi%fo2(:) + candi%fno3(:) + candi%fmno2(:) + candi%ffeoh(:) + candi%fso4(:) + candi%fmet(:)
     !print*, 'simMnFe on'
!This might not be working. Check OPTIONAL COMPONENTS, line 1500.
   ELSE
     candi%rox(:) = candi%fo2(:) + candi%fno3(:) + candi%fso4(:) + candi%fmet(:)
    print*, 'simMnFe off'
   END IF
!Need to make a few specific reaction rates for the OM Model 3
IF (candi%param%OMModel ==3) THEN
        DO i = 1,candi%npt ! Search "Caterpillar"
        candi%rox(:) = candi%fo2(:) + candi%fno3(:) + candi%fmno2(:) + candi%ffeoh(:) + candi%fso4(:) + candi%fmet(:)
 ! END IF ! End if simMnFe
  END DO ! Search "Caterpillar"
END IF ! End if OMModel == 3
!End: Need to make a few specific reaction rates for the OM Model 3


   RETURN
 END SUBROUTINE RATES
!------------------------------------------------------------------------------!






!------------------------------------------------------------------------------!
!rl   ****************************************************************
!rl   Early diagenese model C. CANDI from R. Luff (1996-2003), based on
!rl   CANDI from  B.P. Boudreau.
!rl   ****************************************************************
!rl   VERSION 3.0 Roger Luff (rluff@gmx.de) GEOMAR, Kiel, Germany
!rl   ****************************************************************
!rl   FILE: DIFCOEF2.f
!rl   AIM : Calculates molecular and ionic diffusion coefficients,
!           D(i), for a set of species listed below at given
!           conditions of salinity, S, temperature, T, and pressure, P
!           Also calculated is the shear viscosity of this solution.
!           Diffusion coefficients are:

!           D(1) = H2O
!           D(2) = O2
!           D(3) = CO2
!           D(4) = NH3
!           D(5) = H2S
!           D(6) = H3PO4
!           D(7) = B(OH)3
!           D(8) = HCO3-
!           D(9) = CO3=
!           D(10) = NH4+
!           D(11) = HS-
!           D(12) = NO3-
!           D(13) = H2PO4-
!           D(14) = HPO4=
!           D(15) = PO4(---)
!           D(16) = B(OH)4-
!           D(17) = H+
!           D(18) = OH-
!           D(19) = Ca++
!           D(20) = Mg++
!           D(21) = Fe++
!           D(22) = Mn++
!           D(23) = SO4=
!           D(24) = H4SiO4
!           D(25) = CH4
!           D(26) = Na+
!           D(27) = Cl-
!           D(28) = Br-
!
!           Note: 1) enter S in ppt, T in deg. C and P in atm.
!                 2) diffusion coefficients are in units of cm**2/s
!                 3) H2O viscosity is in unit of centipoise
!------------------------------------------------------------------------------!
 SUBROUTINE difcoef2(v,d,s,t,p)
   !-- Incoming
   INTEGER, PARAMETER :: n=28
   REAL(SEDP), INTENT(INOUT)          :: v
   REAL(SEDP), INTENT(IN)             :: s
   REAL(SEDP), INTENT(IN)             :: t
   REAL(SEDP), INTENT(IN)             :: p
   !-- Outgoing
   REAL(SEDP), INTENT(OUT)            :: d(n)
   !-- Local
   REAL(SEDP) :: tk,ts,r25,rho,v0,a,b,t0,phi,mw,vm,fac
   REAL(SEDP) :: ss,v25,vtk
                                                                  ! SUBROUTINE DIFFUSION COEFFICIENTS


   tk = t + 273.15

   !  Calculate density of pure water at 25 deg C and sample temperature

   ts = 25.0
   r25 = rho_sw(ts)
   rho = rho_sw(t)

   !  Calculate the viscosity for the true sample conditions.

   v = visco(s,t,p)

   !  Start calculations of diffusion coefficients in pure water
   !  at sample temperature.

   v0 = visco(zero,t,one)

   !  Water : from Cohen and Turnbull (1959) and Krynicki et al. (1978)

   a = 12.5D-09*EXP(-5.22D-04*p)
   b = 925.0*EXP(-2.6D-04*p)
   t0 = 95.0 + 2.61D-02*p
   d(1) = a*SQRT(tk)*EXP(-b/(tk-t0))*1.0D+04

                                                                  ! SUBROUTINE DIFFUSION COEFFICIENTS
   !  Dissolved gases : from Wilke and Chang (1955)
   !    note: 1) MW = molecular weight of water
   !          2) VM = molar volumes (cm**3/mol) (Sherwood et al., 1975)

   !  The factor PHI is reduced to 2.26 as suggested by Hayduk and
   !  Laudie (1974).


   phi = 2.26
   mw = 18.0
   a = SQRT(phi*mw)*tk/v0

   !  Oxygen

   vm = 25.6
   d(2) = 7.4D-08*(a/(vm**0.6))

   !  CO2

   vm = 34.0
   d(3) = 7.4D-08*(a/(vm**0.6))

   !  NH3

   vm = 25.8
   d(4) = 7.4D-08*(a/(vm**0.6))

   !  H2S

   vm = 32.9
   d(5) = 7.4D-08*(a/(vm**0.6))

   !  CH4

   vm = 37.7
   d(25) = 7.4D-08*(a/(vm**0.6))

                                                                  ! SUBROUTINE DIFFUSION COEFFICIENTS
   !  The coefficients in pure water for the following species are
   !  calculated by the linear functions of temperature (deg C)
   !  found in Boudreau (in press).

   !  i.e. NO3-,HS-,H2PO4-,CO3=,SO4=,Ca++,Mg++,Mn++,Fe++,NH4+,H+ & OH-,
   !         HCO3-, HPO4=, PO4(3-)


   fac = 1.0D-06
   d(8) = (5.06 + 0.275*t)*fac
   d(9) = (4.33 + 0.199*t)*fac
   d(10) = (9.5 + 0.413*t)*fac
   d(11) = (10.4 + 0.273*t)*fac
   d(12) = (9.50 + 0.388*t)*fac
   d(13) = (4.02 + 0.223*t)*fac
   d(14) = (3.26 + 0.177*t)*fac
   d(15) = (2.62 + 0.143*t)*fac
   d(17) = (54.4 + 1.555*t)*fac
   d(18) = (25.9 + 1.094*t)*fac
   d(19) = (3.60 + 0.179*t)*fac
   d(20) = (3.43 + 0.144*t)*fac
   d(21) = (3.31 + 0.150*t)*fac
   d(22) = (3.18 + 0.155*t)*fac
   d(23) = (4.88 + 0.232*t)*fac
   d(26) = (6.06 + 0.297*t)*fac
   d(27) = (9.60 + 0.438*t)*fac
   d(28) = (10.0 + 0.441*t)*fac

   ! H3PO4 : Least (1984) determined D(H3PO4) at 25 deg C and 0 ppt S.
   !         Assume that this value can be scaled by the Stokes-Einstein
   !         relationship to any other temperature.

   d(6) = 0.87D-05
   ts = 25.0
   ss = 0.0
   v25 = visco(ss,ts,one)
   vtk = visco(ss,t,one)
   d(6) = d(6)*v25/298.15*tk/vtk

   !  B(OH)3 : Mackin (1986) determined D(B(OH)3) at 25 deg C and
   !           about 29.2 ppt S.
   !           Assume that this value can be scaled by the Stokes-Einstein
   !           relationship to any other temperature.

   d(7) = 1.12D-05
   ts = 25.0
   ss = 29.2
   v25 = visco(ss,ts,one)
   vtk = visco(zero,t,one)
   d(7) = d(7)*v25/298.15*tk/vtk

   !  B(OH)4 : No information on this species whatsoever! Boudreau and
   !           Canfield (1988) assume it is 12.5% smaller than B(OH)3.

   d(16) = 0.875*d(7)

   !  H4SiO4 : Wollast and Garrels (1971) found D(H4SiO4) at 25 deg C
   !           and 36.1 ppt S.
   !           Assume that this value can be scaled by the Stokes-Einstein
   !           relationship to any other temperature.

   d(24) = 1.0D-05
   ts = 25.0
   ss = 36.1
   v25 = visco(ss,ts,one)
   vtk = visco(zero,t,one)
   d(24) = d(24)*v25/298.15*tk/vtk

   !  To correct for salinity, the Stokes-Einstein relationship is used.
   !  This is not quite accurate, but is at least consistent.

   fac = v0/v
   d(:) = d(:)*fac
                                                                 ! SUBROUTINE DIFFUSION COEFFICIENTS
   RETURN

   !---------------------------------------------------------------------------!
 CONTAINS

   !---------------------------------------------------------------------------!
   !rl   Calculates the density of pure water (Sweet Water) using the
   !     equation given by Bigg (1967). (See also Millero and
   !     Poisson, 1981)
   !rl   ****************************************************************
   !---------------------------------------------------------------------------!
   REAL(SEDP) FUNCTION rho_sw(t)

   REAL(SEDP), INTENT(IN)             :: t

   rho_sw =  (999.842594 + t*(6.793952D-02 + t*(-9.09529D-03  &
       + t*(1.001685D-04 + t*(-1.120083D-06 + t*6.536336D-09)))))/1000.0
   RETURN
   END FUNCTION rho_sw
   !---------------------------------------------------------------------------!
   !rl   ****************************************************************
   !rl   AIM : Calculates the shear viscosity of water using the equation
   !           given by Kukulka et al. (1987).
   !           Calculated viscosity is in centipoise.
   !           Valid for 0<T<30 and 0<S<36.
   !rl   ****************************************************************
   !---------------------------------------------------------------------------!
   REAL(SEDP) FUNCTION visco(s,t,p)

   REAL(SEDP), INTENT(IN)             :: s
   REAL(SEDP), INTENT(IN)             :: t
   REAL(SEDP), INTENT(IN)             :: p

   visco =  1.7910 - t*(6.144D-02 - t*(1.4510D-03 - t*1.6826D-05))  &
          - 1.5290D-04*p + 8.3885D-08*p*p + 2.4727D-03*s  &
          + (6.0574D-06*p - 2.6760D-09*p*p)*t + (t*(4.8429D-05  &
          - t*(4.7172D-06 - t*7.5986D-08)))*s
   RETURN
   END FUNCTION visco
   !---------------------------------------------------------------------------!

 END SUBROUTINE difcoef2
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!     GET MOLECULAR DIFFUSION COEFFICIENTS FOR SOLUTES FROM THE
!     PROGRAM DIFCOEF2.f.  USE A 50/50 WEIGHTING SCHEME FOR THE
!     AMMONIA AND SULFIDE SPECIES TO GET DTHN4 AND DTH2S.
!     USE AN 85% HPO4, 14% PO4 AND 1% H2PO4 DISTRIBUTION TO
!     GET THE DTPO4 VALUE.  DTCO2 IS SET TO THE VALUE FOR HCO3-.
!     Diffusion coefficients are calculated in units of cm^2/yr
!     ("YEAR" is num of secs /yr, so this assumes df = cm2/s) !##MATT
!------------------------------------------------------------------------------!
 SUBROUTINE setdifcoef(candi,df)
   !-- Incoming
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   REAL(SEDP), INTENT(IN) :: df(:)

!  DIFFC(:,:)     = 0.0
!  DIFFC(:,:)   = df(27)*year   !  MATT: Set default to Cl df value
   candi%DIFFC        = df(27)*year   !  MATT: Set default to Cl df value

   candi%DIFFC(candi%o2y,:)   = df(2)*year  !O2
   candi%DIFFC(candi%no3y,:)  = df(12)*year !NO3
   candi%DIFFC(candi%so4y,:)  = df(23)*year !SO4
   candi%DIFFC(candi%po4ly,:) = (0.85D+00*df(14)+0.14D+00*df(15)+0.01D+00*df(13))*year !TPO4
   candi%DIFFC(candi%nh4y,:)  = ((0.5*df(4)+0.5*df(10)))*year !TNH4
   IF(candi%param%simMnFe) THEN
     candi%DIFFC(candi%mniiy,:)=df(22)*year !MnII
     candi%DIFFC(candi%feiiy,:)=df(21)*year !FeII
   END IF
   candi%DIFFC(candi%ch4y,:)  = df(25)*year    !CH4

   candi%DIFFC(candi%hsy,:)   = (0.5*df(5)+0.5*df(11))*year !TH2S
   candi%DIFFC(candi%hco3y,:) = df(8)*year                  !TCO2


   RETURN
 END SUBROUTINE setdifcoef
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!rl   ****************************************************************
!rl   Early diagenese model C. CANDI from R. Luff (1996-2003), based on
!rl   CANDI from  B.P. Boudreau.
!rl   ****************************************************************
!rl   VERSION 3.0 Roger Luff (rluff@gmx.de) GEOMAR, Kiel, Germany
!rl   ****************************************************************
!rl   FILE: SED.f
!rl   AIM : A subroutine for the sediment properties.
!           Calculates the value of the depth-dependent porosity (P),
!           derivative of the porosity with depth (DPDX), the
!           value of the bioturbation coefficient and the value of the
!           bioirrigation coefficient.
!           Poros= porosity at depth X
!           DPDX = spatial derivative of porosity at depth X
!rl   ****************************************************************
!------------------------------------------------------------------------------!
 SUBROUTINE SedProperties(candi)
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   !-- Local
   REAL(SEDP) :: pe
   INTEGER  :: i

   !---------------------------------------------------------------------------!
   !-- Set Bioturbation Profile and depth derivatives of bioturb(X)
   IF (candi%param%imix == 0) THEN
   ! Gaussian decrease (Christensen, 1982) JGR
     candi%bioturb(:)=candi%param%db0*EXP(-(candi%rpar(:)*candi%rpar(:))/(2.0*candi%param%xs*candi%param%xs))
   END IF

   IF (candi%param%imix == 1) THEN
     WHERE (candi%rpar <= candi%param%x1)
       candi%bioturb = candi%param%db0
     ENDWHERE
     WHERE (candi%rpar > candi%param%x1.AND.candi%rpar <= candi%param%x2)
       candi%bioturb = candi%param%db0*(candi%param%x2-candi%rpar)/(candi%param%x2-candi%param%x1)
     ENDWHERE
     WHERE (candi%rpar > candi%param%x2)
       candi%bioturb = 0.0
     ENDWHERE
   END IF
   IF (candi%param%imix == 2) THEN
     candi%bioturb(:)=candi%param%db0*EXP(-candi%rpar(:)/candi%param%xs)
   END IF

   !     in case of Gaussian decrease and exponential decreases set
   !     the bioturbation coefficient to zero in the last three
   !     layers to make sure that no organic matter will be transported
   !     out of the sediment column by bioturbation. This is at least
   !     important to check the mass balance for C.
   IF (candi%param%imix == 0.OR.candi%param%imix == 2) THEN
     candi%param%x1=candi%rpar(candi%npt-6)
     candi%param%x2=candi%rpar(candi%npt-3)
     DO i=candi%npt-5,candi%npt
       IF(candi%rpar(i) <= candi%param%x2) THEN
         candi%bioturb(i) =  candi%bioturb(candi%npt-6)*(candi%param%x2-candi%rpar(i))/(candi%param%x2-candi%param%x1)
       ELSE
         candi%bioturb(i) = 0.0
       END IF
     END DO
   END IF

   !---------------------------------------------------------------------------!
   !-- Set depth derivatives of Porosity(X),
   IF (candi%param%bp == 0.0.OR.(candi%param%p0-candi%param%p00) == 0.0) THEN
     candi%poros(:)= candi%param%p0
   ELSE
     candi%poros(:)= (candi%param%p0-candi%param%p00)*EXP(-candi%param%bp*candi%rpar(:)) + candi%param%p00
   END IF
   !#ifdef porosvar
   !do i=1,npt
   !  poros_bg(i) = poros(i)
   !enddo
   !
   !#endif

   !---------------------------------------------------------------------------!
   !-- Set and dependent process variables of porosity like w00,
   !-- u00, t2 and thier derivates
   CALL porosi(candi)

   !---------------------------------------------------------------------------!
   !-- Calculate the bioirrigation coefficient
   candi%cirrig(:) = candi%param%alpha0
   WHERE (candi%rpar > candi%param%xirrig) candi%cirrig=candi%param%alpha0*EXP(-1.5*(candi%rpar-candi%param%xirrig))
   !-- ensure not irrigation in the last layer
   candi%cirrig(candi%npt) = 0.0

   RETURN
 END SUBROUTINE SedProperties
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!rl   ****************************************************************
!rl   Early diagenese model C. CANDI from R. Luff (1996-2003), based on
!rl   CANDI from  B.P. Boudreau.
!rl   ****************************************************************
!rl   VERSION 3.0 Roger Luff (rluff@gmx.de) GEOMAR, Kiel, Germany
!rl   ****************************************************************
!rl   FILE: POROSI.f
!rl   AIM : Calculate process variables depending on porosity like
!rl         w00, u00, t2 and thier derivates
!rl         Porosity init will be done in SED, if Poros is a function
!rl         of time, the setting will be done in candi.f
!rl   ****************************************************************
!------------------------------------------------------------------------------!
 SUBROUTINE porosi(candi)
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   INTEGER  :: i
   REAL(SEDP) :: ventflow1

   !rl   add diffusive sublayer
   !     if (ibc1.eq.2) then
   !       do i=1,npt
   !         if (rpar(i) .le. 0.1) then
   !           subdepth = i
   !           poros(i) = 0.95
   !         endif
   !       enddo
   !     else
   !       subdepth = 1
   !     endif

   !rl   alternative numerical approximation of dpdx
   DO i=2,candi%npt-1
     candi%dpdx(i) = (-candi%dh(i+1)*candi%poros(i-1)/(candi%dh(i)*(candi%dh(i)+candi%dh(i+1)))&
           +(candi%dh(i+1)-candi%dh(i))*candi%poros(i)/(candi%dh(i)*candi%dh(i+1))&
           +candi%dh(i)*candi%poros(i+1)/(candi%dh(i+1)*(candi%dh(i)+candi%dh(i+1))))
     candi%dbdx(i) = (-candi%dh(i+1)*candi%bioturb(i-1)/(candi%dh(i)*(candi%dh(i)+candi%dh(i+1)))&
           +(candi%dh(i+1)-candi%dh(i))*candi%bioturb(i)/(candi%dh(i)*candi%dh(i+1))&
           +candi%dh(i)*candi%bioturb(i+1)/(candi%dh(i+1)*(candi%dh(i)+candi%dh(i+1))))
   END DO
   candi%dpdx(1)=(candi%poros(2)-candi%poros(1))/candi%dh(1)
   candi%dpdx(candi%npt)=(candi%poros(candi%npt)-candi%poros(candi%npt-1))/candi%dh(candi%npt-1)
   candi%dbdx(1)=(candi%bioturb(2)-candi%bioturb(1))/candi%dh(1)
   candi%dbdx(candi%npt)=(candi%bioturb(candi%npt)-candi%bioturb(candi%npt-1))/candi%dh(candi%npt-1)
   !rl   use the smallest porosity value to calculate the size of ventflow
   candi%porosmin = MINVAL(candi%poros)

   IF (candi%porosmin > 0.8) THEN
     WRITE (*,*) "Attention ventflow calculation may be wrong !"
     WRITE (*,*) "The minimal value for porosity should not be greater than 0.8 "
     WRITE (*,*) " - but it is: ", candi%porosmin
   END IF
   !rl   Ventflow at the basis of the core
   ventflow1 = candi%param%ventflow * (candi%porosmin**3)/((1.0-candi%porosmin)**2)*((1.0-candi%param%p00)**2)/(candi%param%p00**3)
   candi%ps(:)  = one - candi%poros(:)
   candi%psp(:) = candi%ps(:)/candi%poros(:)
   candi%pps(:) = candi%poros(:)/candi%ps(:)
   !       Burial velocity of solids
   candi%wvel(:) = candi%param%w00*(one - candi%param%p00)/(candi%ps(:))
   !       Porewater velocity du to burial of sediment in case of const.
   !       porosity it has the same direction and amount of wvel if the
   !       porosity is not const. it is always smaller the Wvel but still
   !       same direction.
   candi%uvel(:) = (candi%param%w00-ventflow1)*candi%param%p00/candi%poros(:)
   DO i=2,candi%npt-1
     candi%dudx(i) = (-candi%dh(i+1)*candi%uvel(i-1)/(candi%dh(i)*(candi%dh(i)+candi%dh(i+1)))&
           +(candi%dh(i+1)-candi%dh(i))*candi%uvel(i)/(candi%dh(i)*candi%dh(i+1))&
           +candi%dh(i)*candi%uvel(i+1)/(candi%dh(i+1)*(candi%dh(i)+candi%dh(i+1))))
   END DO
   candi%dudx(candi%npt)=(candi%uvel(candi%npt)-candi%uvel(candi%npt-1))/(candi%dh(candi%npt-1))
   !     TORT2  Calculates the tortuosity at a given depth and time.

   !     T2 = the square of the tortuosity at depth X (as it appears
   !          in Berners(1980) version of the diagenetic equations.
   !     DT2DX = spatial derivative of T2 at depth X.
   !     P     = porosity

   !rl   tortuosity Eq. after Boudreau 1996 "The diffusive tortuosity of fine..."

   !rl   Eq. 6 Archies Law
   IF (candi%param%torteq == 1) THEN
     candi%t2(:) = candi%poros(:)**(1-candi%param%an)
   END IF
   !rl   Eq. 7 Burger Friecke
   IF (candi%param%torteq == 2) THEN
     candi%t2(:) = candi%poros(:) + candi%param%aa*(one-candi%poros(:))
   END IF
   !rl   Eq. 9 Weissberg
   IF (candi%param%torteq == 3) THEN
     candi%t2(:) = 1.0 - candi%param%ab*LOG(candi%poros(:))
   END IF
   !rl
   candi%t2(:) = 1.0-2.0*LOG(candi%poros(:))

   ! DT2DX = spatial derivative of T2 at depth X.
   DO i=2,candi%npt-1
     candi%dt2dx(i) = (-candi%dh(i+1)*candi%t2(i-1)/(candi%dh(i)*(candi%dh(i)+candi%dh(i+1)))&
           +(candi%dh(i+1)-candi%dh(i))*candi%t2(i)/(candi%dh(i)*candi%dh(i+1))&
           +candi%dh(i)*candi%t2(i+1)/(candi%dh(i+1)*(candi%dh(i)+candi%dh(i+1))))
   END DO
   candi%dt2dx(1)=(candi%t2(2)-candi%t2(1))/candi%dh(1)
   candi%dt2dx(candi%npt)=(candi%t2(candi%npt)-candi%t2(candi%npt-1))/candi%dh(candi%npt)

   candi%dff(:) = (candi%t2(:)*candi%dpdx(:)/candi%poros(:)-candi%dt2dx(:))/(candi%t2(:)**2)

   RETURN
 END SUBROUTINE porosi
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!rl   AIM : Calculates the fluxes between sediment layers and between
!rl         sediment and bottom water.
!rl   ****************************************************************
!rl   calculate the fluxes of speci between two layers (i=2: flux from
!rl   sediment to bottom water
!------------------------------------------------------------------------------!
 REAL(SEDP) FUNCTION dfluxes(candi,layer,speci)
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   INTEGER, INTENT(IN)                      :: layer
   INTEGER, INTENT(IN)                      :: speci
   INTEGER  :: i
   REAL(SEDP) :: surfconz

   IF (layer == 2) THEN
     surfconz = candi%top_bound(speci)
   ELSE
     surfconz = candi%y(speci,layer-1)
   ENDIF

   !--   calculate the flux out of the sediment (-) or into it (+)
   !--   see BB (1997) p. 324, 360 for forth order accurate, flux formula
   !--   checked by Boudreau, pers comm.

   candi%fluxirr = 0.0
   IF (layer < candi%npt-2) THEN
     !rl calculate the flux from irrigartion
     IF (candi%param%irrg == 1) THEN
       candi%fluxirr = (firrig(candi,speci,layer))/2.0*candi%dh(layer)*candi%poros(layer)
       DO i=layer+1,candi%npt-1
         candi%fluxirr = candi%fluxirr + (firrig(candi,speci,i))*candi%dh(i)*candi%poros(i)
       END DO ! End do i
     END IF !End if irrig = 1
     candi%fluxadv = candi%poros(layer)*candi%uvel(layer)*candi%y(speci,layer)
     candi%fluxdiff = -candi%DIFFC(speci,layer)/candi%t2(layer)*candi%poros(layer) *(candi%y(speci,layer+1)-surfconz)/(2.0*candi%dh(layer))
                  !cm2/y          /1        *1              *mmol/L                      / cm

   ELSE
     candi%fluxadv = -candi%poros(candi%npt-1)*candi%uvel(candi%npt-1)*candi%y(speci,candi%npt-1)
     candi%fluxdiff = candi%DIFFC(speci,candi%npt-1)/candi%t2(candi%npt-1)*candi%poros(candi%npt-1)*(candi%y(speci,candi%npt)-candi%y(speci,candi%npt-2))/(2.0*candi%dh(candi%npt-1))
   END IF
   dfluxes = candi%fluxadv + candi%fluxdiff + candi%fluxirr

   IF (layer == 2) THEN
     dfluxes = dfluxes * candi%BC%flux_scale
   END IF
   RETURN
 END FUNCTION dfluxes
!------------------------------------------------------------------------------!






!------------------------------------------------------------------------------!
!--   2D to 1D to solve y / ydot with vode (was yaufyfex)
!------------------------------------------------------------------------------!
 SUBROUTINE Copy2Dinto1D(neq,yfex,ydotfex,nSPECIES,y,ydot)
   !-- Incoming
   INTEGER,    INTENT(IN)  :: neq
   REAL(SEDP), INTENT(OUT) :: yfex(neq)
   REAL(SEDP), INTENT(OUT) :: ydotfex(neq)
   INTEGER,    INTENT(IN)  :: nSPECIES
   REAL(SEDP), INTENT(IN) :: y(:,:),ydot(:,:)
   !-- Local
   INTEGER  :: ii,i,j,ii1


   DO ii = 1,neq
     ii1 = ii-1
     i = MOD(ii1,nSPECIES)
     j = ((ii1)/nSPECIES)+1
     yfex(ii) = y(i,j)
     ydotfex(ii) = ydot(i,j)
   END DO

   RETURN
 END SUBROUTINE Copy2Dinto1D
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!--   1D to 2D to solve y / ydot with vode
!------------------------------------------------------------------------------!
 SUBROUTINE Copy1Dinto2D(neq,yfex,nSPECIES,y)
   !-- Incoming
   INTEGER,    INTENT(IN)  :: neq
   REAL(SEDP), INTENT(IN)  :: yfex(neq)
   INTEGER,    INTENT(IN)  :: nSPECIES
   REAL(SEDP), INTENT(OUT) :: y(:,:)
   !-- Local
   INTEGER  :: ii,i,j,ii1

   DO ii = 1,neq
     ii1 = ii-1
     i = MOD(ii1,nSPECIES)
     j = ((ii1)/nSPECIES)+1
     y(i,j) = yfex(ii)
   END DO

   RETURN
END SUBROUTINE Copy1Dinto2D
!------------------------------------------------------------------------------!







!------------------------------------------------------------------------------!
! SetGrid()
!
! Subroutine to grid the sediment given the total sediment depth
! The grid may be non-uniform
!------------------------------------------------------------------------------!
 SUBROUTINE SetGrid(candi,SedDepth)
   !-- Incoming
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   REAL(SEDP), INTENT(IN)  :: SedDepth
   !-- Local
   INTEGER  :: i
   REAL(SEDP)              :: s, dx


   IF (candi%job == 0) THEN
     !-- Equal vertical grid
     WRITE(*,"(8X,'CANDI being set with an equal vertical grid',/)")
     s = SedDepth/REAL((candi%npt-1),SEDP)
     DO i=1,candi%npt
        candi%dh(i)   = s
     END DO
   ELSE
     !-- Exponential increase of grid size with depth
     WRITE(*,"(8X,'CANDI set with an exponential increase of grid size with depth',/)")
     s  = 0.0
     dx = 2 !10.0
     DO i = 1,candi%npt
       s = s+1.0/(1.0+EXP((REAL(i,SEDP)-REAL(candi%npt,SEDP)/2.0)/dx))
     END DO

     s = SedDepth/(REAL(candi%npt,SEDP)-s)
     DO i=1,candi%npt
       candi%dh(i)=s-s/(1.0+exp((REAL(i,SEDP)-REAL(candi%npt,SEDP)/2.0)/dx))
     END DO
     candi%dh(2) = candi%dh(1)

   END IF

   !-- Manual grid definition:
   !--     change here the vertical spacing for staggered grid
   !DO i=1,npt
   !  IF (i <= 20)            dh(i)   = 0.05
   !  IF (i > 20.AND.i <= 40) dh(i)   = 0.10
   !  IF (i > 40.AND.i <= 60) dh(i)   = 0.15
   !  IF (i > 60.AND.i <= 80) dh(i)   = 0.20
   !  IF (i > 80)             dh(i)   = 0.25
   !END DO

   !-- Calculate necessary depth variables:
   !--   dh = distance between two layers
   !--   dh2= dh*dh for the central differences
   !--   dhf= dh for the foreward differences
   !--   dhr= dh for the backward differences
   candi%dh2(1)=0.25*(candi%dh(1)+candi%dh(1))*(candi%dh(1)+candi%dh(2))
   candi%dhf(1)=0.5 *(candi%dh(2)+candi%dh(1))
   candi%dhr(1)=candi%dh(1)
   DO i=2,candi%npt-1
     candi%dh2(i)=0.25*(candi%dh(i-1)+candi%dh(i))*(candi%dh(i)+candi%dh(i+1))
     candi%dhf(i)=0.5 *(candi%dh(i+1)+candi%dh(i))
     candi%dhr(i)=0.5 *(candi%dh(i-1)+candi%dh(i))
   END DO
   candi%dh2(candi%npt)=0.25*(candi%dh(candi%npt-1)+candi%dh(candi%npt))*(candi%dh(candi%npt)+candi%dh(candi%npt))
   candi%dhf(candi%npt)=candi%dh(candi%npt)
   candi%dhr(candi%npt)=0.5*(candi%dh(candi%npt-1)+candi%dh(candi%npt))
   !--  2 times dh, only for the first layer = (candi%dh/2+candi%dh2/2)*2
   candi%tdh = candi%dh(1)+candi%dh(2)
   candi%rpar(1)=0.0
   DO i=2,candi%npt
     candi%rpar(i)=candi%rpar(i-1)+candi%dh(i)
   END DO


 END SUBROUTINE SetGrid
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
SUBROUTINE myerror(candi,errorcode, layer, text)
!rl   ****************************************************************
!rl   Early diagenese model C. CANDI from R. Luff (1996-2003), based on
!rl   CANDI from  B.P. Boudreau.
!rl   ****************************************************************
!rl   VERSION 3.0 Roger Luff (rluff@gmx.de) GEOMAR, Kiel, Germany
!rl   ****************************************************************
!rl   FILE: MYERROR.f
!rl   AIM : Errorhandling routine. Print out the errror number and an
!rl         explanation and stop the simulation. If the rouitine is
!rl         called, a major error occured, no output will be written
!rl         into a file.
!rl   ****************************************************************
!USE CONSTANTS
!USE PHYSICS
!USE STEUERPARA

!IMPLICIT NONE

   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
INTEGER, INTENT(IN)                      :: errorcode
INTEGER, INTENT(IN)                      :: layer
CHARACTER (LEN=*), INTENT(IN)       :: text


IF (errorcode < 0) THEN
  PRINT *, 'ERROR occured during execution, number:',errorcode
  PRINT *, 'Program aborded, reason : '
  IF (errorcode == -1) THEN
    PRINT *, 'VODE excess work done on this CALL. Perhaps '
    PRINT *, 'wrong MF.'
  END IF
  IF (errorcode == -2) THEN
    PRINT *, 'VODE excess accuracy requested.'
    PRINT *, ' Tolerances too small.'
  END IF
  IF (errorcode == -3) THEN
    PRINT *, 'VODE illegal input detected. '
    PRINT *, '(See printed message.)'
  END IF
  IF (errorcode == -4) THEN
    PRINT *, 'VODE repeated error test failures.'
    PRINT *, ' (Check all input.)'
  END IF
  IF (errorcode == -5) THEN
    PRINT *, 'VODE repeated convergence failures.'
    PRINT *, '(Perhaps bad Jacobian supplied or wrong choice '
    PRINT *, 'of MF or tolerances.)'
  END IF
  IF (errorcode == -6) THEN
    PRINT *, 'VODE error weight became zero during problem.'
    PRINT *, '(Solution component i vanished, and ATOL or '
    PRINT *, 'ATOL(i) = 0.0'
  END IF
  IF (errorcode == -10) THEN
    PRINT *, 'LOG error while calculating pH'
    PRINT *, 'Argument le 0.0 '
  END IF
  IF (errorcode == -11) THEN
    PRINT *, 'Incorrect setup file, number of layers are not'
    PRINT *, 'the same in STEUER.DAT and AUFSETZ.DAT '
  END IF
  IF (errorcode == -12) THEN
    PRINT *, 'Incorrect program version, number of species are'
    PRINT *, 'not the same in this compitation of C. CANDI and '
    PRINT *, 'AUFSETZ.DAT'
  END IF
  IF (errorcode == -13) THEN
    PRINT *, 'Incorrect program version, type of equilibrium'
    PRINT *, 'calculation are not the same in this compitation'
    PRINT *, ' of C. CANDI and AUFSETZ.DAT '
  END IF
  IF (errorcode == -14) THEN
    PRINT *, 'Incorrect makefile setup, no one or more than one'
    PRINT *, 'equilibrium-module has been choosen.'
    PRINT *, 'You need to compile the model again with the '
    PRINT *, 'correct settings.'
  END IF
  IF (errorcode == -15) THEN
    PRINT *, 'CONFIGURATION ERROR'
    PRINT *, 'Modul C12 can only run in combination with the'
    PRINT *, 'advancement module.'
    PRINT *, 'You need to compile the model again with the '
    PRINT *, 'correct settings.'
  END IF
  IF (errorcode == -16) THEN
    PRINT *, 'CONFIGURATION ERROR'
    PRINT *, 'ARRAY for fluxdata to small, increase MAXCORG'
    PRINT *, 'ERROR OCCURED IN SUBROUTINE STEUER'
    PRINT *, candi%corginp, ' > ' , maxcorg
  END IF
  IF (errorcode == -17) THEN
    PRINT *, 'CONFIGURATION ERROR'
    PRINT *, 'Too much layers. Reduce the number of points to'
    PRINT *, candi%maxnpts,' or increase maxnpt in candi.inc. In '
    PRINT *, 'latter case latter do not forget to delete all'
    PRINT *, '*.o files and  compile candi again.'
  END IF
  IF (errorcode == -18) THEN
    PRINT *, 'CONFIGURATION ERROR'
    PRINT *, 'Length of area with double precision (XLdouble) >='
    PRINT *, 'length of total simulation (xl):'
    PRINT *, candi%xldouble, ' >= ' , candi%param%xl
  END IF

  IF (errorcode == -20) THEN
    PRINT *, 'U-Velocity < 0, this may cause to problems with'
    PRINT *, 'the upstream advection equation in FEX.'
  END IF
  IF (errorcode == -21) THEN
    PRINT *, 'ERROR READING DATAFILE: AUFSETZ.DAT'
    PRINT *, 'ERROR OCCURED IN SUBROUTINE READAUF'
  END IF
!rl      some Newton Raphson ERROR messages
  IF (errorcode == -30) THEN
    PRINT *, 'N.R. ERROR ',text,': To much iterations'
  END IF
  IF (errorcode == -31) THEN
    PRINT *, 'N.R. ERROR CONCENTRATION of ',TRIM(text),' < 0'
  END IF
!rl      some equilibrium factor ERROR messages
  IF (errorcode == -35) THEN
    PRINT *, 'ERROR IN EQFAC: SUM OF ',TRIM(text),' FACTOR NOT 1'
  END IF

!rl      some Module ERROR messages
  IF (errorcode == -41) THEN
    PRINT *, 'Error reading old initial simulation candi:'
    PRINT *, 'Modul FeII was not activated in the setup file'
  END IF
  IF (errorcode == -42) THEN
    PRINT *, 'Error reading old initial simulation candi:'
    PRINT *, 'Modul Tracer was not activated in the setup file'
  END IF
  IF (errorcode == -43) THEN
    PRINT *, 'Error reading old initial simulation candi:'
    PRINT *, 'Modul FeS was not activated in the setup file'
  END IF
  IF (errorcode == -44) THEN
    PRINT *, 'Error reading old initial simulation candi:'
    PRINT *, 'Modul CaXCO3 was not activated in the setup file'
  END IF
  IF (errorcode == -45) THEN
    PRINT *, 'Error reading old initial simulation candi:'
    PRINT *, 'Modul CaCO3 was not activated in the setup file'
  END IF
  IF (errorcode == -46) THEN
    PRINT *, 'Error reading old initial simulation candi:'
    PRINT *, 'Modul C12 was not activated in the setup file'
  END IF
  IF (errorcode == -47) THEN
    PRINT *, 'Error reading old initial simulation candi:'
    PRINT *, 'Modul CaCO3C12 was not activated in the setup file'
  END IF
  IF (errorcode == -48) THEN
    PRINT *, 'Error reading old initial simulation candi:'
    PRINT *, 'Modul MNFE was not activated in the setup file'
  END IF
  IF (errorcode == -49) THEN
    PRINT *, 'Error reading old initial simulation candi:'
    PRINT *, 'Modul PO4Solid was not activated in the setup file'
  END IF
  IF (layer /= 0) THEN
    PRINT *, 'Error in Layer: ', layer
  END IF
  stop
END IF
RETURN
END SUBROUTINE myerror
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
!rl   ****************************************************************
!rl   Early diagenese model C. CANDI from R. Luff (1996-2003), based on
!rl   CANDI from  B.P. Boudreau.
!rl   ****************************************************************
!rl   VERSION 3.0 Roger Luff (rluff@gmx.de) GEOMAR, Kiel, Germany
!rl   ****************************************************************
!rl   FILE: MASS.f
!rl   AIM : Prints out the fluxes at both (upper - lower) boundaries to
!rl         calculate a mass ballance in case of a steady state and
!rl         calculates vertical integrated rates.
!rl   ****************************************************************
!------------------------------------------------------------------------------!
 SUBROUTINE mass(candi,steadiness,fid)
   !-- Incoming
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   REAL(SEDP), INTENT(IN)         :: steadiness
   INTEGER,    INTENT(IN)         :: fid
   !-- Local
   INTEGER  :: i,j
   REAL(SEDP) :: cflux(7)
   REAL(SEDP) :: reacflux, fluxout !, dfluxes

   CALL vertrates(candi,cflux)
   reacflux = 0.0
   DO j = 1,6
     reacflux=reacflux+cflux(j)
   END DO

!   fluxout = (y(cor1y,npt)+y(cor2y,npt)+y(cor0y,npt))*wvel(npt)*ps(npt)
   fluxout = (candi%y(candi%DOCLy,candi%npt)+candi%y(candi%POCLy,candi%npt)+candi%y(candi%DOCRy,candi%npt)+candi%y(candi%POCRy,candi%npt))*candi%wvel(candi%npt)*candi%ps(candi%npt)
   ! note: should split to wvel and uvel for solids and liquids


   WRITE(fid,'(a50)') '============================================='
   WRITE(fid,'(a50)') 'Results from Corg-mass conservation check'
   WRITE(fid,'(a25,f20.10)') 'Flux in (fg0+fg1+fg2) : ',(candi%fg0+candi%fg1+candi%fg2)
   WRITE(fid,'(a25,f20.10)') 'Flux out              : ',fluxout
   WRITE(fid,'(a25,f20.10)') 'Flux through reac     : ',reacflux
   WRITE(fid,'(a25,f20.10)') 'O2-POC-flux           : ',cflux(1)
   WRITE(fid,'(a25,f20.10)') 'NO-POC-flux           : ',cflux(2)
   IF(candi%param%simMnFe) THEN
   WRITE(fid,'(a25,f20.10)') 'Mn-POC-flux           : ',cflux(3)
   WRITE(fid,'(a25,f20.10)') 'Fe-POC-flux           : ',cflux(4)
   END IF
   WRITE(fid,'(a25,f20.10)') 'SO-POC-flux           : ',cflux(5)
   WRITE(fid,'(a25,f20.10)') 'CH-POC-flux           : ',cflux(6)
   WRITE(fid,'(a25,f20.10)') 'Corgrate              : ',cflux(7)
   WRITE(fid,'(a50)') '============================================='
   WRITE(fid,'(a60)') 'Definition of the fluxes: out of sediment (-) into it (+)'
   WRITE(fid,'(a25)') 'Units: [umol/cm**2/a]'
   WRITE(fid,'(a20,f20.10)') 'Flux O2         :', dfluxes(candi,2,candi%o2y)
   WRITE(fid,'(a20,f20.10)') 'Flux O2 Theory  :',((candi%sc+2.0*candi%sn)/candi%sc)* cflux(1)+2.0*candi%srsox+candi%srmnox+candi%srfeox &
     +2.0*candi%srch4ox
   WRITE(fid,'(a20,f20.10)') 'Flux NO3        :',dfluxes(candi,2,candi%no3y)
   WRITE(fid,'(a20,f20.10)') 'Flux SO4        :',dfluxes(candi,2,candi%so4y)
   WRITE(fid,'(a20,f20.10)') 'Flux PO4        :',dfluxes(candi,2,candi%po4ly)
   WRITE(fid,'(a20,f20.10)') 'Flux PO4 Theory :',-cflux(7)*candi%sp/candi%sc
   WRITE(fid,'(a20,f20.10)') 'Flux NH4        :',dfluxes(candi,2,candi%nh4y)
   IF(candi%param%simMnFe) THEN
   WRITE(fid,'(a20,f20.10)') 'Flux MnII       :',dfluxes(candi,2,candi%mniiy)
   WRITE(fid,'(a20,f20.10)') 'Flux FeII       :',dfluxes(candi,2,candi%feiiy)
   END IF
   WRITE(fid,'(a20,f20.10)') 'Flux CH4        :',dfluxes(candi,2,candi%ch4y)

   WRITE(fid,'(a50)') '============================================='
   WRITE(fid,'(a50)') 'Fluxes at the bottom:                        '
   WRITE(fid,'(a20,f20.10)') 'Flux O2         :',dfluxes(candi,candi%npt,candi%o2y)
   WRITE(fid,'(a20,f20.10)') 'Flux NO3        :',dfluxes(candi,candi%npt,candi%no3y)
   WRITE(fid,'(a20,f20.10)') 'Flux SO4        :',dfluxes(candi,candi%npt,candi%so4y)
   WRITE(fid,'(a20,f20.10)') 'Flux PO4        :',dfluxes(candi,candi%npt,candi%po4ly)
   WRITE(fid,'(a20,f20.10)') 'Flux NH4        :',dfluxes(candi,candi%npt,candi%nh4y)
   IF(candi%param%simMnFe) THEN
   WRITE(fid,'(a20,f20.10)') 'Flux MnII       :',dfluxes(candi,candi%npt,candi%mniiy)
   WRITE(fid,'(a20,f20.10)') 'Flux FeII       :',dfluxes(candi,candi%npt,candi%feiiy)
   END IF
   WRITE(fid,'(a20,f20.10)') 'Flux CH4        :',dfluxes(candi,candi%npt,candi%ch4y)

   WRITE(fid,'(a20,f20.10)') 'Flux OM         :', -(candi%y(candi%DOCLy,candi%npt)+candi%y(candi%POCLy,candi%npt)+candi%y(candi%POCRy,candi%npt)+candi%y(candi%DOCRy,candi%npt) &
       )*candi%ps(candi%npt)*candi%wvel(candi%npt)

   WRITE(fid,'(a50)') '============================================='
   WRITE(fid,'(a50)') 'Vertical integrated secondary redox-rates'
   WRITE(fid,'(a25,f20.10)') 'RNHOX                 : ',candi%srnhox
   WRITE(fid,'(a25,f20.10)') 'RSOX                  : ',candi%srsox
   WRITE(fid,'(a25,f20.10)') 'RMNOX                 : ',candi%srmnox
   WRITE(fid,'(a25,f20.10)') 'RFEOX                 : ',candi%srfeox
   WRITE(fid,'(a25,f20.10)') 'RMNFE                 : ',candi%srmnfe
   WRITE(fid,'(a25,f20.10)') 'RFENO3                : ',candi%srfeno3
   WRITE(fid,'(a25,f20.10)') 'RMNO2TS               : ',candi%srmno2ts
   WRITE(fid,'(a25,f20.10)') 'RCH4OX                : ',candi%srch4ox
   WRITE(fid,'(a25,f20.10)') 'RCH4SO4               : ',candi%srch4so4
   WRITE(fid,'(a25,f20.10)') 'RNO3TS                : ',candi%srno3ts
   WRITE(fid,'(a25,f20.10)') 'RTSFE3                : ',candi%srtsfe3
   IF(candi%param%simFeS) THEN
   WRITE(fid,'(a25,f20.10)') 'RFESOX                : ',candi%srfesox
   WRITE(fid,'(a25,f20.10)') 'RFESFE3               : ',candi%srfesfe3
   WRITE(fid,'(a25,f20.10)') 'RFESMN4               : ',candi%srfesmn4
   END IF
   WRITE(fid,'(a25,f20.10)') 'RMnNO3                : ',candi%srmnno3
   IF(candi%simFeII) THEN
   WRITE(fid,'(a25,f20.10)') 'RFe1OX                : ',candi%srfe1ox
   WRITE(fid,'(a25,f20.10)') 'RFe2OX                : ',candi%srfe2ox
   WRITE(fid,'(a25,f20.10)') 'RFe1NO3               : ',candi%srfe1no3
   WRITE(fid,'(a25,f20.10)') 'RFe2NO3               : ',candi%srfe2no3
   END IF
   WRITE(fid,'(a50)') '============================================='
   WRITE(fid,'(a50)') 'Vertical integrated equilibrium rates'
   WRITE(fid,'(a25,f20.10)') 'RDCO2                : ',candi%srdco2
   WRITE(fid,'(a25,f20.10)') 'RDHCO3               : ',candi%srdhco3
   WRITE(fid,'(a25,f20.10)') 'RDCO3                : ',candi%srdco3
   WRITE(fid,'(a50)') '============================================='

   RETURN
 END SUBROUTINE mass
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!rl   ****************************************************************
!rl   Early diagenese model C. CANDI from R. Luff (1996-2003), based on
!rl   CANDI from  B.P. Boudreau.
!rl   ****************************************************************
!rl   VERSION 3.0 Roger Luff (rluff@gmx.de) GEOMAR, Kiel, Germany
!rl   ****************************************************************
!rl   FILE: VERTRATES.f
!rl   AIM : Calcultes the vertical integrated reaction rates.
!rl         incl for corg
!rl         cflux(1):organic matter degradation by oxygen
!rl         cflux(2): "        "       "        "  nitrate
!rl          ....
!rl         cflux(7): corg degradation
!rl         Performs simple numerical integration using the (compound)
!rl         Traperoidal Rule.
!rl   ****************************************************************
!------------------------------------------------------------------------------!
 SUBROUTINE vertrates(candi,cflux)
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   !-- Outgoing
   REAL(SEDP), INTENT(OUT)            :: cflux(7)
   !-- Local
   REAL(SEDP) :: r(7)
   INTEGER  :: i,j

   CALL RATES(candi)

   candi%srnhox   = candi%rnh4ox(2) *candi%dh(2)*candi%poros(2)* HALF
   candi%srsox    = candi%rtsox(2)  *candi%dh(2)*candi%poros(2)* HALF
   candi%srmnox   = candi%rmnox(2)  *candi%dh(2)*candi%poros(2)* HALF
   candi%srfeox   = candi%rfeox(2)  *candi%dh(2)*candi%poros(2)* HALF
   candi%srmnfe   = candi%rfemnA(2) *candi%dh(2)*candi%ps(2)   * HALF
   candi%srfeno3  = candi%rfeno3(2) *candi%dh(2)*candi%ps(2)   * HALF
   candi%srmno2ts = candi%rtsmnA(2) *candi%dh(2)*candi%ps(2)   * HALF
   candi%srch4ox  = candi%rch4ox(2) *candi%dh(2)*candi%poros(2)* HALF
   candi%srch4so4 = candi%rch4so4(2)*candi%dh(2)*candi%poros(2)* HALF
   candi%srno3ts  = candi%rtsno3(2) *candi%dh(2)*candi%poros(2)* HALF
   candi%srtsfe3  = candi%rtsfeA(2) *candi%dh(2)*candi%poros(2)* HALF
   IF(candi%param%simFeS) THEN
     candi%srfesox  = candi%rfesox(2) *candi%dh(2)*candi%ps(2)  * HALF
     candi%srfesfe3 = candi%rfesfeA(2)*candi%dh(2)*candi%ps(2)  * HALF
     candi%srfesmn4 = candi%rfesmnA(2)*candi%dh(2)*candi%ps(2)  * HALF
   END IF
   candi%srmnno3  = candi%rmnno3(2)*candi%dh(2)*candi%poros(2) * HALF
     !                                        SUBROUTINE vertrates(cflux)

   r(1) =candi%rgC(2)*candi%fo2(2)
   r(2) =candi%rgC(2)*candi%fno3(2)
   r(3) =candi%rgC(2)*candi%fmno2(2)
   r(4) =candi%rgC(2)*candi%ffeoh(2)
   r(5) =candi%rgC(2)*candi%fso4(2)
   r(6) =candi%rgC(2)*candi%fmet(2)
   r(7) =candi%rgC(2)*candi%rox(2)
   DO j = 1,7
     cflux(j) = r(j)*candi%ps(2)*candi%dh(2)* HALF
   END DO

   DO i = 3,(candi%npt-1)
     candi%srnhox   = candi%srnhox   + candi%rnh4ox(i) *candi%dh(i)*candi%poros(i)
     candi%srsox    = candi%srsox    + candi%rtsox(i)  *candi%dh(i)*candi%poros(i)
     candi%srmnox   = candi%srmnox   + candi%rmnox(i)  *candi%dh(i)*candi%poros(i)
     candi%srfeox   = candi%srfeox   + candi%rfeox(i)  *candi%dh(i)*candi%poros(i)
     candi%srmnfe   = candi%srmnfe   + candi%rfemnA(i) *candi%dh(i)*candi%ps(i)
     candi%srfeno3  = candi%srfeno3  + candi%rfeno3(i) *candi%dh(i)*candi%ps(i)
     candi%srmno2ts = candi%srmno2ts + candi%rtsmnA(i) *candi%dh(i)*candi%ps(i)
     candi%srch4ox  = candi%srch4ox  + candi%rch4ox(i) *candi%dh(i)*candi%poros(i)
     candi%srch4so4 = candi%srch4so4 + candi%rch4so4(i)*candi%dh(i)*candi%poros(i)
     candi%srno3ts  = candi%srno3ts  + candi%rtsno3(i) *candi%dh(i)*candi%poros(i)
     candi%srtsfe3  = candi%srtsfe3  + candi%rtsfeA(i) *candi%dh(i)*candi%ps(i)
     IF(candi%param%simFeS) THEN
       candi%srfesox  = candi%srfesox  + candi%rfesox(i) *candi%dh(i)*candi%ps(i)
       candi%srfesfe3 = candi%srfesfe3 + candi%rfesfeA(i)*candi%dh(i)*candi%ps(i)
       candi%srfesmn4 = candi%srfesmn4 + candi%rfesmnA(i)*candi%dh(i)*candi%ps(i)
     END IF
     candi%srmnno3  = candi%srmnno3  + candi%rmnno3(i)*candi%dh(i)*candi%poros(i)
     !                                        SUBROUTINE vertrates(cflux)

     r(1) = candi%rgC(i)*candi%fo2(i)
     r(2) = candi%rgC(i)*candi%fno3(i)
     r(3) = candi%rgC(i)*candi%fmno2(i)
     r(4) = candi%rgC(i)*candi%ffeoh(i)
     r(5) = candi%rgC(i)*candi%fso4(i)
     r(6) = candi%rgC(i)*candi%fmet(i)
     r(7) = candi%rgC(i)*candi%rox(i)

     DO j=1,7
       cflux(j)=cflux(j)+r(j)*candi%ps(i)*candi%dh(i)
     END DO

   END DO
     !                                        SUBROUTINE vertrates(cflux)

   candi%srnhox   = candi%srnhox   + candi%rnh4ox(candi%npt) *candi%dh(candi%npt)*candi%poros(candi%npt)* HALF
   candi%srsox    = candi%srsox    + candi%rtsox(candi%npt)  *candi%dh(candi%npt)*candi%poros(candi%npt)* HALF
   candi%srmnox   = candi%srmnox   + candi%rmnox(candi%npt)  *candi%dh(candi%npt)*candi%poros(candi%npt)* HALF
   candi%srfeox   = candi%srfeox   + candi%rfeox(candi%npt)  *candi%dh(candi%npt)*candi%poros(candi%npt)* HALF
   candi%srmnfe   = candi%srmnfe   + candi%rfemnA(candi%npt) *candi%dh(candi%npt)*candi%ps(candi%npt)* HALF
   candi%srfeno3  = candi%srfeno3  + candi%rfeno3(candi%npt) *candi%dh(candi%npt)*candi%ps(candi%npt)* HALF
   candi%srmno2ts = candi%srmno2ts + candi%rtsmnA(candi%npt) *candi%dh(candi%npt)*candi%ps(candi%npt)* HALF
   candi%srch4ox  = candi%srch4ox  + candi%rch4ox(candi%npt) *candi%dh(candi%npt)*candi%poros(candi%npt)* HALF
   candi%srch4so4 = candi%srch4so4 + candi%rch4so4(candi%npt)*candi%dh(candi%npt)*candi%poros(candi%npt)* HALF
   candi%srno3ts  = candi%srno3ts  + candi%rtsno3(candi%npt) *candi%dh(candi%npt)*candi%poros(candi%npt)* HALF
   candi%srtsfe3  = candi%srtsfe3  + candi%rtsfeA(candi%npt) *candi%dh(candi%npt)*candi%ps(candi%npt)* HALF
   IF(candi%param%simFeS) THEN
     candi%srfesox  = candi%srfesox  + candi%rfesox(candi%npt) *candi%dh(candi%npt)*candi%ps(candi%npt)* HALF
     candi%srfesfe3 = candi%srfesfe3 + candi%rfesfeA(candi%npt)*candi%dh(candi%npt)*candi%ps(candi%npt)* HALF
     candi%srfesmn4 = candi%srfesmn4 + candi%rfesmnA(candi%npt)*candi%dh(candi%npt)*candi%ps(candi%npt)* HALF
   END IF
   candi%srmnno3  = candi%srmnno3  + candi%rmnno3(candi%npt)*candi%dh(candi%npt)*candi%poros(candi%npt)* HALF

!                                        SUBROUTINE vertrates(cflux)

   r(1) = candi%rgC(candi%npt)*candi%fo2(candi%npt)
   r(2) = candi%rgC(candi%npt)*candi%fno3(candi%npt)
   r(3) = candi%rgC(candi%npt)*candi%fmno2(candi%npt)
   r(4) = candi%rgC(candi%npt)*candi%ffeoh(candi%npt)
   r(5) = candi%rgC(candi%npt)*candi%fso4(candi%npt)
   r(6) = candi%rgC(candi%npt)*candi%fmet(candi%npt)
   r(7) = candi%rgC(candi%npt)*candi%rox(candi%npt)

   DO j = 1,7
     cflux(j)=cflux(j)+r(j)*candi%ps(candi%npt)*candi%dh(candi%npt)*0.5
   END DO

   RETURN
 END SUBROUTINE vertrates
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
!rl   ****************************************************************
!rl   Early diagenese model C. CANDI from R. Luff (1996-2003), based on
!rl   CANDI from  B.P. Boudreau.
!rl   ****************************************************************
!rl   VERSION 3.0 Roger Luff (rluff@gmx.de) GEOMAR, Kiel, Germany
!rl   ****************************************************************
!rl   FILE: STATUS.f
!rl   AIM : calculate status of steady state and print runtime information
!rl         on the screen
!rl   ****************************************************************
!------------------------------------------------------------------------------!
SUBROUTINE CheckSteadyStatus(candi,steadiness)
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   REAL(SEDP), INTENT(OUT)            :: steadiness
   REAL(SEDP)                         :: MAX
   INTEGER  :: j,i,idx

   MAX = 0.0
   DO i = 0,candi%nSPECIES-1
     steadiness = 0.0
     DO j = 1, candi%npt
       steadiness = steadiness + ABS(candi%ydot(i,j))
     END DO
     IF (steadiness > MAX) THEN
       MAX = steadiness
       idx = i
     END IF
   END DO
   steadiness = MAX

   RETURN
 END SUBROUTINE CheckSteadyStatus
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
!rl   ****************************************************************
!rl   Early diagenese model C. CANDI from R. Luff (1996-2003), based on
!rl   CANDI from  B.P. Boudreau.
!rl   ****************************************************************
!rl   VERSION 3.0 Roger Luff (rluff@gmx.de) GEOMAR, Kiel, Germany
!rl   ****************************************************************
!rl   FILE: WRITERATES.f
!rl   AIM : Stores the reaction rates into a file. The
!rl   name of the output file has to be defined in STEUER.DAT.
!rl   ****************************************************************
!------------------------------------------------------------------------------!
 SUBROUTINE WriteRates (candi,fid)
   !-- Incoming
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   INTEGER, INTENT(IN) :: fid
   !-- Local
   REAL(SEDP) :: mdummy
   INTEGER  :: i

   !rl   Calculate and print the mass-balance results in MASS.f
   !rl   aviod warnings:
   mdummy = 0.0
   CALL mass(candi,mdummy,fid)

   1012 FORMAT (100F30.15)
   !rl   Get the rates from the last timestep and write it to the file
   OPEN(16,FILE=candi%filerat,STATUS='UNKNOWN')
   CALL RATES(candi)
   DO i=1,candi%npt
     WRITE(16,1012) &
                   candi%rpar(i), &
                   candi%rgC(i)*candi%rox(i),&
                   candi%rgC(i)*candi%fo2(i),&
                   candi%rgC(i)*candi%fno3(i),&
                   candi%rgC(i)*candi%fmno2(i),&
                   candi%rgC(i)*candi%ffeoh(i),&
                   candi%rgC(i)*candi%fso4(i),&
                   candi%rgC(i)*candi%fmet(i),&
                   candi%RNH4OX(i),&
                   candi%RMnOX(i),&
                   candi%RFeOX(i),&
                   candi%RTSOX(i),&
                   candi%RCH4OX(i),&
                   candi%RFeSOX(i),&
                   candi%RFeS2OX(i),&
                   candi%RNH4NO2(i),&
                   candi%RMnNO3(i),&
                   candi%RFeNO3(i),&
                   candi%RTSNO3(i),&
                   candi%RFeMnA(i),&
                   candi%RFeMnB(i),&
                   candi%RTSMnA(i),&
                   candi%RTSMnB(i),&
                   candi%RFeSMnA(i),&
                   candi%RFeSMnB(i),&
                   candi%RTSFeA(i),&
                   candi%RTSFeB(i),&
                   candi%RFeSFeA(i),&
                   candi%RFeSFeB(i),&
                   candi%RCH4SO4(i),&
                   candi%RMnAge(i),&
                   candi%RFeOHAppt(i),&
                   candi%RFeOHBppt(i),&
                   candi%RFeAge(i),&
                   candi%RFeSppt(i),&
                   candi%RPyrite(i),&
                   candi%RXSppt(i),&
                   candi%RSidppt(i),&
                   candi%RRodppt(i),&
                   candi%RCalppt(i),&
                   candi%RMnO2Appt(i),&
                   candi%RMnO2Bppt(i),&
                   candi%RPO4ads(i),&
                   candi%RNH4ads(i),&
                   candi%RMnO2Bppt(i)
   END DO

   CLOSE(16)

   RETURN
 END SUBROUTINE WriteRates
!------------------------------------------------------------------------------!


!New file by Dan
!------------------------------------------------------------------------------!
!rl   ****************************************************************
!rl   ****************************************************************
!rl   ****************************************************************
!rl   ****************************************************************
!------------------------------------------------------------------------------!
 SUBROUTINE WriteRates2 (candi,fid)
   !-- Incoming
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   INTEGER, INTENT(IN) :: fid
   !-- Local
   REAL(SEDP) :: mdummy
   INTEGER  :: i
   !rl   Calculate and print the mass-balance results in MASS.f
   !rl   aviod warnings:
   mdummy = 0.0
   CALL mass(candi,mdummy,fid)
   1012 FORMAT (100F30.15)
   !rl   Get the rates from the last timestep and write it to the file
   OPEN(16,FILE=candi%filerat2,STATUS='UNKNOWN')
   CALL RATES(candi)
   DO i=1,candi%npt
     WRITE(16,1012) &
                   candi%rpar(i), &
                   candi%RDenO2DHyd(i), &
                   candi%RDenNO3DHyd(i), &
                   candi%RFerDHyd(i), &
                   candi%RAerOAc(i)
   END DO
   CLOSE(16)
   RETURN
 END SUBROUTINE WriteRates2
!------------------------------------------------------------------------------!

!End Dan

!New file by Dan
!------------------------------------------------------------------------------!
!rl   ****************************************************************
!rl   ****************************************************************
!rl   ****************************************************************
!rl   ****************************************************************
!------------------------------------------------------------------------------!
 SUBROUTINE WriteOAcRates (candi,fid)
   !-- Incoming
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   INTEGER, INTENT(IN) :: fid
   !-- Local
   REAL(SEDP) :: mdummy
   INTEGER  :: i
   !rl   Calculate and print the mass-balance results in MASS.f
   !rl   aviod warnings:
   mdummy = 0.0
   CALL mass(candi,mdummy,fid)
   1012 FORMAT (100F30.15)
   !rl   Get the rates from the last timestep and write it to the file
   OPEN(16,FILE=candi%fileratOAc,STATUS='UNKNOWN')
   CALL RATES(candi)
   DO i=1,candi%npt
     WRITE(16,1012) &
                   candi%rpar(i), &
                   candi%RAerOAc(i), &
                   candi%RDenOAc(i), &
                   candi%RManOAc(i), &
                   candi%RIroOAc(i), &
                   candi%RSulOAc(i), &
                   candi%RMetOAc(i), &
                   candi%ROAc(i)
   END DO
   CLOSE(16)
   RETURN
 END SUBROUTINE WriteOAcRates
!------------------------------------------------------------------------------!

!End Dan
!------------------------------------------------------------------------------!
FUNCTION GetAEDColNum(candi,iden) RESULT(column)
   !-- Incoming                                                         !
   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   CHARACTER (LEN=*) :: iden                                            !
   !-- Outgoing                                                         !
   INTEGER  :: column,i
   !-- Local


   column = -1
   DO i = 0,SIZE(candi%CANDIVarNames)-1
     IF(candi%CANDIVarNames(i) == iden ) THEN
       column = i
       EXIT
     END IF
   END DO
   print *,'candi col ',iden,column
   !print *,'CANDIVarNames',   candi%CANDIVarNames(i)





END FUNCTION GetAEDColNum


!------------------------------------------------------------------------------!
END MODULE aed2_sedcandi
