ptf #
*** Files *************************************************
LacDeNeuchatel\InitialConditions_2003_02_20.dat
LacDeNeuchatel\Grid_Neuchatel.dat
LacDeNeuchatel\Morph_Neuchatel.dat
LacDeNeuchatel\SimForceStationComb_Neuchatel_1994_2014_WFILT.dat3
LacDeNeuchatel\Abs_Secchi.dat
LacDeNeuchatel_Results_wfilt\
LacDeNeuchatel\Output_depth.dat
LacDeNeuchatel\Output_time.dat
LacDeNeuchatel\Qinp0.dat
LacDeNeuchatel\Qout0.dat
LacDeNeuchatel\Tinp0.dat
LacDeNeuchatel\Sinp0.dat
\
\
*** Model setup *******************************************
300          Timestep dt [s]
32873        Start time [d]
40543        End time [d]
*** Model, conditions and output selection ****************
1            Turbulence model (1:k-epsilon, 2:MY)
0            Biogeochemistry model (0/default:off, 1:on)
1            Stability function (1:constant, 2:quasi-equilibrium)
1            Flux condition (0:Dirichlet condition, 1:no-flux)
3            Forcing (1:Wind+Temp+SolRad, 2:(1)+Vap, 3:(2)+Cloud, 4:Wind+HeatFlux+SolRad)
1			 Use filtered wind to compute seiche energy (0/default:off, 1:on) (if 1:on, one more column is needed in forcing file)
2            Seiche normalization (1:max N^2, 2:integral)
3            Wind drag model (1/default:constant, 2:ocean (increasing), 3:lake (Wüest and Lorke 2003))
0			 Inflow placement (0/default:manual, 1:density-driven)
0            Pressure gradients (0:off, 1:Svensson 1978, 2:?)
1            Enable salinity transport (0:off, 1/default:on)
0            Display simulation (0:off, 1:when data is saved, 2:at each iteration, 3:extra display)
0            Display diagnose (0:off, 1:standard display, 2:extra display)
10           Averaging data
*** Model parameters **************************************
#lat       # Lat [°]		Latitude for Coriolis parameter
#pair      # p_air [mbar]	Air pressure
#alpha     # a_seiche [-]	Fraction of wind energy to seiche energy
#qNN       # q_NN			Fit parameter for distribution of seiche energy
#fwind     # f_wind	[-]		Fraction of forcing wind to wind at 10m (W10/Wf)
#C10       # C10 [-]		Wind drag coefficient (used if wind drag model is 1:lazy)
#CD        # CD [-]			Bottom friction coefficient
#fgeo      # fgeo [W/m2]	Geothermal heat flux
#kmin      # k_min [J/kg]	Minimal value for TKE
#p1        # p1				Fit parameter for absorption of IR radiation from sky
#p2        # p2				Fit parameter for convective and latent heat fluxes
#beta      # beta [-]		Fraction of short-wave radiation directly absorbed as heat?
#albsw     # albsw [-]		Albedo for reflection of short-wave radiation
