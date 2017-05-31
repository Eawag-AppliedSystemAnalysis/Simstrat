module Forcing
  use SimstratModel
  implicit none
  save

  real(dp), parameter :: r_a = 0.03_dp           ! Ratio of reflected to total long-wave iradiance
  real(dp), parameter :: B0 = 0.61_dp            ! Bowen constant
  real(dp), parameter :: sig = 5.67e-8_dp        ! Stefan-Boltzmann constant [W/m2/K4]

  private
  public ReadForcing, doForcing

contains
  !Read forcing file to get values A_cur at given datum
  !AG 2014: revision + correction
  !####################################################################
  subroutine ReadForcing(datum,A_s,A_e,A_cur,nval,idx)
  !####################################################################

        implicit none

        ! Global variables
        integer, intent(in) :: idx,nval
        real(dp), intent(in) :: datum  !Required date
        real(dp), intent(inout) :: A_cur(1:nval), A_s(1:nval), A_e(1:nval) !output values, Start and end values

        ! Local variables
        real(dp) :: tb_start, tb_end !Start and end time
        integer :: eof,i

        save tb_start, tb_end
        save eof

        if (idx==1) then
           if(disp_dgn/=0) write(6,*) 'Forcing : ',trim(ForcingName)
           if(disp_dgn==2) write(6,*) 'Starting to read forcing file...'
           open(20,status='old',file=ForcingName)
           eof = 0
           !Read first values
           read(20,*,end=9)
           read(20,*,end=9) tb_start, (A_s(i),i=1,nval)
           if(datum<tb_start) write(6,*) 'Warning: first forcing date after simulation start time.'
           read(20,*,end=7) tb_end, (A_e(i),i=1,nval)
        end if

        if (datum<=tb_start .or. eof==1) then !If datum before first date or end of file reached
            goto 8
        else
            do while (datum>tb_end)         !Move to appropriate interval to get correct values
                tb_start = tb_end
                A_s(1:nval) = A_e(1:nval)
                !Read next value
                read(20,*,end=7) tb_end, (A_e(i),i=1,nval)
            end do
            !Linearly interpolate values at correct datum
            A_cur(1:nval) = A_s(1:nval) + (datum-tb_start) * (A_e(1:nval)-A_s(1:nval))/(tb_end-tb_start)
        end if
        return

  7     eof = 1
        if(datum>tb_start) write(6,*) 'Warning: last forcing date before simulation end time.'
  8     A_cur(1:nval) = A_s(1:nval)       !Take first value of current interval
        return

  9     write(6,*) 'Error reading forcing file (no data found).'
        stop

  end subroutine ReadForcing



  !Compute appropriate forcing parameters at given datum
  !AG 2014: revision
  !######################################################################
  subroutine doForcing(datum,T_surf,idx,tx,ty,u_taus,rad0,heat,SST,u10,v10,wind,Wf)
  !######################################################################
        implicit none

        ! Global iables
        integer, intent(in) :: idx
        real(dp), intent(in) :: datum, T_surf
        real(dp), intent(inout) :: tx,ty,u_taus, rad0,heat,SST,u10,v10,wind,Wf


        ! Local iables
        real(dp) :: tau
        real(dp) :: A_s(7), A_e(7), A_cur(7)
        real(dp) :: fu, Vap_wat, heat0
        real(dp) :: T_atm, F_glob, Vap_atm, Cloud
        real(dp) :: H_A, H_K, H_V, H_W

        save A_s, A_e

        if (NBC==1) then
            call ReadForcing(datum,A_s,A_e,A_cur,4+WindFilt,idx)
            u10 = A_cur(1)*f_wind       !MS 2014: added f_wind
            v10 = A_cur(2)*f_wind       !MS 2014: added f_wind
            wind = sqrt(u10**2+v10**2)  !AG 2014
            SST = A_cur(3) !Sea surface temperature
            rad0 = A_cur(4)*(1-albsw)*(1-beta_sol) ! MS: added beta_sol and albsw
            heat = 0.0_dp
            if(WindFilt==1) Wf = A_cur(5)  !AG 2014
        else if (NBC>=2) then
            if (NBC==2) then            ! date, U,V,Tatm,Hsol,Vap
                call ReadForcing(datum,A_s,A_e,A_cur,5+WindFilt,idx)
                u10 = A_cur(1)*f_wind      !MS 2014: added f_wind
                v10 = A_cur(2)*f_wind      !MS 2014: added f_wind
                T_atm = A_cur(3)
                F_glob = A_cur(4)*(1-albsw)
                Vap_atm = A_cur(5)
                Cloud = 0.5_dp
                if(WindFilt==1) Wf = A_cur(6) !AG 2014
            else if (NBC==3) then       ! date,U10,V10,Tatm,Hsol,Vap,Clouds
                call ReadForcing(datum,A_s,A_e,A_cur,6+WindFilt,idx)
                u10 = A_cur(1)*f_wind      !MS 2014: added f_wind
                v10 = A_cur(2)*f_wind      !MS 2014: added f_wind
                T_atm = A_cur(3)
                F_glob = A_cur(4)*(1-albsw)
                Vap_atm = A_cur(5)
                Cloud = A_cur(6)
                if (Cloud<0 .or. Cloud>1) then
                    write(6,*) 'Cloudiness should always be between 0 and 1.'
                    stop
                end if
                if(WindFilt==1) Wf = A_cur(7) !AG 2014
            else if (NBC==4) then       ! date,U10,V10,Hnet,Hsol
                call ReadForcing(datum,A_s,A_e,A_cur,4+WindFilt,idx)
                u10 = A_cur(1)*f_wind      !MS 2014: added f_wind
                v10 = A_cur(2)*f_wind      !MS 2014: added f_wind
                heat0 = A_cur(3)           !MS 2014
                F_glob = A_cur(4)*(1-albsw)
                if(WindFilt==1) Wf = A_cur(5) !AG 2014
            else
                write(6,*) 'Error: wrong forcing type (must be 1, 2, 3 or 4).'
                stop
            end if
            wind = sqrt(u10**2+v10**2)  !AG 2014


            if (NBC/=4) then ! in the water column

                ! Wind function (Livingstone & Imboden 1989)
                fu = 4.40_dp+1.82_dp*wind+0.26_dp*(T_surf-T_atm)
                !fu = 5.44+2.19*wind+0.24*(T_surf-T_atm)
                fu = fu*p_windf    ! Provided fitting factor p_windf (~1)
                ! Water vapor saturation pressure in air at water temperature (Gill 1992) [millibar]
                Vap_wat = 10**((0.7859_dp+0.03477_dp*T_surf)/(1+0.00412_dp*T_surf))
                Vap_wat = Vap_wat*(1+1e-6_dp*p_air*(4.5_dp+0.00006_dp*T_surf**2))
                ! Solar short-wave radiation absorbed
                rad0 = F_glob*(1-beta_sol) !MS: added beta_sol

                ! Long-wave radiation from sky (Livingstone & Imboden 1989)
                ! H_A = 1.24*sig*(1-r_a)*(1+0.17*Cloud**2)*(Vap_atm/(273.15+T_atm))**(1./7)*(273.15+T_atm)**4
                ! Long-wave radiation according to Dilley and O'Brien
                ! see Flerchinger et al. (2009)
                H_A = (1-r_a) * ((1-0.84_dp*Cloud) * (59.38_dp + 113.7_dp*((T_atm+273.15_dp)/273.16_dp)**6&
                &   + 96.96_dp*sqrt(465*Vap_atm/(T_atm+273.15_dp)*0.04_dp)) / 5.67e-8_dp / &
                &   (T_atm+273.15_dp)**4 + 0.84_dp * Cloud) *5.67e-8_dp*(T_atm+273.15_dp)**4
                H_A = H_A*p_radin    ! Provided fitting factor p_radin (~1)
                ! Long-wave radiation from water body (black body)
                H_W = -0.97_dp*sig*(T_surf+273.15_dp)**4
                ! Flux of sensible heat (convection)
                H_K = -B0*fu*(T_surf-T_atm)
                ! Flux of latent heat (evaporation, condensation)
                H_V = -fu*(Vap_wat-Vap_atm)
                ! Global heat flux (positive: air to water, negative: water to air)
                heat = H_A + H_W + H_K + H_V + F_glob*beta_sol !MS: added term with beta_sol
            else
                heat = heat0 + F_glob*beta_sol !MS: added term with beta_sol
            end if
            if ((T_surf<0).and.(heat<0)) heat=0.
        end if

        !Drag coefficient as a function of wind speed (AG 2014)
        if (ModC10==2) then !Ocean model
            C10 = -0.000000712_dp*wind**2+0.00007387_dp*wind+0.0006605_dp
        else if (ModC10==3) then !Lake model (Wüest and Lorke 2003)
            if (wind<=0.1) then
                C10 = 0.06215_dp
            else if (wind<=3.85_dp) then
                C10 = 0.0044_dp*wind**(-1.15_dp)
            else !Polynomial approximation of Charnock's law
                C10 = -0.000000712_dp*wind**2+0.00007387_dp*wind+0.0006605_dp
                !C10 = -0.000000385341*wind**2+0.0000656519*wind+0.000703768
                !C10 = 0.0000000216952*wind**3-0.00000148692*wind**2+0.0000820705*wind+0.000636251
            end if
  !         if (abs(u10)<=4.7 .and. abs(u10)>0.5) then
  !             C10 = 1.78e-3*wind**(-1.3561)
  !         elseif (abs(u10)>4.7) then
  !             C10 = 2.13e-4
  !         else
  !             C10 = 0.0046
  !         end if
        end if

        tau = C10*rho_air/rho_0*wind**2
        u_taus = sqrt(tau)

        tx = C10*rho_air/rho_0*wind*u10
        ty = C10*rho_air/rho_0*wind*v10

        !tx = tau*u10/wind
        !ty = tau*v10/wind

  !     if (u10==0.) then
  !         tx = 0.
  !     else if (abs(u10)<4.7 .and. abs(u10)>0.1) then
  !         tx = u10*u10*u10/abs(u10)*rho_air/rho_0*2.38e-3*abs(u10)**-1.5
  !     else
  !         tx = u10*u10*u10/abs(u10)*rho_air/rho_0*2.34e-4
  !     end if
  !     if (v10==0.) then
  !         ty = 0.
  !     else if (abs(v10)<4.7 .and. abs(v10)>0.1) then
  !         ty = v10*v10*v10/abs(v10)*rho_air/rho_0*2.38e-3*abs(v10)**-1.5
  !     else
  !         ty = v10*v10*v10/abs(v10)*rho_air/rho_0*2.34e-4
  !     end if

        return
  end subroutine doForcing


end module Forcing
