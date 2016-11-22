!Interpolation of yi on grid zi (based on the given y on grid z)
!####################################################################
subroutine Interp(z,y,num_z,zi,yi,num_zi)
!####################################################################

      double precision z(0:num_z),y(0:num_z),zi(0:num_zi),yi(0:num_zi)
      integer num_z,num_zi

      integer posk1,posk2,posi,i

      !Assign closest value if out of given grid
      posk1 = 0
      do while (zi(posk1)<=z(0))
         yi(posk1) = y(0)
         posk1 = posk1+1
      end do
      posk2 = num_zi
      do while (zi(posk2)>=z(num_z))
         yi(posk2) = y(num_z)
         posk2 = posk2-1
      end do

      !Linear interpolation
      posi = 0
      do i=posk1,posk2
         do while (zi(i)>z(posi+1))
            posi = posi+1
         end do
         yi(i)=y(posi)+(zi(i)-z(posi))*(y(posi+1)-y(posi))/(z(posi+1)-z(posi))
      end do

      return
end


!Interpolation of yi on grid zi (based on the given y on grid z)
!####################################################################
subroutine Interp_nan(z,y,num_z,zi,yi,num_zi)
!####################################################################

      double precision z(0:num_z),y(0:num_z),zi(0:num_zi),yi(0:num_zi)
      integer num_z,num_zi

      integer posk1,posk2,posi,i

      !Assign NaN if out of given grid
      posk1 = 0
      do while (zi(posk1)<z(0))
         yi(posk1) = 0.
         yi(posk1) = 0./yi(posk1)   ! NaN
         posk1 = posk1+1
      end do
      posk2 = num_zi
      do while (zi(posk2)>z(num_z))
         yi(posk1) = 0.
         yi(posk2) = 0./yi(posk1)   ! NaN
         posk2 = posk2-1
      end do

      !Linear interpolation
      posi = 0
      do i=posk1,posk2
         do while (zi(i)>z(posi+1))
            posi = posi+1
         end do
         yi(i)=y(posi)+(zi(i)-z(posi))*(y(posi+1)-y(posi))/(z(posi+1)-z(posi))
      end do

      return
end


!Integrate discrete function y[x] using the trapezoidal rule
!####################################################################
subroutine Integrate(x,y,inty,num)
!####################################################################

      implicit none
      include 'common_parameters.i'

      double precision x(1:num), y(1:num), inty(1:num)
      integer num, i

      inty(1) = 0
      do i=2,num
         inty(i) = inty(i-1) + (x(i)-x(i-1))*(y(i)+y(i-1))/2
      end do

      return
end


!Compute appropriate forcing parameters at given datum
!AG 2014: revision
!######################################################################
subroutine Forcing(datum,T,idx,tx,ty,u_taus,rad0,heat,SST,u10,v10,wind,Wf)
!######################################################################

      implicit none
      include 'common_parameters.i'

      ! Global iables
      double precision datum,T(0:xl),tx,ty,u_taus
      double precision rad0,heat,SST,u10,v10,wind,Wf
      integer idx

      ! Local iables
      double precision tau
      double precision A_s(7), A_e(7), A_cur(7)
      double precision r_a, B0, sig, fu, Vap_wat, heat0
      double precision T_atm, F_glob, Vap_atm, Cloud
      double precision H_A, H_K, H_V, H_W

      save A_s, A_e

      if (NBC==1) then
          call ReadForcing(datum,A_s,A_e,A_cur,4+WindFilt,idx)
          u10 = A_cur(1)*f_wind       !MS 2014: added f_wind
          v10 = A_cur(2)*f_wind       !MS 2014: added f_wind
          wind = sqrt(u10**2+v10**2)  !AG 2014
          SST = A_cur(3) !Sea surface temperature
          rad0 = A_cur(4)*(1-albsw)*(1-beta_sol) ! MS: added beta_sol and albsw
          heat = 0.
          if(WindFilt==1) Wf = A_cur(5)  !AG 2014
      else if (NBC>=2) then
          if (NBC==2) then            ! date, U,V,Tatm,Hsol,Vap
              call ReadForcing(datum,A_s,A_e,A_cur,5+WindFilt,idx)
              u10 = A_cur(1)*f_wind      !MS 2014: added f_wind
              v10 = A_cur(2)*f_wind      !MS 2014: added f_wind
              T_atm = A_cur(3)
              F_glob = A_cur(4)*(1-albsw)
              Vap_atm = A_cur(5)
              Cloud = 0.5
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
              r_a = 0.03           ! Ratio of reflected to total long-wave iradiance
              !albsw = 0.20        ! Fraction of reflected short-wave radiation (MS: made this a parameter)
              B0 = 0.61            ! Bowen constant
              sig = 5.67e-8        ! Stefan-Boltzmann constant [W/m2/K4]

              ! Wind function (Livingstone & Imboden 1989)
              fu = 4.40+1.82*wind+0.26*(T(xl)-T_atm)
              !fu = 5.44+2.19*wind+0.24*(T(xl)-T_atm)
              fu = fu*p_windf    ! Provided fitting factor p_windf (~1)
              ! Water vapor saturation pressure in air at water temperature (Gill 1992) [millibar]
              Vap_wat = 10**((0.7859+0.03477*T(xl))/(1+0.00412*T(xl)))
              Vap_wat = Vap_wat*(1+1e-6*p_air*(4.5+0.00006*T(xl)**2))
              ! Solar short-wave radiation absorbed
              rad0 = F_glob*(1-beta_sol) !MS: added beta_sol

              ! Long-wave radiation from sky (Livingstone & Imboden 1989)
              H_A = 1.24*sig*(1-r_a)*(1+0.17*Cloud**2)*(Vap_atm/(273.15+T_atm))**(1./7)*(273.15+T_atm)**4
              H_A = H_A*p_radin    ! Provided fitting factor p_radin (~1)
              ! Long-wave radiation from water body (black body)
              H_W = -0.97*sig*(T(xl)+273.15)**4
              ! Flux of sensible heat (convection)
              H_K = -B0*fu*(T(xl)-T_atm)
              ! Flux of latent heat (evaporation, condensation)
              H_V = -fu*(Vap_wat-Vap_atm)
              ! Global heat flux (positive: air to water, negative: water to air)
              heat = H_A + H_W + H_K + H_V + F_glob*beta_sol !MS: added term with beta_sol
          else
              heat = heat0 + F_glob*beta_sol !MS: added term with beta_sol
          end if
          if ((T(xl)<0).and.(heat<0)) heat=0.
      end if

      !Drag coefficient as a function of wind speed (AG 2014)
      if (ModC10==2) then !Ocean model
          C10 = -0.000000712*wind**2+0.00007387*wind+0.0006605
      else if (ModC10==3) then !Lake model (Wüest and Lorke 2003)
          if (wind<=0.1) then
              C10 = 0.06215
          else if (wind<=3.85) then
              C10 = 0.0044*wind**(-1.15)
          else !Polynomial approximation of Charnock's law
              C10 = -0.000000712*wind**2+0.00007387*wind+0.0006605
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
end


!Read forcing file to get values A_cur at given datum
!AG 2014: revision + correction
!####################################################################
subroutine ReadForcing(datum,A_s,A_e,A_cur,nval,idx)
!####################################################################

      implicit none
      include 'common_parameters.i'

      double precision datum, A_cur(1:nval) !Required date, output values
      double precision tb_s, tb_e !Start and end time
      double precision A_s(1:nval), A_e(1:nval) !Start and end values
      integer idx,eof,i,nval

      save tb_s, tb_e
      save eof

      if (idx==1) then
         if(disp_dgn/=0) write(6,*) 'Forcing : ',trim(ForcingName)
         if(disp_dgn==2) write(6,*) 'Starting to read forcing file...'
         open(20,status='old',file=ForcingName)
         eof = 0
         !Read first values
         read(20,*,end=9)
         read(20,*,end=9) tb_s, (A_s(i),i=1,nval)
         if(datum<tb_s) write(6,*) 'Warning: first forcing date after simulation start time.'
         read(20,*,end=7) tb_e, (A_e(i),i=1,nval)
      end if

      if (datum<=tb_s .or. eof==1) then !If datum before first date or end of file reached
          goto 8
      else
          do while (datum>tb_e)         !Move to appropriate interval to get correct values
              tb_s = tb_e
              A_s(1:nval) = A_e(1:nval)
              !Read next value
              read(20,*,end=7) tb_e, (A_e(i),i=1,nval)
          end do
          !Linearly interpolate values at correct datum
          A_cur(1:nval) = A_s(1:nval) + (datum-tb_s) * (A_e(1:nval)-A_s(1:nval))/(tb_e-tb_s)
      end if
      return

7     eof = 1
      if(datum>tb_s) write(6,*) 'Warning: last forcing date before simulation end time.'
8     A_cur(1:nval) = A_s(1:nval)       !Take first value of current interval
      return

9     write(6,*) 'Error reading forcing file (no data found).'
      stop

end


!AG 2014: revision + correction
!###############################################################
subroutine Absorption(datum,ga1,zk,idx)
!###############################################################

      implicit none
      include 'common_parameters.i'

      double precision zk(0:xl) !Required depths
      double precision datum, ga1(0:xl) !Required date, output values
      double precision tb_s, tb_e !Start and end time
      double precision z_ga1(0:mxl), dummy !Read depths
      double precision ga1_rs(0:mxl), ga1_re(0:mxl) !Read start and end values
      double precision ga1_s(0:mxl), ga1_e(0:mxl) !Interpolated start and end values
      integer idx,eof,i,num_z

      save tb_s, tb_e              ! MS 2014
      save z_ga1, ga1_s, ga1_e     ! MS 2014
      save eof, num_z              ! MS 2014

      if (idx==1) then
          if(disp_dgn/=0) write(6,*) 'Light attenuation : ',trim(AbsorpName)
          if(disp_dgn==2) write(6,*) 'Starting to read light attenuation file...'
          open(30,status='old',file=AbsorpName)
          eof = 0
          !Read depths
          read(30,*,end=9)
          read(30,*,end=9) num_z
          read(30,*,end=9) dummy, (z_ga1(i),i=0,num_z-1)
          !Make depths positives
          do i=0,num_z-1
             z_ga1(i) = abs(z_ga1(i))
          end do
          !Read first values
          read(30,*,end=9) tb_s, (ga1_rs(i),i=0,num_z-1)
          if(datum<tb_s) write(6,*) 'Warning: first light attenuation date after simulation start time.'
          call Interp(z_ga1, ga1_rs, num_z-1, zk, ga1_s, xl)
          read(30,*,end=7) tb_e, (ga1_re(i),i=0,num_z-1)
          call Interp(z_ga1, ga1_re, num_z-1, zk, ga1_e, xl)
      end if

      if (datum<=tb_s .or. eof==1) then !If datum before first date or end of file reached
          goto 8
      else
          do while (datum>tb_e)         !Move to appropriate interval to get correct value
              tb_s = tb_e
              ga1_s(0:xl) = ga1_e(0:xl)
              !Read next value
              read(30,*,end=7) tb_e, (ga1_re(i),i=0,num_z-1)
              call Interp(z_ga1, ga1_re, num_z-1, zk, ga1_e, xl)
          end do
          !Linearly interpolate value at correct datum (for all depths)
          ga1(0:xl) = ga1_s(0:xl) + (datum-tb_s) * (ga1_e(0:xl)-ga1_s(0:xl))/(tb_e-tb_s)
      end if
      return

7     eof = 1
      if(datum>tb_s) write(6,*) 'Warning: last light attenuation date before simulation end time.'
8     ga1(0:xl) = ga1_s(0:xl)           !Take first value of current interval
      return

9     write(6,*) 'Error reading light attenuation file (no data found).'
      stop
end


!Read and compute inflow/ouflow parameters
!AG 2014: revision + correction
!####################################################################
subroutine Lateral(datum,idx,zu,zk,z_zero,Qvert,Q_inp)
!####################################################################

      implicit none
      include 'common_parameters.i'

      ! Global Declarations
      double precision datum,zu(0:xl),zk(0:xl),z_zero,Qvert(0:xl),Q_inp(1:4,0:xl)
      integer idx

      ! Local Declarations
      integer i, j, xlp, num_z(1:4), fnum(1:4), eof(1:4)
      character*20 fname(1:4)
      double precision z_Inp(1:4,mxl), dummy
      double precision Inp_rs(1:4,mxl), Inp_re(1:4,mxl)
      double precision Q_s(1:4,mxl), Q_e(1:4,mxl)
      double precision Q_rs(1:4,mxl), Q_re(1:4,mxl)
      double precision tb_s(1:4), tb_e(1:4)

      save Q_s, Q_e         ! MS 2014
      save xlp, num_z, eof  ! MS 2014
      save tb_s, tb_e       ! MS 2014

      fname = ['inflow           ','outflow          ','input temperature','input salinity   ']
      fnum = [41,42,43,44]
      do i=1,4
         if (idx==1) then
            if(disp_dgn==2) write(6,*) 'Starting to read '//trim(fname(i))//' file...'
            eof(i) = 0
            read(fnum(i),*,end=9)      ! Skip first row: description of columns
            !Read input depths and convert coordinate system
            read(fnum(i),*,end=9) num_z(i)
            read(fnum(i),*,end=9) dummy, (z_Inp(i,j),j=1,num_z(i))
            z_Inp(i,1:num_z(i)) = z_zero + z_Inp(i,1:num_z(i))
            !Read first values
            read(fnum(i),*,end=9) tb_s(i),(Inp_rs(i,j),j=1,num_z(i))
            call Integrate(z_Inp(i,1:num_z(i)),Inp_rs(i,1:num_z(i)),Q_rs(i,1:num_z(i)),num_z(i))
            call Interp(z_Inp(i,1:num_z(i)),Q_rs(i,1:num_z(i)),num_z(i)-1,zk(1:xl),Q_s(i,1:xl),xl-1)

            read(fnum(i),*,end=7) tb_e(i),(Inp_re(i,j),j=1,num_z(i))
            call Integrate(z_Inp(i,1:num_z(i)),Inp_re(i,1:num_z(i)),Q_re(i,1:num_z(i)),num_z(i))
            call Interp(z_Inp(i,1:num_z(i)),Q_re(i,1:num_z(i)),num_z(i)-1,zk(1:xl),Q_e(i,1:xl),xl-1)
            xlp = xl
         end if

         if (xlp<xl) then
             Q_s(i,xlp+1:xl) = Q_s(i,xlp)
             Q_e(i,xlp+1:xl) = Q_e(i,xlp)
         end if

         if ((datum<=tb_s(i)).or.(eof(i)==1)) then       ! if datum before first date or end of file reached
             goto 8
         else
             do while (.not.((datum>=tb_s(i)).and.(datum<=tb_e(i)))) ! do until datum between dates
                 tb_s(i) = tb_e(i)             ! move one step in time
                 Q_s(i,1:xl) = Q_e(i,1:xl)
                 read(fnum(i),*,end=7) tb_e(i),(Inp_re(i,j),j=1,num_z(i))
                 call Integrate(z_Inp(i,1:num_z(i)),Inp_re(i,1:num_z(i)),Q_re(i,1:num_z(i)),num_z(i))
                 call Interp(z_Inp(i,1:num_z(i)),Q_re(i,1:num_z(i)),num_z(i)-1,zk(1:xl),Q_e(i,1:xl),xl-1)
             end do
             !do j=xl+1,mxl !(AG 2014:commented) useful???? (+previously not done at first timestep)
                 !Q_e(i,j) = Q_e(i,xl)
                 !Q_s(i,j) = Q_s(i,xl)
             !end do
             !Linearly interpolate value at correct datum (for all depths)
             if(tb_e(i)<=tb_s(i)) then
                write(6,*) 'Error: dates in ',trim(fname(i)),' file must always be increasing.'
                stop
             end if
             do j=1,xl
                 Q_inp(i,j) = Q_s(i,j) + (datum-tb_s(i)) * (Q_e(i,j)-Q_s(i,j))/(tb_e(i)-tb_s(i))
             end do
         end if
         goto 11

7        eof(i) = 1
8        Q_inp(i,1:xl) = Q_s(i,1:xl)              ! Set to closest available value
         goto 11

9        write(6,*) 'No data found in ',trim(fname(i)),' file. Check number of depths. Values set to zero.'
         eof(i) = 1
         Q_inp(i,0:xl) = 0.
         Q_s(i,1:xl) = 0.

11       continue
      end do      ! end do i=1,4

      Qvert(1:xl) = Q_inp(1,1:xl)+Q_inp(2,1:xl)
      !Set all Q to the differences (from the integrals)
      do i=1,4
          do j=1,xl-1
              Q_inp(i,xl-j+1) = Q_inp(i,xl-j+1)-Q_inp(i,xl-j)
          end do
      end do
      xlp = xl

      !AG 2014: commented this, as the shift causes violation of salinity and temperature balances
      !do i=1,xl-1                               ! Outflow always at bottom, otherwise
          !Q_inp(2,i) = Q_inp(2,i+1)             ! vertical advection and outflow not
      !end do                                    ! consistent (preliminary version!)

      return
end

!Read first column of inflow/outflow files and compute parameters
!Assuming inflow will plunge according to its density, entraining water from other layers
!####################################################################
subroutine Lateral_rho(datum,idx,zu,zk,z_zero,h,T,S,rho,Qvert,Q_inp)
!####################################################################

      implicit none
      include 'common_parameters.i'

      ! Global Declarations
      double precision datum,zu(0:xl),zk(0:xl),z_zero,h(0:xl),T(0:xl),S(0:xl),rho(0:xl)
      double precision Qvert(0:xl),Q_inp(1:4,0:xl)
      integer idx

      ! Local Declarations
      integer i, j, k, i1, i2, num_z(1:4), fnum(1:4), eof(1:4)
      character*20 fname(1:4)
      double precision z_Inp(1:4,mxl), dummy
      double precision Inp_rs(1:4,mxl), Inp_re(1:4,mxl), Inp(1:4)
      double precision Q_s(1:4,mxl), Q_e(1:4,mxl)
      double precision Q_rs(1:4,mxl), Q_re(1:4,mxl)
      double precision tb_s(1:4), tb_e(1:4)
      double precision T_in, S_in, rho_in
      double precision slope, hang, CD_in, Ri, g_red, E
      double precision h_in(0:xl), Q_in(0:xl)

      save Inp_rs, Inp_re, Q_s, Q_e
      save num_z, eof
      save tb_s, tb_e

      fname = ['inflow           ','outflow          ','input temperature','input salinity   ']
      fnum = [41,42,43,44]
      do i=1,4
         if (idx==1) then
            if(disp_dgn==2) write(6,*) 'Starting to read '//trim(fname(i))//' file...'
            eof(i) = 0
            !Skip first three rows (except for Qout)
            read(fnum(i),*,end=9)
            read(fnum(i),*,end=9) num_z(i)
            read(fnum(i),*,end=9) dummy, (z_Inp(i,j),j=1,num_z(i))
            if(i/=2) num_z(i)=1
            if(i==2) z_Inp(i,1:num_z(i)) = z_zero + z_Inp(i,1:num_z(i))
            !Read first values
            read(fnum(i),*,end=9) tb_s(i),(Inp_rs(i,j),j=1,num_z(i))
            if(i==2) call Integrate(z_Inp(i,1:num_z(i)),Inp_rs(i,1:num_z(i)),Q_rs(i,1:num_z(i)),num_z(i))
            if(i==2) call Interp(z_Inp(i,1:num_z(i)),Q_rs(i,1:num_z(i)),num_z(i)-1,zk(1:xl),Q_s(i,1:xl),xl-1)
            read(fnum(i),*,end=7) tb_e(i),(Inp_re(i,j),j=1,num_z(i))
            if(i==2) call Integrate(z_Inp(i,1:num_z(i)),Inp_re(i,1:num_z(i)),Q_re(i,1:num_z(i)),num_z(i))
            if(i==2) call Interp(z_Inp(i,1:num_z(i)),Q_re(i,1:num_z(i)),num_z(i)-1,zk(1:xl),Q_e(i,1:xl),xl-1)
         end if

         if ((datum<=tb_s(i)).or.(eof(i)==1)) then       ! if datum before first date or end of file reached
             goto 8
         else
             do while (.not.((datum>=tb_s(i)).and.(datum<=tb_e(i)))) ! do until datum between dates
                 tb_s(i) = tb_e(i)             ! move one step in time
                 if(i/=2) Inp_rs(i,1) = Inp_re(i,1)
                 if(i==2) Q_s(i,1:xl) = Q_e(i,1:xl)
                 read(fnum(i),*,end=7) tb_e(i),(Inp_re(i,j),j=1,num_z(i))
                 if(i==2) call Integrate(z_Inp(i,1:num_z(i)),Inp_re(i,1:num_z(i)),Q_re(i,1:num_z(i)),num_z(i))
                 if(i==2) call Interp(z_Inp(i,1:num_z(i)),Q_re(i,1:num_z(i)),num_z(i)-1,zk(1:xl),Q_e(i,1:xl),xl-1)
             end do
             !Linearly interpolate value at correct datum
             if (i/=2) then
                 Inp(i) = Inp_rs(i,1) + (datum-tb_s(i)) * (Inp_re(i,1)-Inp_rs(i,1))/(tb_e(i)-tb_s(i))
             else
                 if(tb_e(i)<=tb_s(i)) then
                    write(6,*) 'Error: dates in ',trim(fname(i)),' file must always be increasing.'
                    stop
                 end if
                 do j=1,xl
                     Q_inp(i,j) = Q_s(i,j) + (datum-tb_s(i)) * (Q_e(i,j)-Q_s(i,j))/(tb_e(i)-tb_s(i))
                 end do
             end if
         end if
         goto 11

7        eof(i) = 1
8        if(i/=2) Inp(i) = Inp_rs(i,1)            ! Set to closest available value
         if(i==2) Q_inp(i,1:xl) = Q_s(i,1:xl)     ! Set to closest available value
         goto 11

9        write(6,*) 'No data found in ',trim(fname(i)),' file. Check number of depths. Values set to zero.'
         eof(i) = 1
         if(i/=2) Inp(i) = 0.
         if(i/=2) Inp_rs(i,1) = 0.
         if(i==2) Q_inp(i,0:xl) = 0.
         if(i==2) Q_s(i,1:xl) = 0.

11       continue
      end do      ! end do i=1,4

      !Set Qout to the differences (from the integrals)
      do j=1,xl-1
          Q_inp(2,xl-j+1) = Q_inp(2,xl-j+1)-Q_inp(2,xl-j)
      end do

      if (Inp(1)>1E-15) then
          slope = pi/72 !Slope of inflow
          hang = pi/3 !Stream half-angle
          Q_in(xl) = Inp(1) !Inflow flow rate [m3/s]
          T_in = Inp(3) !Inflow temperature [°C*m3/s]
          S_in = Inp(4) !Inflow salinity [‰*m3/s]
          rho_in = rho_0*(0.9998395+T_in*(6.7914e-5+T_in*(-9.0894e-6+T_in*&
                   (1.0171e-7+T_in*(-1.2846e-9+T_in*(1.1592e-11+T_in*(-5.0125e-14))))))+&
                   (8.181e-4+T_in*(-3.85e-6+T_in*(4.96e-8)))*S_in) !Inflow density [kg/m3]
          g_red = g*(rho_in-rho(xl))/rho_in !Reduced gravity [m/s2]
          if (g_red>0) then
              CD_in = CD*10 !Inflow drag coefficient
              Ri = CD_in*(1+0.21*CD_in*sin(hang))/(sin(hang)*tan(slope)) !Richardson number
              E = 1.6*CD_in**1.5/Ri !Entrainment coefficient
              h_in(xl) = (2*Inp(1)**2*Ri*tan(slope)**2/g_red)**0.2 !Inflow thickness [m]
              do k=xl,1,-1
                  if(rho_in<=rho(k)) exit
                  h_in(k-1) = 1.2*E*(zu(k)-zu(k-1))/sin(slope) + h_in(k)
                  Q_in(k-1) = Q_in(k)*(h_in(k-1)/h_in(k))**(5./3.)
                  Q_inp(2,k) = Q_inp(2,k) - (Q_in(k-1)-Q_in(k))
                  T_in = (T_in*Q_in(k)+T(k)*(Q_in(k-1)-Q_in(k)))/Q_in(k-1)
                  S_in = (S_in*Q_in(k)+S(k)*(Q_in(k-1)-Q_in(k)))/Q_in(k-1)
                  rho_in = (rho_in*Q_in(k)+rho(k)*(Q_in(k-1)-Q_in(k)))/Q_in(k-1)
              end do
          else
              k=xl
          end if

          if (k==xl) then !surface flow
              i1=xl
              i2=xl-1
          else
              do i1=k,xl !extend upwards
                  if(zu(i1)>zu(k)+h_in(k)/2) exit
              end do
              do i2=k,1,-1 !extend downwards
                  if(zu(i2)<zu(k)-h_in(k)/2) exit
              end do
          end if
          i1 = i1-1
          i2 = i2-1

          Q_inp(1,:) = 0
          Q_inp(3:4,:) = 0
          do i=i1,i2+1,-1
              Q_inp(1,i) = Q_in(k)/(zk(i1)-zk(i2))*h(i)
              Q_inp(3,i) = T_in*Q_inp(1,i)
              Q_inp(4,i) = S_in*Q_inp(1,i)
          end do
      end if

      Qvert(1) = Q_inp(1,1)+Q_inp(2,1)
      do i=2,xl
          Qvert(i) = Qvert(i-1)+Q_inp(1,i)+Q_inp(2,i)
      end do

      return
end
