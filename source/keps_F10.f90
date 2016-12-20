!     +---------------------------------------------------------------+
!     |  k-epsilon model for simulation of                            |
!     |  vertical transport in reservoirs                             |  
!     +---------------------------------------------------------------+

      program keps

!     +---------------------------------------------------------------+
!     |                                                               |
!     |  Declarations of constants and variables                      |
!     |                                                               |
!     +---------------------------------------------------------------+

      implicit none

!     +---------------------------------------------------------------+
!     |  Constant Declarations                                        |
!     +---------------------------------------------------------------+

      include 'const_parameter.i'

      double precision h(0:mxl)
      double precision Az(0:mxl),dAdz(0:mxl)
      double precision form_1(0:mxl),form_2(0:mxl)
      double precision form_k1(0:mxl),form_k2(0:mxl),form_beps(0:mxl)
      double precision volume

      common /morph_variables1/ h,Az,dAdz
      common /morph_variables3/ volume
      common /form_variables  / form_1,form_2
      common /form_variablesk/form_k1,form_k2,form_beps
      common /form_variables3/meanint

!     +---------------------------------------------------------------+
!     |  Variable Declarations                                        |
!     +---------------------------------------------------------------+

      double precision M(0:mxl),meanint(0:mxl)
	  integer i
!     +---------------------------------------------------------------+
!     |  Constant specification                                       |
!     +---------------------------------------------------------------+

      call Initialization
	  call Form

  
!      write(6,*)
!      write(6,*) " ------------------------ "
!      write(6,*) " INITIALIZATION SUCCESSFUL"
!      write(6,*) " ------------------------ "
!      write(6,*)
       call keps_simulation(M)
	  
	  stop
	  end

!     +---------------------------------------------------------------+
!     |   Main loop                                                   |
!     +---------------------------------------------------------------+

      subroutine keps_simulation(M)
	  
      implicit none

      include 'const_parameter.i'


      double precision Sini(0:mxl),Tini(0:mxl),uini(0:mxl),vini(0:mxl)
      double precision kini(0:mxl),epsini(0:mxl),Lini(0:mxl)
      double precision numini(0:mxl),nuhini(0:mxl) 
      double precision E_Seicheini,datumini
      double precision dragini,txini,tyini

      double precision zu(0:mxl),zk(0:mxl),h(0:mxl), z_zero
      double precision Az(0:mxl),dAdz(0:mxl)
      double precision form_1(0:mxl),form_2(0:mxl)
      double precision form_k1(0:mxl),form_k2(0:mxl),form_beps(0:mxl)
      double precision meanint(0:mxl)
      double precision volume

      double precision tout_ctr1(0:9000),tout_ctr2(0:9000)
      integer write_tout

      common /ini_values1/ Sini,Tini
      common /ini_values2/ uini,vini
      common /ini_values3/ epsini,kini,Lini
      common /ini_values4/ numini,nuhini
      common /ini_values5/ dragini,txini,tyini,E_Seicheini,datumini


      common /morph_variables1/ h,Az,dAdz
      common /morph_varaibles2/ zu,zk,z_zero
      common /morph_variables3/ volume
      common /form_variables  / form_1,form_2
      common /form_variablesk/form_k1,form_k2,form_beps
      common /form_variables3/meanint

      common /savet1  /   write_tout
      common /savet2  /   tout_ctr1
      common /savet3  /   tout_ctr2

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    LOCAL VARIABLES
     
      integer step,itera,i,std

      double precision P(0:mxl),B(0:mxl),NN(0:mxl)
      double precision S(0:mxl),T(0:mxl),u(0:mxl),v(0:mxl)

	  double precision Qvert(0:mxl), Q_inp(1:4,0:mxl)
	  
      double precision cmue1(0:mxl),cmue2(0:mxl)
      double precision k(0:mxl),ko(0:mxl)
      double precision eps(0:mxl),L(0:mxl)
      double precision num(0:mxl),nuh(0:mxl) 
      double precision ga1(0:mxl),u10,v10
      double precision E_Seiche, P_Seiche(0:mxl)
      
      double precision gamma

      double precision SST,datum
      double precision u_taub,drag,u_taus
      double precision tx,ty,I_0,heat

      double precision M(0:mxl)

      double precision uav(0:mxl),vav(0:mxl)
      double precision Tav(0:mxl),Sav(0:mxl)
      double precision kav(0:mxl)
      double precision epsav(0:mxl)
      double precision numav(0:mxl),nuhav(0:mxl)
      double precision Bav(0:mxl),Pav(0:mxl),NNav(0:mxl)
      double precision P_Seicheav(0:mxl),E_Seicheav


      integer iav


      integer date_time (8)
      character (len = 12) real_clock (3)

!     +---------------------------------------------------------------+
!     |  Constant specification                                       |
!     +---------------------------------------------------------------+

      iav=0
!    Set initial values (calculated in Initialization)
      do i=0,xl
		  u(i)=uini(i)
		  v(i)=vini(i)
		  k(i)=kini(i)
		  eps(i)=epsini(i)
		  L(i)=Lini(i)
		  T(i)=Tini(i)
		  S(i)=Sini(i)
		  num(i)=numini(i)
		  nuh(i)=nuhini(i)
      end do
      drag=dragini
      tx=txini
      ty=tyini
      E_Seiche=E_Seicheini
      datum=datumini
	  
      gamma=Az(xl)/(volume**1.5)/ sqrt(rho_0) *CDeff

      std=0
      step=0  
      itera=0

      if (write_tout.eq.0) then
             if (tout_ctr2(0).eq.0) then           ! first output = initial time

	        if (disp_simulation .eq. 1)  write(6,990) datum,T(xl),T(xl-5)
                call write_out (datum,u,v,T,S,k,eps,num,nuh,
     &                        B,P,NN,P_Seiche,E_Seiche,zu,zk,M)

                itera=0
                step=1
             end if
      end if

	  call write_out_new (datum,std, u,v,T,S,k,eps,nuh,					! Write initial conditions
     &                        B,P,NN,P_Seiche,E_Seiche,zu,zk,M,Qvert)	! preliminary version!!!
	  
900   std=std+1

      itera=itera+1

      if (write_tout.eq.0) then
              if (itera.eq.1) then  
                   dt=tout_ctr2(step)
	      end if
              if ( itera.eq.int(tout_ctr1(step)) ) then
                   itera=0
                   step=step+1
              end if
      end if

       datum = datum + dt/86400.0
         call  Forcing(datum,T,std,tx,ty,u_taus,I_0,heat,SST,u10,v10)
         call  Absorption(datum, ga1, zk, std)
		 if (adv == 1) then
		    call Lateral(datum,std,zu,z_zero,h,Qvert,Q_inp)
			call Advection(Qvert,Q_inp,u,v,T,S,k,eps,num,nuh,zu,zk,h,Az)
		 end if
         u_taub=sqrt(drag*(u(1)*u(1)+v(1)*v(1)))

         call StabilityFunctions(u,v,k,eps,T,S,meanint,NN,
     &                       cmue1,cmue2)
         call Coriolis(u,v)
         call uvEquation(u,v,num,h,Az,drag,tx,ty,dAdz,
     &                   form_1,form_2)

         call Temperature(nuh,datum,I_0,h,Az,T,heat,SST,ga1,
     &                    dAdz,form_1,form_2)

         call Salinity(S,h,Az,nuh,form_1,form_2)

         call Production(u,v,NN,meanint,num,nuh,P,B)
         call Seiche(E_Seiche,P_Seiche,u10,v10,u,v,Az,dAdz,
     &                NN,h,tx,ty,gamma)
         call TKE(num,P,B,eps,L,h,Az,u_taus,u_taub,k,ko,
     &                P_Seiche,form_k1,form_k2)
         call Dissipation(cmue1,cmue2,P,B,k,ko,h,Az,eps,L,num,nuh,NN,
     &                u_taus,u_taub,P_Seiche,form_k1,form_k2,form_beps)
	     !call date_and_time (real_clock (1), real_clock (2), real_clock (3), date_time)
         !write(*,*) real_clock(2)

!
! check if averaging and average data
!
        if (igoal.lt.0) then
		   iav=iav+1
                 call  avstate(iav,u,v,T,S,k,eps,num,nuh,
     &                        B,P,NN,P_Seiche,E_Seiche,
     &                        uav,vav,Tav,Sav,kav,epsav,numav,nuhav,
     &                        Bav,Pav,NNav,P_Seicheav,E_Seicheav)
         end if
!
! write results to file
!
        if ((itera.eq.write_tout).or.(datum.ge.t_end)) then                         
           
           itera=0
	   
           if (igoal.lt.0) then
                  call  prep_outav(iav,   
     &                        uav,vav,Tav,Sav,kav,epsav,numav,nuhav,
     &                        Bav,Pav,NNav,P_Seicheav,E_Seicheav)
                  call write_out_new (datum,std,uav,vav,Tav,Sav,kav,epsav,
     &                        nuhav,Bav,Pav,NNav,P_Seicheav,
     &                        E_Seicheav,zu,zk,M,Qvert)

	      if (disp_simulation .eq. 1) write(6,991) datum,T(xl),T(xl-5),
     &                                   Tav(xl),Tav(xl-5)
              iav=0

	   else
	      if (disp_simulation .eq. 1) write(6,990) datum,T(xl),T(xl-1),zk(xl)

              call write_out_new (datum,std, u,v,T,S,k,eps,nuh,
     &                        B,P,NN,P_Seiche,E_Seiche,zu,zk,M,Qvert)
           end if
        end if


  	if (disp_simulation .eq. 2)  then
                 write(6,990) datum,T(xl),T(xl-1)
        elseif (disp_simulation .eq. 3) then
              write(6,999) datum,T(xl),T(xl-2),u(xl),u(xl-2),
     &                           eps(xl),eps(xl-1),num(xl),num(xl-2)
        end if
 
      if (.not.((datum).ge.(t_end))) goto 900		! Besser: do while Schlaufe

      close(20)
      close(80)

 990     format (F12.3, F12.5, F12.5, F12.5)
 991     format (F12.3, F12.5, F12.5, F12.5, F12.5)
 999     format (F12.5,F12.5,F12.5,F12.5,F12.5,E12.5,E12.5,F12.5,F12.5)
 700     format(F5.2, F5.2)  

      stop
      end
 

!     +---------------------------------------------------------------+
!     |  Subroutines                                                  |
!     +---------------------------------------------------------------+


!     ###############################################################
      subroutine Coriolis(u,v)
!     ###############################################################

      implicit none

      include 'const_parameter.i'

!     +-------------------------------------------------------------+
!     |   Global variables                                          |
!     +-------------------------------------------------------------+
      double precision u(0:xl),v(0:xl)

!     +-------------------------------------------------------------+
!     |   Local variables                                           |
!     +-------------------------------------------------------------+
      integer i
      double precision    ua

!     +-------------------------------------------------------------+
!     |    Calculation                                              |
!     +-------------------------------------------------------------+

	  do i=1,xl
         ua=u(i)
         u(i)=u(i)*cos(Cori*dt)+v(i)*sin(Cori*dt)
         v(i)=-ua*sin(Cori*dt)+v(i)*cos(Cori*dt)
      end do

      return
      end

!     ###############################################################
      subroutine uvEquation(u,v,num,h,Az,drag,tx,ty,dAdz,
     &                   form_1,form_2)
!     ###############################################################

      implicit none

      include 'const_parameter.i'

!     +-------------------------------------------------------------+
!     |   Global variables                                          |
!     +-------------------------------------------------------------+
      double precision u(0:xl),v(0:xl),h(0:xl),Az(0:xl),dAdz(0:xl)
      double precision num(0:xl)
      double precision drag,tx,ty
      double precision form_1(0:xl),form_2(0:xl)
!     +-------------------------------------------------------------+
!     |   Local variables                                           |
!     +-------------------------------------------------------------+
      double precision au(0:mxl),bu(0:mxl),cu(0:mxl),du(0:mxl)
      integer i
      double precision intu, intv, length
      
	  if (pgrad == 1) then
	     intu = 0
		 intv = 0
	     do i=1,xl
            intu = intu+u(i)
			intv = intv+v(i)
         end do
		 length = sqrt(Az(xl))
      end if

!     +-------------------------------------------------------------+
!     |    Calculation of u-equation                                |
!     +-------------------------------------------------------------+
      do i=2,xl-1
         cu(i)=dt*num(i)*form_2(i)
         au(i)=dt*num(i-1)*form_1(i)
         bu(i)=1.-au(i)-cu(i)
		 if (pgrad == 1) then
		    du(i) = u(i)  - 3.14*31.4*intu/xl*dt/86400*1000*9.81*depth/length/length
         elseif (pgrad == 2) then
            du(i) = u(i) - drag*u(i)*sqrt(u(i)*u(i)+v(i)*v(i))
     &           *dAdz(i)/Az(i)*dt
         else
            du(i) = u(i)
         end if
      end do

      cu(1)=dt*num(1)*form_2(1)
      bu(1)=1.-cu(1)+drag*dt/h(1)*sqrt(u(1)*u(1)+v(1)*v(1))
      du(1) = u(1)

      au(xl)=dt*num(xl-1)*form_1(xl)
      bu(xl)=1.-au(xl)
      du(xl) = u(xl) + tx*dt/h(xl)

      call Tridiagonal(1,xl,au,bu,cu,du,u)

!     +-------------------------------------------------------------+
!     |    Calculation of v-equation                                |
!     +-------------------------------------------------------------+
      do i=2,xl-1
		 if (pgrad == 1) then
		    du(i) = v(i)  - 3.14*31.4*intv/xl*dt/86400*1000*9.81*depth/length/length
         elseif (pgrad.eq.2) then
            du(i) = v(i) - drag*v(i)*sqrt(u(i)*u(i)+v(i)*v(i))
     &           *dAdz(i)/Az(i)*dt
         else
            du(i) = v(i)
         end if
      end do
      du(1) = v(1)
      du(xl) = v(xl) + ty*dt/h(xl)

      call Tridiagonal(1,xl,au,bu,cu,du,v)

      return
      end

!     ###############################################################
      subroutine Tridiagonal(fi,lt,au,bu,cu,du,value)
!     ###############################################################

      implicit none

      include 'const_parameter.i'


!     +----------------------------------------------------------------+
!     | Global variables                                               |
!     +----------------------------------------------------------------+
      double precision au(0:xl),bu(0:xl),cu(0:xl),du(0:xl),
     &                 value(0:xl)
      integer fi,lt

!     +----------------------------------------------------------------+
!     | Local variables                                                |
!     +----------------------------------------------------------------+
      double precision ru(0:mxl),qu(0:mxl)
      integer i

      ru(lt)=au(lt)/bu(lt)
      qu(lt)=du(lt)/bu(lt)

      do i=lt-1,fi+1,-1
         ru(i)=au(i)/(bu(i)-cu(i)*ru(i+1))
         qu(i)=(du(i)-cu(i)*qu(i+1))/(bu(i)-cu(i)*ru(i+1))
      end do

      qu(fi)=(du(fi)-cu(fi)*qu(fi+1))/(bu(fi)-cu(fi)*ru(fi+1))

      value(fi)=qu(fi)
      do i=fi+1,lt
         value(i)=qu(i)-ru(i)*value(i-1)
      end do

!     +----------------------------------------------------------------+
!     | End of subroutine Tridiagonal                                  |
!     +----------------------------------------------------------------+
      return
      end



!     +----------------------------------------------------------------+
!     | Equation of state                                              |
!     +----------------------------------------------------------------+
!     ###############################################################
      subroutine Buoyancy(T,S,meanint,NN)
!     ###############################################################

      implicit none

      include 'const_parameter.i'

!     +----------------------------------------------------------------+
!     | Global variables                                               |
!     +----------------------------------------------------------------+

      double precision T(0:xl), S(0:xl), NN(0:xl),meanint(0:xl)

!     +----------------------------------------------------------------+
!     | Local variables                                                |
!     +----------------------------------------------------------------+
      double precision a(0:mxl),rho(0:mxl),buoy(0:mxl)
      double precision rho0t(0:mxl),rho0st(0:mxl)
 
      integer i

       if  (salctr.eq.0) then 				! salinity is zero everywhere
        do i=1,xl-1
           a(i)= -68.0+T(i)*(18.2091+T(i)*(-0.30866+T(i)*
     &                  (5.3445e-3+T(i)*(-6.0721e-5+T(i)*(3.1441e-7)))))
!  	if (press.ne.0) then   ! ignore this pressure thing for alpha in first approximation
!		a(i)= a(i) + 0.3682 +T(i)*(-1.520e-2 +T(i)*(1.91e-4))).*p(i)
!	end if
           a(i)= 1.0e-6*a(i)
       	   NN(i)=g*a(i)*(meanint(i)*(t(i)-t(i+1))+(t(i)+273.15)/cp)
       end do

       else
	if (delsal.eq.0) then			! salinity gradient is zero everywhere

            do i=1,xl-1
	      a(i)= -68.0+T(i)*(18.2091+T(i)*(-0.30866+T(i)*(5.3445e-3+T(i)*
     &		  (-6.0721e-5+T(i)*(3.1441e-7)))))
	      a(i)= a(i) + (4.599 + T(i)*(-0.1999 +T(i)*(2.79e-3)))*S(i)
!             if (press.ne.0) then   ! ignore this pressure thing for alpha in first approximation
!		 a(i) = a(i) + (0.3682 +T(i)*(-1.520e-2 +T(i)*(1.91e-4))-4.613e-3.*S(i)).*p(i)
!	      end if
	      a(i) = 1.0e-6*a(i)
	      NN(i)=g*a(i)*(meanint(i)*(t(i)-t(i+1))+(t(i)+273.15)/cp)
            end do
         else
!	  if (press.ne.0) then
!		rho0t=(0.9998395+t.*(6.7914e-5 +t.*(-9.0894e-6+t.*(1.0171e-7+t.*(-1.2846e-9 +t.*(1.1592e-11
!                       +t.*(-5.0125e-14)))))));
!		kbart = (19652.17 +t.*(148.113 +t.*(-2.293 +t.*(1.256e-2+t.*(-4.18e-5)))));
!		kbarpt = ((3.2726 +t.*(-2.147e-4 +t.*(1.128e-4))).*p);
!		kbar = kbart+kbarpt
!		rho0st = ((8.181e-4 +t.*(-3.85e-6 +t.*(4.96e-8))).*s);
!		kbarspt = ((53.238 -0.313.*t +5.728e-3.*p).*s);
!		kbar=kbar+kbarspt;        
!		rho = 1000.0*(rho0t + rho0st)./(1.0 - p./kbar);
!			if ~isempty(f)
!				rhoS0 = rho0t./(1.0 - p./(kbart+kbarpt));
!  				rho=rhoS0.*(1-fc)+fc.*rho;
!			end
!		rho = 1000.0*rho0t./(1.0 - p./kbar);
!          else	
             do i=0,xl
		rho0t(i)=(0.9998395+T(i)*(6.7914e-5 +T(i)*(-9.0894e-6+T(i)*
     &		(1.0171e-7+T(i)*
     &           (-1.2846e-9 +T(i)*(1.1592e-11 +T(i)*(-5.0125e-14)))))))
		rho0st(i) = ((8.181e-4 +T(i)*(-3.85e-6 +T(i)*(4.96e-8)))*S(i))

		rho(i)=1000*(rho0t(i)+rho0st(i))
!		if (fc.ne.0) then
!  			rho(i)=rho0t(i)*(1-fc)+fc*rho(i);
!		end
                buoy(i)=-g*(rho(i)-rho_0)/rho_0
            end do

            do i=1,xl-1
              NN(i)=meanint(i)*(buoy(i+1)-buoy(i))
            end do


	  end if


         end if
         NN(0)=NN(1)
         NN(xl)=NN(xl-1)

         return
         end


!     ###############################################################
      subroutine Temperature(nuh,datum,I_0,h,Az,T,heat,SST,ga1,
     &                dAdz,form_1,form_2)
!     ###############################################################

      implicit none
 
      include 'const_parameter.i'


!     +----------------------------------------------------------------+
!     | Global variables                                               |
!     +----------------------------------------------------------------+

      double precision nuh(0:xl),datum,I_0,
     &     h(0:xl),Az(0:xl),T(0:xl),heat,SST, ga1(0:xl),dAdz(0:xl)
      double precision form_1(0:xl),form_2(0:xl)

!     +----------------------------------------------------------------+
!     | Local variables                                                |
!     +----------------------------------------------------------------+

      double precision rad(0:mxl)
      double precision au(0:mxl),bu(0:mxl),cu(0:mxl),du(0:mxl)
      integer i

!     +----------------------------------------------------------------+
!     | Calculation                                                    |
!     +----------------------------------------------------------------+

      rad(xl)=I_0/rho_0/cp
      do i=xl-1,0,-1
         rad(i) = rad(i+1)*exp(-h(i)*ga1(xl-i))
      end do

      do i=2,xl-1
         cu(i)=dt*nuh(i)*form_2(i)
         au(i)=dt*nuh(i-1)*form_1(i)
         bu(i)=1.-au(i)-cu(i)
         du(i)=T(i)+(rad(i)-rad(i-1))/h(i)*dt              
      end do

      cu(1)=dt*nuh(1)*form_2(1)
      bu(1)=1.-cu(1)
      du(1)=T(1)+(rad(1)-rad(0))/h(1)*dt

      au(xl)=dt*nuh(xl-1)*form_1(xl)
      bu(xl)=1.-au(xl)
      du(xl)=T(xl)+(rad(xl)-rad(xl-1))/h(xl)*dt
     &         +heat/rho_0/cp*dt/h(xl)

      if (NBC .eq. 1) then
         bu(xl)=1.
         au(xl)=0
         cu(xl)=0
         du(xl)=SST
      end if


      if (fgeo.ne.0) then
         do i=1,xl
            du(i)=du(i)+fgeo_add(i)*dt
         end do
      end if

      call Tridiagonal(1,xl,au,bu,cu,du,T)    

	T(0)=T(1)     ! set value to boundary value at zu(1)

      return
      end


!     ###############################################################
      subroutine Salinity(S,h,Az,nuh,form_1,form_2)
!     ###############################################################

      implicit none

      include 'const_parameter.i'

      double precision S(0:xl),h(0:xl),Az(0:xl),nuh(0:xl)
      double precision form_1(0:xl),form_2(0:xl)
      double precision au(0:mxl),bu(0:mxl),cu(0:mxl),du(0:mxl)
      integer i

      if (ModSal.ne.0.0) then
         do i=2,xl-1
            cu(i)=dt*nuh(i)*form_2(i)
            au(i)=dt*nuh(i-1)*form_1(i)
!            bu(i)=1.+dt/SalRel-au(i)-cu(i)
	      bu(i)=1.-au(i)-cu(i)
!            du(i)=S(i)+dt/SalRel
		  du(i) = S(i)
         end do

         cu(1)=-4.*dt*nuh(1)*Az(1)/(h(1)+h(2))/h(1)/(Az(1)+Az(0))
!         bu(1)=1.+dt/SalRel-cu(1)
	   bu(1)=1.-cu(1)
!         du(1)=S(1)+dt/SalRel
	   du(1)=S(1)

         au(xl)=-4.*dt*nuh(xl-1)*Az(xl-1)/(h(xl)+h(xl-1))
     &        /h(xl)/(Az(xl)+Az(xl-1))
!         bu(xl)=1.+dt/SalRel-au(xl)
!         du(xl)=S(xl)+dt/SalRel
	   bu(xl)=1.-au(xl)
         du(xl)=S(xl)

!	   if (fsed.ne.0) then		! priliminary test for mineralisation
!	      do i=1,xl
!		     du(i)=du(i)+fsed_add(i)*dt
!		  end do
!	   end if

         call Tridiagonal(1,xl,au,bu,cu,du,S)

	S(0)=S(1)     ! set value to boundary value at zu(1)

      end if

      return
      end


!     +----------------------------------------------------------------+
!     | Calculation of cmue                                            |
!     +----------------------------------------------------------------+

!     ###############################################################
      subroutine cmue_cn(beta,cmue1,cmue2)
!     ###############################################################

      implicit none

      include 'const_parameter.i'

!     +----------------------------------------------------------------+
!     | Global variables                                               |
!     +----------------------------------------------------------------+

      double precision beta,cmue1,cmue2

!     +----------------------------------------------------------------+
!     | Local variables                                                |
!     +----------------------------------------------------------------+
!      cmue1=cde*cm0                 !Burchard Version
!      cmue2=cde*cm0/Prndtl          !  ,,        ,,
      cmue1 = cmue                   ! Standard version of k-eps model
      cmue2 = cmue/Prndtl            ! cmue = 0.09

      return
      end


!     ###############################################################
      subroutine cmue_qe(beta,cmue1,cmue2)
!     ###############################################################

      implicit none
 
      include 'const_parameter.i'
 
      double precision beta,cmue1,cmue2
 
      double precision gh,sm,sh
 
      gh=-cde*cde*0.5*beta
      
      if (gh.gt.0.02) gh=gh-(gh-0.02)**2/(gh+0.0233-2*0.02) 
      if (gh.lt.-0.28) gh = -0.28
      sm=1-3*c1-6*a1/b1-3*a2*gh*((b2-3*a2)*(1-6*a1/b1)-3*c1*(b2+6*a1))
      sm=a1*sm/((1-3*a2*gh*(6*a1+b2))*(1-9*a1*a2*gh))
      sh=a2*(1-6*a1/b1)/(1-3*a2*gh*(6*a1+b2))
      cmue1=sqrt(2.)*cde*sm
      cmue2=sqrt(2.)*cde*sh

      return
      end


!     ###############################################################
      subroutine StabilityFunctions(u,v,k,eps,T,S,meanint,NN,
     &                       cmue1,cmue2)
!     ###############################################################
 
      implicit none
 
      include 'const_parameter.i'
!     +----------------------------------------------------------------+
!     | Global variables                                               |
!     +----------------------------------------------------------------+ 
      double precision u(0:xl),v(0:xl),k(0:xl)
      double precision eps(0:xl),NN(0:xl),cmue1(0:xl)
      double precision cmue2(0:xl),T(0:xl),S(0:xl)
      double precision meanint(0:xl)
	 
!     +----------------------------------------------------------------+
!     | Local variables                                                |
!     +----------------------------------------------------------------+
      integer i
      double precision kkee
      double precision beta(0:mxl)

      call Buoyancy(T,S,meanint,NN)
 
      do i=1,xl-1
         kkee=k(i)*k(i)/eps(i)/eps(i)
         beta(i)=kkee*NN(i)
 
         if (stab.eq.1) call cmue_cn(beta(i),cmue1(i),cmue2(i))
         if (stab.eq.2) call cmue_qe(beta(i),cmue1(i),cmue2(i))
 
       end do

      cmue1(0)=cmue1(1)
      cmue1(xl)=cmue1(xl-1)
      cmue2(0)=cmue2(1)
      cmue2(xl)=cmue2(xl-1)
 
      return
      end


!     ###############################################################
      subroutine Seiche(E_Seiche,P_Seiche,u10,v10,u,v,Az,dAdz,
     &                 NN,h,tx,ty,gamma)
!     ###############################################################

      implicit none

      include 'const_parameter.i'

!     +----------------------------------------------------------------+
!     | Global variables                                               |
!     +----------------------------------------------------------------+
      double precision E_Seiche,P_Seiche(0:xl),u10,v10
      double precision u(0:xl), v(0:xl), NN(0:xl), h(0:xl)
      double precision Az(0:xl), dAdz(0:xl)
      double precision tx, ty
      double precision gamma, minNN
!     +----------------------------------------------------------------+
!     | Local variables                                                |
!     +----------------------------------------------------------------+
      double precision W10,PS,PW,f_norm
      double precision Distrib(0:mxl)
      integer i

	minNN = 0
      if (alpha_seiche .ne. 0) then
         do i=1,xl-1
		  Distrib(i) = max(NN(i)**q_NN,minNN) / Az(i)*dAdz(i)
         end do
         if (norm_ctr.eq.1) then
              f_norm=0
              do i=1,xl-1
                  if (NN(i).gt.f_norm) f_norm=NN(i)
              end do
              f_norm=(f_norm**q_NN) * Az(xl) *rho_0     ! to create normalization
                                                        ! factor c if max NN

         else if (norm_ctr.eq.2) then
            f_norm=0
            do i=1, xl-1
                f_norm = f_norm+Distrib(i)*Az(i)*h(i)
            end do
            f_norm=f_norm*rho_0
         end if

         if ( f_norm .eq. 0.0 ) then
            do i=1, xl-1
	         Distrib(i) = 1/h(i)
!               Distrib(i) = 1/h(i)*max(NN(i)**q_NN,minNN)
!	         f_norm = f_norm + max(NN(i)**q_NN,minNN)
            end do
            f_norm=Az(xl)*rho_0
!		  f_norm=f_norm*Az(xl)*rho_0
         end if

!        Seiche energy per surface are changed to total energy of Seiche     
         W10 = sqrt(u10**2+v10**2)
         PW=alpha_Seiche*Az(xl)*rho_air*C10*W10*W10*W10
         PS=E_Seiche**(1.5)*gamma

         E_Seiche = E_Seiche + (PW - PS)*dt 

!        must limit so that E_seiche does not become negative
         if (E_seiche .lt. 0) then
            PS = (PS*dt+E_Seiche)/dt
            E_Seiche=0
         end if

!         write(6,*) E_Seiche,f_norm
         do i=1,xl-1
            P_Seiche(i) = 1/f_norm*Distrib(i)*PS
            P_Seiche(i) = P_Seiche(i)*(1-10*sqrt(CDeff))
         end do 
         Distrib(0)=0
         Distrib(xl)=0

      else			! if alpha =0
         do i=0,xl
            P_Seiche(i) = 0
         end do
      end if
      return 
      end


!     ###############################################################
      subroutine Production(u,v,NN,meanint,num,nuh,P,B)
!     ###############################################################

      implicit none

      include 'const_parameter.i'

!     +----------------------------------------------------------------+
!     | Global variables                                               |
!     +----------------------------------------------------------------+
      double precision u(0:xl),v(0:xl),meanint(0:xl),num(0:xl),
     &     nuh(0:xl),P(0:xl),B(0:xl),NN(0:xl)

!     +----------------------------------------------------------------+
!     | Local variables                                                |
!     +----------------------------------------------------------------+
      integer i

      do i=1,xl-1

!     +----------------------------------------------------------------+
!     | Shear Production of Turbulence                                 |
!     +----------------------------------------------------------------+
         P(i)=(u(i+1)-u(i))*(u(i+1)-u(i))
         P(i)=(P(i)+(v(i+1)-v(i))*(v(i+1)-v(i)))*meanint(i)*meanint(i)
         P(i)=P(i)*num(i)
!     +----------------------------------------------------------------+
!     | Buoyancy Production of Turbulence                              |
!     +----------------------------------------------------------------

         B(i)=-nuh(i)*NN(i)

      end do

      return
      end


!     ###############################################################
      subroutine TKE(num,P,B,eps,L,h,Az,u_taus,u_taub,k,ko,P_Seiche,
     &               form_k1,form_k2)
!     ###############################################################
 
      implicit none
 
      include 'const_parameter.i'
 
!     +----------------------------------------------------------------+
!     | Global variables                                               |
!     +----------------------------------------------------------------+
      double precision num(0:xl),P(0:xl),B(0:xl),k(0:xl),eps(0:xl),
     &        h(0:xl),Az(0:xl),ko(0:xl),P_Seiche(0:xl),L(0:xl)
      double precision u_taus,u_taub
      double precision form_k1(0:xl),form_k2(0:xl)
 
!     +----------------------------------------------------------------+
!     | Local variables                                                |
!     +----------------------------------------------------------------+
      double precision avh(0:mxl),au(0:mxl),bu(0:mxl),cu(0:mxl),
     &     du(0:mxl)
      double precision pminus(0:mxl),pplus(0:mxl),Prod,Buoy,Diss                
      integer i
 
      do i=0,xl
          ko(i)=k(i)                               ! ko = TKE at old time step
      end do
 
      do i=2,xl-1                                  ! average num for TKE
         avh(i)=0.5/sig_k*(num(i-1)+num(i))
      end do

      if ((fluxcond .eq. 1) .and. (Mod .eq. 1)) then
         avh(1) = 0
         avh(xl) = 0
      else
         avh(1)=u_taub**4*2/(eps(0)+eps(1))        ! = 0 for no shear stress 
         avh(xl)=u_taus**4*2/(eps(xl)+eps(xl-1))   ! = 0 for no shear stress
      end if
 
      do i=1,xl-1
         Prod = P(i)+P_Seiche(i)                   ! Add seiche energy
         Buoy=B(i)
         Diss=eps(i)
         if (Prod+Buoy.gt.0) then
            pplus(i)=Prod+Buoy
            pminus(i)=Diss
         else
            pplus(i)=Prod
            pminus(i)=Diss-Buoy
         end if
      end do
 
      do i=1,xl-1
         au(i)=dt*avh(i)*form_k1(i)
         cu(i)=dt*avh(i+1)*form_k2(i)
         bu(i)=1.-au(i)-cu(i)+pminus(i)*dt/k(i)
         du(i)=(1.+pplus(i)*dt/k(i))*k(i)
      end do
 
      if ((fluxcond .eq. 1) .and. (Mod .eq. 1)) then
         call Tridiagonal(1,xl-1,au,bu,cu,du,k)
         k(0)  = k(1)                              ! Define TKE at boundary;
         k(xl) = k(xl-1)                           ! no-flux condition
      else
         cu(0)=0
         bu(0)=1.
         du(0)=u_taub*u_taub/sqrt(cm0*cde)
 
         bu(xl)=1.
         au(xl)=0
         cu(xl)=0
         du(xl)=u_taus*u_taus/sqrt(cm0*cde)
 
         call Tridiagonal(0,xl,au,bu,cu,du,k)
      end if

      do i=0,xl
         if (k(i).lt.k_min) k(i)=k_min             ! Lower limit of TKE
      end do

      return
      end
 

!     ###############################################################
      subroutine Dissipation(cmue1,cmue2,P,B,k,ko,h,Az,eps,L,num,nuh,
     &           NN,u_taus,u_taub,P_Seiche,form_k1,form_k2,form_beps)
!     ###############################################################
 
      implicit none
 
      include 'const_parameter.i'
 
!     +----------------------------------------------------------------+
!     | Global variables                                               |
!     +----------------------------------------------------------------+ 
      double precision cmue1(0:xl),P(0:xl),B(0:xl),k(0:xl),
     &                 eps(0:xl),h(0:xl),Az(0:xl)
      double precision num(0:xl),nuh(0:xl),cmue2(0:xl),ko(0:xl)
      double precision L(0:xl),NN(0:xl),P_Seiche(0:xl)
      double precision u_taus,u_taub
      double precision form_k1(0:xl),form_k2(0:xl),form_beps(0:xl)
 
!     +----------------------------------------------------------------+
!     | Local variables                                                |
!     +----------------------------------------------------------------+
      double precision avh(0:mxl),au(0:mxl),bu(0:mxl)
      double precision cu(0:mxl),du(0:mxl)
      double precision flux(0:mxl)
      double precision pminus(0:mxl),pplus(0:mxl),Prod,Buoy,Diss,cee3
      integer i
      double precision epslim 

 
      do i=1,xl                                          ! Average num for Diss
         avh(i)=0.5/sig_e*(num(i-1)+num(i))
      end do

      if ((fluxcond .eq. 1) .and. (Mod .eq. 1)) then
         flux(0)    = avh(1 )*(cde*((ko(1   ))**1.5)
     &                /(kappa*(K_s+0.5*h(1 )))**2.)
         flux(xl) = avh(xl)*(cde*((ko(xl-1))**1.5)
     &                /(kappa*(z0 +0.5*h(xl)))**2.)
         do i=1,xl-1
           flux(i) = num(i)/sig_e*(cde*((ko(i))**1.5)       
     &                /(kappa*(z0 +0.25*(h(i)+h(i+1))))**2.)
         end do

!         do i=1,xl-1
!            flux(i) = -flux(i)*eps(i)/ko(i)
!         end do

         avh(1 ) = 0
         avh(xl) = 0

      else
         avh(1)=u_taub**4*2/sig_e/(eps(0)+eps(1))        ! = 0 for no shear stress
         avh(xl)=u_taus**4*2/sig_e/(eps(xl)+eps(xl-1))   ! = 0 for no shear stress
      end if

      do i=1,xl-1
         if (B(i).gt.0) then                             
            cee3=1.
         else
            cee3=ce3
         end if

         Prod=ce1*eps(i)/ko(i)*(P(i)+P_Seiche(i))        ! New code plus seiche
         Buoy=cee3*eps(i)/ko(i)*B(i)                     
         Diss=ce2*eps(i)*eps(i)/ko(i)                    
         if (Prod+Buoy.gt.0) then
            pplus(i)=Prod+Buoy
            pminus(i)=Diss
         else
            pplus(i)=Prod
            pminus(i)=Diss-Buoy
         end if
      end do
 
      do i=1,xl-1
         au(i)=dt*avh(i)*form_k1(i)
         cu(i)=dt*avh(i+1)*form_k2(i)
         bu(i)=1.-au(i)-cu(i)+pminus(i)*dt/eps(i)
         du(i)=(1.+pplus(i)*dt/eps(i))*eps(i)
      end do
 
      if ((fluxcond .eq. 1) .and. (Mod .eq. 1)) then
         do i=1,xl-1
            du(i)=du(i)+flux(i)*dt*form_beps(i)      ! form_beps = 1/A * dA/dz             
         end do                                      ! (at epsilon posions)
         if (Az(0) .ne. 0) then                      ! Flux from bottom only!
            du(1)=du(1)+flux(0)*dt/
     &           (Az(0)+Az(1))/(Az(1)*(h(1)+h(2)))
         end if
         du(xl-1)=du(xl-1)+flux(xl)*dt/
     &        (Az(xl)+Az(xl-1))/(Az(xl-1)*(h(xl)+h(xl-1)))
         call Tridiagonal(1,xl-1,au,bu,cu,du,eps)
         eps(0) = eps(1)    + (cde*((ko(1   ))**1.5) ! Define eps at boundaries
     &                /(kappa*(K_s+1.0*h(1 )))**2.)*h(1)
         eps(xl)= eps(xl-1) + (cde*((ko(xl-1))**1.5)
     &                /(kappa*(z0 +1.0*h(xl)))**2.)*h(xl)
      else
         cu(0)=0
         bu(0)=1.
         du(0)=cde*sqrt(k(0)*k(0)*k(0))/kappa/K_s
 
         bu(xl)=1.
         au(xl)=0
         du(xl)=cde*sqrt(k(xl)*k(xl)*k(xl))/kappa/z0      

         call Tridiagonal(0,xl,au,bu,cu,du,eps)
      end if

      do i=0,xl
         if (NN(i).gt.0) then
            epslim=0.212*k(i)*sqrt(NN(i)) 
         else
            epslim=epsmin 
         end if 
         if (eps(i).lt.epslim) eps(i)=epslim
         if (eps(i) .lt. 0) then
            write(*,*) 'Dissipation negative'
         end if

         num(i)=cmue1(i)*k(i)*k(i)/eps(i)+1.5e-6
         nuh(i)=cmue2(i)*k(i)*k(i)/eps(i)+1.5e-7
         L(i)=cde*sqrt(ko(i)*ko(i)*ko(i))/eps(i)

      end do

      num(0 ) = kappa*u_taub*K_s+avhmin
      num(xl) = kappa*u_taus*z0+avhmin
      return
      end


!     ###############################################################
      subroutine Advection(Qvert,Q_inp,u,v,T,S,k,eps,num,nuh,zu,zk,h,Az)
!     ###############################################################

	  implicit none


      include 'const_parameter.i'

      double precision Qvert(0:mxl)	! Vertical advection [m/s]
      double precision Q_inp(1:4,0:mxl)	! Q_inp(1):inflow [m3/s]
										! Q_inp(2):outflow [m3/s]
										! Q_inp(3):T-input [Tm3/s]
										! Q_inp(4):S-input [‰m3/s]
      double precision u(0:mxl), v(0:mxl), T(0:mxl), S(0:mxl)
      double precision k(0:mxl), eps(0:mxl), num(0:mxl), nuh(0:mxl)
      double precision  h(0:mxl), Az(0:mxl), zu(0:mxl), zk(0:mxl)

      integer i
      double precision du(0:mxl), dv(0:mxl), dTemp(0:mxl), dS(0:mxl)
      double precision dk(0:mxl), deps(0:mxl), dnum(0:mxl), dnuh(0:mxl)
      double precision form_adv1(1:mxl), form_adv2(1:mxl)
      double precision dh, dh1, dh2		! depth difference
      double precision dt1, dt2			! first and second time step
	  double precision top

      dh = Qvert(xl)/Az(xl)*dt				! New depth difference
	  if (dh == 0) then						! If volume does not change
		 dt1 = dt							! One normal time step
         dt2 = 0							! second step not necessary

	  else if ((dh+zk(xl)) >= depth) then	! If surface level reached
	     dt1 = (depth - zk(xl))/dh*dt
         dt2 = dt-dt1
      
	  else if (((dh+h(xl)) .gt. h(xl-1)/2) .and.	! If top box > 0.5*lower box
     &			((dh+h(xl)) .lt. 2*h(xl-1))) then	! and top box < 2*lower box
         dt1 = dt							! One normal time step
         dt2 = 0							! second step not necessary	  


      else if ((dh+h(xl)) .le. 
     &		h(xl-1)/2) then					! If top box <= 0.5*lower box
		 dt1 = abs((h(xl)-h(xl-1)/2)/dh*dt)	! Step till top box = lower box / 2
         dt2 = dt-dt1						! Rest
   

      else									! If top box >= 2*lower box
         dt1 = (2*h(xl-1)-h(xl))/dh*dt		! Step till top box = 2*lower box
         dt2 = dt-dt1						! Rest
      end if 

      do i=1,xl								! form factor
         form_adv1(i) = dt1/(Az(i)*h(i))
      end do
      dh1 = dh*dt1/dt						! depth difference for dt1

      do i=1,xl								! Advection out of box i, always negative	
	     if ((i == xl) .and. (Qvert(i) > 0)) then
		    top = 0
		 else
			top = 1
	     end if
         du(i)    = -top*abs(form_adv1(i)*Qvert(i))*u(i)
         dv(i)    = -top*abs(form_adv1(i)*Qvert(i))*v(i)
         dTemp(i) = -top*abs(form_adv1(i)*Qvert(i))*T(i)
         dS(i)    = -top*abs(form_adv1(i)*Qvert(i))*S(i)
         if (i > 1) then
		    if (Qvert(i-1) > 0) then
			   du(i)    = du(i) + form_adv1(i)*Qvert(i-1)*u(i-1)
			   dv(i)    = dv(i) + form_adv1(i)*Qvert(i-1)*v(i-1)
			   dTemp(i) = dTemp(i) + form_adv1(i)*Qvert(i-1)*T(i-1)
			   dS(i)    = dS(i) + form_adv1(i)*Qvert(i-1)*S(i-1)
			end if
         end if
         if (i < xl) then
		    if (Qvert(i+1) < 0) then
               du(i)    = du(i) - form_adv1(i)*Qvert(i+1)*u(i+1)
			   dv(i)    = dv(i) - form_adv1(i)*Qvert(i+1)*v(i+1)
			   dTemp(i) = dTemp(i) - form_adv1(i)*Qvert(i+1)*T(i+1)
			   dS(i)    = dS(i) - form_adv1(i)*Qvert(i+1)*S(i+1)
			end if
         end if
      end do

      do i=1,xl				! Inflow and outflow
         dTemp(i) = dTemp(i) + form_adv1(i)*(Q_inp(3,i)+Q_inp(2,i)*T(i))
         dS(i)    = dS(i) + form_adv1(i)*(Q_inp(4,i)+Q_inp(2,i)*S(i))
      end do

      do i=1,xl				! Take first time step
         u(i) = u(i) + du(i)
         v(i) = v(i) + dv(i)
         T(i) = T(i) + dTemp(i)
         S(i) = S(i) + dS(i)
      end do
		
      u(xl) = u(xl) - u(xl) *dh1/(h(xl)+dh1)	! Consider variation of variables
      v(xl) = v(xl) - v(xl) *dh1/(h(xl)+dh1)	! due to change in volume
      T(xl) = T(xl) - T(xl) *dh1/(h(xl)+dh1)
      S(xl) = S(xl) - S(xl) *dh1/(h(xl)+dh1)

	  if (dh == 0) then
         dh = 0							! dummy

	  else if ((dh+zk(xl)) >= depth) then	! If surface level reached
         h(xl) = h(xl) + dh1			! new magnitude of top box
         zu(xl) = zu(xl) + dh1/2		! new centre coordinate of top box
		 zk(xl) = zk(xl) + dh1			! new surface coordinate of top box
		 dh2 = 0						! No change in volume

	  else if (((dh+h(xl)) .gt. h(xl-1)/2) .and. 
     &		((dh+h(xl)) .lt. 2*h(xl-1))) then
         h(xl) = h(xl) + dh1			! new magnitude of top box
         zu(xl) = zu(xl) + dh1/2		! new centre coordinate of top box
		 zk(xl) = zk(xl) + dh1			! new surface coordinate of top box
		 return

	  else if ((dh+h(xl)) .le.			! If top box <= lower box / 2
     &			h(xl-1)/2) then	  
	     zk(xl-1) = zk(xl)
		 zu(xl-1) = (zk(xl-1)+zk(xl-2))/2
		 h(xl-1)  = (zk(xl-1)-zk(xl-2))	
		 U(xl-1) = (0.5*U(xl)*Az(xl)+U(xl-1)*Az(xl-1))/(0.5*Az(xl)+Az(xl-1))
		 V(xl-1) = (0.5*V(xl)*Az(xl)+V(xl-1)*Az(xl-1))/(0.5*Az(xl)+Az(xl-1))
		 T(xl-1) = (0.5*T(xl)*Az(xl)+T(xl-1)*Az(xl-1))/(0.5*Az(xl)+Az(xl-1))
		 S(xl-1) = (0.5*S(xl)*Az(xl)+S(xl-1)*Az(xl-1))/(0.5*Az(xl)+Az(xl-1))
		 k(xl-1) = (0.5*k(xl)*Az(xl)+k(xl-1)*Az(xl-1))/(0.5*Az(xl)+Az(xl-1))
		 eps(xl-1) = (0.5*eps(xl)*Az(xl)+eps(xl-1)*Az(xl-1))/(0.5*Az(xl)+Az(xl-1))
		 Qvert(xl-1) = (0.5*Qvert(xl)*Az(xl)+Qvert(xl-1)*Az(xl-1))/(0.5*Az(xl)+Az(xl-1))
         xl = xl-1						! Reduce number of boxes
		 dh2 = (Qvert(xl))/Az(xl)*dt2
		 call Form
      else	
	     h(xl+1)   = h(xl)/2
		 h(xl)     = h(xl)/2
		 zk(xl+1)  = zk(xl)
		 zk(xl)    = zk(xl) - h(xl)/2
		 zu(xl+1)  = zu(xl) + h(xl)/4
		 zu(xl)    = zu(xl) - h(xl)/4 
		 u(xl+1)   = u(xl)
		 v(xl+1)   = v(xl)
		 T(xl+1)   = T(xl)
		 S(xl+1)   = S(xl)
		 k(xl+1)   = k(xl)
		 eps(xl+1) = eps(xl)
	     Qvert(xl+1) = Qvert(xl)			! Vertical discharge of new box
		 xl        = xl + 1				! Increase number of boxes
		 !dh2 = dh*dt2/dt
		 dh2 = (Qvert(xl))/Az(xl)*dt2
	     call Form
      end if 

      do i=1,xl							! form factor
         form_adv2(i) = dt2/(Az(i)*h(i))
      end do     
		
	  if (dt2 > 0.0) then							
      do i=1,xl							! Second time step (if necessary)
	     if (((i == xl) .and. (Qvert(i) > 0)) .and. (zk(xl) < depth)) then
		    top = 0
		 else
			top = 1
	     end if
         du(i)    = -top*abs(form_adv2(i)*Qvert(i))*u(i)
         dv(i)    = -top*abs(form_adv2(i)*Qvert(i))*v(i)
         dTemp(i) = -top*abs(form_adv2(i)*Qvert(i))*T(i)
         dS(i)    = -top*abs(form_adv2(i)*Qvert(i))*S(i)
         if (i .gt. 1) then
		    if (Qvert(i-1) > 0) then
               du(i)    = du(i) + form_adv2(i)*Qvert(i-1)*u(i-1)
			   dv(i)    = dv(i) + form_adv2(i)*Qvert(i-1)*v(i-1)
			   dTemp(i) = dTemp(i) + form_adv2(i)*Qvert(i-1)*T(i-1)
			   dS(i)    = dS(i) + form_adv2(i)*Qvert(i-1)*S(i-1)
			end if
         end if
         if (i .lt. xl) then
		    if (Qvert(i+1) < 0) then
               du(i)    = du(i) - form_adv2(i)*Qvert(i+1)*u(i+1)
			   dv(i)    = dv(i) - form_adv2(i)*Qvert(i+1)*v(i+1)
			   dTemp(i) = dTemp(i) - form_adv2(i)*Qvert(i+1)*T(i+1)
			   dS(i)    = dS(i) - form_adv2(i)*Qvert(i+1)*S(i+1)
		    end if
         end if
      end do

      do i=1,xl				! Inflow and outflow
         dTemp(i) = dTemp(i) + form_adv2(i)*(Q_inp(3,i)+Q_inp(2,i)*T(i))
         dS(i)    = dS(i) + form_adv2(i)*(Q_inp(4,i)+Q_inp(2,i)*S(i))
      end do

      do i=1,xl				! Take second time step
         u(i) = u(i) + du(i)
         v(i) = v(i) + dv(i)
         T(i) = T(i) + dTemp(i)
         S(i) = S(i) + dS(i)
      end do

      u(xl) = u(xl) - u(xl) *dh2/(h(xl)+dh2)	! Consider variation of variables
      v(xl) = v(xl) - v(xl) *dh2/(h(xl)+dh2)	! due to change in volume
      T(xl) = T(xl) - T(xl) *dh2/(h(xl)+dh2)
      S(xl) = S(xl) - S(xl) *dh2/(h(xl)+dh2)

      end if
      return 
      end


!     ###############################################################
      subroutine Form
!     ###############################################################

      implicit none
 
      include 'const_parameter.i'

!     +---------------------------------------------------------------+
!     |  Global Declarations                                        |
!     +---------------------------------------------------------------+
      double precision h(0:mxl)
      double precision Az(0:mxl),dAdz(0:mxl)
      double precision form_1(0:mxl),form_2(0:mxl)
      double precision form_k1(0:mxl),form_k2(0:mxl),form_beps(0:mxl)
      double precision volume, meanint(0:mxl)

      common /morph_variables1/ h,Az,dAdz
      common /morph_variables3/ volume
      common /form_variables  / form_1,form_2
      common /form_variablesk/form_k1,form_k2,form_beps
      common /form_variables3/meanint
!     +---------------------------------------------------------------+
!     |  Local Declarations                                        |
!     +---------------------------------------------------------------+
	  integer i

	  do i=1,xl-1
        form_1(i) =-4.*Az(i-1)/(h(i)+h(i-1))/h(i)/(Az(i)+Az(i-1)) 
        form_2(i) =-4.*Az(i)/(h(i)+h(i+1))/h(i)/(Az(i)+Az(i-1))
        form_k1(i)=-(Az(i)+Az(i+1))/(h(i)+h(i+1))/h(i+1)/Az(i)
        form_k2(i)=-(Az(i)+Az(i-1))/(h(i)+h(i+1))/h(i)/Az(i)
        form_beps(i)=( (Az(i)-Az(i-1))/h(i)+(Az(i+1)-Az(i))/h(i+1) )
     &                            /2/Az(i)
      end do 
      form_1(xl) =-4.*Az(xl-1)/(h(xl)+h(xl-1))/h(xl)/(Az(xl)+Az(xl-1)) 
      form_2(xl) =-4.*Az(xl)/(h(xl)+h(xl+1))/h(xl)/(Az(xl)+Az(xl-1))

      volume=0
      do i=0,xl-1
         meanint(i)=2.0/(h(i)+h(i+1))
         volume=volume+h(i+1)*(Az(i)+Az(i+1))/2
      end do  

	  return
	  end

 
!     ###############################################################
      subroutine  avstate(iav,u,v,T,S,k,eps,num,nuh,
     &                        B,P,NN,P_Seiche,E_Seiche,
     &                        uav,vav,Tav,Sav,kav,epsav,numav,nuhav,
     &                        Bav,Pav,NNav,P_Seicheav,E_Seicheav)

!     ###############################################################

      implicit none

      include 'const_parameter.i'

!     +---------------------------------------------------------------+
!     |  Global Declarations                                        |
!     +---------------------------------------------------------------+

      double precision u(0:xl),v(0:xl)
      double precision T(0:xl),S(0:xl)
      double precision k(0:xl)
      double precision eps(0:xl)
      double precision num(0:xl),nuh(0:xl)
      double precision B(0:xl),P(0:xl),NN(0:xl),P_Seiche(0:xl),E_Seiche

      double precision uav(0:xl),vav(0:xl)
      double precision Tav(0:xl),Sav(0:xl)
      double precision kav(0:xl)
      double precision epsav(0:xl)
      double precision numav(0:xl),nuhav(0:xl)
      double precision Bav(0:xl),Pav(0:xl),NNav(0:xl)
      double precision P_Seicheav(0:xl),E_Seicheav
	  integer iav
!      double precision h(0:xl)
!      double precision SST

!     +---------------------------------------------------------------+
!     |  Local Declarations                                        |
!     +---------------------------------------------------------------+
 
	  integer i



        if (igoal_s(1).eq.1) then
               if (iav.eq.1) then
		    do i=0,xl
		        uav(i)=u(i)
                    end do
               else
		    do i=0,xl
		        uav(i)=uav(i)+u(i)
                    end do
               end if
        end if
        if (igoal_s(2).eq.1) then
               if (iav.eq.1) then
		    do i=0,xl
		        vav(i)=v(i)
                    end do
               else
		    do i=0,xl
		        vav(i)=vav(i)+v(i)
                    end do
               end if
        end if
        if (igoal_s(3).eq.1) then
               if (iav.eq.1) then
		    do i=0,xl
		        Tav(i)=T(i)
                    end do
               else
		    do i=0,xl
		        Tav(i)=Tav(i)+T(i)
                    end do
               end if
        end if        
        if (igoal_s(4).eq.1) then
               if (iav.eq.1) then
		    do i=0,xl
		        Sav(i)=S(i)
                    end do
               else
		    do i=0,xl
		        Sav(i)=Sav(i)+S(i)
                    end do
               end if
        end if        
        if (igoal_s(5).eq.1) then
               if (iav.eq.1) then
		    do i=0,xl
		        kav(i)=k(i)
                    end do
               else
		    do i=0,xl
		        kav(i)=kav(i)+k(i)
                    end do
               end if
        end if
        if (igoal_s(6).eq.1) then
               if (iav.eq.1) then
		    do i=0,xl
		        epsav(i)=eps(i)
                    end do
               else
		    do i=0,xl
		        epsav(i)=epsav(i)+eps(i)
                    end do
               end if
        end if
        if (igoal_s(7).eq.1) then
               if (iav.eq.1) then
		    do i=0,xl
		        numav(i)=num(i)
                    end do
               else
		    do i=0,xl
		        numav(i)=numav(i)+num(i)
                    end do
               end if
        end if
        if (igoal_s(8).eq.1) then
               if (iav.eq.1) then
		    do i=0,xl
		        nuhav(i)=nuh(i)
                    end do
               else
		    do i=0,xl
		        nuhav(i)=nuhav(i)+nuh(i)
                    end do
               end if
        end if
        if (igoal_s(9).eq.1) then
               if (iav.eq.1) then
		    do i=0,xl
		        Bav(i)=B(i)
                    end do
               else
		    do i=0,xl
		        Bav(i)=Bav(i)+B(i)
                    end do
               end if
        end if
        if (igoal_s(10).eq.1) then
               if (iav.eq.1) then
		    do i=0,xl
		        Pav(i)=P(i)
                    end do
               else
		    do i=0,xl
		        Pav(i)=Pav(i)+k(i)
                    end do
               end if
        end if
        if (igoal_s(11).eq.1) then
               if (iav.eq.1) then
		    do i=0,xl
		        P_Seicheav(i)=P_Seiche(i)
                    end do
               else
		    do i=0,xl
		        P_Seicheav(i)=P_Seicheav(i)+P_Seiche(i)
                    end do
               end if
        end if
        if (igoal_s(12).eq.1) then
               if (iav.eq.1) then
		    do i=0,xl
		        NNav(i)=NN(i)
                    end do
               else
                    do i=0,xl
		        NNav(i)=NNav(i)+NN(i)
                    end do
               end if
        end if
        if (igoal_s(13).eq.1) then
               if (iav.eq.1) then
		    E_Seicheav=E_Seiche
               else
		    E_Seicheav=E_Seicheav+E_Seiche
               end if
        end if



	  
          return
	  end
   
!     ###############################################################
      subroutine write_out (datum,u,v,T,S,k,eps,num,nuh,
     &                       B,P,NN,P_Seiche,E_Seiche,zu,zk,M)
!     ###############################################################
      implicit none
 
      include 'const_parameter.i'

!     +---------------------------------------------------------------+
!     |  Global Declarations                                        |
!     +---------------------------------------------------------------+

      double precision u(0:xl),v(0:xl)
      double precision T(0:xl),S(0:xl)
      double precision k(0:xl)
      double precision eps(0:xl)
      double precision num(0:xl),nuh(0:xl)
      double precision datum
      double precision B(0:xl),P(0:xl),NN(0:xl),P_Seiche(0:xl),E_Seiche
	  double precision M(0:xl)
      double precision zu(0:xl),zk(0:xl)
 
!      double precision h(0:xl)
!      double precision SST

!     +---------------------------------------------------------------+
!     |  Local Declarations                                        |
!     +---------------------------------------------------------------+
      double precision xs(0:mxl)	 
	  integer i, istart

  

      write(80) datum

      if (save_vctr.eq.0) then


           if (igoal_s(1).eq.1) then
             write(80)(u(indexu_save(i)), i=1,xlu)
           end if
           if (igoal_s(2).eq.1) then
             write(80)(v(indexu_save(i)), i=1,xlu)
           end if
           if (igoal_s(3).eq.1) then
             write(80)(T(indexu_save(i)), i=1,xlu)
           end if          
           if (igoal_s(4).eq.1) then
             write(80)(S(indexu_save(i)), i=1,xlu)
           end if
           if (igoal_s(5).eq.1) then
             write(80)(k(indexu_save(i)), i=1,xlk)
           end if             
           if (igoal_s(6).eq.1) then
             write(80)(eps(indexu_save(i)), i=1,xlk)
           end if             
           if (igoal_s(7).eq.1) then
             write(80)(num(indexu_save(i)), i=1,xlk)
           end if             
           if (igoal_s(8).eq.1) then
             write(80)(nuh(indexu_save(i)), i=1,xlk)
           end if             
           if (igoal_s(9).eq.1) then
             write(80)(B(indexu_save(i)), i=1,xlk)
           end if             
           if (igoal_s(10).eq.1) then
             write(80)(P(indexu_save(i)), i=1,xlk)
           end if             
           if (igoal_s(11).eq.1) then
             write(80)(P_Seiche(indexu_save(i)), i=1,xlk)
           end if             
           if (igoal_s(12).eq.1) then
             write(80)(NN(indexu_save(i)), i=1,xlk)
           end if 
           if (igoal_s(13).eq.1) then
             write(80) E_Seiche
           end if 
	  
	else

             if (ifit.eq.0) then
	        if (igoal_s(1).eq.1) then
		         call Interp(zu,u,xl,zs,xs,save_vctr)
		         write(80) (xs(i), i=save_vctr, 0, -1)      
                end if                   
	        if (igoal_s(2).eq.1) then
		         call Interp(zu,v,xl,zs,xs,save_vctr)
		         write(80) (xs(i), i=save_vctr, 0, -1)      
                end if                   
	        if (igoal_s(3).eq.1) then
		         call Interp(zu,T,xl,zs,xs,save_vctr)
		         write(80) (xs(i), i=save_vctr, 0, -1)      
                end if                   
	        if (igoal_s(4).eq.1) then
		         call Interp(zu,S,xl,zs,xs,save_vctr)
		         write(80) (xs(i), i=save_vctr, 0, -1)      
                end if                   
	        if (igoal_s(5).eq.1) then
		         call Interp(zu,k,xl,zs,xs,save_vctr)
		         write(80) (xs(i), i=save_vctr, 0, -1)      
                end if                   
	        if (igoal_s(6).eq.1) then
		         call Interp(zu,eps,xl,zs,xs,save_vctr)
		         write(80) (xs(i), i=save_vctr, 0, -1)      
                end if                   
	        if (igoal_s(7).eq.1) then
		         call Interp(zu,num,xl,zs,xs,save_vctr)
		         write(80) (xs(i), i=save_vctr, 0, -1)      
                end if                   
	        if (igoal_s(8).eq.1) then
		         call Interp(zu,nuh,xl,zs,xs,save_vctr)
		         write(80) (xs(i), i=save_vctr, 0, -1)      
                end if                   
	        if (igoal_s(9).eq.1) then
		         call Interp(zu,B,xl,zs,xs,save_vctr)
		         write(80) (xs(i), i=save_vctr, 0, -1)      
                end if                   
	        if (igoal_s(10).eq.1) then
		         call Interp(zu,P,xl,zs,xs,save_vctr)
		         write(80) (xs(i), i=save_vctr, 0, -1)      
                end if                   
	        if (igoal_s(11).eq.1) then
		         call Interp(zu,P_Seiche,xl,zs,xs,save_vctr)
		         write(80) (xs(i), i=save_vctr, 0, -1)      
                end if                   
	        if (igoal_s(12).eq.1) then
		         call Interp(zu,NN,xl,zs,xs,save_vctr)
		         write(80) (xs(i), i=save_vctr, 0, -1)      
                end if                   
	        if (igoal_s(13).eq.1) then
		         write(80) E_Seiche      
                end if                   

             else
                istart=0
		if (igoal_s(1).eq.1) then
		         call Interp(zu,u,xl,zs,xs,save_vctr)
                         do i=istart+save_vctr, istart, -1
		             M(i)=xs(i)
	                 end do
                         istart=istart+save_vctr+1               
                end if                   
	        if (igoal_s(2).eq.1) then
		         call Interp(zu,v,xl,zs,xs,save_vctr)
                         do i=istart+save_vctr, istart, -1
		             M(i)=xs(i)
	                 end do
                         istart=istart+save_vctr+1               
                end if                   
	        if (igoal_s(3).eq.1) then
		         call Interp(zu,T,xl,zs,xs,save_vctr)
                         do i=istart+save_vctr, istart, -1
		             M(i)=xs(i)
	                 end do
                         istart=istart+save_vctr+1               
                end if                   
	        if (igoal_s(4).eq.1) then
		         call Interp(zu,S,xl,zs,xs,save_vctr)
                         do i=istart+save_vctr, istart, -1
		             M(i)=xs(i)
	                 end do
                         istart=istart+save_vctr+1               
                end if                   
	        if (igoal_s(5).eq.1) then
		         call Interp(zu,k,xl,zs,xs,save_vctr)
                         do i=istart+save_vctr, istart, -1
		             M(i)=xs(i)
	                 end do
                         istart=istart+save_vctr+1               
                end if                   
	        if (igoal_s(6).eq.1) then
		         call Interp(zu,eps,xl,zs,xs,save_vctr)
                         do i=istart+save_vctr, istart, -1
		             M(i)=xs(i)
	                 end do
                         istart=istart+save_vctr+1               
                end if                   
	        if (igoal_s(7).eq.1) then
		         call Interp(zu,num,xl,zs,xs,save_vctr)
                         do i=istart+save_vctr, istart, -1
		             M(i)=xs(i)
	                 end do
                         istart=istart+save_vctr+1               
                end if                   
	        if (igoal_s(8).eq.1) then
		         call Interp(zu,nuh,xl,zs,xs,save_vctr)
                         do i=istart+save_vctr, istart, -1
		             M(i)=xs(i)
	                 end do
                         istart=istart+save_vctr+1               
                end if                   
	        if (igoal_s(9).eq.1) then
		         call Interp(zu,B,xl,zs,xs,save_vctr)
                         do i=istart+save_vctr, istart, -1
		             M(i)=xs(i)
	                 end do
                         istart=istart+save_vctr+1               
                end if                   
	        if (igoal_s(10).eq.1) then
		         call Interp(zu,P,xl,zs,xs,save_vctr)
                         do i=istart+save_vctr, istart, -1
		             M(i)=xs(i)
	                 end do
                         istart=istart+save_vctr+1               
                end if                   
	        if (igoal_s(11).eq.1) then
		         call Interp(zu,P_Seiche,xl,zs,xs,save_vctr)
                         do i=istart+save_vctr, istart, -1
		             M(i)=xs(i)
	                 end do
                         istart=istart+save_vctr+1               
                end if                   
	        if (igoal_s(12).eq.1) then
		         call Interp(zu,NN,xl,zs,xs,save_vctr)
                         do i=istart+save_vctr, istart, -1
		             M(i)=xs(i)
	                 end do
                         istart=istart+save_vctr+1               
                end if                  
	        if (igoal_s(13).eq.1) then
		         M(istart)= E_Seiche      
                end if                   

	     endif   ! ifit
         endif    ! save -ctr


          return
	  end
  
!     ###############################################################
      subroutine prep_outav(iav,uav,vav,Tav,Sav,kav,epsav,numav,nuhav,
     &                     Bav,Pav,NNav,P_Seicheav,E_Seicheav)
!     ###############################################################
      implicit none

      include 'const_parameter.i'

!     +---------------------------------------------------------------+
!     |  Global Declarations                                        |
!     +---------------------------------------------------------------+


      double precision uav(0:xl),vav(0:xl)
      double precision Tav(0:xl),Sav(0:xl)
      double precision kav(0:xl)
      double precision epsav(0:xl)
      double precision numav(0:xl),nuhav(0:xl)
      double precision Bav(0:xl),Pav(0:xl),NNav(0:xl)
      double precision P_Seicheav(0:xl),E_Seicheav
      integer iav
!      double precision h(0:xl)
!      double precision SST

!     +---------------------------------------------------------------+
!     |  Local Declarations                                        |
!     +---------------------------------------------------------------+
	 
	  integer i

  
!            print 999,datum,T(xl),SST,num(xl),E_Seiche,heat,I_0,ga1(1)
          
                 !
        if (igoal_s(1).eq.1) then
		    do i=0,xl
		        uav(i)=uav(i)/iav
                    end do
        end if
        if (igoal_s(2).eq.1) then
		    do i=0,xl
		        vav(i)=vav(i)/iav
                    end do
        end if
        if (igoal_s(3).eq.1) then
		    do i=0,xl
		        Tav(i)=Tav(i)/iav
                    end do
        end if
        if (igoal_s(4).eq.1) then
		    do i=0,xl
		        Sav(i)=Sav(i)/iav
                    end do
        end if
        if (igoal_s(5).eq.1) then
		    do i=0,xl
		        kav(i)=kav(i)/iav
                    end do
        end if
        if (igoal_s(6).eq.1) then
		    do i=0,xl
		        epsav(i)=epsav(i)/iav
                    end do
        end if
        if (igoal_s(7).eq.1) then
		    do i=0,xl
		        numav(i)=numav(i)/iav
                    end do
        end if
        if (igoal_s(8).eq.1) then
		    do i=0,xl
		        nuhav(i)=nuhav(i)/iav
                    end do
        end if
        if (igoal_s(9).eq.1) then
		    do i=0,xl
		        Bav(i)=Bav(i)/iav
                    end do
        end if
        if (igoal_s(10).eq.1) then
		    do i=0,xl
		        Pav(i)=Pav(i)/iav
                    end do
        end if
        if (igoal_s(11).eq.1) then
		    do i=0,xl
		        P_Seicheav(i)=P_Seicheav(i)/iav
                    end do
        end if
        if (igoal_s(12).eq.1) then
		    do i=0,xl
		        NNav(i)=NNav(i)/iav
                    end do
        end if
        if (igoal_s(13).eq.1) then
		        E_Seicheav=E_Seicheav/iav
        end if
        

          return
	  end


!     ###############################################################
      subroutine write_out_new (datum,std,u,v,T,S,k,eps,nuh,
     &                       B,P,NN,P_Seiche,E_Seiche,zu,zk,M,Qvert)
!     ###############################################################
      implicit none
 
      include 'const_parameter.i'

!     +---------------------------------------------------------------+
!     |  Global Declarations                                        |
!     +---------------------------------------------------------------+

	  double precision datum
	  integer std
      double precision u(0:xl),v(0:xl),T(0:xl),S(0:xl),k(0:xl),eps(0:xl),nuh(0:xl)
      double precision B(0:xl),P(0:xl),NN(0:xl),P_Seiche(0:xl),E_Seiche
	  double precision zu(0:xl),zk(0:xl),M(0:xl),Qvert(0:xl)
	  integer write_tout
      double precision tout_ctr1(0:9000),tout_ctr2(0:9000)
	  double precision E_Seicheini,datumini
      double precision dragini,txini,tyini


      common /ini_values5/ dragini,txini,tyini,E_Seicheini,datumini
	  common /savet1  /   write_tout
	  common /savet2  /   tout_ctr1
      common /savet3  /   tout_ctr2

!     +---------------------------------------------------------------+
!     |  Local Declarations                                        |
!     +---------------------------------------------------------------+
      double precision xs(0:mxl)
	  double precision zext(0:mxl)
	  real xt(0:mxl)
	  integer i, istart
	  
      zext(0) = zk(0)
	  do i=2,xl-1
	     zext(i-1) = zu(i)
	  end do
	  zext(xl-1) = zk(xl)

	  write(80) datum
	  call Interp_results(zext,u,xl-1,zs, xs,save_vctr)
	  write(80) (xs(i), i=save_vctr, 0, -1)
	  call Interp_results(zext,v,xl-1,zs,xs,save_vctr)
	  write(80) (xs(i), i=save_vctr, 0, -1)      
	  call Interp_results(zext,T,xl-1,zs,xs,save_vctr)
	  write(80) (xs(i), i=save_vctr, 0, -1)      
	  call Interp_results(zext,S,xl-1,zs,xs,save_vctr)
	  write(80) (xs(i), i=save_vctr, 0, -1)
	  
	  call Interp_results(zext,Qvert,xl-1,zs,xs,save_vctr)
	  write(80) (xs(i), i=save_vctr, 0, -1)
	        
      call Interp_results(zk,k,xl,zs,xs,save_vctr)
	  write(80) (xs(i), i=save_vctr, 0, -1)      
	  call Interp_results(zk,eps,xl,zs,xs,save_vctr)
	  write(80) (xs(i), i=save_vctr, 0, -1)      
	  call Interp_results(zk,nuh,xl,zs,xs,save_vctr)
	  write(80) (xs(i), i=save_vctr, 0, -1)      
	  call Interp_results(zk,B,xl,zs,xs,save_vctr)
	  write(80) (xs(i), i=save_vctr, 0, -1)      
	  call Interp_results(zk,P,xl,zs,xs,save_vctr)
	  write(80) (xs(i), i=save_vctr, 0, -1)      
	  call Interp_results(zk,P_Seiche,xl,zs,xs,save_vctr)
	  write(80) (xs(i), i=save_vctr, 0, -1)      
	  call Interp_results(zk,NN,xl,zs,xs,save_vctr)
	  write(80) (xs(i), i=save_vctr, 0, -1)      

	  write(80) E_Seiche, zk(xl), zu(xl), u(xl), v(xl), T(xl),S(xl), 
     &			k(xl),eps(xl),nuh(xl), B(xl), P(xl),P_Seiche(xl),NN(xl), Qvert(xl)

	  return	 
	  end
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	  include 'keps_utilities_F10.f90'
      include 'keps_initialization_F10.f90'	 
	 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! 	Dokumentation of changes to original program
!	
!	24.7.98: 
!	include 'utilities.f'  to  include 'keps_utilities.f'
!   This affects input parameter files:
!	Now two files required: 1: keps_parameter.dat
!							2: file_name.dat
!	keps_parameter.dat: contains parameters for 
!								- k-eps and MY equations
!								- freshwater
!   file_name.data: contains name of file with lake and model specific parameters
!					see Parameterdefault.dat
