!     ###############################################################
      subroutine Initialization
!     ###############################################################

      implicit none

      include 'const_parameter.i'



      double precision Sini(0:mxl),Tini(0:mxl),uini(0:mxl),vini(0:mxl) 
      double precision kini(0:mxl),epsini(0:mxl),Lini(0:mxl)
      double precision numini(0:mxl),nuhini(0:mxl) 
      double precision E_Seicheini,datumini
      double precision dragini,txini,tyini
      double precision h(0:mxl)
      double precision Az(0:mxl),dAdz(0:mxl)
 

      common /ini_values1/ Sini,Tini
      common /ini_values2/ uini,vini
      common /ini_values3/ epsini,kini,Lini
      common /ini_values4/ numini,nuhini
      common /ini_values5/ dragini,txini,tyini,E_Seicheini,datumini


      common /morph_variables1/ h,Az,dAdz


!     +-------------------------------------------------------------+
!     |   Local variables                                           |
!     +-------------------------------------------------------------+
      integer i

	  
!     +-------------------------------------------------------------+
!     |    Calculation                                              |
!     +-------------------------------------------------------------+

      call ParameterList
      call Morph
      call InitCond

      if (stab.eq.1) cm0=0.5625
      if (stab.eq.2) cm0=0.556170995512353
      cde=cm0*cm0*cm0
      sig_e=kappa*kappa*cm0/(ce2-ce1)/cde 
      L_min  = 0.2  

      txini=0
      tyini=0

      do i=0,xl
         Lini(i)  = L_min
         numini(i)=0
         nuhini(i)=0
      end do

      dragini=kappa/(log(30/K_s*(K_s/30+h(1)/2)))
      dragini=dragini*dragini 

      open (20,status='unknown',file=ForcingName)
      open (30, status='unknown', file=AbsorpName)

	  ! open (41, status='unknown', file=QinpName)
	  ! open (42, status='unknown', file=QoutName)
	  ! open (43, status='unknown', file=TinpName)
	  ! open (44, status='unknown', file=SinpName)

	  call check_advection

      open (80, access = 'SEQUENTIAL', status='unknown',
     &     FORM='unformatted', file=OutName)

      datumini = t_start

      call save_ini
 
      E_Seicheini = 0
     
      if (fgeo.ne.0) then
           do i=1,xl
               fgeo_add(i)=fgeo/rho_0/cp*dAdz(i)/Az(i) ! calculation per kg 
           end do
           if (Az(0) .ne. 0) then
              fgeo_add(1) = fgeo_add(1) + 
     &             fgeo/rho_0/cp*2*Az(0)/((Az(0)+Az(1))*h(1))
           end if
      end if

	fsed = 2.5e-9

	if (fsed.ne.0) then
           do i=1,xl
               fsed_add(i)=fsed*dAdz(i)/Az(i) ! calculation per kg 
           end do
           if (Az(0) .ne. 0) then
              fsed_add(1) = fsed_add(1) + 
     &             fsed*2*Az(0)/((Az(0)+Az(1))*h(1))
           end if
      end if


!  salinity controle for buoyancy functions
      if (ModSal.eq.0) then
          do i=0,xl
            if (Sini(i).ne.0) salctr=1
          end do
          if (salctr.eq.1) then          
             do i=1,xl
                if (Sini(i)-Sini(i-1).ne.0) delsal=1
             end do
          end if
      else
         salctr=1
         delsal=1
      end if
 
      return
      end

!     ###############################################################
      subroutine  Grid
!     ###############################################################

      implicit none

      include 'const_parameter.i'

!     +----------------------------------------------------------------+
!     | Global variables                                               |
!     +----------------------------------------------------------------+
      double precision zu(0:mxl),zk(0:mxl),h(0:mxl),z_zero
      double precision Az(0:mxl),dAdz(0:mxl)
 
      common /morph_variables1/ h,Az,dAdz
      common /morph_varaibles2/ zu,zk,z_zero

!     +----------------------------------------------------------------+
!     | local variables                                                |
!     +----------------------------------------------------------------+
      integer ictr,i
      double precision gr(0:mxl)

      open (15,status='old',file=GridName)
      read(15,*)
      do ictr=0, mxl
         read(15, *, END = 9) gr(ictr)
      end do
 9    if (ictr .eq. mxl) then
         write(*,*) 'Only first ',mxl,' values of file read'
      end if
      close(15)

      if (ictr.eq.1) then                 ! if constant spacing
        xl=gr(0)
        do i=1,xl
           h(i)=depth/xl
        end do

      else                                ! variable spacing
        xl=ictr-1

        if (gr(0).ne.0) then              ! include top value if not included
          xl=xl+1
          do i=xl,1,-1
            gr(i)=gr(i-1)
          end do
          gr(0)=0
        end if

        if (gr(xl).gt.depth) then         ! if maxdepth grid larger than morphology
    
            do while ((gr(xl).gt.depth).and.(xl.gt.0))
                xl=xl-1
            end do
        end if
         
        if (gr(xl).lt.depth) then         ! include bottom value if not included
          xl=xl+1
          gr(xl)=depth
        end if

        do i=1,xl                         ! calculate h
          h(1+xl-i)=gr(i)-gr(i-1)
        end do

      end if

  
      zu(0)=0
      zk(0)=0
      zu(1)=h(1)/2
      zk(1)=h(1)

      do i=2,xl 
         zu(i)=zu(i-1)+0.5*(h(i-1)+h(i))
         zk(i)=zk(i-1)+h(i)
      end do
	  xl_ini = xl
  
      if ((disp_diagnose.eq.1).or.(disp_diagnose.eq.2)) then
	     write(6,*) " Grid File      : "//GridName
	     write(6,*) " No grid Points : ", xl
         write(6,*) " Grid from bottom (0m) to top : "
         write(6,1120) (zk(i), i=0,xl)
         write(6,*)

      end if  
 1120 Format(F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2)
      return
      end

!     ###############################################################
      subroutine morph
!     ###############################################################

      implicit none

      include 'const_parameter.i'
!     +----------------------------------------------------------------+
!     | Global variables                                               |
!     +----------------------------------------------------------------+
      double precision zu(0:mxl),zk(0:mxl),h(0:mxl),z_zero
      double precision Az(0:mxl),dAdz(0:mxl)
 
      common /morph_variables1/ h,Az,dAdz
      common /morph_varaibles2/ zu,zk,z_zero

!     +----------------------------------------------------------------+
!     | Local variables                                                |
!     +----------------------------------------------------------------+
      double precision z(0:mxl), A(0:mxl), zr(0:mxl), Ar(0:mxl)
      integer num, i

      open (15,status='old',file=MorphName)
      read(15,*)							! Read header of depth and area
      do i=0, 1000							! Read depth and area (first 1000 values)
         read(15, *, END = 9)zr(i), Ar(i)
      end do
 9    if (i .eq. 1000) then
         write(*,*) 'Only first 1000 values of file read'
      end if 
      close (15)
	  
	  
      num = i								! Number of area values
      do i=0,num-1							! Reverse order of values
         z(i) = -zr(num-i-1)
         A(i) = Ar(num-i-1)
      end do

      depth = z(0) - z(num-1)				! depth = max. - min. depth
	  z_zero = z(0)							! zero point of z
      do i=0,num-1							! z-coordinate is positive upwards,
         z(i) = z_zero - z(i)				! zero point is at reservoir bottom
      end do

!     create Grid for constant spacing
      call Grid

      call Interp(z, A, num-1, zk, Az, xl)

      do i=1,xl
         dAdz(i) = (Az(i)-Az(i-1))/(zk(i)-zk(i-1))
      end do


      if ((disp_diagnose.eq.1).or.(disp_diagnose.eq.2)) then
         write(6,*) " Morphology File     : "//MorphName
         write(6,*) " No of depth points : ", num
         write(6,*) " Data from File : depth   Area "
         do i=0,num-1
            write(6,'(F8.2, G20.8)') zr(i),Ar(i)
         end do
         write(6,*)
      end if 
	  
      return
      end

!     ###############################################################
      subroutine ParameterList
!     ###############################################################

      implicit none

      include 'const_parameter.i'

      integer i
      double precision Lat
	  
	  
      include 'incl_set_fixparameter.i'       ! this sets the fixed parameters 
	                                        ! for k_epsilon model etc.
	  

!      File with Name and Directory of User Defined Parameters
!        WRITE(*,111)
! 111    FORMAT(' Name of Parameter File ? ' $)

       ! Original version with variable name of parameter file 
       ! read(5,'(A)') ParName
	   ! Version with fixed parameter file
	   ParName = 'kepsilon.par'

!	write(*,*) ParName

!	ParName = 'kepsilon.par'



!     File with Name and Directory of User Defined Parameters
!      open (17,status='old',file='parameter_filename.dat')
!      read(17,'(A)') ParName
!      close(17)
      

!     Reading User defined parameters from file
	  open (16,status='old',file=ParName)  
      read(16,*) 
      read(16,'(A)') InitName
      read(16,'(A)') GridName
      read(16,'(A)') MorphName
      read(16,'(A)') ForcingName
      read(16,'(A)') AbsorpName
      read(16,'(A)') OutName
      read(16,'(A)') zoutName
      read(16,'(A)') toutName
	  read(16,'(A)') QinpName	  
	  read(16,'(A)') QoutName
	  read(16,'(A)') TinpName
	  read(16,'(A)') SinpName
      read(16,*) 
      read(16,*) dt
      read(16,*) t_start
      read(16,*) t_end
      read(16,*) 
      read(16,*) Mod
      read(16,*) stab
      read(16,*) fluxcond
      read(16,*) NBC
      read(16,*) pgrad
      read(16,*) disp_simulation
      read(16,*) disp_diagnose
      read(16,*) 
      read(16,*) alpha_Seiche
      read(16,*) q_NN
      read(16,*) CDeff
      read(16,*) 
      read(16,*) Lat
      read(16,*) CD
      read(16,*) fgeo
      read(16,*) k_min
      read(16,*) air_press
      read(16,*) norm_ctr
      read(16,*) p_radin
      read(16,*) p_windfun  
      read(16,*) igoal
	  read(16,*) ModSal

      close(16)


!     Calculate coriolis parameter from latitude
      Cori=2*7.29e-5*sin(Lat*2*pi/360)
!     change igoal to a different type of information
      if (abs(igoal).eq.10) then
           do i=1,13
              igoal_s(i)=1
           end do
      elseif  (abs(igoal).eq.11) then
           do i=1,3
              igoal_s(i)=1
           end do
      elseif (abs(igoal).eq.12) then
           igoal_s(3)=1
           igoal_s(8)=1
	
      else
           igoal_s(abs(igoal))=1
      end if

      return
      end


!     ###############################################################
      subroutine InitCond
!     ###############################################################

      implicit none

      include 'const_parameter.i'

      double precision Sini(0:mxl),Tini(0:mxl),uini(0:mxl),vini(0:mxl)
      double precision kini(0:mxl),epsini(0:mxl),Lini(0:mxl)
      double precision numini(0:mxl),nuhini(0:mxl) 

      common /ini_values1/ Sini,Tini
      common /ini_values2/ uini,vini
      common /ini_values3/ epsini,kini,Lini
      common /ini_values4/ numini,nuhini

      double precision zu(0:mxl),zk(0:mxl),h(0:mxl),z_zero
	  double precision Az(0:mxl),dAdz(0:mxl)
      common /morph_variables1/ h,Az,dAdz
      common /morph_varaibles2/ zu,zk,z_zero


      double precision z_temp(0:mxl),
     &                 u_temp(0:mxl),v_temp(0:mxl),
     &                 T_temp(0:mxl),S_temp(0:mxl),
     &                 k_temp(0:mxl),eps_temp(0:mxl)
      double precision zini(0:mxl)
	  double precision zini_depth, zmax
      
      integer i,num
      open (2,status='old',file=InitName)		! Opens file with initial conditions
      ! read(2,*)									! Read header of inital depth
      ! read(2,*) zini_depth						! Read initial depth
	  read(2,*)									! Read header of initial u,v, T, etc
      do i=0,mxl								! Read initial u,v,T, etc
         read(2, *, END=9) zini(i), uini(i), vini(i), Tini(i), 
     &                     Sini(i), kini(i), epsini(i)
      end do

 9    num = i-1									! Number of values
      do i=0,num
	     zini(i) = -zini(i)
      end do
	  zini_depth = zini(0)						! Initial depth = first depth of 
												! initial condition list

	  do i=0,xl
	     if (zk(i) > (z_zero-zini_depth)) then	! If above initial water level
		    zmax = zk(i)
		    zk(i) = z_zero-zini_depth
			zu(i) = (zk(i)+zk(i-1))/2
			h(i)  = zk(i) - zk(i-1)
			! Az(i) = Az(i-1) + h(i)/(zmax-zk(i-1))*(Az(i)-Az(i-1))
			xl = i
			if (h(xl) .le. 0.5*h(xl-1)) then	! If top box too small
			   zk(xl-1) = zk(xl)				! Combine two upper boxes
			   zu(xl-1) = (zk(xl-1)+zk(xl-2))/2
			   h(xl-1)  = h(xl)+h(xl-1)
			   xl = xl -1						! reduce number of boxes
			end if
			exit      
		 end if
	  end do      

      do i=0,num
         z_temp(num-i) = z_zero - zini(i)
         u_temp(num-i) = uini(i)
         v_temp(num-i) = vini(i)
         T_temp(num-i) = Tini(i)
         S_temp(num-i) = Sini(i)
         k_temp(num-i) = kini(i)
         eps_temp(num-i) = epsini(i)
      end do

      if (num .eq. 0) then
         write (*,*) 'Only one value!'
         do i=0,xl
            uini(i) = u_temp(0)
            vini(i) = v_temp(0)
            Tini(i) = T_temp(0)
            Sini(i) = S_temp(0)
            kini(i) = k_temp(0)
            epsini(i)= eps_temp(0)
         end do
      else
         call Interp(z_temp, u_temp, num, zu, uini, xl_ini)
         call Interp(z_temp, v_temp, num, zu, vini, xl_ini)
         call Interp(z_temp, T_temp, num, zu, Tini, xl_ini)
         call Interp(z_temp, S_temp, num, zu, Sini, xl_ini)
         call Interp(z_temp, k_temp, num, zk, kini, xl_ini)
         call Interp(z_temp, eps_temp, num, zk, epsini, xl_ini)
      end if
 
      close(2)

	  
	  if ((disp_diagnose.eq.1).or.(disp_diagnose.eq.2)) then
		  write(6,*) " Initial Conditions File : "//InitName
		  do i=num,0,-1
		     write(6,'(F8.2,F10.3,F10.3,F10.3,G14.3,G14.3,G14.3)') 
     &			 depth-z_temp(i),u_temp(i),v_temp(i),
     &                       T_temp(i),k_temp(i),eps_temp(i)
	          end do
                  write(6,*)
	  end if
      return 
      end


 
!     ###############################################################
      subroutine save_ini
!     ###############################################################
      implicit none

      include 'const_parameter.i'

      double precision zu(0:mxl),zk(0:mxl),h(0:mxl), z_zero
      double precision Az(0:mxl),dAdz(0:mxl)
      double precision tout_ctr1(0:9000),tout_ctr2(0:9000)
      integer write_tout
 
      common /morph_variables1/ h,Az,dAdz
      common /morph_varaibles2/ zu,zk,z_zero

      common /savet1  /   write_tout
      common /savet2  /   tout_ctr1
      common /savet3  /   tout_ctr2

!     +---------------------------------------------------------------+
!     |  Global Declarations                                        |
!     +---------------------------------------------------------------+


!     +---------------------------------------------------------------+
!     |  Local Declarations                                        |
!     +---------------------------------------------------------------+

      double precision t_out(0:9000),test(0:mxl)
      integer i,j,tctr

!     reads z-vector for vertical output	  
      open (16,status='old',file=zoutName)
      read(16,*)
      do i=0, 1000
         read(16, *, END = 59) zs(i)
      end do
 59    if (i .eq. 1000) then
         write(*,*) 'Only first 1000 values of file read'
!      else
!         write(*,*) 'Morphology ok'
      end if 
      close (16)
      save_vctr=i-1



      write (80) igoal
      write (80) save_vctr

      if (save_vctr.eq.0) then                      ! takes evry nth index for output   
        j=0
        depth_save=int(zs(0))
        do i=xl,0,-depth_save
          j=j+1
          indexk_save(j)=i
          indexu_save(j)=i
        end do
        xlk=j
        xlu=j
      
        if (indexk_save(xlk).gt.0) then
           xlk=xlk+1
           indexk_save(xlk)=0
        end if
        if (indexu_save(xlu).gt.1) then
           xlu=xlu+1
           indexu_save(xlu)=1
        end if

!       writes vertical vector to file
	  
        if (abs(igoal).gt.9) then
          write (80) xlu
          write (80) (zk(xl)-zu(indexu_save(i)), i=1,xlu)
          write (80) xlk
          write (80) (zk(xl)-zk(indexk_save(i)), i=1,xlk)
        elseif (abs(igoal).lt.5) then
          write (80) xlu
          write (80) (zk(xl)-zu(indexu_save(i)), i=1,xlu)
        elseif (abs(igoal).le.8) then
          write (80) xlk
          write (80) (zk(xl)-zk(indexk_save(i)), i=1,xlk)
        else 
          write(6,*) "wrong choice of index to goal function"
          write(6,*) "        PROGRAM STOPS   "
          return
	end if
	
      else				     ! if depth vector is given by user
          write (80) save_vctr+1
          write (80) (zs(i), i=0,save_vctr)	

          ! do i=0,save_vctr
          !   test(i)=zs(save_vctr-i)
          ! end do
          ! do i=0,save_vctr
          !   !zs(i)=zk(xl)-test(i)
		  !  	zs(i) = z_zero - test(i)
          ! end do
		do i=0, save_vctr
		   test(i) = zs(i)
		end do
		do i=0, save_vctr
		   zs(save_vctr-i) = z_zero + test(i)
	    end do
      end if



!     reads t-vector for temporal output	  
      open (16,status='old',file=toutName)
      read(16,*) 
      do i=0, 9000
         read(16, *, END = 199) t_out(i)
      end do
 199   if (i .eq. 9000) then
         write(*,*) 'Only first 9000 values of file read'
      end if 
      close (16)
      tctr=i-1

       if (tctr.ne.0) then	                 ! indices to save

	   if (t_start.gt.t_out(0)) then
			 write(6,*) " Starting time is larger than first Output time"
			 write(6,*) " Program Stops"
			 goto 9999
           end if	

            
    	   tout_ctr1(0)=(t_out(0)-t_start)*86400.0/dt
           if ( int(tout_ctr1(0)).eq.0 ) then
               tout_ctr2(0)=(t_out(0)-t_start)*86400.0
               test(0)=t_start+tout_ctr2(0)
           else
	       tout_ctr2(0)=((t_out(0)-t_start)*86400.0)/int(tout_ctr1(0))
               test(0)=t_start           
               do j=1,int(tout_ctr1(0))
                 test(0)=test(0)+tout_ctr2(0)/86400.0
               end do
           end if


           do i=1,tctr
	      tout_ctr1(i)=(t_out(i)-test(i-1))*86400.0/dt
              if ( int(tout_ctr1(i)).eq.0 ) then
                 tout_ctr2(i)=(t_out(i)-test(i-1))*86400.0
                 test(i)=test(i-1)+tout_ctr2(i)
                 write(6,*) " Warning: output delta t is smaller than dt
     &                                  for iteration "
              else
	         tout_ctr2(i)=((t_out(i)-test(i-1))*86400.0)/int(tout_ctr1(i))
                 test(i)=test(i-1)           
                 do j=1,int(tout_ctr1(i))
                    test(i)=test(i)+tout_ctr2(i)/86400.0
                 end do
              end if
           end do

	   t_end=test(tctr)
	   write_tout=0
      else                                        ! intervall givenfor output
		 write_tout=int(t_out(0))
      end if


	  
      if ((disp_diagnose.eq.1).or.(disp_diagnose.eq.2)) then
             write(6,*)
	     write(6,*) " Depths of Output : "//zoutName
	     write(6,'(7F8.2)') (depth-zs(i), i=save_vctr,0,-1)

             write(6,*)
             write(6,*) "Time of Output : "//toutName
             if (write_tout.eq.0) then


             write(6,'(6F10.2)') (test(i), i=0,tctr)

             else
                write(6,*) " Intervall steps : ",write_tout 
                write(6,*) " Intervall time  : ",write_tout*dt/86400
             end if

       end if
	 
      
9999  return
      end
          	   


!     ###############################################################
      subroutine check_advection
!     ###############################################################

	  implicit none

      include 'const_parameter.i'

	  !integer adv
	  integer i, j, num_z, filenum(1:4)
	  double precision z_Inp(0:mxl), dummy

	  open (41, status='unknown', file=QinpName)
	  open (42, status='unknown', file=QoutName)
	  open (43, status='unknown', file=TinpName)
	  open (44, status='unknown', file=SinpName)	  

	  filenum = [41,42,43,44] 
      adv = 0

	  do i=1,4
	     read(filenum(i),*)						! Read first row: description of columns			
	     read(filenum(i),*)num_z				! Read number of input depths (static)
	     read(filenum(i),*)dummy, (z_Inp(j), j=1,num_z) ! Read input depths
		 if (eof(filenum(i)) == .TRUE.) then
		    adv = adv + 1
		 end if
		 close(filenum(i))	  
	  end do

	  open (41, status='unknown', file=QinpName)
	  open (42, status='unknown', file=QoutName)
	  open (43, status='unknown', file=TinpName)
	  open (44, status='unknown', file=SinpName)


      if (adv == 4) then
	     adv = 0
	  else
	     adv = 1
	  end if

	  return
	  end
