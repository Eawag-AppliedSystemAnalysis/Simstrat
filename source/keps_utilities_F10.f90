!     ###############################################################
      subroutine Interp(z,y,num_z,zi,yi,num_zi)
!     ###############################################################

      integer num_z, num_zi 
      integer posk1, posk2, posi,i
      double precision z(0:num_z),y(0:num_z),zi(0:num_zi),yi(0:num_zi)
      
      posk1 = 0
      posk2 = num_zi
      posi = 0

      do while (zi(posk1) .le. z(0))
         yi(posk1) = y(0)
         posk1 = posk1+1
      end do 
      
      do while (zi(posk2) .ge. z(num_z))
         yi(posk2) = y(num_z)
         posk2 = posk2-1
      end do

      do i=posk1,posk2 
         do while (zi(i) .gt. z(posi+1))
            posi = posi + 1
         end do
            yi(i)=y(posi)+(zi(i)-z(posi))
     &                              *(y(posi+1)-y(posi))
     &                              /(z(posi+1)-z(posi))
      end do
      
      return
      end


!     ###############################################################
      subroutine Interp_results(z,y,num_z,zi,yi,num_zi)
!     ###############################################################

      integer num_z, num_zi 
      integer posk1, posk2, posi,i
      double precision z(0:num_z),y(0:num_z),zi(0:num_zi),yi(0:num_zi)
	  double precision c1,c2
      
      posk1 = 0
      posk2 = num_zi
      posi = 0

      do while (zi(posk1) .lt. z(0))
	     c1 = 0
		 c2 = 0
         yi(posk1) = c1/c2						! NaN
         posk1 = posk1+1
      end do 
      
      do while (zi(posk2) .gt. z(num_z))
	     c1 = 0
		 c2 = 0
         yi(posk2) = c1/c2						! NaN
         posk2 = posk2-1
      end do

      do i=posk1,posk2 
         do while (zi(i) .gt. z(posi+1))
            posi = posi + 1
         end do
            yi(i)=y(posi)+(zi(i)-z(posi))
     &                              *(y(posi+1)-y(posi))
     &                              /(z(posi+1)-z(posi))
      end do
      
      return
      end

!     #################################################################
      subroutine Forcing(datum,T,first,tx,ty,u_taus,I_0,heat
     &                              ,SST,u10,v10) 
!     #################################################################

      implicit none

      include 'const_parameter.i'

!     +----------------------------------------------------------------+
!     | Global variables                                               |
!     +----------------------------------------------------------------+

      double precision datum,T(0:xl)
      double precision tx,ty,u_taus,I_0,heat,SST,u10,v10
      integer first

!     +----------------------------------------------------------------+
!     | Local variables                                                |
!     +----------------------------------------------------------------+
      double precision tau
      double precision tb_s, tb_e
      double precision A_start(7), A_end(7), A_cur(7)
      double precision r_a, r_s, B0, fu, Vap_wat
      double precision T_atm, F_glob, Vap_atm, Cloud

      static tb_s, tb_e
      static A_start, A_end, A_cur

      integer eof
      static eof

      if (NBC .eq. 1) then
         call ReadForcing (20, datum, tb_s, tb_e, 
     &             A_start, A_end, A_cur, 4, eof, first)
         u10 = A_cur(1)
         v10 = A_cur(2)
         SST = A_cur(3)
         I_0 = A_cur(4)
         heat = 0
      else if (NBC .ge. 2) then

  
         if (NBC.eq.2) then			! date, U,V,Tatm,Hsol,Vap
              call ReadForcing (20, datum, tb_s, tb_e, 
     &             A_start, A_end, A_cur, 5, eof, first)
            u10 = A_cur(1)
            v10 = A_cur(2)
            T_atm = A_cur(3)
            F_glob = A_cur(4)
            Vap_atm = A_cur(5)
            Cloud = 0.5             
         else if (NBC.eq.3) then		! date,U10,V10,Tatm,Hsol,Vap,Clouds	
            call ReadForcing (20, datum, tb_s, tb_e, 
     &             A_start, A_end, A_cur, 6, eof, first)
            u10 = A_cur(1)
            v10 = A_cur(2)
            T_atm = A_cur(3)
            F_glob = A_cur(4)
            Vap_atm = A_cur(5)
            Cloud = A_cur(6)     
	   else if (NBC.eq.4) then		! date,U10,V10,Hnet,Hsol 
		  call ReadForcing (20, datum, tb_s, tb_e, 
     &             A_start, A_end, A_cur, 4, eof, first)
	      u10 = A_cur(1)
            v10 = A_cur(2)
            heat = A_cur(3)
            F_glob = A_cur(4)
         end if
  
         r_a = 0.03            ! Ratio of reflected to total 
                               ! long-wave iradiance
         r_s = 0.20            ! Fraction of reflected short-wave
                               ! radiation
         B0 = 0.61             ! Bowen constant

                               ! Long-wave radiation from sky
                               ! Livingstone and Imboden 1989
         H_A =   (1+0.17*Cloud*Cloud)
     &         *  1.24*(Vap_atm/(T_atm+273.15))**(1.0/7.0)
     &         *  5.67e-8*(1-r_a)*(273.15+T_atm)**4
         H_A = H_A*p_radin     ! fitting factor p_radin (e.g. 1.09 for Aegerisee)

                               ! Wind function according to Livingstone and Imboden (1989)
         fu = 4.4+1.82*sqrt(u10*u10+v10*v10)+0.26*(T(xl)-T_atm)
!         fu = 5.44+2.19*sqrt(u10*u10+v10*v10)+0.24*(T(xl)-T_atm)
         fu=fu*p_windfun       ! fitting function (eg. 1.09 for Aegerisee)

                               ! Flux of sensible heat (convection)
         H_K = -B0*fu*(T(xl)-T_atm)
                               ! Water vapor saturation pressure in air at
                               ! temperature of water (nach Gill 1992) (pressure in milibar)
         Vap_wat =  10**( (0.7859+0.03477*T(xl)) / (1+0.00412*T(xl)) )
         if (air_press.ne.0) then
           Vap_wat=Vap_wat
     &         *(1+1.0e-6*air_press*(4.5+0.00006*T(xl)*T(xl)))
         end if

                               ! Flux of latent heat (evaporation,
                               ! condensation)
         H_V = -fu*(Vap_wat-Vap_atm)
                               ! Long-wave radiation from water body
         H_W = -0.97*5.67e-8*(T(xl)+273.15)**4

         I_0 = (1-r_s)*F_glob      ! Solar short-wave radiation absorbed
                               ! in the water column
	   if (NBC.ne.4) then
            heat = H_A + H_W + H_V + H_K
	   end if
         if ( (T(xl).lt.0).and.(heat.lt.0) )  heat=0	

      end if


	if (CD .lt. 0.0) then
		if (abs(u10) .le. 4.7 .and. abs(u10) .gt. 0.5) then
			C10 = 1.78e-3*(sqrt(u10*u10+v10*v10))**(-1.3561)
		elseif (abs(u10) .gt. 4.7) then
			C10 = 2.13e-4
		else
			C10 = 0.0046
		end if
	else
		C10 = CD
	end if

      tau = (u10**2+v10**2)*rho_air/rho_0*C10
      if (u10 .eq. 0.0) then
	 tx = 0
         ty = tau
      else if (v10 .eq. 0.0) then
         ty = 0
         tx = tau
      else
         tx = v10/abs(v10)*tau/sqrt(1+(u10/v10)**2)
         ty = u10/abs(u10)*tau/sqrt(1+(v10/u10)**2)
      end if

      u_taus = sqrt(tau)

!	 if (u10 .eq. 0.0) then
!	    tx = 0
!         else if (abs(u10) .lt. 4.7 .and. abs(u10) .gt. 0.1) then
!	 tx = u10*u10*u10/abs(u10)*1.2/1000*2.38e-3*abs(u10)**-1.5
!	 else
!	 tx = u10*u10*u10/abs(u10)*1.2/1000*2.34e-4
!	 end if
	 
!	 if (v10 .eq. 0.0) then
!	    ty = 0
!	 else if (abs(v10) .lt. 4.7 .and. abs(v10) .gt. 0.1) then
!	 ty = v10*v10*v10/abs(v10)*1.2/1000*2.38e-3*abs(v10)**-1.5
!	 else
!	 ty = v10*v10*v10/abs(v10)*1.2/1000*2.34e-4
!	 end if

      return
      end


!     ###############################################################
      subroutine Absorption(datum, ga1, zk, first)
!     ###############################################################

      implicit none

      include 'const_parameter.i'

      double precision zk(0:xl)


      double precision datum, ga1(0:xl)
      double precision tb_s, tb_e, ratio
      double precision z_ga1(0:mxl), dummy
      double precision ga1_rs(0:mxl), ga1_re(0:mxl)
      double precision ga1_s(0:mxl), ga2_s(0:mxl)
      double precision ga1_e(0:mxl), ga2_e(0:mxl)
      integer first, eof, i,num_z

      static tb_s, tb_e
      static ga1_s, ga1_e, ga2_s, ga2_e, z_ga1
      static eof, num_z

       if (first .eq. 1) then
         eof=0 
         read(30,*, end=10) 
		 
         read(30, *)num_z
         read(30, *)dummy, (z_ga1(i),i=0,num_z-1)
		 do i=0,num_z-1
		    z_ga1(i) = -z_ga1(i)
		 end do

         read(30,*, end=8) tb_s, (ga1_rs(i),i=0,num_z-1)
         read(30,*, end=9) tb_e, (ga1_re(i),i=0,num_z-1)
 9       continue
         call Interp(z_ga1, ga1_rs, num_z-1, zk, ga1_s, xl)
         call Interp(z_ga1, ga1_re, num_z-1, zk, ga1_e, xl)

         if(datum .le. tb_s) then          ! if datum before first date
		  write(6,*) " First datum in Absorption is 
     &		  larger than that required !!"

            do i=0, xl
               ga1(i) = ga1_s(i)           ! take first value
            end do
         else if (datum .ge. tb_s) then    ! if datum after last date
            do i=0, xl
               ga1(i) = ga1_e(i)           ! take last value
            end do
         else
            do while (datum .gt. tb_e)     ! until correct interval
			   tb_s = tb_e				   
               do i=0, xl
                  ga1_s(i) = ga1_e(i)
               end do
               read(30,*, end=10) tb_e, (ga1_re(i),i=0,num_z-1)
               call Interp(z_ga1, ga1_re, num_z-1, zk, ga1_e, xl)
            end do
            ratio = (datum-tb_s)/(tb_e-tb_s)
            do i=0, xl
               ga1(i) = ga1_s(i)+ratio*(ga1_e(i) - ga1_s(i))
            end do

         end if
                  
         if ((disp_diagnose.eq.1).or.(disp_diagnose.eq.2)) then
           write(6,*) " Absorption : "//AbsorpName
          do i=xl,xl-4,-1
              write(6,'(F12.3,G12.4)')depth-zk(i),ga1(i)
           end do
           write(6,*) " ... "
           write(6,*) " ... "
          do i=5,1,-1
             write(6,'(F12.3,G12.4)') depth-zk(i),ga1(i)
           end do
           write(6,*)
           write(6,*) " --------------------------"
           write(6,*) "   SIMULATION IN PROGRESS"
           write(6,*) " --------------------------"
           write(6,*)

         end if

      elseif (datum .gt. tb_e .and. eof .eq. 1) then
         do i=0, xl
            ga1(i) = ga1_e(i)
         end do

      else
         do while (datum .gt. tb_e)
            tb_s = tb_e
            do i=0, xl
                ga1_s(i) = ga1_e(i)
            end do
            read(30,*, end=10) tb_e, (ga1_re(i),i=0,num_z-1)
            call Interp(z_ga1, ga1_re, num_z-1, zk, ga1_e, xl)

         end do
         ratio = (datum-tb_s)/(tb_e-tb_s)
         do i=0, xl
            ga1(i) = ga1_s(i)+ratio*(ga1_e(i) - ga1_s(i))
         end do

      end if

      return

 8    call Interp(z_ga1, ga1_rs, num_z-1, zk, ga1, xl)   

 10   do i=0,xl
         ga1(i) = ga1_e(i)
      end do

      return
      end


!     ###############################################################
      subroutine ReadForcing (FileNumber, datum, tb_s, tb_e, 
     &     A_start, A_end, A_cur, size, eof, first)
!     ###############################################################

	implicit none

      include 'const_parameter.i'

      integer FileNumber, size, eof, first
      double precision datum, tb_s, tb_e
      double precision A_start(size), A_end(size), A_cur(size)
      double precision ratio
!      character*500 string
      integer i

      if (first .eq. 1) then
         eof=0
         rewind(Filenumber)
         read(FileNumber,*, end=10)
         read(FileNumber,*, end=10) tb_s, (A_end(i),i=1,size) 

         if ((disp_diagnose.eq.1).or.(disp_diagnose.eq.2)) then
           write(6,*) " Forcing : "//ForcingName
           write(6,'(4G14.5)') tb_s, (A_end(i),i=1,size)
           write(6,*)
         end if
         read(FileNumber,*, end=10) tb_e,   (A_end(i),i=1,size)

         if(datum .lt. tb_s) then
			write(6,*) " First datum in Forcing is larger than that required !!"
			
!            do i=1,size
!               A_cur(i) = A_start(i)
!            end do
         else
            do while (datum .gt. tb_e)
               tb_s = tb_e
               do i=1,size
                  A_start(i) = A_end(i)
               end do
               read(FileNumber,*, end=10) tb_e, (A_end(i),i=1,size) 
            end do
            ratio = (datum-tb_s)/(tb_e-tb_s)
            do i=1,size
               A_cur(i) = A_start(i)+ratio*(A_end(i) - A_start(i))
            end do
         end if
      else if (datum .gt. tb_e .and. eof .eq. 1) then
         do i=1,size
            A_cur(i) = A_end(i)
         end do
      else
         do while (datum .gt. tb_e)
            tb_s = tb_e
            do i=1,size
               A_start(i) = A_end(i)
            end do
            read(FileNumber,*, end=10) tb_e,   (A_end(i),i=1,size)
         end do
         ratio = (datum-tb_s)/(tb_e-tb_s)
         do i=1,size
            A_cur(i) = A_start(i)+ratio*(A_end(i) - A_start(i))
         end do    
         if (disp_diagnose.eq.2) then
           write(6,'(" Forcing : ",4G14.5)') tb_s, (A_end(i),i=1,size)
           write(6,*)
         end if

      end if

      return

 10   do i=1,size
         A_cur(i) = A_start(i)
      end do
      eof = 1


      return
      end


!     ###############################################################
      subroutine Lateral(datum,first,zu,z_zero,h,Qvert,Q_inp)
!     ###############################################################
  
	  implicit none

      include 'const_parameter.i'

!     +---------------------------------------------------------------+
!     |  Global Declarations                                        |
!     +---------------------------------------------------------------+
	  integer first
	  double precision datum, zu(0:xl), h(0:xl), z_zero, Qvert(0:xl), Q_inp(1:4,0:mxl)

!     +---------------------------------------------------------------+
!     |  Local Declarations                                        |
!     +---------------------------------------------------------------+
	  integer i, j, num_z(1:4), filenum(1:4)
	  double precision z_Inp(1:4,0:mxl), dummy
	  double precision Inp_rs(1:4,0:mxl), Inp_re(1:4,0:mxl)
	  double precision Q_s(1:4,0:mxl), Q_e(1:4,0:mxl)
	  double precision Q_rs(1:4,0:mxl), Q_re(1:4,0:mxl)
	  double precision tb_s(1:4), tb_e(1:4), ratio(1:4)

	  static Q_s, Q_e
	  static num_z

	  filenum = [41,42,43,44]							
	  do i=1,4
		 if (first .eq. 1) then
		    if (eof(filenum(i)) .eq. .FALSE.) then
		       read(filenum(i),*)						! Read first row: description of columns			
			   read(filenum(i),*)num_z(i)				! Read number of input depths (static)
			   read(filenum(i),*)dummy, (z_Inp(i,j), j=1,num_z(i)) ! Read input depths 
			   !read(filenum(i),*)dummy, (z_Inp(i,j), j=num_z(i),1,-1) ! Read input depths 

			   do j=1,num_z(i)							! Transfer coordinate system
			      !z_Inp(i,j) = z_zero - z_Inp(i,j)
				  z_Inp(i,j) = z_zero + z_Inp(i,j)
			   end do

			   if (eof(filenum(i)) .eq. .TRUE.) then	! If empty -> Inp(i,all) = 0
			      do j=0,xl
				     Q_inp(i,j) = 0
			      end do
			   else										! else read start time & inputs
			      read(filenum(i),*)tb_s(i), (Inp_rs(i,j),j=1,num_z(i))	
				  !read(filenum(i),*)tb_s(i), (Inp_rs(i,j),j=num_z(i),1,-1)	

				  call integrate(z_Inp(i,1:num_z(i)), Inp_rs(i,1:num_z(i)), 
     &								Q_rs(i, 1:num_z(i)), num_z(i))
			      if (eof(filenum(i)) .eq. .FALSE.) then		! if eof == FALSE: read end time
			         read(filenum(i),*)tb_e(i), (Inp_re(i,j),j=1,num_z(i))	
					 !read(filenum(i),*)tb_e(i), (Inp_re(i,j),j=num_z(i),1,-1)	

					 call integrate(z_Inp(i,1:num_z(i)), Inp_re(i,1:num_z(i)), 
     &								Q_re(i, 1:num_z(i)), num_z(i))		   
					 call Interp(z_Inp(i,1:num_z(i)), Q_rs(i,1:num_z(i)), 
     &								num_z(i)-1, zu, Q_s(i,0:xl), xl)
				     call Interp(z_Inp(i,1:num_z(i)), Q_re(i,1:num_z(i)), 
     &								num_z(i)-1, zu, Q_e(i,0:xl), xl)

			         if(datum .le. tb_s(i)) then				! if datum before first date
			            do j=1,xl								! take start Input
						   Q_inp(i,j) = Q_s(i,j)
				        end do
					 else
					    do while ( .not. ((datum .ge. tb_s(i)) 	! do until datum between dates
     &							.and. (datum .le. tb_e(i))))
						   tb_s(i) = tb_e(i)					! move one step in time
						   do j=1, xl
							  Q_s(i,j) = Q_e(i,j)
                           end do
						   if (eof(filenum(i)) .ne. .TRUE.) then
                              read(filenum(i),*) tb_e(i), (Inp_re(i,j),j=1,num_z(i))
							  !read(filenum(i),*) tb_e(i), (Inp_re(i,j),j=num_z(i),1,-1)
							  call integrate(z_Inp(i,1:num_z(i)), Inp_re(i,1:num_z(i)), 
     &								Q_re(i, 1:num_z(i)), num_z(i))		   

	 						  call Interp(z_Inp(i,1:num_z(i)), Q_re(i,1:num_z(i)), 
     &										num_z(i)-1, zu, Q_e(i,0:xl), xl)
						   else							! if eof: exit do while loop
						      exit
					       end if
						end do

                        ratio(i) = (datum-tb_s(i))/(tb_e(i)-tb_s(i))
                        do j=1, xl
						   Q_inp(i,j) = Q_s(i,j)+ratio(i)*(Q_e(i,j) - Q_s(i,j))
                        end do
			         end if

				  else									! Interpolate Q_rs -> Q_s
					 call Interp(z_Inp(i,1:num_z(i)), Q_rs(i,1:num_z(i)), num_z(i)-1, zu, 
     &								Q_s(i,0:xl), xl)

					 do j=1,xl							! Q_inp = Q_s
						Q_inp(i,j) = Q_s(i,j)
		             end do
			      end if

			   end if

			end if								! if (eof(filenum(i)) .eq. .FALSE.) then

		 else									! If not first time:
		    do while ( .not. ((datum .ge. tb_s(i)) 	! do until datum between dates
     &							.and. (datum .le. tb_e(i))))
			   do j=1, xl
				  Q_s(i,j) = Q_e(i,j)
               end do
			   if (eof(filenum(i)) .ne. .TRUE.) then
                  read(filenum(i),*) tb_e(i), (Inp_re(i,j),j=1,num_z(i))
				  !read(filenum(i),*) tb_e(i), (Inp_re(i,j),j=num_z(i),1,-1)

				  call integrate(z_Inp(i,1:num_z(i)), Inp_re(i,1:num_z(i)), 
     &								Q_re(i, 1:num_z(i)), num_z(i))		   
	 			  call Interp(z_Inp(i,1:num_z(i)), Q_re(i,1:num_z(i)), 
     &										num_z(i)-1, zu, Q_e(i,0:xl), xl)
			   else								! if eof: exit do while loop
			      exit
			   end if
			end do

            ratio(i) = (datum-tb_s(i))/(tb_e(i)-tb_s(i))
		    do j=xl+1,mxl
			   Q_e(i,j) = Q_e(i,xl)
			   Q_s(i,j) = Q_s(i,xl)
			end do
            do j=1, xl
			   Q_inp(i,j) = Q_s(i,j)+ratio(i)*(Q_e(i,j) - Q_s(i,j))
            end do
	     end if									! end if first .eq. 1

      end do									! do i=1,4

	  do i=1,xl
	     Qvert(i) = Q_inp(1,i)+Q_inp(2,i)
	  end do
	  do i=1,4
	     do j=1,xl-1
		    Q_inp(i,xl-j+1) = (Q_inp(i,xl-j+1)-Q_inp(i,xl-j)) 
		 end do
	  end do

	  do i=1,xl-1								! Outflow always at bottom; otherwise 
	     Q_inp(2,i) = Q_inp(2,i+1)				! vertical advection and outflow not
	  end do									! consistent (preliminary version)!
	  
	  return
	  end
	 

!     ###############################################################
      subroutine Integrate(x,y,inty, size)
!     ###############################################################

	  implicit none

      include 'const_parameter.i'

	  double precision x(1:mxl), y(1:mxl), inty(1:mxl)
	  integer size, i

	  inty(1) = 0
	  do i=2,size
	     inty(i) = inty(i-1) + (x(i)-x(i-1))*(y(i)+y(i-1))/2
	  end do

	  return
	  end
