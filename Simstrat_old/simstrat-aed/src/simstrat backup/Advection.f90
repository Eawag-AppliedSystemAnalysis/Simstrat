module Advection
  use SimstratModel
  use utilities
  implicit none
  save

  real(dp), parameter :: slope = 0.05_dp
  real(dp), parameter :: phi = pi/4
  real(dp), parameter :: sinatan_slope = sin(atan(slope))
  real(dp), parameter :: tan_phi = tan(phi)
  

  private
  public doAreaFactors, doAdvection, doLateral, doLateral_rho

contains

  !####################################################################
  subroutine doAreaFactors()
  !####################################################################
      
      implicit none
      

      integer :: i

      AreaFactor_1(1:nz) = -4*Az(0:nz-1)/(h(1:nz)+h(0:nz-1))/h(1:nz)/(Az(1:nz)+Az(0:nz-1))
      AreaFactor_2(1:nz) = -4*Az(1:nz)/(h(1:nz)+h(2:nz+1))/h(1:nz)/(Az(1:nz)+Az(0:nz-1))
      AreaFactor_k1(1:nz-1) = -(Az(1:nz-1)+Az(2:nz))/(h(1:nz-1)+h(2:nz))/h(2:nz)/Az(1:nz-1)
      AreaFactor_k2(1:nz-1) = -(Az(1:nz-1)+Az(0:nz-2))/(h(1:nz-1)+h(2:nz))/h(1:nz-1)/Az(1:nz-1)
      AreaFactor_eps(1:nz-1) = 0.5_dp*((Az(1:nz-1)-Az(0:nz-2))/h(1:nz-1)+(Az(2:nz)-Az(1:nz-1))/h(2:nz))/Az(1:nz-1)

      meanint(0:nz-1) = 2.0_dp/(h(0:nz-1)+h(1:nz))

      volume=0
      do i=0,nz-1
          volume = volume + 0.5_dp*h(i+1)*(Az(i)+Az(i+1))
      end do

      return
  end subroutine doAreaFactors

  !####################################################################
  subroutine doAdvection(Q_vert,Q_inp,U,V,T,S,k,eps)
  !####################################################################
    
    implicit none

    !Global variables
    real(dp), intent(inout) :: U(0:), V(0:), T(0:), S(0:)
    real(dp), intent(inout) :: k(0:), eps(0:)
    real(dp), intent(inout) :: Q_vert(0:) ! Vertical advection [m3/s] (exchange between boxes)
    real(dp), intent(in) :: Q_inp(1:,0:) ! Inflow [m3/s], Outflow [m3/s], T-input [°C*m3/s], S-input [‰*m3/s]

    ! Local variables
    real(dp) :: dh, dh_i(1:2), h_div_2, h_mult_2      ! depth differences
    real(dp) :: dt_i(1:2)          ! first and second time step
    real(dp) :: dU(0:nz_grid), dV(0:nz_grid), dTemp(0:nz_grid), dS(0:nz_grid)
    real(dp) :: AreaFactor_adv(1:nz_grid)
    real(dp) :: top

    integer :: i, t_i
    
    lake_level_old = z_upp(nz)
    !Depth difference compared to previous timestep
    dh = Q_vert(nz)/Az(nz)*dt
    h_div_2 = 0.5_dp*h(nz-1)    ! Take second highest box since the top box might not be at the full height
    h_mult_2 = 2_dp*h(nz-1)

    !Split timestep depending on situation
    if (dh==0.) then                          ! If volume does not change, take one normal time step
        dt_i(1) = dt
    else if ((dh+z_upp(nz))>=depth) then        ! If surface level reached, take a step until surface
        dt_i(1) = (depth - z_upp(nz))/dh*dt
    else if (((dh+h(nz))>h_div_2) .and.&   ! If top box>0.5*lower box and <2*lower box, take one time step
            ((dh+h(nz))<h_mult_2)) then
        dt_i(1) = dt
    else if ((dh+h(nz))<=h_div_2) then     ! If top box<=0.5*lower box, first step until top box=0.5*lower box
        dt_i(1) = abs((h(nz)-h_div_2)/dh)*dt
    else                                     ! If top box>=2*lower box, first step until top box = 2*lower box
        dt_i(1) = abs((2*h(nz-1)-h(nz))/dh)*dt
    end if
    dt_i(2) = dt-dt_i(1)                       ! Rest of timestep

    ! FB 2016: Revision
    do t_i=1,2 !First and (if needed) second timestep
        AreaFactor_adv(1:nz) = dt_i(t_i)/(Az(1:nz)*h(1:nz))     ! Area factor for dt(t_i)
        dh_i(t_i) = dh*dt_i(t_i)/dt                         ! Depth difference for dt(t_i)

        do i=1,nz
            if (i==nz .and. Q_vert(i)>0) then
                top = 0
            else
                top = 1
            end if
             ! Advective flow out of box i, always negative
            dU(i)    = -top*abs(Q_vert(i))*U(i)
            dV(i)    = -top*abs(Q_vert(i))*V(i)
            dTemp(i) = -top*abs(Q_vert(i))*T(i)
            dS(i)    = -top*abs(Q_vert(i))*S(i)
            if (i>1 .and. Q_vert(i-1)>0) then       ! Advective flow into box i, from above
                dU(i)    = dU(i) + Q_vert(i-1)*U(i-1)
                dV(i)    = dV(i) + Q_vert(i-1)*V(i-1)
                dTemp(i) = dTemp(i) + Q_vert(i-1)*T(i-1)
                dS(i)    = dS(i) + Q_vert(i-1)*S(i-1)
            end if
            if (i<nz .and. Q_vert(i+1)<0) then      ! Advective flow into box i, from below
                dU(i)    = dU(i) - Q_vert(i+1)*U(i+1)
                dV(i)    = dV(i) - Q_vert(i+1)*V(i+1)
                dTemp(i) = dTemp(i) - Q_vert(i+1)*T(i+1)
                dS(i)    = dS(i) - Q_vert(i+1)*S(i+1)
            end if
        end do

        ! dT = dT(vertical advection) + dT(inflow) + dT(negative outflow), units: °C*m^3/s
        dTemp(1:nz) = dTemp(1:nz) + (Q_inp(3,1:nz)+Q_inp(2,1:nz)*T(1:nz))
        ! dS = dS(vertical advection) + dS(inflow) + dT(negative outflow), units: ‰*m^3/s
        dS(1:nz) = dS(1:nz) + (Q_inp(4,1:nz)+Q_inp(2,1:nz)*S(1:nz))

        ! Add change to the state variable
        U(1:nz) = U(1:nz) + AreaFactor_adv(1:nz)*dU(1:nz)
        V(1:nz) = V(1:nz) + AreaFactor_adv(1:nz)*dV(1:nz)
        T(1:nz) = T(1:nz) + AreaFactor_adv(1:nz)*dTemp(1:nz)
        S(1:nz) = S(1:nz) + AreaFactor_adv(1:nz)*dS(1:nz)

        ! Variation of variables due to change in volume
        U(nz) = U(nz)*h(nz)/(h(nz)+dh_i(t_i))
        V(nz) = V(nz)*h(nz)/(h(nz)+dh_i(t_i))
        T(nz) = T(nz)*h(nz)/(h(nz)+dh_i(t_i))
        S(nz) = S(nz)*h(nz)/(h(nz)+dh_i(t_i))

        if (dh==0) then                          ! If volume does not change, return
          return
        else if ((dh+z_upp(nz))>=depth) then        ! If surface level reached
          if (.not. (z_upp(nz)==depth)) then
            h(nz)       = h(nz) + dh_i(t_i)                      ! New size of top box
            z_cent(nz)  = z_cent(nz) + 0.5_dp*dh_i(t_i)                  ! New centre coordinate of top box
            z_upp(nz)   = z_upp(nz) + dh_i(t_i)                    ! New surface coordinate of top box
          else
            return
          end if

        else if (((dh_i(t_i) + h(nz))>h_div_2) .and.& 
                ((dh_i(t_i) + h(nz))<(h_mult_2))) then    ! and top box<2*lower box
            h(nz)       = h(nz) + dh_i(t_i)                      ! New size of top box
            z_cent(nz)  = z_cent(nz) + 0.5_dp*dh_i(t_i)                  ! New centre coordinate of top box
            z_upp(nz)   = z_upp(nz) + dh_i(t_i)                    ! New surface coordinate of top box

        else if (t_i==1 .and. (dh+h(nz))<=h_div_2) then     ! If top box<=0.5*lower box, merge 2 boxes
            z_upp(nz-1)   = z_upp(nz) + dh_i(t_i)
            z_cent(nz-1)  = 0.5_dp*(z_upp(nz-1) + z_upp(nz-2))
            h(nz-1)       = (z_upp(nz-1) - z_upp(nz-2))

            ! New values of the state variables are weighted averages
            U(nz-1) = (0.5_dp*U(nz)*Az(nz)+U(nz-1)*Az(nz-1))/(0.5_dp*Az(nz)+Az(nz-1))
            V(nz-1) = (0.5_dp*V(nz)*Az(nz)+V(nz-1)*Az(nz-1))/(0.5_dp*Az(nz)+Az(nz-1))
            T(nz-1) = (0.5_dp*T(nz)*Az(nz)+T(nz-1)*Az(nz-1))/(0.5_dp*Az(nz)+Az(nz-1))
            S(nz-1) = (0.5_dp*S(nz)*Az(nz)+S(nz-1)*Az(nz-1))/(0.5_dp*Az(nz)+Az(nz-1))
            k(nz-1) = (0.5_dp*k(nz)*Az(nz)+k(nz-1)*Az(nz-1))/(0.5_dp*Az(nz)+Az(nz-1))
            eps(nz-1) = (0.5_dp*eps(nz)*Az(nz)+eps(nz-1)*Az(nz-1))/(0.5_dp*Az(nz)+Az(nz-1))
            Q_vert(nz-1) = (0.5_dp*Q_vert(nz)*Az(nz)+Q_vert(nz-1)*Az(nz-1))/(0.5_dp*Az(nz)+Az(nz-1))
            nz = nz-1                     ! Reduce number of boxes
            call doAreaFactors()

        else if (t_i==1 .and. (dh+h(nz))>=h_mult_2) then    ! If top box>=2*lower box, add one box
            z_upp(nz+1)   = z_upp(nz) + dh_i(t_i)
            z_upp(nz)     = z_upp(nz) - (h(nz-1)-dh_i(t_i))
            z_cent(nz+1)  = z_cent(nz) + 0.5_dp*dh_i(t_i)
            z_cent(nz)    = z_cent(nz) - 0.5_dp*(h(nz-1)-dh_i(t_i))
            h(nz+1)       = z_upp(nz+1) - z_upp(nz)
            h(nz)         = z_upp(nz) - z_upp(nz-1)

            ! State variables in the newly createdbox are the same as in the box below
            U(nz+1)   = U(nz)
            V(nz+1)   = V(nz)
            T(nz+1)   = T(nz)
            S(nz+1)   = S(nz)
            k(nz+1)   = k(nz)
            eps(nz+1) = eps(nz)
            Q_vert(nz+1) = Q_vert(nz)       ! Vertical discharge of new box
            nz = nz + 1                   ! Increase number of boxes
            call doAreaFactors()

        end if ! dh==0
    end do      !end do t_i=1,2

    return
  end subroutine doAdvection

  !Read and compute inflow/ouflow parameters
  !AG 2014: revision + correction
  !FB 2016: revision + correction
  !####################################################################
  subroutine doLateral(datum,idx,Q_vert,Q_inp)
  !####################################################################
        implicit none
        

        ! Global Declarations
        real(dp), intent(in) :: datum
        ! Q_inp is the input at each depth for each time step, Q_vert is the integrated net water input
        real(dp), intent(inout) :: Q_inp(1:,0:), Q_vert(0:)    ! Zero is not used
        integer, intent(in) :: idx

        ! Local Declarations
        real(dp) :: tb_start(1:4), tb_end(1:4)        ! Input depths, start time, end time
        real(dp) :: Q_read_start(1:4,0:nz_input), Q_read_end(1:4,0:nz_input)      ! Integrated input
        real(dp) :: dummy

        integer :: nval(1:4), nval_deep(1:4), nval_surface(1:4)      ! Number of values
        integer :: i, j, nz_old, eof(1:4)   ! indexes, nz from previous time step, end of file
        integer :: fnum(1:4)  ! File number
        character*20 :: fname(1:4)

        ! Variables that need to keep their values when returning to the main program (MS 2014)
        save tb_start, tb_end, nz_old, nval, eof
        save nval_surface, nval_deep
        ! z_Inp, Inp_read_start,Inp_read_end, Q_start and Q_end are saved automatically as they are
        ! global dynamic arrays (FB 2016)

        fname = ['inflow           ','outflow          ','input temperature','input salinity   ']
        fnum = [41,42,43,44]

        ! FB 2016: Major revision to include surface inflow
        ! Do this for inflow, outflow, temperature and salinity
        do i=1,4
          if (idx==1) then   ! First iteration

            ! Allocate arrays if not already done, this saves memory compared to declaring with nz_max
            if(.not.allocated(z_Inp)) allocate(z_Inp(1:4,0:nz_input))     ! Input depths  
            if(.not.allocated(Inp_read_start)) allocate(Inp_read_start(1:4,0:nz_input))
            if(.not.allocated(Inp_read_end)) allocate(Inp_read_end(1:4,0:nz_input))
            if(.not.allocated(Q_start)) allocate(Q_start(1:4,1:nz_grid))  ! Input interpolated on grid
            if(.not.allocated(Q_end)) allocate(Q_end(1:4,1:nz_grid))      ! Input interpolated on grid

            ! Open file and start to read
            if(disp_dgn==2) write(6,*) 'Starting to read '//trim(fname(i))//' file...'
            eof(i) = 0
            read(fnum(i),*,end=9)      ! Skip first row: description of columns

            ! Read number of deep and surface inflows
            read(fnum(i),*,end=9) nval_deep(i), nval_surface(i)

            ! Total number of values to read
            nval(i) = nval_deep(i) + nval_surface(i)

            ! Read input depths and convert coordinate system
            read(fnum(i),*,end=9) dummy, (z_Inp(i,j),j=0,nval(i)-1)

            z_Inp(i,0:nval_deep(i)-1) = z_zero + z_Inp(i,0:nval_deep(i)-1)

            ! Allocate array for depth of the surface inflow (only done for i=1)
            if(.not.allocated(depth_surfaceFlow)) allocate(depth_surfaceFlow(1:4,1:(nval_surface(i)+2)))

            ! Read first line
            read(fnum(i),*,end=9) tb_start(i),(Inp_read_start(i,j),j=0,nval(i)-1)

            ! Integrate the inflow (direct interpolation of inflow is not correct)
            call Integrate(z_Inp(i,:),Inp_read_start(i,:),Q_read_start(i,:),nval_deep(i))

            ! Very important: once the inflowing quantitiy is integrated, it necessarily has to be
            ! interpolated on the z_upp grid starting with index 1
            call Interp(z_Inp(i,:),Q_read_start(i,:),nval_deep(i)-1,z_upp(1:nz),Q_start(i,1:nz),nz)

            ! Read second line and treatment of deep inflow
            read(fnum(i),*,end=7) tb_end(i),(Inp_read_end(i,j),j=0,nval(i)-1)
            call Integrate(z_Inp(i,:),Inp_read_end(i,:),Q_read_end(i,:),nval_deep(i))
            call Interp(z_Inp(i,:),Q_read_end(i,:),nval_deep(i)-1,z_upp(1:nz),Q_end(i,1:nz),nz)

            ! Add surface flow for both in- and outflow
            if(nval_surface(i)>0) then
              do j=1,nval_surface(i)
                depth_surfaceFlow(i,j) = z_Inp(i,nval_deep(i)-1+j)
                call SurfaceFlow(Inp_read_start(i,nval_deep(i)-1+j),Q_start(i,:),depth_surfaceFlow(i,j),i)
                call SurfaceFlow(Inp_read_end(i,nval_deep(i)-1+j),Q_end(i,:),depth_surfaceFlow(i,j),i)
              end do
            end if

            nz_old = nz
          end if  ! idx==1

          ! If lake level changes and if there is surface inflow, adjust inflow depth to keep them at the surface
          if((.not. lake_level_old==z_upp(nz)).and. (nval_surface(i)>0)) then

            ! Recalculate Q_start from deep inflows
            call Integrate(z_Inp(i,:),Inp_read_start(i,:),Q_read_start(i,:),nval_deep(i))
            call Interp(z_Inp(i,:),Q_read_start(i,:),nval_deep(i)-1,z_upp(1:nz),Q_start(i,1:nz),nz)

            ! Recalculate Q_end from deep inflows
            call Integrate(z_Inp(i,:),Inp_read_end(i,:),Q_read_end(i,:),nval_deep(i))
            call Interp(z_Inp(i,:),Q_read_end(i,:),nval_deep(i)-1,z_upp(1:nz),Q_end(i,1:nz),nz)

            ! Add surface flow
            if(nval_surface(i)>0) then
              do j=1,nval_surface(i)
                call SurfaceFlow(Inp_read_start(i,nval_deep(i)-1+j),Q_start(i,:),depth_surfaceFlow(i,j),i)
                call SurfaceFlow(Inp_read_end(i,nval_deep(i)-1+j),Q_end(i,:),depth_surfaceFlow(i,j),i)
              end do
            end if
          end if ! end if not lake_level_old...

          ! Temporal treatment of inflow
          if ((datum<=tb_start(i)).or.(eof(i)==1)) then ! if datum before first date or end of file reached
            goto 8
          else
            do while (.not.((datum>=tb_start(i)).and.(datum<=tb_end(i)))) ! Do until datum between dates
              tb_start(i) = tb_end(i)             ! Move one step in time
              Q_start(i,1:nz) = Q_end(i,1:nz)
              read(fnum(i),*,end=7) tb_end(i),(Inp_read_end(i,j),j=0,nval(i)-1)
              ! Treat deep inflow
              call Integrate(z_Inp(i,:),Inp_read_end(i,:),Q_read_end(i,:),nval_deep(i))
              call Interp(z_Inp(i,:),Q_read_end(i,:),nval_deep(i)-1,z_upp(1:nz),Q_end(i,1:nz),nz)

              ! Add surface flow
              if(nval_surface(i)>0) then
                do j=1,nval_surface(i)
                  call SurfaceFlow(Inp_read_end(i,nval_deep(i)-1+j),Q_end(i,:),depth_surfaceFlow(i,j),i)
                end do
              end if

            end do  ! end do while
          end if

          !Linearly interpolate value at correct datum (for all depths)
          do j=1,nz
            Q_inp(i,j) = Q_start(i,j) + (datum-tb_start(i)) * (Q_end(i,j)-Q_start(i,j))/(tb_end(i)-tb_start(i))
          end do

          goto 11

7         eof(i) = 1
8         Q_inp(i,1:nz) = Q_start(i,1:nz)              ! Set to closest available value
          goto 11

9         write(6,*) 'No data found in ',trim(fname(i)),' file. Check number of depths. Values set to zero.'
          eof(i) = 1
          Q_inp(i,0:nz) = 0.0_dp
          Q_start(i,1:nz) = 0.0_dp

11        continue
        end do      ! end do i=1,4

        ! Q_vert is the integrated difference between in- and outflow (starting at the lake bottom)
        ! Q_vert is located on the upper grid, m^3/s
        Q_vert(1:nz) = Q_inp(1,1:nz)+Q_inp(2,1:nz)

        ! Set all Q to the differences (from the integrals)
        ! Q_inp is located on the central grid, element 1 remains unchanged since element 0 is 0
        do i=1,4
            do j=1,nz-1
                Q_inp(i,nz-j+1) = Q_inp(i,nz-j+1)-Q_inp(i,nz-j)
            end do
        end do

        nz_old = nz

        return
  end subroutine doLateral


  !Read first column of inflow/outflow files and compute parameters
  !Assuming inflow will plunge according to its density, entraining water from other layers
  !####################################################################
  subroutine doLateral_rho(datum,idx,T,S,rho,Q_vert,Q_inp)
  !####################################################################
        implicit none
        

        ! Global Declarations
        real(dp), intent(in) :: datum,T(0:),S(0:),rho(0:)
        real(dp), intent(inout) :: Q_vert(0:),Q_inp(1:,0:)
        integer, intent(in) :: idx

        ! Local Declarations
        integer :: i, j, k, i1, i2, nval(1:4), fnum(1:4), eof(1:4)
        character*20 :: fname(1:4)
        real(dp) :: z_Inp(1:4,1:nz_max), dummy
        real(dp) :: Inp_read_start(1:4,1:nz_max), Inp_read_end(1:4,nz_max), Inp(1:4)
        real(dp) :: Q_start(1:4,nz_max), Q_end(1:4,nz_max)
        real(dp) :: Q_read_start(1:4,nz_max), Q_read_end(1:4,nz_max)
        real(dp) :: tb_start(1:4), tb_end(1:4)
        real(dp) :: T_in, S_in, rho_in
        real(dp) :: CD_in, Ri, g_red, E
        real(dp) :: h_in(0:nz), Q_in(0:nz)

        save Inp_read_start, Inp_read_end
        save nval, eof
        save tb_start, tb_end

        fname = ['inflow           ','outflow          ','input temperature','input salinity   ']
        fnum = [41,42,43,44]
        do i=1,4
           if (idx==1) then
              if(disp_dgn==2) write(6,*) 'Starting to read '//trim(fname(i))//' file...'
              eof(i) = 0
              !Skip first three rows (except for Qout)
              read(fnum(i),*,end=9)
              read(fnum(i),*,end=9) nval(i)
              read(fnum(i),*,end=9) dummy, (z_Inp(i,j),j=1,nval(i))
              if(i/=2) nval(i)=1
              if(i==2) z_Inp(i,1:nval(i)) = z_zero + z_Inp(i,1:nval(i))
              !Read first values
              read(fnum(i),*,end=9) tb_start(i),(Inp_read_start(i,j),j=1,nval(i))
              if(i==2) call Integrate(z_Inp(i,:),Inp_read_start(i,:),Q_read_start(i,:),nval(i))
              if(i==2) call Interp(z_Inp(i,:),Q_read_start(i,:),nval(i)-1,z_upp(1:nz),Q_start(i,:),nz-1)
              read(fnum(i),*,end=7) tb_end(i),(Inp_read_end(i,j),j=1,nval(i))
              if(i==2) call Integrate(z_Inp(i,:),Inp_read_end(i,:),Q_read_end(i,:),nval(i))
              if(i==2) call Interp(z_Inp(i,:),Q_read_end(i,:),nval(i)-1,z_upp(1:nz),Q_end(i,:),nz-1)
           end if

           if ((datum<=tb_start(i)).or.(eof(i)==1)) then       ! if datum before first date or end of file reached
               goto 8
           else
               do while (.not.((datum>=tb_start(i)).and.(datum<=tb_end(i)))) ! do until datum between dates
                   tb_start(i) = tb_end(i)             ! move one step in time
                   if(i/=2) Inp_read_start(i,1) = Inp_read_end(i,1)
                   if(i==2) Q_start(i,1:nz) = Q_end(i,1:nz)
                   read(fnum(i),*,end=7) tb_end(i),(Inp_read_end(i,j),j=1,nval(i))
                   if(i==2) call Integrate(z_Inp(i,:),Inp_read_end(i,:),Q_read_end(i,:),nval(i))
                   if(i==2) call Interp(z_Inp(i,:),Q_read_end(i,:),nval(i)-1,z_upp(1:nz),Q_end(i,:),nz-1)
               end do
               !Linearly interpolate value at correct datum
               if (i/=2) then
                   Inp(i) = Inp_read_start(i,1) + (datum-tb_start(i)) * &
     &(Inp_read_end(i,1)-Inp_read_start(i,1))/(tb_end(i)-tb_start(i))
               else
                   do j=1,nz
                       Q_inp(i,j) = Q_start(i,j) + (datum-tb_start(i)) * &
     &(Q_end(i,j)-Q_start(i,j))/(tb_end(i)-tb_start(i))
                   end do
               end if
           end if
           goto 11

  7        eof(i) = 1
  8        if(i/=2) Inp(i) = Inp_read_start(i,1)            ! Set to closest available value
           if(i==2) Q_inp(i,1:nz) = Q_start(i,1:nz)     ! Set to closest available value
           goto 11

  9        write(6,*) 'No data found in ',trim(fname(i)),' file. Check number of depths. Values set to zero.'
           eof(i) = 1
           if(i/=2) Inp(i) = 0.0_dp
           if(i/=2) Inp_read_start(i,1) = 0.0_dp
           if(i==2) Q_inp(i,0:nz) = 0.0_dp
           if(i==2) Q_start(i,1:nz) = 0.0_dp

  11       continue
        end do      ! end do i=1,4

        !Set Qout to the differences (from the integrals)
        do j=1,nz-1
            Q_inp(2,nz-j+1) = Q_inp(2,nz-j+1)-Q_inp(2,nz-j)
        end do

        if (Inp(1)>1E-15) then
            
            Q_in(nz) = Inp(1) !Inflow flow rate [m3/s]
            T_in = Inp(3) !Inflow temperature [°C*m3/s]
            S_in = Inp(4) !Inflow salinity [‰*m3/s]
            rho_in = rho_0*(0.9998395_dp+T_in*(6.7914e-5_dp+T_in*(-9.0894e-6_dp+T_in*&
                     (1.0171e-7_dp+T_in*(-1.2846e-9_dp+T_in*(1.1592e-11_dp+T_in*(-5.0125e-14_dp))))))+&
                     (8.181e-4_dp+T_in*(-3.85e-6_dp+T_in*(4.96e-8_dp)))*S_in) !Inflow density [kg/m3]
            g_red = g*(rho_in-rho(nz))/rho_in !Reduced gravity [m/s2]
            if (g_red>0) then
                CD_in = CD*10 !Inflow drag coefficient
                Ri = CD_in*(1+0.21_dp*CD_in*sinatan_slope)/(sinatan_slope*tan_phi) !Richardson number
                E = 1.6_dp*CD_in**1.5_dp/Ri !Entrainment coefficient
                h_in(nz) = (2*Inp(1)**2*Ri*slope**2/g_red)**0.2_dp !Inflow thickness [m]
                do k=nz,1,-1
                    if(rho_in<=rho(k)) exit
                    h_in(k-1) = 1.2_dp*E*(z_cent(k)-z_cent(k-1))/slope + h_in(k)
                    Q_in(k-1) = Q_in(k)*(h_in(k-1)/h_in(k))**(5.0_dp/3.0_dp)
                    Q_inp(2,k) = Q_inp(2,k) - (Q_in(k-1)-Q_in(k))
                    T_in = (T_in*Q_in(k)+T(k)*(Q_in(k-1)-Q_in(k)))/Q_in(k-1)
                    S_in = (S_in*Q_in(k)+S(k)*(Q_in(k-1)-Q_in(k)))/Q_in(k-1)
                    rho_in = (rho_in*Q_in(k)+rho(k)*(Q_in(k-1)-Q_in(k)))/Q_in(k-1)
                end do
            else
                k=nz
            end if

            if (k==nz) then !surface flow
                i1=nz
                i2=nz-1
            else
                do i1=k,nz !extend upwards
                    if(z_cent(i1)>z_cent(k)+0.5_dp*h_in(k)) exit
                end do
                do i2=k,1,-1 !extend downwards
                    if(z_cent(i2)<z_cent(k)-0.5_dp*h_in(k)) exit
                end do
            end if
            i1 = i1-1
            i2 = i2-1

            Q_inp(1,:) = 0.0_dp
            Q_inp(3:4,:) = 0.0_dp
            do i=i1,i2+1,-1
                Q_inp(1,i) = Q_in(k)/(z_upp(i1)-z_upp(i2))*h(i)
                Q_inp(3,i) = T_in*Q_inp(1,i)
                Q_inp(4,i) = S_in*Q_inp(1,i)
            end do
        end if

        Q_vert(1) = Q_inp(1,1)+Q_inp(2,1)  ! Mass balance in bottom
        do i=2,nz                         ! Mass balance in all cells
            Q_vert(i) = Q_vert(i-1)+Q_inp(1,i)+Q_inp(2,i)
        end do

        return
  end subroutine doLateral_rho

  !####################################################################
  subroutine SurfaceFlow(Q_surface, Q_total, depth_surface,i)
  !####################################################################
        implicit none

        ! Global declarations
        real(dp), intent(in) :: Q_surface, depth_surface
        real(dp), intent(inout) :: Q_total(1:)
        integer, intent(in) :: i

        ! Local variables
        real(dp) :: Q_surf_integral, add_surf
        integer :: j

        Q_surf_integral = Q_surface

        do j=1,nz

          if (((Q_surf_integral>0) .and. (.not. i==2)) .or. ((Q_surf_integral<0) .and. (i==2))) then
            Q_total(nz+1-j) = Q_total(nz+1-j) + Q_surf_integral
          else
            exit
          end if

          add_surf = Q_surface*h(nz+1-j)/depth_surface
          Q_surf_integral = Q_surf_integral - add_surf
        end do
        return

  end subroutine
end module Advection