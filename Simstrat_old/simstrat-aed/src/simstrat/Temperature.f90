module Temperature
    use SimstratModel
    use utilities
    implicit none

    private
    public doTemperature

contains
    !####################################################################
    subroutine doTemperature(nuh,rad0,rad,T,heat,SST,absorb)
    !####################################################################

        implicit none

        ! Global variables
        real(dp), intent(inout) :: T(0:), rad(0:)
        real(dp), intent(in) :: nuh(0:), absorb(0:), rad0, heat, SST

        ! Local variables
        real(dp) :: au(0:nz),bu(0:nz),cu(0:nz),du(0:nz)
        integer :: i

        ! Calculation
        !Radiation reaching each layer
        rad(nz) = rad0/rho_0/cp ![°C*m/s]
        do i=nz-1,0,-1
            rad(i) = rad(i+1)*exp(-h(i)*absorb(nz-i)) !Attenuated by absorption
        end do

        !Build the three diagonals
        au(1) = 0.0_dp
        au(2:nz) = dt*nuh(1:nz-1)*AreaFactor_1(2:nz)
        cu(1:nz-1) = dt*nuh(1:nz-1)*AreaFactor_2(1:nz-1)
        cu(nz) = 0.0_dp
        bu(1:nz) = 1-au(1:nz)-cu(1:nz)
        du(1:nz) = T(1:nz)+(rad(1:nz)-rad(0:nz-1))/h(1:nz)*dt
        du(nz) = du(nz) + heat/rho_0/cp*dt/h(nz)

        if (NBC==1) then
            bu(nz) = 1.0_dp
            au(nz) = 0.0_dp
            cu(nz) = 0.0_dp
            du(nz) = SST
            write(*,*) "NBC=1"
        end if

        !Add geothermal heat flux
        if(fgeo/=0) du(1:nz)=du(1:nz)+fgeo_add(1:nz)*dt

        !Solve
        call Tridiagonal(1,nz,au,bu,cu,du,T)
        T(0) = T(1)     ! set value to boundary value at z_cent(1)

        return
    end subroutine doTemperature

end module Temperature
