module UVEquation
    use SimstratModel
    use utilities
    implicit none

    private
    public doUVEquation

contains

    !####################################################################
    subroutine doUVEquation(U,V,num,drag,tx,ty)
    !####################################################################

        implicit none

        ! Global variables
        real(dp), intent(inout) :: U(0:),V(0:)
        real(dp), intent(in) :: num(0:), drag,tx,ty

        ! Local variables
        real(dp) :: au(0:nz),bu(0:nz),cu(0:nz),du(0:nz)
        real(dp) :: intU, intV, length

        if (Pgrad==1) then
            intU = sum(U(1:nz))
            intV = sum(V(1:nz))
            length = sqrt(Az(nz))
        end if

        ! Build the diagonals
        cu(1:nz-1) = dt*num(1:nz-1)*AreaFactor_2(1:nz-1)
        cu(nz) = 0.0_dp
        au(1) = 0.0_dp
        au(2:nz) = dt*num(1:nz-1)*AreaFactor_1(2:nz)
        bu(1:nz) = 1-au(1:nz)-cu(1:nz)
        bu(1) = bu(1) + drag*sqrt(U(1)**2+V(1)**2)/h(1)*dt

        ! Calculation of U-equation
        du(1) = U(1)
        if (Pgrad==1) then !Svensson 1978
            !du(2:nz-1) = U(2:nz-1) - 10*pi**2*intU/nz*rho_0*g*depth/length**2*dt/86400
            du(2:nz-1) = U(2:nz-1) - pi**2*rho_0*g*intU/nz*depth/length**2
        elseif (Pgrad==2) then !???
            du(2:nz-1) = U(2:nz-1) - drag*U(2:nz-1)*sqrt(U(2:nz-1)**2+V(2:nz-1)**2)*dAdz(2:nz-1)/Az(2:nz-1)
        else
            du(2:nz-1) = U(2:nz-1)
        end if
        du(nz) = U(nz) + tx*dt/h(nz)
        ! Solve
        call Tridiagonal(1,nz,au,bu,cu,du,U)
        ! Calculation of V-equation
        du(1) = V(1)
        if (Pgrad==1) then !Svensson 1978
            !du(2:nz-1) = V(2:nz-1) - 10*pi**2*intV/nz*rho_0*g*depth/length**2*dt/86400
            du(2:nz-1) = V(2:nz-1) - pi**2*rho_0*g*intV/nz*depth/length**2*dt
        elseif (Pgrad==2) then !???
            du(2:nz-1) = V(2:nz-1) - drag*V(2:nz-1)*sqrt(U(2:nz-1)**2+V(2:nz-1)**2)*dAdz(2:nz-1)/Az(2:nz-1)*dt
        else
            du(2:nz-1) = V(2:nz-1)
        end if
        du(nz) = V(nz) + ty*dt/h(nz)
        ! Solve
        call Tridiagonal(1,nz,au,bu,cu,du,V)

        return
    end subroutine doUVEquation
end module UVEquation
