module Transport
    use SimstratModel
    use utilities
    implicit none

    private
    public doTransportEquation

contains
    !####################################################################
    subroutine doTransportEquation(C,dC,nuh)
    !####################################################################
        
        implicit none

        ! Global variables
        real(dp), intent(inout) :: C(0:)
        real(dp), intent(in) :: dC(0:),nuh(0:)

        ! Local variables
        real(dp) :: au(0:nz),bu(0:nz),cu(0:nz),du(0:nz)
        !real(dp) SalRel
        !SalRel = 172800.

        !Build the three diagonals
        cu(1:nz-1) = dt*nuh(1:nz-1)*AreaFactor_2(1:nz-1)
        cu(nz) = 0.0_dp
        au(1) = 0.0_dp
        au(2:nz) = dt*nuh(1:nz-1)*AreaFactor_1(2:nz)
        bu(1:nz) = 1.0_dp-au(1:nz)-cu(1:nz)
        du(1:nz) = C(1:nz) + dC(1:nz)*dt !AG 2014: added dC*dt term for source/sink

        !Solve
        call Tridiagonal(1,nz,au,bu,cu,du,C)
        C(0) = C(1)

        return
    end subroutine doTransportEquation

end module Transport
