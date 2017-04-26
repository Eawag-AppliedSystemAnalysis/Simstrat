module simstrat_stability_module
  use simstrat_kinds
  use simstrat_model_module
  use simstrat_model_constants
  implicit none

  private
  public water_column_stability

contains
  subroutine water_column_stability(stab, salctr, delsal, h_centres, k, eps, T, S, cde, rho, NN, cmue1, cmue2, nz_mfq, nz_tq)
    implicit none
    ! Global variables
    real(RK), dimension(:), intent(in) :: h_centres
    real(RK), dimension(:), intent(in) :: T,S
    real(RK), dimension(:), intent(in) :: k,eps
    real(RK), intent(in) :: cde
    real(RK), dimension(:), intent(inout) :: cmue1,cmue2,NN,rho
    integer, intent(in) :: stab, salctr, delsal, nz_mfq, nz_tq

    ! Local variables
    real(RK), dimension(nz_mfq+1) :: beta
    integer :: i

    !NN on volume centres
    call buoyancy(T,S,rho,NN,h_centres,salctr,delsal,nz_mfq)

    if (stab==1) then
        call cmue_cn(cmue1,cmue2)
    else if (stab==2) then
        beta = NN(1:nz_mfq+1)*(k(1:nz_mfq+1)/eps(1:nz_mfq+1))**2
        call cmue_qe(beta,cde,cmue1,cmue2)
    end if

  end subroutine

  !####################################################################
  subroutine buoyancy(T,S,rho,NN,h_centres,salctr,delsal,nz_mfq)
  !####################################################################
      implicit none

      ! Global variables
      real(RK), dimension(:), intent(in) :: T,S,h_centres
      real(RK), dimension(:), intent(inout) :: NN,rho
      integer, intent(in) :: salctr, delsal, nz_mfq

      ! Local variables
      real(RK), dimension(nz_mfq+1) :: T_faces, rho_faces
      real(RK), dimension(nz_mfq) :: a
      real(RK), dimension(nz_mfq) :: rho0t,rho0st
      real(RK), dimension(nz_mfq) :: dTdz
      integer :: i

      ! Approxiamte T at the volume faces by shifting index
      T_faces(1) = T(1)
      T_faces(2:nz_mfq+1) = T(1:nz_mfq)

      ! T gradient at volume centres
      dTdz(1:nz_mfq) = (T_faces(1:nz_mfq)-T_faces(2:nz_mfq+1))/h_centres(1:nz_mfq)
      !dTdz(1) = dTdz(2) !continue gradient at the bottom
      

      if (salctr==0) then                ! salinity is zero everywhere
        a = -68.0_RK+T(1:nz_mfq)*(18.2091_RK+T(1:nz_mfq)*(-0.30866_RK+T(1:nz_mfq)*&
                   (5.3445e-3_RK+T(1:nz_mfq)*(-6.0721e-5_RK+T(1:nz_mfq)*(3.1441e-7_RK)))))
        a = 1.0e-6_RK*a
        NN(2:nz_mfq+1) = g*a*dTdz+(T(1:nz_mfq)+273.15_RK)/cp
      else
          if (delsal==0) then        ! salinity gradient is zero everywhere
            a = -68.0_RK+T(1:nz_mfq)*(18.2091_RK+T(1:nz_mfq)*(-0.30866_RK+T(1:nz_mfq)*&
                  (5.3445e-3_RK+T(1:nz_mfq)*(-6.0721e-5_RK+T(1:nz_mfq)*(3.1441e-7_RK)))))
            a = a + (4.599_RK+T(1:nz_mfq)*(-0.1999_RK+T(1:nz_mfq)*(2.79e-3_RK)))*S(i)
            a = 1.0e-6_RK*a
            NN(2:nz_mfq+1) = g*a*dTdz+(T(1:nz_mfq)+273.15_RK)/cp
          else
            !on volume centres
            rho0t= 0.9998395_RK+T(1:nz_mfq)*(6.7914e-5_RK+T(1:nz_mfq)*(-9.0894e-6_RK+T(1:nz_mfq)*&
                (1.0171e-7_RK+T(1:nz_mfq)*(-1.2846e-9_RK+T(1:nz_mfq)*(1.1592e-11_RK+T(1:nz_mfq)*(-5.0125e-14_RK))))))
            rho0st= (8.181e-4_RK+T(1:nz_mfq)*(-3.85e-6_RK+T(1:nz_mfq)*(4.96e-8_RK)))*S(1:nz_mfq)
            rho(1:nz_mfq) = rho_0*(rho0t+rho0st)
            !interpolate rho on volume faces by shifting
            rho_faces(1) = rho(1)
            rho_faces(2:nz_mfq+1) = rho(1:nz_mfq)

            !interpolate on volume faces
            NN(2:nz_mfq+1) = g*(rho_faces(1:nz_mfq)-rho_faces(2:nz_mfq+1))/h_centres(1:nz_mfq)
          end if
      end if

      !NN is calculated at the volume centres and approximated at the volume faces by shifting one index up
      NN(1) = NN(2)

      return
  end subroutine

  ! Calculation of cmue
  !####################################################################
  subroutine cmue_cn(cmue1,cmue2)
  !####################################################################
    implicit none

    ! Global variables
    real(RK), dimension(:), intent(inout) :: cmue1, cmue2

    !Standard version of k-eps model
    cmue1 = cmue
    cmue2 = cmue/Prndtl
    !Burchard Version
    !cmue1=cde*cm0
    !cmue2=cde*cm0/Prndtl

    return
  end subroutine

  ! Calculation of cmue
  !####################################################################
  subroutine cmue_qe(beta,cde,cmue1,cmue2)
  !####################################################################
    implicit none

    ! Global variables
    real(RK), dimension(:), intent(in) :: beta
    real(RK), intent(in) :: cde
    real(RK), dimension(:), intent(inout) :: cmue1,cmue2

    ! Local variables
    integer :: n
    real(RK), dimension(size(beta)) :: gh, sm, sh

    integer :: i

    n = size(beta)
    

    gh = -cde**2*0.5_RK*beta
    where(gh > 0.02_RK) gh = gh-(gh-0.02_RK)**2/(gh+0.0233_RK-2*0.02_RK)
    where(gh < -0.28_RK) gh = -0.28_RK
    
    sm = 1.0_RK-3*c1-6*a1/b1-3*a2*gh*((b2-3*a2)*(1.0_RK-6*a1/b1)-3*c1*(b2+6*a1))
    sm = a1*sm/((1.0_RK-3*a2*gh*(6*a1+b2))*(1.0_RK-9*a1*a2*gh))
    sh = a2*(1.0_RK-6*a1/b1)/(1.0_RK-3*a2*gh*(6*a1+b2))

    cmue1(1:n) = sqrt(2.0_RK)*cde*sm
    cmue2(1:n) = sqrt(2.0_RK)*cde*sh

    return
  end subroutine
end module
