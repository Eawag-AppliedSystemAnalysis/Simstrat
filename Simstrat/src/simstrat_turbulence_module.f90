module simstrat_turbulence_module
  use simstrat_kinds
  use simstrat_model_module
  use simstrat_model_constants
  use simstrat_finite_volume_implementation
  use utilities
  implicit none

  private
  public shear_buoyancy_production, seiche_production, k_terms, eps_terms

contains
  subroutine k_terms(datum, model, boundary_cond, fluxes, sources)
    implicit none
    class(SimstratModel), intent(in) :: model
    real(RK), dimension(:), intent(inout) :: boundary_cond, fluxes, sources
    real(RK), intent(in) :: datum

    ! Local variables
    real(RK), dimension(model%discretization%nz_tq) :: pminus, pplus
    real(RK) :: Prod, Buoy, Diss
    integer :: i

    associate(nz_tq => model%discretization%nz_tq,&
              dt=>model%dt,&
              k=>model%k,&
              P_Seiche=>model%P_Seiche, P=>model%P, B=>model%B, eps=>model%eps)

      do i=1,nz_tq
        Prod = P(i)+P_Seiche(i)                   ! Add seiche energy
        Buoy = B(i)
        Diss = eps(i)
        !if (Prod+Buoy>0) then
        !    pplus(i)=Prod+Buoy
        !    pminus(i)=Diss
        !else
        !    pplus(i)=Prod
        !    pminus(i)=Diss-Buoy
        !end if
        pminus(i) = Diss
        pplus(i) = Prod+Buoy
      end do
      boundary_cond(1:nz_tq) = pminus(1:nz_tq)/k(1:nz_tq)
      sources(1:nz_tq) = pplus(1:nz_tq)
    end associate
  end subroutine

  subroutine eps_terms(datum, model, boundary_cond, fluxes, sources)
    implicit none
    class(SimstratModel), intent(in) :: model
    real(RK), dimension(:), intent(inout) :: boundary_cond, fluxes, sources
    real(RK), intent(in) :: datum

    ! Local variables
    real(RK), dimension(model%discretization%nz_tq) :: pminus, pplus
    real(RK) :: Prod, Buoy, Diss, cee3

    integer :: i

    select type(discretization=>model%discretization)
      class is (StaggeredFiniteVolumeDiscretization)

        associate(nz_tq => discretization%nz_tq,&
                  k=>model%k,&
                  P_Seiche=>model%P_Seiche, P=>model%P, B=>model%B, eps=>model%eps, ko=>model%ko, ModFlux=>model%ModFlux, Mod=>model%Mod, h_faces=>discretization%h_faces, a_faces=>discretization%a_faces, cde=>model%cde, num_eps=>model%num_eps)
          if(ModFlux == 1 .and. Mod == 1) then
            sources(1:nz_tq) = num_eps(1:nz_tq)/sig_e*(cde*((ko(1:nz_tq))**1.5_RK)/(kappa*(z0+0.25_RK*(h_faces(1:nz_tq)+h_faces(1:nz_tq))))**2)*discretization%form_eps(1:nz_tq)
            sources(1) = sources(1) + num_eps(1)*(cde*((ko(1))**1.5_RK)/(kappa*(K_s+0.5_RK*h_faces(1)))**2)*(a_faces(1)+a_faces(2))/(a_faces(1)*(h_faces(1)+h_faces(2)))
            sources(nz_tq) = sources(nz_tq) + num_eps(nz_tq)*(cde*((ko(nz_tq))**1.5_RK)/(kappa*(z0+0.5_RK*h_faces(nz_tq)))**2)*(a_faces(nz_tq)+a_faces(nz_tq-1))/(a_faces(nz_tq)*(h_faces(nz_tq)+h_faces(nz_tq-1)))
          end if

          do i=1,nz_tq
            if (B(i)>0) then
                cee3=1.0_RK
            else
                cee3=ce3
            end if
            Prod=ce1*eps(i)/ko(i)*(P(i)+P_Seiche(i))        ! New code plus seiche
            Buoy=cee3*eps(i)/ko(i)*B(i)
            Diss=ce2*eps(i)*eps(i)/ko(i)
            if (Prod+Buoy>0) then
                pplus(i)=Prod+Buoy
                pminus(i)=Diss
            else
                pplus(i)=Prod
                pminus(i)=Diss-Buoy
            end if
          end do
          !write(*,*) pminus(nz_tq-1)
          !write(*,*) sources(nz_tq-1)
          boundary_cond(1:nz_tq) = pminus(1:nz_tq)/eps(1:nz_tq)
          sources(1:nz_tq) = pplus(1:nz_tq)
          !boundary_cond(2:nz_tq-1) = pminus(2:nz_tq-1)/eps(2:nz_tq-1)
          !sources(2:nz_tq-1) = sources(2:nz_tq-1) - pminus(2:nz_tq-1) + pplus(2:nz_tq-1)
          !sources(1:nz_tq) = sources(1:nz_tq) -pminus(1:nz_tq)+pplus(1:nz_tq)
        end associate
      class default
        call error('eps_terms not implemented for this discretization scheme')
    end select
  end subroutine

  !####################################################################
  pure subroutine shear_buoyancy_production(U,V,NN,num,nuh,h_faces,P,B,nz_mfq)
  !####################################################################
    implicit none

    ! Global variables
    integer, intent(in) :: nz_mfq
    real(RK), dimension(nz_mfq), intent(in) :: U,V
    real(RK), dimension(nz_mfq), intent(in) :: h_faces
    real(RK), dimension(nz_mfq+1), intent(in) :: NN, num, nuh
    real(RK), dimension(:), intent(inout) :: P,B

    real(RK), dimension(nz_mfq+1) :: dUdz, dVdz

    !dUdz and dVdz at volume faces
    dUdz(2:nz_mfq+1)=(U(2:nz_mfq)-U(1:nz_mfq-1))/h_faces(1:nz_mfq)
    dUdz(1) = dUdz(2) !could also be zero
    dVdz(2:nz_mfq+1)=(V(2:nz_mfq)-V(1:nz_mfq-1))/h_faces(1:nz_mfq)
    dVdz(1) = dVdz(2) !could also be zero

    ! Equation 5 (left) of Goudsmit, 2002
    ! On volume faces
    P(1:nz_mfq+1) = (dUdz**2+dVdZ**2)*num

    ! Equation 5 (right) of Goudsmit, 2002
    ! On volume faces
    B(1:nz_mfq+1) = -nuh(1:nz_mfq+1)*NN(1:nz_mfq+1)

    return
  end subroutine

  !####################################################################
  subroutine seiche_production(E_Seiche,P_Seiche,a_seiche,u10,v10,CD,C10,WindFilt,Wf,NN,gamma,q_NN,ModSNorm,h_faces,a_faces,dAdz_faces,dt,nz_mfq)
  !####################################################################

        implicit none

        ! Global variables
        integer, intent(in) :: nz_mfq
        real(RK), dimension(nz_mfq+1), intent(in) :: NN, a_faces, dAdz_faces, h_faces
        real(RK), intent(in) :: a_seiche, u10, v10, CD, C10, Wf, gamma, q_NN,dt
        logical, intent(in) :: WindFilt
        integer, intent(in) :: ModSNorm
        real(RK), dimension(:), intent(inout) :: P_Seiche
        real(RK), intent(inout) :: E_Seiche

        ! Local variables
        real(RK) :: W10, PS, PW, f_norm, minNN
        real(RK), dimension(nz_mfq+1) :: Distrib
        integer :: i

        Distrib = 0.0_RK
        minNN = 0.0_RK
        if (a_seiche/=0.0_RK) then      ! a_seiche is defined in the par-file
          do i=1,nz_mfq+1
            Distrib(i) = max(NN(i)**q_NN,minNN) / a_faces(i)*dAdz_faces(i)
          end do
            !Seiche normalization factor
            f_norm = 0.0_RK
            if (ModSNorm==1) then !max NN
                do i=1,nz_mfq+1
                    if(NN(i)>f_norm) f_norm = NN(i)
                end do
                f_norm = (f_norm**q_NN)*a_faces(nz_mfq+1)*rho_0
            else if (ModSNorm==2) then !integral
                do i=1,nz_mfq+1
                    f_norm = f_norm+Distrib(i)*a_faces(i)*h_faces(i)
                end do
                f_norm = f_norm*rho_0
            end if

            if (f_norm==0.0_RK) then
              Distrib(1:nz_mfq+1)=1.0_RK/h_faces(1:nz_mfq+1)
              f_norm= a_faces(nz_mfq+1)*rho_0
            end if

            !Seiche energy per surface converted to total seiche energy (equation 10 in Goudsmit, 2002)
            if (WindFilt) then !use filtered wind (AG 2014)
                PW = a_seiche*a_faces(nz_mfq+1)*rho_air*C10*Wf**3
            else !use real wind
                W10 = sqrt(u10**2+v10**2)
                PW = a_seiche*a_faces(nz_mfq+1)*rho_air*C10*W10**3
            end if
            PS = E_Seiche**(1.5_RK)*gamma

            E_Seiche = E_Seiche + (PW - PS)*dt

            !Limit so that E_Seiche does not become negative
            if (E_Seiche<0.) then
                PS = (PS*dt+E_Seiche)/dt
                E_Seiche = 0.0_RK
            end if

            ! Equation 24 in Goudsmit, 2002
            P_Seiche(1:nz_mfq+1) = 1.0_RK/f_norm*Distrib*PS*(1.0_RK-10*sqrt(CD))
        else          !if alpha==0
            P_seiche = 0.0_RK
        end if

        return
    end subroutine
end module
