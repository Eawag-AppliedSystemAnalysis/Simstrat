module simstrat_wind_shear_module
  use simstrat_kinds
  use simstrat_model_module
  use simstrat_model_constants
  use simstrat_finite_volume_implementation
  use utilities
  implicit none

  private
  public coriolis, u_terms, v_terms

contains
  !####################################################################
  pure subroutine coriolis(U,V,Cori,dt,nz)
  !####################################################################
    implicit none

    ! Global variables
    real(RK), dimension(:), intent(inout) :: U, V
    real(RK), intent(in) :: Cori,dt
    integer, intent(in) :: nz

    ! Local variables
    real(RK), dimension(nz) :: U_old

    ! Calculation
    U_old(1:nz) = U(1:nz)
    U(1:nz) =  U(1:nz)*cos(Cori*dt) + V(1:nz)*sin(Cori*dt)
    V(1:nz) =-U_old(1:nz)*sin(Cori*dt) + V(1:nz)*cos(Cori*dt)

    return
  end subroutine
  
  subroutine u_terms(datum, model, boundary_cond, fluxes, sources)
    implicit none
    class(SimstratModel), intent(in) :: model
    real(RK), dimension(:), intent(inout) :: boundary_cond, fluxes, sources
    real(RK), intent(in) :: datum

    real(RK) :: intU, length

    select type(discretization=>model%discretization)
      class is (StaggeredFiniteVolumeDiscretization)
        associate(nz_mfq => model%discretization%nz_mfq,&
                  U=>model%u, V=>model%v,&
                  drag=>model%drag,&
                  tx=>model%tx,&
                  Pgrad=>model%Pgrad,&
                  depth=>model%depth, h_centres=>discretization%h_centres, a_centres=>discretization%a_centres, dAdz=>discretization%dAdz)

          if (Pgrad == 1) then
            intU = sum(U(1:nz_mfq))
            length = sqrt(a_centres(nz_mfq))
          end if

          boundary_cond(1) = drag*sqrt(U(1)**2+V(1)**2)/h_centres(1)
          if (Pgrad==1) then !Svensson 1978
              sources(2:nz_mfq-1) = -pi**2*rho_0*g*intU/nz_mfq*depth/length**2
          elseif (Pgrad==2) then !???
              sources(2:nz_mfq-1) = -drag*U(2:nz_mfq-1)*sqrt(U(2:nz_mfq-1)**2+V(2:nz_mfq-1)**2)*dAdz(2:nz_mfq-1)/a_centres(2:nz_mfq-1)
          end if
          sources(nz_mfq) = tx/h_centres(nz_mfq)
        end associate

      class default
        call error('u_terms not implemented for this discretization scheme')
    end select

  end subroutine

  subroutine v_terms(datum, model, boundary_cond, fluxes, sources)
    implicit none
    class(SimstratModel), intent(in) :: model
    real(RK), dimension(:), intent(inout) :: boundary_cond, fluxes, sources
    real(RK), intent(in) :: datum

    real(RK) :: intV, length

    select type(discretization=>model%discretization)
      class is (StaggeredFiniteVolumeDiscretization)
        associate(nz_mfq => model%discretization%nz_mfq,&
                  U=>model%u, V=>model%v,&
                  drag=>model%drag,&
                  ty=>model%ty,&
                  Pgrad=>model%Pgrad,&
                  depth=>model%depth, h_centres=>discretization%h_centres, a_centres=>discretization%a_centres, dAdz=>discretization%dAdz)

          if (Pgrad == 1) then
            intV = sum(V(1:nz_mfq))
            length = sqrt(a_centres(nz_mfq))
          end if

          boundary_cond(1) = drag*sqrt(U(1)**2+V(1)**2)/h_centres(1)
          if (Pgrad==1) then !Svensson 1978
              sources(2:nz_mfq-1) = -pi**2*rho_0*g*intV/nz_mfq*depth/length**2
          elseif (Pgrad==2) then !???
              sources(2:nz_mfq-1) = -drag*V(2:nz_mfq-1)*sqrt(U(2:nz_mfq-1)**2+V(2:nz_mfq-1)**2)*dAdz(2:nz_mfq-1)/a_centres(2:nz_mfq-1)
          end if
          sources(nz_mfq) = ty/h_centres(nz_mfq)
        end associate

      class default
        call error('u_terms not implemented for this discretization scheme')
    end select

  end subroutine

end module
