!module that collects the subroutins used by the simulation chain
!some of the source terms are distributed to other modules
module simstrat_lake_physics_module
  use simstrat_kinds
  use simstrat_model_module
  use simstrat_model_constants
  use simstrat_finite_volume_implementation
  use simstrat_stability_module
  use simstrat_turbulence_module
  use simstrat_wind_shear_module
  use utilities

  implicit none

  private
  !source terms
  public temperature_terms, salinity_terms, u_terms, v_terms, k_terms, eps_terms
  !additionsl subroutines
  public coriolis, water_column_stability
  public shear_buoyancy_production, seiche_production

contains
  !###################################################################
  ! PDE terms
  !###################################################################

  subroutine temperature_terms(datum, model, boundary_cond, fluxes, sources)
    implicit none
    class(SimstratModel), intent(in) :: model
    real(RK), dimension(:), intent(inout) :: boundary_cond, fluxes, sources
    real(RK), intent(in) :: datum

    integer :: i

    real(RK), dimension(model%discretization%nz_mfq) :: rad, dRaddz
    rad = 0.0_RK
    dRaddz = 0.0_RK

    select type(discretization=>model%discretization)
      class is (StaggeredFiniteVolumeDiscretization)

        associate(nz_mfq=>model%discretization%nz_mfq, fgeo=>model%fgeo, h_centres=>discretization%h_centres, h_faces=>discretization%h_faces, heat=>model%heat, SST=>model%SST, ForcingType=>model%ForcingType, rad0=>model%rad0, ga1=>model%ga1)
          ! Calculation
          !Radiation reaching each layer
          rad(nz_mfq) = rad0/rho_0/cp ![°C*m/s]
          do i=nz_mfq-1,1,-1
              rad(i) = rad(i+1)*exp(-h_centres(i)*ga1(i)) !Attenuated by absorption
          end do

          dRaddz(2:nz_mfq) = (rad(2:nz_mfq)-rad(1:nz_mfq-1))/h_faces(2:nz_mfq)
          dRaddz(1) = dRaddz(2)

          !Build the three diagonals
          sources(1:nz_mfq) = dRaddz(1:nz_mfq)
          sources(nz_mfq) = sources(nz_mfq) + heat/rho_0/cp/h_centres(nz_mfq)

          !if (ForcingType==1) then
          !    md(nz_mfq) = 1.
          !    ld(nz_mfq) = 0.
          !    ud(nz_mfq) = 0.
          !    sources(nz_mfq) = SST
          !end if

          !Add geothermal heat flux
          !if(fgeo/=0) du(1:xl)=du(1:xl)+fgeo_add(1:xl)*dt
        end associate
    end select


  end subroutine

  subroutine salinity_terms(datum, model, boundary_cond, fluxes, sources)
    implicit none
    class(SimstratModel), intent(in) :: model
    real(RK), dimension(:), intent(inout) :: boundary_cond, fluxes, sources
    real(RK), intent(in) :: datum

  end subroutine
end module
