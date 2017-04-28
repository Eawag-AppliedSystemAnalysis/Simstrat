module strat_grid
  use strat_kinds
  use utilites
  implicit none
  private

type, public :: GridConfig
  real(RK) :: nz_grid_max
  real(RK), dimension(:), allocatable :: grid_read ! Grid definition
  real(RK), dimension(:), allocatable :: A_read ! Area definition
  real(RK), dimension(:), allocatable :: z_A_read ! Height of area definitions
  real(RK) :: depth
  logical :: equidistant_grid
end type

type, public :: StaggeredGrid
  real(RK), dimension(:), allocatable :: h
  real(RK), dimension(:), allocatable :: z_cent
  real(RK), dimension(:), allocatable :: z_upp
  real(RK), dimension(:), allocatable :: Az
  real(RK), dimension(:), allocatable :: dAz
  real(RK), dimension(:), allocatable :: meanint
  real(RK) :: volume

  real(RK), dimension(:), allocatable :: AreaFactor_1
  real(RK), dimension(:), allocatable :: AreaFactor_2
  real(RK), dimension(:), allocatable :: AreaFactor_k1
  real(RK), dimension(:), allocatable :: AreaFactor_k2
  real(RK), dimension(:), allocatable :: AreaFactor_eps
  real(RK) :: nz_grid
  real(RK): z_zero

  type(GridConfig), allocatable, pointer :: config
contains
    procedure, pass :: init => grid_init
    procedure, pass :: init_grid_points => grid_init_grid_points
    procedure, pass :: init_z_axes => grid_init_z_axes
    procedure, pass :: init_morphology => grid_init_morphology
    procedure, pass :: init_areas => grid_init_areas
    procedure, pass :: update_area_factors => grid_update_area_factors
end type

contains
  subroutine grid_init(self, config)
    implicit none
    class(StaggeredGrid), intent(inout) :: self
    class(GridConfig), intent(in) :: config

    ! Assign config
    self%config => config

    ! Use read config to determine grid size
    call self%init_morphology()
    call self%init_grid_points()


    ! Allocate arrays according to size
    call self%memory_init()

    ! Call init functions
    call self%init_z_axes()

    call self%init_z_axes()
    call self%init_areas()
    call self%calc_area_factors()
  end subroutine grid_init

  subroutine grid_memory_init(self)
    implicit none
    class(StaggeredGrid), intent(inout) :: self

    associate(nz_grid => self%nz_grid)

    ! Allocate memory for initialization using nz_grid
    allocate(self%h(0:nz_grid))         ! Depth axis with center of boxes
    allocate(self%z_cent(0:nz_grid))         ! Depth axis with center of boxes
    allocate(self%z_upp(0:nz_grid))          ! Depth axis with upper border of boxes
    allocate(self%Az(0:nz_grid))             ! Az is defined on z_upp
    allocate(self%dAdz(1:nz_grid))           ! dAz is the difference between Az and thus on z_cent

    ! Area factors used in calculations
    allocate(self%AreaFactor_1(1:nz_grid))
    allocate(self%AreaFactor_2(1:nz_grid))
    allocate(self%AreaFactor_k1(1:nz_grid))
    allocate(self%AreaFactor_k2(1:nz_grid))
    allocate(self%AreaFactor_eps(1:nz_grid))

    allocate(self%meanint(0:nz_grid))        ! Inverse ratio of mean height of two adjacent boxes

  end subroutine grid_memory_init

  ! Calculates h and definite number of nz_grid
  ! Depending on configuration
  subroutine grid_init_grid_points(self)
    implicit none
    integer i
    class(StaggeredGrid), intent(inout) :: self

    if (self%config%equidistant_grid) then     ! Constant spacing
      nz_grid=int(self%config%grid_rea(0))
    else ! Variable spacing according to config

      ! Grid size given through number of grid_rea
      self%nz_grid = size(self%config%grid_rea)

      !Include top value if not included
      if (self%config%grid_rea(0)/=0.) then
        self%nz_grid=self%nz_grid+1
        do i=self%nz_grid,1,-1
            self%config%grid_rea(i)=self%config%grid_read(i-1)
        end do
        self%config%grid_rea(0)=0.0_dp
      end if

      !If maxdepth grid larger than morphology
      if (self%config%grid_rea(self%nz_grid)>self%config%depth) then
          do while ((self%config%grid_rea(self%nz_grid)>self%config%depth).and.(self%nz_grid>0.))
              self%nz_grid=self%nz_grid-1
          end do
      end if

      !Include bottom value if not included
      if (self%config%grid_rea(self%nz_grid)< self%config%depth) then
          self%nz_grid=self%nz_grid+1
          self%config%grid_rea(self%nz_grid)=self%config%depth
      end if
    end if

    !Construct H
    allocate(self%h(0:nz_grid))
    h(0) = 0                ! Note that h(0) has no physical meaning but helps with some calculations

    if (self%config%equidistant_grid) then
      ! Equidistant grid
      self%h(1:self%nz_grid) = self%config%depth/self%nz_grid
    else
       ! Set up h according to configuraiton
      do i=1,self%nz_grid
         self%h(1+self%nz_grid-i)=self%config%grid_rea(i)-self%config%grid_rea(i-1)
      end do
    end if
  end subroutine grid_init_grid_points

  ! Initializes z_cent and z_upp
  subroutine grid_init_z_axes(self)
    implicit none
    class(StaggeredGrid), intent(inout) :: self
    integer :: i

     !Compute position of layer center and top
     self%z_cent(0)=0.0_dp
     self%z_upp(0)=0.0_dp
     do i=1,self%nz_grid
         self%z_cent(i)=self%z_cent(i-1)+0.5_dp*(self%h(i-1)+self%h(i))
         self%z_upp(i)=self%z_upp(i-1)+self%h(i)
     end do
     do i=1,self%nz_grid
         self%z_cent(i)=nint(1e6_dp*self%z_cent(i))/1e6_dp
         self%z_upp(i)=nint(1e6_dp*self%z_upp(i))/1e6_dp
     end do
     ! needed?
     !lake_level_old = z_upp(nz)
  end subroutine grid_init_z_axes


  subroutine grid_init_morphology(self)
    implicit none
    class(StaggeredGrid), intent(inout) :: self
    real(RK), dimension(:), allocatable z_tmp, A_tmp
    integer :: num_A_read
    associate(cfg => self%config)

    z_tmp = cfg%z_read
    A_tmp = cfg%A_read
    num_read = size(self%config%A_read)

    do i=0,num_read-1                          ! Reverse order of values
         cfg%z_read(i) = -z_tmp(num_read-i-1)
         cfg%A_read(i) = A_tmp(num_read-i-1)
    end do

    self%depth = cfg%z_read(0) - cfg%z_read(num_read-1)              ! depth = max - min depth
    self%z_zero = cfg%z_read(0)
    ! lake_level_old = z_zero %needed?
    cfg%z_read(0:num_read-1) = z_zero - cfg%z_read(0:num_read-1)      ! z-coordinate is positive upwards, zero point is at reservoir bottom
  end subroutine grid_init_morphology

  subroutine grid_init_areas(self)
    implicit none
    class(StaggeredGrid), intent(inout) :: self

    associate(nz_grid  => self%nz_grid)
    associate(dAdz  => self%dAdz)
    associate(z_upp  => self%z_upp)


    ! Interpolate area (A) at all depths (z_upp)
    call Interp(z_read, A_read, num_read-1, z_upp, Az, nz_grid)

   ! Compute area derivative (= projected sediment area over layer thickness)
    dAdz(1:nz_grid) = (Az(1:nz_grid)-dAz(0:nz_grid-1))/(z_upp(1:nz_grid)-z_upp(0:nz_grid-1))
  end subroutine grid_init_areas


  subroutine grid_update_area_factors(self)
    implicit none
    class(StaggeredGrid), intent(inout) :: self
    implicit none

    integer :: i
    associate(Az => self%Az)
    associate(h => self%h)
    associate(nz => nz_grid)

    self%AreaFactor_1(1:nz) = -4*Az(0:nz-1)/(h(1:nz)+h(0:nz-1))/h(1:nz)/(Az(1:nz)+Az(0:nz-1))
    self%AreaFactor_2(1:nz) = -4*Az(1:nz)/(h(1:nz)+h(2:nz+1))/h(1:nz)/(Az(1:nz)+Az(0:nz-1))
    self%AreaFactor_k1(1:nz-1) = -(Az(1:nz-1)+Az(2:nz))/(h(1:nz-1)+h(2:nz))/h(2:nz)/Az(1:nz-1)
    self%AreaFactor_k2(1:nz-1) = -(Az(1:nz-1)+Az(0:nz-2))/(h(1:nz-1)+h(2:nz))/h(1:nz-1)/Az(1:nz-1)
    self%AreaFactor_eps(1:nz-1) = 0.5_dp*((Az(1:nz-1)-Az(0:nz-2))/h(1:nz-1)+(Az(2:nz)-Az(1:nz-1))/h(2:nz))/Az(1:nz-1)

    self%meanint(0:nz-1) = 2.0_dp/(h(0:nz-1)+h(1:nz))

    self%volume=0
    do i=0,nz-1
       self%volume = self%volume + 0.5_dp*h(i+1)*(Az(i)+Az(i+1))
    end do
  end subroutine grid_update_area_factors




end module strat_grid
