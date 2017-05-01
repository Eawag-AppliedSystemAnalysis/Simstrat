module strat_grid
  use strat_kinds
  use utilities
  implicit none
  private

type, public :: GridConfig
  integer :: nz_grid_max
  integer :: nz_grid
  real(RK) :: depth
  real(RK), dimension(:), allocatable :: grid_read ! Grid definition
  real(RK), dimension(:), allocatable :: A_read ! Area definition
  real(RK), dimension(:), allocatable :: z_A_read ! Height of area definitions
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
  integer :: nz_grid
  real(RK) :: z_zero
  real(RK) :: lake_level_old

  type(GridConfig), allocatable :: config
contains
    procedure, pass :: init => grid_init
    procedure, pass :: memory_init => grid_memory_init
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
    self%config = config

    ! Use read config to determine grid size
    call self%init_morphology()
    call self%init_grid_points()

    ! Allocate arrays according to size
    call self%memory_init()

    ! Call init functions
    call self%init_z_axes()
    call self%init_areas()
    call self%update_area_factors()
  end subroutine grid_init

  subroutine grid_memory_init(self)
    implicit none
    class(StaggeredGrid), intent(inout) :: self

    associate(nz_grid => self%config%nz_grid)

    ! Allocate memory for initialization using nz_grid
    ! h is already allocated by grid point init
    allocate(self%z_cent(0:nz_grid))         ! Depth axis with center of boxes
    allocate(self%z_upp(0:nz_grid))          ! Depth axis with upper border of boxes
    allocate(self%Az(0:nz_grid))             ! Az is defined on z_upp
    allocate(self%dAz(1:nz_grid))           ! dAz is the difference between Az and thus on z_cent

    ! Area factors used in calculations
    allocate(self%AreaFactor_1(1:nz_grid))
    allocate(self%AreaFactor_2(1:nz_grid))
    allocate(self%AreaFactor_k1(1:nz_grid))
    allocate(self%AreaFactor_k2(1:nz_grid))
    allocate(self%AreaFactor_eps(1:nz_grid))

    allocate(self%meanint(0:nz_grid))        ! Inverse ratio of mean height of two adjacent boxes

    end associate
  end subroutine grid_memory_init

  ! Calculates h and definite number of nz_grid
  ! Depending on configuration
  subroutine grid_init_grid_points(self)
    implicit none
    integer i
    class(StaggeredGrid), intent(inout) :: self

    self%nz_grid=int(self%config%nz_grid)

    if (.not. self%config%equidistant_grid) then   ! variable spacing

      ! Grid size given through number of grid_read
      self%nz_grid = self%config%nz_grid

      !Include top value if not included
      if (self%config%grid_read(0)/=0.) then
        self%nz_grid=self%nz_grid+1
        do i=self%nz_grid,1,-1
            self%config%grid_read(i)=self%config%grid_read(i-1)
        end do
        self%config%grid_read(0)=0.0_RK
      end if

      !If maxdepth grid larger than morphology
      if (self%config%grid_read(self%nz_grid)>self%config%depth) then
          do while ((self%config%grid_read(self%nz_grid)>self%config%depth).and.(self%nz_grid>0.))
              self%nz_grid=self%nz_grid-1
          end do
      end if

      !Include bottom value if not included
      if (self%config%grid_read(self%nz_grid)< self%config%depth) then
          self%nz_grid=self%nz_grid+1
          self%config%grid_read(self%nz_grid)=self%config%depth
      end if
    end if

    !Construct H
    allocate(self%h(0:self%nz_grid))
    self%h(0) = 0                ! Note that h(0) has no physical meaning but helps with some calculations

    if (self%config%equidistant_grid) then
      ! Equidistant grid
      self%h(1:self%nz_grid) = self%config%depth/self%nz_grid
    else
       ! Set up h according to configuraiton
      do i=1,self%nz_grid
         self%h(1+self%nz_grid-i)=self%config%grid_read(i)-self%config%grid_read(i-1)
      end do
    end if
  end subroutine grid_init_grid_points

  ! Initializes z_cent and z_upp
  subroutine grid_init_z_axes(self)
    implicit none
    class(StaggeredGrid), intent(inout) :: self
    integer :: i

     !Compute position of layer center and top
     self%z_cent(0)=0.0_RK
     self%z_upp(0)=0.0_RK
     do i=1,self%nz_grid
         self%z_cent(i)=self%z_cent(i-1)+0.5_RK*(self%h(i-1)+self%h(i))
         self%z_upp(i)=self%z_upp(i-1)+self%h(i)
     end do
     do i=1,self%nz_grid
         self%z_cent(i)=nint(1e6_RK*self%z_cent(i))/1e6_RK
         self%z_upp(i)=nint(1e6_RK*self%z_upp(i))/1e6_RK
     end do
     ! needed?
     self%lake_level_old = self%z_upp(self%nz_grid)
  end subroutine grid_init_z_axes


  subroutine grid_init_morphology(self)
    implicit none
    class(StaggeredGrid), intent(inout) :: self
    integer :: num_read
    associate(cfg => self%config)

    num_read = size(self%config%A_read)

    self%z_zero = cfg%z_A_read(0)
    self%lake_level_old =self%z_zero !needed?
    cfg%z_A_read(0:num_read-1) = self%z_zero - cfg%z_A_read(0:num_read-1)      ! z-coordinate is positive upwards, zero point is at reservoir bottom

    end associate
  end subroutine grid_init_morphology

  subroutine grid_init_areas(self)
    implicit none
    class(StaggeredGrid), intent(inout) :: self
    integer :: num_read
    associate(cfg => self%config, nz_grid => self%nz_grid, dAz => self%dAz,  z_upp => self%z_upp,Az => self%Az)
    num_read = size(self%config%A_read)

    ! Interpolate area (A) at all depths (z_upp)
    call Interp(cfg%z_A_read, cfg%A_read, num_read-1, self%z_upp, Az, nz_grid)

   ! Compute area derivative (= projected sediment area over layer thickness)
    dAz(1:nz_grid) = (Az(1:nz_grid)-dAz(0:nz_grid-1))/(z_upp(1:nz_grid)-z_upp(0:nz_grid-1))

    end associate
  end subroutine


  subroutine grid_update_area_factors(self)
    implicit none
    class(StaggeredGrid), intent(inout) :: self

    integer :: i
    associate(Az => self%Az, &
              h => self%h, &
              nz => self%nz_grid)


    self%AreaFactor_1(1:nz) = -4*Az(0:nz-1)/(h(1:nz)+h(0:nz-1))/h(1:nz)/(Az(1:nz)+Az(0:nz-1))
    self%AreaFactor_2(1:nz) = -4*Az(1:nz)/(h(1:nz)+h(2:nz+1))/h(1:nz)/(Az(1:nz)+Az(0:nz-1))
    self%AreaFactor_k1(1:nz-1) = -(Az(1:nz-1)+Az(2:nz))/(h(1:nz-1)+h(2:nz))/h(2:nz)/Az(1:nz-1)
    self%AreaFactor_k2(1:nz-1) = -(Az(1:nz-1)+Az(0:nz-2))/(h(1:nz-1)+h(2:nz))/h(1:nz-1)/Az(1:nz-1)
    self%AreaFactor_eps(1:nz-1) = 0.5_RK*((Az(1:nz-1)-Az(0:nz-2))/h(1:nz-1)+(Az(2:nz)-Az(1:nz-1))/h(2:nz))/Az(1:nz-1)

    self%meanint(0:nz-1) = 2.0_RK/(h(0:nz-1)+h(1:nz))

    self%volume=0
    do i=0,nz-1
       self%volume = self%volume + 0.5_RK*h(i+1)*(Az(i)+Az(i+1))
    end do
  end associate
  end subroutine grid_update_area_factors

end module strat_grid
