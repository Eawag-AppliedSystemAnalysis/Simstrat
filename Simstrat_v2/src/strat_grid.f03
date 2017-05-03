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
  integer :: nz_inuse
  integer :: nz_grid_max
  real(RK) :: z_zero
  real(RK) :: lake_level_old
  real(RK) :: depth

contains
    procedure, pass :: init => grid_init
    procedure, pass :: memory_init => grid_memory_init
    procedure, pass :: init_grid_points => grid_init_grid_points
    procedure, pass :: init_z_axes => grid_init_z_axes
    procedure, pass :: init_morphology => grid_init_morphology
    procedure, pass :: init_areas => grid_init_areas
    procedure, pass :: update_area_factors => grid_update_area_factors
    procedure, pass :: update_depth => grid_update_depth
    procedure, pass :: interpolate_to_upp => grid_interpolate_to_upp
    procedure, pass :: interpolate_to_cent => grid_interpolate_to_cent
    procedure, pass :: interpolate_from_upp => grid_interpolate_from_upp
    procedure, pass :: interpolate_from_cent => grid_interpolate_from_cent
end type


contains

  subroutine grid_init(self, config)
    implicit none
    class(StaggeredGrid), intent(inout) :: self
    class(GridConfig), intent(inout) :: config

    ! Assign config
    self%nz_grid = config%nz_grid
    self%nz_grid_max = config%nz_grid_max

    ! Use read config to determine grid size
    call self%init_morphology(config)
    call self%init_grid_points(config)

    ! Allocate arrays according to size
    call self%memory_init()

    ! Call init functions
    call self%init_z_axes()
    call self%init_areas(config)
    call self%update_area_factors()
  end subroutine grid_init

  subroutine grid_memory_init(self)
    implicit none
    class(StaggeredGrid), intent(inout) :: self

    associate(nz_grid => self%nz_grid)

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
  subroutine grid_init_grid_points(self, config)
    implicit none
    integer i
    class(StaggeredGrid), intent(inout) :: self
    class(GridConfig), intent(inout) :: config

    self%nz_grid=int(config%nz_grid)

    if (.not. config%equidistant_grid) then   ! variable spacing

      ! Grid size given through number of grid_read
      self%nz_grid = config%nz_grid

      !Include top value if not included
      if (config%grid_read(0)/=0.) then
        self%nz_grid=self%nz_grid+1
        do i=self%nz_grid,1,-1
            config%grid_read(i)=config%grid_read(i-1)
        end do
        config%grid_read(0)=0.0_RK
      end if

      !If maxdepth grid larger than morphology
      if (config%grid_read(self%nz_grid)>config%depth) then
          do while ((config%grid_read(self%nz_grid)>config%depth).and.(self%nz_grid>0.))
              self%nz_grid=self%nz_grid-1
          end do
      end if

      !Include bottom value if not included
      if (config%grid_read(self%nz_grid)< config%depth) then
          self%nz_grid=self%nz_grid+1
          config%grid_read(self%nz_grid)=config%depth
      end if
    end if

    !Construct H
    allocate(self%h(0:self%nz_grid))
    self%h(0) = 0                ! Note that h(0) has no physical meaning but helps with some calculations

    if (config%equidistant_grid) then
      ! Equidistant grid
      self%h(1:self%nz_grid) = config%depth/self%nz_grid
    else
       ! Set up h according to configuraiton
      do i=1,self%nz_grid
         self%h(1+self%nz_grid-i)=config%grid_read(i)-config%grid_read(i-1)
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
     self%lake_level_old = self%z_upp(self%nz_inuse)
     write(*,*) "Warning, nz_inuse not set yet"
  end subroutine grid_init_z_axes


  subroutine grid_init_morphology(self, config)
    implicit none
    class(StaggeredGrid), intent(inout) :: self
    class(GridConfig), intent(inout) :: config
    integer :: num_read

    num_read = size(config%A_read)

    self%z_zero = config%z_A_read(0)
    self%lake_level_old =self%z_zero !needed?
    config%z_A_read(0:num_read-1) = self%z_zero - config%z_A_read(0:num_read-1)      ! z-coordinate is positive upwards, zero point is at reservoir bottom

  end subroutine grid_init_morphology

  subroutine grid_init_areas(self, config)
    implicit none
    class(StaggeredGrid), intent(inout) :: self
    class(GridConfig), intent(inout) :: config

    integer :: num_read
    associate(nz_grid => self%nz_grid, dAz => self%dAz,  z_upp => self%z_upp,Az => self%Az)
    num_read = size(config%A_read)

    ! Interpolate area (A) at all depths (z_upp)
    call Interp(config%z_A_read, config%A_read, num_read-1, self%z_upp, Az, nz_grid)

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
              nz => self%nz_inuse)


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


  subroutine grid_update_depth(self, new_depth)
    implicit none
    class(StaggeredGrid), intent(inout) :: self
    real(RK), intent(inout) :: new_depth
    integer :: i
    real(RK) ::  zmax
    associate(nz_grid => self%nz_grid, &
              nz_inuse => self%nz_inuse, &
              z_upp => self%z_upp, &
              z_cent => self%z_cent, &
              h => self%h, &
              z_zero => self%z_zero)

    do i=0,nz_grid
        if (z_upp(i) >= (z_zero-new_depth)) then    ! If above initial water level
            zmax = z_upp(i)
            z_upp(i)= z_zero-new_depth
            z_cent(i)= (z_upp(i)+z_upp(i-1))/2
            h(i) = z_upp(i) - z_upp(i-1)
            !Az(i) = Az(i-1) + h(i)/(zmax-z_upp(i-1))*(Az(i)-Az(i-1))
            nz_inuse = i
            if (h(nz_inuse)<=0.5*h(nz_inuse-1)) then         ! If top box is too small
                z_upp(nz_inuse-1) = z_upp(nz_inuse)                ! Combine the two upper boxes
                z_cent(nz_inuse-1) = (z_upp(nz_inuse-1)+z_upp(nz_inuse-2))/2
                h(nz_inuse-1)  = h(nz_inuse)+h(nz_inuse-1)
                nz_inuse = nz_inuse-1                        ! Reduce number of boxes
            end if
            exit
        end if
    end do

    end associate
  end subroutine


  ! Interpolates values of y on grid z onto array yi and grid z_cent
  subroutine grid_interpolate_to_cent(self, z,y,num_z,yi)
    implicit none
    class(StaggeredGrid), intent(in) :: self
    real(RK), dimension(:), intent(in) :: z,y
    real(RK), dimension(:), intent(out) :: yi
    integer, intent(in) :: num_z

    call Interp(z, y, num_z, self%z_cent, yi, self%nz_grid-1)
  end subroutine

  subroutine grid_interpolate_to_upp(self, z,y,num_z,yi)
    implicit none
    class(StaggeredGrid), intent(in) :: self
    real(RK), dimension(:), intent(in) :: z,y
    real(RK), dimension(:), intent(out) :: yi
    integer, intent(in) :: num_z

    call Interp(z, y, num_z, self%z_upp, yi, self%nz_grid)
  end subroutine

  subroutine grid_interpolate_from_cent(self, z, y, num_z, yi)
    class(StaggeredGrid), intent(in) :: self
    real(RK), dimension(:), intent(in) :: z,y
    real(RK), dimension(:), intent(out) :: yi
    integer, intent(in) :: num_z

    ! TO do: Interp or Interp_NAN for grid boundaries??
    call Interp(self%z_cent, y, self%nz_grid-1, z, yi, num_z)
  end subroutine

  subroutine grid_interpolate_from_upp(self, z, y, num_z, yi)
    class(StaggeredGrid), intent(in) :: self
    real(RK), dimension(:), intent(in) :: z,y
    real(RK), dimension(:), intent(out) :: yi
    integer, intent(in) :: num_z

    call Interp(self%z_upp, y, self%nz_grid, z, yi, num_z)
  end subroutine

end module strat_grid
