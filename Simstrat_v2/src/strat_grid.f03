!     +---------------------------------------------------------------+
!     | Grid module
!     |  - Contains class to store, use and modify the grid
!     +---------------------------------------------------------------+


module strat_grid
   use strat_kinds
   use utilities
   implicit none
   private

   ! Class that holds initial configuration for grid setup
   type, public :: GridConfig
      integer :: nz_grid_max
      integer :: nz_grid
      real(RK) :: depth
      real(RK), dimension(:), allocatable :: grid_read ! Grid definition
      real(RK), dimension(:), allocatable :: A_read ! Area definition
      real(RK), dimension(:), allocatable :: z_A_read ! Height of area definitions
      logical :: equidistant_grid
   end type

   ! StaggeredGrid implementation
   type, public :: StaggeredGrid
      real(RK), dimension(:), allocatable :: h        ! Box height
      real(RK), dimension(:), allocatable :: z_face   ! Holds z-values of faces
      real(RK), dimension(:), allocatable :: z_volume ! Holds z-values of volume centers
      real(RK), dimension(:), allocatable :: Az       ! Areas
      real(RK), dimension(:), allocatable :: dAz      ! Difference of areas
      real(RK), dimension(:), allocatable :: meanint  ! ?
      real(RK) :: volume

      !Area factors
      real(RK), dimension(:), allocatable :: AreaFactor_1
      real(RK), dimension(:), allocatable :: AreaFactor_2
      real(RK), dimension(:), allocatable :: AreaFactor_k1
      real(RK), dimension(:), allocatable :: AreaFactor_k2
      real(RK), dimension(:), allocatable :: AreaFactor_eps
      
      integer :: nz_grid      ! Number of allocated grid cells
      integer :: nz_occupied  ! number of grid cells in use as per current lake depth
      integer :: nz_grid_max  ! Hard limit of grid cells for reading files of unknnown length etc

      integer :: ubnd_vol, ubnd_fce, l_vol, l_fce   ! Upper and lenght for volume (vol) and face(fce) grids
      real(RK) :: z_zero
      real(RK) :: lake_level_old
      real(RK) :: depth

   contains
      !Many methods...
      ! Init methods used for setup
      procedure, pass :: init => grid_init                              ! Main initialization method
      procedure, pass :: memory_init => grid_memory_init
      procedure, pass :: init_grid_points => grid_init_grid_points
      procedure, pass :: init_z_axes => grid_init_z_axes
      procedure, pass :: init_morphology => grid_init_morphology
      procedure, pass :: init_areas => grid_init_areas

      ! Update methods, used for internal update of parameters
      procedure, pass :: update_area_factors => grid_update_area_factors  ! Called to recalculate area factors
      procedure, pass :: update_depth => grid_update_depth                ! Called for recalculating depth
      procedure, pass :: update_nz => grid_update_nz                      ! updates lower and upper bounds based on nz_occupied

      ! Interpolation methods
      procedure, pass :: interpolate_to_face => grid_interpolate_to_face    ! Interpolate quantitiy onto face grid
      procedure, pass :: interpolate_to_face_from_second => grid_interpolate_to_face_from_second  ! Interpolate quantity onto face grid, ignoring first value
      procedure, pass :: interpolate_to_vol => grid_interpolate_to_vol      ! Interpolate quantity onto volume grid
      procedure, pass :: interpolate_from_face => grid_interpolate_from_face  ! Interpolate quantity that is stored on face grid onto arbitrary grid
      procedure, pass :: interpolate_from_vol => grid_interpolate_from_vol    !  and so on

      ! Manipulation methods
      procedure, pass :: grow => grid_grow          ! Add a new box
      procedure, pass :: shrink => grid_shrink      ! Shrink grid by one box (= merge uppermost 2 boxes)
      procedure, pass :: modify_top_box => grid_modify_top_box  ! Change size of topmost box to reflect niveau change


   end type

contains

  ! Set up grid at program start
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
   end subroutine grid_init

   ! Allocate necessary memory based on sizes
   subroutine grid_memory_init(self)
      implicit none
      class(StaggeredGrid), intent(inout) :: self

      associate (nz_grid=>self%nz_grid)
         ! Allocate memory for initialization using nz_grid
         ! h is already allocated by grid point init
         allocate (self%z_volume(0:nz_grid)) ! Depth axis with center of boxes
         self%z_volume(0) = 0 ! trick for array access - index not in use

         allocate (self%z_face(nz_grid + 1)) ! Depth axis with faceer border of boxes
         allocate (self%Az(nz_grid + 1)) ! Az is defined on the faces
         allocate (self%dAz(nz_grid)) ! dAz is the difference between Az and thus defined on the volume

         ! Area factors used in calculations
         allocate (self%AreaFactor_1(nz_grid)) ! defined on volumes
         allocate (self%AreaFactor_2(nz_grid)) ! defined on volumes
         allocate (self%AreaFactor_k1(nz_grid)) ! defined on volumes
         allocate (self%AreaFactor_k2(nz_grid)) ! defined on volumes
         allocate (self%AreaFactor_eps(nz_grid)) ! defined on volumes
         allocate (self%meanint(nz_grid)) ! Inverse ratio of mean height of two adjacent boxes

      end associate
   end subroutine grid_memory_init

   ! Calculates h and definite number of nz_grid
   ! Depending on configuration
   subroutine grid_init_grid_points(self, config)
      implicit none
      integer i
      class(StaggeredGrid), intent(inout) :: self
      class(GridConfig), intent(inout) :: config

      self%nz_grid = int(config%nz_grid)

      if (.not. config%equidistant_grid) then ! variable spacing

         ! Grid size given through number of grid_read
         self%nz_grid = config%nz_grid

         !Include top value if not included
         if (config%grid_read(1) /= 0.) then
            write (*, *) "Grid invalid - top value not included"
            !Todo: check indexing of this code
            !  self%nz_grid=self%nz_grid+1
            !  do i=self%nz_grid,1,-1
            !      config%grid_read(i)=config%grid_read(i-1)
            !  end do
            !  config%grid_read(1)=0.0_RK
         end if

         !If maxdepth grid larger than morphology
         if (config%grid_read(self%nz_grid) > config%depth) then
            write (*, *) "Grid invalid - maxdepth of grid larger than morphology"
            !Todo: check indexing of this code
            !do while ((config%grid_read(self%nz_grid)>config%depth).and.(self%nz_grid>0.))
            !    self%nz_grid=self%nz_grid-1
            !end do
         end if

         !Include bottom value if not included
         if (config%grid_read(self%nz_grid) < config%depth) then
            write (*, *) "Grid invalid - Bottom value not included"
            !Todo: check indexing of this code
            !self%nz_grid=self%nz_grid+1
            !config%grid_read(self%nz_grid)=config%depth
         end if
      end if

      !Construct H
      allocate (self%h(0:self%nz_grid + 1))
      self%h(0) = 0 ! Note that h(0) has no physical meaning but helps with some calculations
      self%h(self%nz_grid + 1) = 0 ! Note that h(0) has no physical meaning but helps with some calculations
      if (config%equidistant_grid) then
         ! Equidistant grid
         self%h(1:self%nz_grid) = config%depth/self%nz_grid
      else
         ! Set up h according to configuraiton
         do i = 1, self%nz_grid
            self%h(1 + self%nz_grid - i) = config%grid_read(i) - config%grid_read(i - 1)
         end do
      end if
   end subroutine grid_init_grid_points

   ! Initializes z_volume and z_face
   subroutine grid_init_z_axes(self)
      implicit none
      class(StaggeredGrid), intent(inout) :: self
      integer :: i

      !Compute position of layer center and top
      self%z_volume(1) = 0.0_RK
      self%z_face(1) = 0.0_RK
      do i = 1, self%nz_grid
         self%z_volume(i) = self%z_volume(i - 1) + 0.5_RK*(self%h(i - 1) + self%h(i))
         self%z_volume(i) = nint(1e6_RK*self%z_volume(i))/1e6_RK
      end do

      do i = 2, self%nz_grid + 1
         self%z_face(i) = self%z_face(i - 1) + self%h(i - 1)
         self%z_face(i) = nint(1e6_RK*self%z_face(i))/1e6_RK
      end do

      ! needed?
      !self%lake_level_old = self%z_face(self%nz_occupied)
      write (*, *) "Warning, nz_occupied not set yet"
   end subroutine grid_init_z_axes

   ! Init morphology variables
   subroutine grid_init_morphology(self, config)
      implicit none
      class(StaggeredGrid), intent(inout) :: self
      class(GridConfig), intent(inout) :: config
      integer :: num_read

      num_read = size(config%A_read)

      self%z_zero = config%z_A_read(1) ! z_zero is the uppermost depth (might be above zero)
      self%lake_level_old = self%z_zero !needed?

      config%z_A_read = self%z_zero - config%z_A_read ! z-coordinate is positive upwards, zero point is at reservoir bottom

   end subroutine grid_init_morphology

   ! Init areas...
   subroutine grid_init_areas(self, config)
      implicit none
      class(StaggeredGrid), intent(inout) :: self
      class(GridConfig), intent(inout) :: config

      integer :: num_read
      associate (nz_grid=>self%nz_grid, dAz=>self%dAz, z_face=>self%z_face, Az=>self%Az)
         num_read = size(config%A_read)

         ! Interpolate area (A) at all depths (z_face)
         call Interp(config%z_A_read, config%A_read, num_read, self%z_face, Az, nz_grid + 1)

         ! Compute area derivative (= projected sediment area over layer thickness)
         dAz(1:nz_grid) = (Az(2:nz_grid + 1) - Az(1:nz_grid))/(z_face(2:nz_grid + 1) - z_face(1:nz_grid))

      end associate
   end subroutine

   ! Updates area factors (good method names are self-explanatory)
   subroutine grid_update_area_factors(self)
      implicit none
      class(StaggeredGrid), intent(inout) :: self

      integer :: i
      associate (Az=>self%Az, &
                 h=>self%h, &
                 nz=>self%nz_occupied)

         !todo: Verify array indexes and boundaries (especially h)
         self%AreaFactor_1(1:nz) = -4*Az(1:nz)/(h(1:nz) + h(0:nz - 1))/h(1:nz)/(Az(2:nz + 1) + Az(1:nz))
         self%AreaFactor_2(1:nz) = -4*Az(2:nz + 1)/(h(1:nz) + h(2:nz + 1))/h(1:nz)/(Az(2:nz + 1) + Az(1:nz))
         self%AreaFactor_k1(1:nz - 1) = -(Az(2:nz) + Az(3:nz + 1))/(h(1:nz - 1) + h(2:nz))/h(2:nz)/Az(2:nz)
         self%AreaFactor_k1(nz) = 0
         self%AreaFactor_k2(1:nz - 1) = -(Az(2:nz) + Az(1:nz - 1))/(h(1:nz - 1) + h(2:nz))/h(1:nz - 1)/Az(2:nz)
         self%AreaFactor_k2(nz) = 0
         self%AreaFactor_eps(1:nz - 1) = 0.5_RK*((Az(2:nz) - Az(1:nz - 1))/h(1:nz - 1) + (Az(3:nz + 1) - Az(2:nz))/h(2:nz))/Az(2:nz)
         self%AreaFactor_eps(nz) = 0
         self%meanint(1:nz) = 2.0_RK/(h(1:nz) + h(2:nz + 1))

         self%volume = 0
         do i = 1, nz
            self%volume = self%volume + 0.5_RK*h(i)*(Az(i) + Az(i + 1))
         end do
      end associate
   end subroutine grid_update_area_factors

   subroutine grid_modify_top_box(self, dh)
      implicit none
      class(StaggeredGrid), intent(inout) :: self
      real(RK) :: dh

      self%h(self%ubnd_vol) = self%h(self%ubnd_vol) + dh
      self%z_volume(self%ubnd_vol) = self%z_volume(self%ubnd_vol) + 0.5_RK*dh
      self%z_face(self%ubnd_fce) = self%z_face(self%ubnd_fce) + dh

   end subroutine

   !Shrink grid (mainly used by advection)
   subroutine grid_shrink(self, dh)
      implicit none
      class(StaggeredGrid), intent(inout) :: self
      real(RK) :: dh

      ! update grid
      self%z_face(self%ubnd_fce - 1) = self%z_face(self%ubnd_fce) + dh
      self%z_volume(self%ubnd_vol - 1) = 0.5_RK*(self%z_face(self%ubnd_fce - 1) + &
                                                 self%z_face(self%ubnd_fce - 2))

      self%h(self%ubnd_vol - 1) = (self%z_face(self%ubnd_fce - 1) - &
                                   self%z_face(self%ubnd_fce - 2))

      ! update number of occupied cells
      self%nz_occupied = self%nz_occupied - 1

      ! update boundaries (ubnd_fce and ubnd_vol)
      call self%update_nz()

   end subroutine

   !Grow grid (mainly used by advection)
   subroutine grid_grow(self, dh)
      implicit none
      class(StaggeredGrid), intent(inout) :: self
      real(RK) :: dh

      associate (z_face=>self%z_face, &
                 ubnd_fce=>self%ubnd_fce, &
                 z_volume=>self%z_volume, &
                 ubnd_vol=>self%ubnd_vol, &
                 h=>self%h)

         z_face(ubnd_fce + 1) = z_face(ubnd_fce) + dh
         z_face(ubnd_fce) = z_face(ubnd_fce) - (h(ubnd_vol) - dh)

         z_volume(ubnd_vol + 1) = z_volume(ubnd_vol) + 0.5_RK*dh
         z_volume(ubnd_vol) = z_volume(ubnd_vol) - 0.5_RK*(h(ubnd_vol) - dh)

         h(ubnd_vol + 1) = z_face(ubnd_fce + 1) - z_face(ubnd_fce)
         h(ubnd_vol) = z_face(ubnd_fce) - z_face(ubnd_fce - 1)

         ! update number of occupied cells
         self%nz_occupied = self%nz_occupied - 1

         call self%update_nz()

      end associate
   end subroutine

   ! Recalculate depth
   subroutine grid_update_depth(self, new_depth)
      implicit none
      class(StaggeredGrid), intent(inout) :: self
      real(RK), intent(inout) :: new_depth
      integer :: i
      real(RK) ::  zmax
      associate (nz_grid=>self%nz_grid, &
                 nz_occupied=>self%nz_occupied, &
                 z_face=>self%z_face, &
                 z_volume=>self%z_volume, &
                 h=>self%h, &
                 z_zero=>self%z_zero)

         do i = 1, nz_grid + 1
            if (z_face(i) >= (z_zero - new_depth)) then ! If above initial water level
               zmax = z_face(i)

               ! set top face to new water level
               z_face(i) = z_zero - new_depth

               !Adjust new volume center of cell i-1 (belonging to upper face i)
               z_volume(i - 1) = (z_face(i) + z_face(i - 1))/2
               h(i - 1) = z_face(i) - z_face(i - 1)
               ! todo: Why is Az not updated?
               !Az(i) = Az(i-1) + h(i)/(zmax-z_face(i-1))*(Az(i)-Az(i-1))
               nz_occupied = i - 1
               if (h(nz_occupied) <= 0.5*h(nz_occupied - 1)) then ! If top box is too small
                  z_face(nz_occupied) = z_face(nz_occupied + 1) ! Combine the two upper boxes
                  z_volume(nz_occupied - 1) = (z_face(nz_occupied) + z_face(nz_occupied - 1))/2
                  h(nz_occupied - 1) = h(nz_occupied) + h(nz_occupied - 1)
                  nz_occupied = nz_occupied - 1 ! Reduce number of boxes
               end if
               exit
            end if
         end do

         call self%update_nz()
      end associate
   end subroutine

   ! Interpolates values of y on grid z onto array yi and grid z_volume
   subroutine grid_interpolate_to_vol(self, z, y, num_z, yi)
      implicit none
      class(StaggeredGrid), intent(in) :: self
      real(RK), dimension(:), intent(in) :: z, y
      real(RK), dimension(:), intent(out) :: yi
      integer, intent(in) :: num_z
      real(RK), dimension(size(self%z_volume + 2)) :: z_vol_test
      call Interp(z, y, num_z, self%z_volume(1:self%nz_grid), yi, self%nz_grid)
   end subroutine

   subroutine grid_interpolate_to_face(self, z, y, num_z, yi)
      implicit none
      class(StaggeredGrid), intent(in) :: self
      real(RK), dimension(:), intent(in) :: z, y
      real(RK), dimension(:), intent(out) :: yi

      integer, intent(in) :: num_z
      call Interp(z, y, num_z, self%z_face, yi, self%nz_grid + 1)
   end subroutine

   subroutine grid_interpolate_to_face_from_second(self, z, y, num_z, yi)
      implicit none
      class(StaggeredGrid), intent(in) :: self
      real(RK), dimension(:), intent(in) :: z, y
      real(RK), dimension(:), intent(out) :: yi

      integer, intent(in) :: num_z
      call Interp(z, y, num_z, self%z_face(2:self%nz_grid + 1), yi, self%nz_grid)
   end subroutine

   subroutine grid_interpolate_from_vol(self, z, y, num_z, yi)
      class(StaggeredGrid), intent(in) :: self
      real(RK), dimension(:), intent(in) :: z, y
      real(RK), dimension(:), intent(out) :: yi
      integer, intent(in) :: num_z

      ! TO do: Interp or Interp_NAN for grid boundaries??
      call Interp(self%z_volume(1:self%nz_grid), y, self%nz_grid, z, yi, num_z)
   end subroutine

   subroutine grid_interpolate_from_face(self, z, y, num_z, yi)
      class(StaggeredGrid), intent(in) :: self
      real(RK), dimension(:), intent(in) :: z, y
      real(RK), dimension(:), intent(out) :: yi
      integer, intent(in) :: num_z

      call Interp(self%z_face, y, self%nz_grid + 1, z, yi, num_z)
   end subroutine

   ! Update all upper bounds and lengths
   subroutine grid_update_nz(self)
      implicit none
      class(StaggeredGrid), intent(inout) :: self

      self%ubnd_vol = self%nz_occupied
      self%ubnd_fce = self%nz_occupied + 1
      self%l_vol = self%nz_occupied
      self%l_fce = self%nz_occupied + 1
   end subroutine

   pure function grid_convert2height_above_sed(z, z_zero) result(h)
      implicit none
      real(RK), dimension(:), intent(in) :: z
      real(RK), intent(in) :: z_zero

      real(RK), dimension(size(z)) :: h
      integer :: n

      n = size(z)
      h = -z_zero + z(n:1:-1)
   end function

end module strat_grid
