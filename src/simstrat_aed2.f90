!     +---------------------------------------------------------------+
!     |  Simstrat - AED2 interface
!     +---------------------------------------------------------------+

module simstrat_aed2
   use strat_simdata
   use utilities
   use aed2_common
   use aed2_core

   implicit none
   private

   type, public :: SimstratAED2
      class(AED2Config), pointer :: cfg

   contains
      procedure, pass(self), public :: init
      procedure, pass(self), public :: update
   end type SimstratAED2

contains

   subroutine init(self, aed2_cfg)
      implicit none
      class(SimstratAED2) :: self
      class(AED2Config), target :: aed2_cfg

      ! Local variables
      character(len=80) :: fname
      type(aed2_variable_t),pointer :: tvar

      character(len=64) :: models(64)
      namelist /aed2_models/ models
      integer i, status


      fname = 'aed2.nml'

      if ( aed2_init_core('.') /= 0 ) call error("Initialisation of aed2_core failed")
      call aed2_print_version

      ! Create model tree
      write (6,*) "     Processing aed2_models config from ", trim(fname)
      open(50,file=fname,action='read',status='old',iostat=status)
      if ( status /= 0 ) then
         call error("Cannot open file " // trim(fname))
         stop
      end if

      models = ''
      read(50, nml=aed2_models, iostat=status)
      if ( status /= 0 ) then
         call error("Cannot read namelist entry aed2_models")
         stop
      end if

      do i=1,size(models)
         if (models(i)=='') exit
         call aed2_define_model(models(i), 50)
         write(6,*) models(i)
      end do

      !# should be finished with this file
      close(50)
      write (6,*) "      AED2 file parsing completed."



   end subroutine


   subroutine update(self, aed2_cfg)
      implicit none
      class(SimstratAED2) :: self
      class(AED2Config), target :: aed2_cfg


   end subroutine



end module simstrat_aed2