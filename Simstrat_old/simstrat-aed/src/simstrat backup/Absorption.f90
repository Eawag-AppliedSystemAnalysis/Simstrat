module Absorption
    use SimstratModel
    use utilities
    implicit none

    private
    public doAbsorption

contains
  !AG 2014: revision + correction
  !###############################################################
  subroutine doAbsorption(datum,absorb,idx)
  !###############################################################
        implicit none

        ! Global Variables
        real(dp), intent(in) :: datum
        real(dp), intent(inout) :: absorb(0:)
        integer, intent(in) :: idx

        ! Local Variables
        real(dp) :: tb_start, tb_end !Start and end time
        real(dp) :: z_absorb(0:nz_max), dummy !Read depths
        real(dp) :: absorb_read_start(0:nz_max), absorb_read_end(0:nz_max) !Read start and end values
        real(dp) :: absorb_start(0:nz_max), absorb_end(0:nz_max) !Interpolated start and end values
        integer :: eof,i,nval

        ! Save these variables to use them in the next iteration (MS 2014)
        save tb_start, tb_end
        save z_absorb, absorb_start, absorb_end
        save eof, nval


        if (idx==1) then    ! First iteration
            if(disp_dgn/=0) write(6,*) 'Light attenuation : ',trim(AbsorpName)
            if(disp_dgn==2) write(6,*) 'Starting to read light attenuation file...'
            open(30,status='old',file=AbsorpName)
            eof = 0
            !Read depths: first line are names, second line is number of depths available
            read(30,*,end=9)
            read(30,*,end=9) nval
            read(30,*,end=9) dummy, (z_absorb(i),i=0,nval-1)

            !Make depths positives
            do i=0,nval-1
               z_absorb(i) = abs(z_absorb(i))
            end do
            
            !Read first values
            read(30,*,end=9) tb_start, (absorb_read_start(i),i=0,nval-1)
            if(datum<tb_start) write(6,*) 'Warning: first light attenuation date after simulation start time.'
            call Interp(z_absorb, absorb_read_start, nval-1, z_upp, absorb_start, nz)
            read(30,*,end=7) tb_end, (absorb_read_end(i),i=0,nval-1)
            call Interp(z_absorb, absorb_read_end, nval-1, z_upp, absorb_end, nz)
        end if

        if (datum<=tb_start .or. eof==1) then !If datum before first date or end of file reached
            goto 8
        else
            do while (datum>tb_end)         !Move to appropriate interval to get correct value
                tb_start = tb_end
                absorb_start(0:nz) = absorb_end(0:nz)
                !Read next value
                read(30,*,end=7) tb_end, (absorb_read_end(i),i=0,nval-1)
                call Interp(z_absorb, absorb_read_end, nval-1, z_upp, absorb_end, nz)
            end do
            !Linearly interpolate value at correct datum (for all depths)
            absorb(0:nz) = absorb_start(0:nz) + (datum-tb_start) * (absorb_end(0:nz)-absorb_start(0:nz))/(tb_end-tb_start)
        end if
        return

  7     eof = 1
        if(datum>tb_start) write(6,*) 'Warning: last light attenuation date before simulation end time.'
  8     absorb(0:nz) = absorb_start(0:nz)           !Take first value of current interval
        return

  9     write(6,*) 'Error reading light attenuation file (no data found).'
        stop
  end subroutine doAbsorption

end module Absorption
