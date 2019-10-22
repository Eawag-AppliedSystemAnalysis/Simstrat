! Subroutines used by simstrat_aed2.f90 (mainly in the update function)
! Contains absorption and mobility algorithm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Light absorption feedback by AED2 variables

subroutine absorption_updateAED2(self, state)

   ! Arguments
   class(SimstratAED2) :: self
   class(ModelState) :: state

   ! Local variables
   integer :: i
   real(RK) :: bio_extinction
   real(RK), dimension(self%grid%nz_occupied) :: extc_coef

   do i=self%grid%nz_occupied, 1, -1
      bio_extinction = 0.0_RK
      call aed2_light_extinction(self%column, i, bio_extinction)
      extc_coef(i) = self%aed2_cfg%background_extinction + bio_extinction

   end do

   ! Interpolate to faces to be compatible with Simstrat temperature module
   call self%grid%interpolate_to_face(self%grid%z_volume, extc_coef, self%grid%nz_occupied, state%absorb)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The mobility algorithm is taken from the GLM code (http://aed.see.uwa.edu.au/research/models/GLM/)
!
! Assumptions:
! 1) movement direction has at most one change down the layers              *
! 2) sides of the lake slope inward (ie bottom is narrower than top)

subroutine mobility(self, state, min_C, settling_v, conc)
   ! Arguments
   class(SimstratAED2) :: self
   class(ModelState) :: state
   real(RK), intent(in) :: min_C
   real(RK), dimension(:), intent(in) :: settling_v
   real(RK), dimension(:), intent(inout) :: conc

   ! Local variables
   real(RK) :: dtMax, tdt, tmp
   real(RK), dimension(self%grid%ubnd_vol) :: mins, vols, Y
   integer :: dirChng, signum, i, count


   ! determine mobility timestep i.e. maximum time step that particles
   ! will not pass through more than one layer in a time step 
   ! (probably not used in Simstrat as the timestep needs to be small)

   dtMax = state%dt
   dirChng = 0   ! this represents the layer at which direction switches from sinking to rising or visa versa
   signum = sign(1.0_RK,settling_v(1))  ! positive for rising, negative for sinking

   do i = 1,self%grid%ubnd_vol
      ! for convenience
      vols(i) = self%grid%Az_vol(i)*self%grid%h(i)
      mins(i) = min_C*vols(i)
      Y(i) = conc(i)*vols(i)

      !look for the change of direction
      if (signum .ne. sign(1.0_RK,settling_v(i))) then
         signum = -signum
         dirChng = i-1
      end if

      ! check if all movement can be from within this cell
      if (abs(settling_v(i)*state%dt) > self%grid%h(i)) then
         tdt = self%grid%h(i)/abs(settling_v(i))
         if (tdt < dtMax) dtMax = tdt
      end if

      ! check if movement can all be into the next cell.
      ! if movement is settling, next is below, otherwise next is above.
      ! check also for top or bottom in case of oopsies

      if (settling_v(i) > 0.) then
         if ((i < self%grid%ubnd_vol) .and. (abs(settling_v(i))*state%dt) > self%grid%h(i+1)) then
             tdt = self%grid%h(i+1)/abs(settling_v(i))
             if(tdt < dtMax) dtMax = tdt
         end if
      else if ((i > 1) .and. (abs(settling_v(i))*state%dt) > self%grid%h(i-1)) then
         tdt = self%grid%h(i-1)/abs(settling_v(i))
         if(tdt < dtMax) dtMax = tdt
      end if
   end do ! end find maximum time step dtMax
   if (dirChng == 0 .and. settling_v(1) > 0. ) dirChng = self%grid%ubnd_vol ! all rising
   if (dirChng == 0 .and. settling_v(self%grid%ubnd_vol) < 0.) dirChng = self%grid%ubnd_vol ! all sinking

   tdt = dtMax
   count = 0
   do
      ! do this in steps of dtMax, but at least once
      ! each time tdt is dtMax, except, possibly, the last which is whatever was left.
      count = count + 1    ! counter
      if (count*dtMax > state%dt) tdt = state%dt - (count - 1)*dtMax ! Last timestep

      ! 2 possibilities
      ! 1) lower levels rising, upper levels sinking
      ! 2) lower levels sinking, upper levels rising
      if (settling_v(1) > 0. ) then ! lower levels rising
         if (settling_v(self%grid%ubnd_vol) < 0.) then !top levels are sinking
             call Sinking(self, Y, conc, settling_v, vols, mins, tdt, self%grid%ubnd_vol, dirChng, tmp)
             Y(dirChng) = Y(dirChng) + tmp
             conc(dirChng) = Y(dirChng)/vols(dirChng)
         end if
         call Rising(self, Y, conc, settling_v, vols, mins, tdt, 1, dirChng)
      else ! lower levels sinking
         call Sinking(self, Y, conc, settling_v, vols, mins, tdt, dirChng, 1, tmp)
         if ( settling_v(self%grid%ubnd_vol) > 0.) then !top levels are rising
            call Rising(self, Y, conc, settling_v, vols, mins, tdt, dirChng, self%grid%ubnd_vol)
         end if
         if (state%dt > 0.) exit
      end if
   end do

end subroutine

! Rising is the easier on the two since the slope means we dont need to look
! at relative areas (the cell above will always be >= to the current cell)
! all matter is moved to the next cell      
! for each cell except the top :  
! 1) calculate how much is going to move  
! 2) subtract amount that must now move 
! 3) add the amount moved from previous cell to current cell
! 4) fix concentration 
! for the top cell :
! 1) add the amount moved from previous cell to current cell
! 2) fix concentration

subroutine Rising(self, Y, conc, settling_v, vols, mins, dt, start_i, end_i)
   ! Arguments
   class(SimstratAED2) :: self
   real(RK), dimension(:), intent(inout) :: conc, Y
   real(RK), dimension(:), intent(in) :: settling_v, vols, mins
   real(RK) :: dt
   integer :: start_i, end_i

   ! Local variables
   real(RK) :: mov, moved
   integer i

   mov = 0.
   moved = 0.

   do i = start_i,end_i
      ! speed times time (=h) time area * concen = mass to move
      mov = (settling_v(i) * dt) * self%grid%Az_vol(i)*conc(i)
      ! if removing that much would bring it below min conc
      if ((Y(i) + moved - mov) < mins(i) ) mov = Y(i) + moved - mins(i)

      Y(i) = Y(i) + moved - mov;
      conc(i) = Y(i) / vols(i) ! return it to a concentration
      moved = mov ! for the next step
   end do
   ! nothing rises out of the end cell, but we still add that which came from below
   Y(end_i) = Y(end_i) + moved
   conc(end_i) = Y(end_i) / vols(end_i)

end subroutine


! for each cell :
! 1) calculate how much is going to move
! 2) subtract amount that must now move
! 3) add the amount moved from previous cell to current cell
! 4) fix concentration
! 5) compute the amount that will go to the next cell

subroutine Sinking(self, Y, conc, settling_v, vols, mins, dt, start_i, end_i, moved)
   ! Arguments
   class(SimstratAED2) :: self
   real(RK), dimension(:), intent(inout) :: conc, Y
   real(RK), dimension(:), intent(in) :: settling_v, vols, mins
   real(RK) :: dt
   real(RK), intent(out) :: moved
   integer :: start_i, end_i

   ! Local variables
   real(RK) :: mov
   integer :: i

   mov = 0.
   moved = 0.

   do i = start_i,end_i, -1
      ! speed times time (=h) time area * concen = mass to move
      mov = (abs(settling_v(i)) * dt) * self%grid%Az_vol(i) * conc(i)
      ! if removing that much would bring it below min conc
      if ((Y(i) + moved - mov) < mins(i)) mov = Y(i) + moved - mins(i)
      Y(i) = Y(i) + moved - mov
      conc(i) = Y(i) / vols(i) ! return it to a concentration

      ! so now mov has how much has moved out of the cell, but not all
      ! of that will go into the next cell (FB: part of it sediments on the flancs)
      if ( i > 1 )  then
         moved = mov * (self%grid%Az_vol(i-1) / self%grid%Az_vol(i)) ! for the next step
      else
         moved = mov ! we are about to exit anyway.
      end if
   end do
end subroutine