test_suite simstrat_solver


setup
 !nothing
end setup

teardown
  ! This code runs immediately after each test
end teardown

! Example test using all six assertions
test funit_assertions
  real(8), dimension(7) :: rhs = (/ 4, 8, 12, 16, 20, 24, 20 /)
  real(8), dimension(7) :: x = (/ 1, 2, 3, 4, 5, 6, 7/)
  real(8), dimension(6) :: upper_diag = (/ 1,1,1,1,1,1 /)
  real(8), dimension(7) :: middle_diag = 2
  real(8), dimension(6) :: lower_diag = (/1,1,1,1,1,1/)
  real(8), dimension(7) :: calculated_x = 0
  integer :: i

  call solve_tridiag_thomas(lower_diag, middle_diag, upper_diag, rhs, calculated_x, 7)

  do i=1,7
    assert_real_equal(calculated_x(i),x(i))
  end do

end test

end test_suite
