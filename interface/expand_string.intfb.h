interface
 SUBROUTINE EXPAND_STRING( &
 & myproc, &
 & nproc, &
 & timestep, &
 & max_timestep, &
 & s, &
 & n)
 use parkind_wave, only:&
 & jwim
 integer(KIND=JWIM), intent(in) :: myproc, nproc
 integer(KIND=JWIM), intent(in) :: timestep, max_timestep
 integer(KIND=JWIM), intent(in) :: n
 character(len=*), intent(inout) :: s(n)
 END SUBROUTINE expand_string
end interface
