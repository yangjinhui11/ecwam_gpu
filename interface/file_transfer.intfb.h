interface
 SUBROUTINE FILE_TRANSFER (IU06, IUNIT, FILENA, CDOPT)
 use parkind_wave, only:&
 & jwim
 INTEGER(KIND=JWIM) :: IU06, IUNIT
 CHARACTER(LEN=1) :: MODE, CDOPT
 CHARACTER(LEN=*) :: FILENA
 END SUBROUTINE FILE_TRANSFER
end interface
