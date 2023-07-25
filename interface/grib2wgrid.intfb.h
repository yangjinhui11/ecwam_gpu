interface
SUBROUTINE GRIB2WGRID (IU06, KPROMA,                                &
     &                 KGRIB_HANDLE, KGRIB, ISIZE,                  &
     &                 LLUNSTR,                                     &
     &                 NGY, KRGG, KLONRGG_LOC,                      &
     &                 NXS, NXE, NYS, NYE,                          &
     &                 XLON, YLAT,                                  &
     &                 PMISS, PPREC, PPEPS,                         &
     &                 CDATE, IFORP, IPARAM, KZLEV, KKK, MMM, FIELD)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
INTEGER(KIND=JWIM), INTENT(IN) :: IU06, KPROMA, KGRIB_HANDLE, ISIZE
      INTEGER(KIND=JWIM), INTENT(IN) :: NGY, KRGG
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      INTEGER(KIND=JWIM), INTENT(IN) :: KGRIB(ISIZE)
      INTEGER(KIND=JWIM), INTENT(IN) :: KLONRGG_LOC(NGY)
      INTEGER(KIND=JWIM), INTENT(OUT) :: IFORP, IPARAM, KZLEV, KKK, MMM

      REAL(KIND=JWRB), INTENT(IN) :: PMISS, PPREC, PPEPS
      REAL(KIND=JWRB) ,DIMENSION(NXS:NXE, NYS:NYE), INTENT(IN) :: XLON, YLAT
      REAL(KIND=JWRB) ,DIMENSION(NXS:NXE, NYS:NYE), INTENT(OUT) :: FIELD

      CHARACTER(LEN=14), INTENT(OUT) :: CDATE
      LOGICAL, INTENT(IN) :: LLUNSTR
END SUBROUTINE GRIB2WGRID
end interface
