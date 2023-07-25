interface
SUBROUTINE WAMWND (KIJS, KIJL,                 &
     &                   IFROMIJ, JFROMIJ,           &
     &                   NXS, NXE, NYS, NYE, FIELDG, &
     &                   UCUR, VCUR,                 &
     &                   U10, US,                    &
     &                   THW, ADS,                   &
     &                   WSTAR, CITH,                &
     &                   LWCUR, ICODE_WND)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 USE YOWDRVTYPE , ONLY : FORCING_FIELDS
INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: IFROMIJ  ,JFROMIJ
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(IN) :: FIELDG
      INTEGER(KIND=JWIM), INTENT(IN) :: ICODE_WND
      REAL(KIND=JWRB), DIMENSION (KIJS:KIJL), INTENT(IN) :: UCUR, VCUR
      REAL(KIND=JWRB), DIMENSION (KIJS:KIJL), INTENT(INOUT) :: U10, US
      REAL(KIND=JWRB), DIMENSION (KIJS:KIJL), INTENT(OUT) :: THW, ADS, WSTAR, CITH
      LOGICAL, INTENT(IN) :: LWCUR
 END SUBROUTINE WAMWND
end interface
