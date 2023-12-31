interface
SUBROUTINE FLDINTER (IU06, NGPTOTG, NC, NR, NFIELDS,FIELDS,       &
     &                     NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,     &
     &                     IFROMIJ, JFROMIJ, KIJS, KIJL, KIJLMAX,       &
     &                     AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,   &
     &                     ILONRGG, IJBLOCK, PMISS,                     &
     &                     LADEN, ROAIR, LGUST, WSTAR0, LLKC, LWCUR,    &
     &                     LLINTERPOL,                                  &
     &                     DJ1M, DII1M, DIIP1M, JJM, IIM, IIPM, MASK_IN,&
     &                     NXS, NXE, NYS, NYE, FIELDG)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 USE YOWDRVTYPE , ONLY : FORCING_FIELDS
INTEGER(KIND=JWIM), INTENT(IN) :: IU06, NGPTOTG, NC, NR, NFIELDS
      INTEGER(KIND=JWIM), INTENT(IN) :: NGX, NGY, KRGG, KIJS, KIJL, KIJLMAX
      INTEGER(KIND=JWIM), INTENT(IN) :: IPERIODIC
      INTEGER(KIND=JWIM), DIMENSION(NGY), INTENT(IN) :: KLONRGG, JJM
      INTEGER(KIND=JWIM), DIMENSION(NR), INTENT(IN) :: ILONRGG
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: IFROMIJ, JFROMIJ
      INTEGER(KIND=JWIM), DIMENSION(0:NC+1,NR), INTENT(IN) :: IJBLOCK
      INTEGER(KIND=JWIM), DIMENSION(NGX,NGY), INTENT(IN) :: IIM, IIPM
      INTEGER(KIND=JWIM), DIMENSION(NGPTOTG), INTENT(INOUT) :: MASK_IN
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FIELDG


      REAL(KIND=JWRB), INTENT(IN) :: XDELLA, AMOWEP, AMOSOP, AMOEAP, AMONOP, PMISS
      REAL(KIND=JWRB), INTENT(IN) :: ROAIR, WSTAR0
      REAL(KIND=JWRB), DIMENSION(NGY), INTENT(IN) :: ZDELLO, DJ1M
      REAL(KIND=JWRB), DIMENSION(NGPTOTG,NFIELDS), INTENT(IN) :: FIELDS
      REAL(KIND=JWRB), DIMENSION(NGX,NGY), INTENT(IN) :: DII1M, DIIP1M

      LOGICAL, INTENT(IN):: LADEN, LGUST, LWCUR, LLKC, LLINTERPOL
 END SUBROUTINE FLDINTER
end interface
