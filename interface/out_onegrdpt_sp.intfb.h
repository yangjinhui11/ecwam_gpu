interface
 SUBROUTINE OUT_ONEGRDPT_SP(SPEC,USTAR,CDTPRO)
 use parkind_wave, only:&
 & jwrb
 USE YOWPARAM , ONLY : NANG ,NFRE
 REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN) :: SPEC
 REAL(KIND=JWRB), INTENT(IN) :: USTAR
 CHARACTER(LEN=12), INTENT(IN) :: CDTPRO
 END SUBROUTINE OUT_ONEGRDPT_SP
end interface
