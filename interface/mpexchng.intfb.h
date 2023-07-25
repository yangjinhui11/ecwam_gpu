interface
 SUBROUTINE MPEXCHNG(FLD, NDIM2, ND3S, ND3E)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 use yowmpp , only:&
 & ninf,&
 & nsup
 INTEGER(KIND=JWIM), INTENT(IN) :: NDIM2, ND3S, ND3E
 REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NDIM2, ND3S:ND3E), INTENT(INOUT) :: FLD
 END SUBROUTINE MPEXCHNG
end interface
