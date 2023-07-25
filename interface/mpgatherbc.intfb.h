interface
SUBROUTINE MPGATHERBC(IRECV, NSCFLD, &
 & FL1, EMEAN, THQ, FMEAN, &
 & FLPTS, EMPTS, TQPTS, FMPTS)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 use yowcpbo , only:&
 & nbounc
 use yowgrid , only:&
 & nproma_wam,&
 & nchnk
 use yowparam , only:&
 & nang,&
 & nfre
 INTEGER(KIND=JWIM) :: IRECV, NSCFLD
 REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(IN) :: FL1
 REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NCHNK), INTENT(IN) :: EMEAN, FMEAN, THQ
 REAL(KIND=JWRB), DIMENSION(NBOUNC), INTENT(INOUT) :: EMPTS, TQPTS, FMPTS
 REAL(KIND=JWRB), DIMENSION(NBOUNC, NANG, NFRE), INTENT(INOUT) :: FLPTS
END SUBROUTINE MPGATHERBC
end interface
