interface
SUBROUTINE OUTWINT(BOUT)
 use parkind_wave, only:&
 & jwrb
 use parkind_wave, only:&
 & jwrb
 use yowcout , only:&
 & niprmout
 use yowgrid , only:&
 & nproma_wam,&
 & nchnk
 REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NIPRMOUT, NCHNK), INTENT(IN) :: BOUT
 END SUBROUTINE OUTWINT
end interface
