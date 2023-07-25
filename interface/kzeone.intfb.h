interface
 SUBROUTINE KZEONE(X, Y, RE0, IM0, RE1, IM1)
 use parkind_wave, only:&
 & jwru
 REAL(KIND=JWRU) :: X, Y, X2, Y2, RE0, IM0, RE1, IM1, &
 & R1, R2, T1, T2, P1, P2, RTERM, ITERM, L
 END SUBROUTINE KZEONE
end interface
