interface
 REAL(KIND=JWRB) FUNCTION A(XI, XJ, THI, THJ)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB) :: RI, RJ, RK, XI, XJ, THI, THJ, THK, OI, OJ, &
 & OK, FI, FJ, FK, VABS, VDIR, OMEG, A1, A3
 END FUNCTION A
 REAL(KIND=JWRB) FUNCTION B(XI, XJ, THI, THJ)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB) :: DEL, RI, RJ, RK, XI, XJ, THI, THJ, THK, OI, &
 & OJ, OK, FI, FJ, FK, VABS, VDIR, OMEG, A2
 END FUNCTION B
 REAL(KIND=JWRB) FUNCTION C_QL(XK0, XK1, TH0, TH1)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB) :: XK0, XK1, TH0, TH1, OM1, F1, OMEG, B2, B3
 END FUNCTION C_QL
 REAL(KIND=JWRB) FUNCTION U(XI, XJ, XK, XL, THI, THJ, THK, THL)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB) :: XI, XJ, XK, XL, THI, THJ, THK, THL, &
 & OI, OJ, OK, OL, XIK, XJK, XIL, XJL, &
 & OIK, OJK, OIL, OJL, QI, QJ, QIK, QJK, &
 & QIL, QJL, SQIJKL, ZCONST, VABS, OMEG
 END FUNCTION U
 REAL(KIND=JWRB) FUNCTION W2(XI, XJ, XK, XL, THI, THJ, THK, THL)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB) :: XI, XJ, XK, XL, THI, THJ, THK, THL, U
 END FUNCTION W2
 REAL(KIND=JWRB) FUNCTION V2(XI, XJ, XK, XL, THI, THJ, THK, THL)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB) :: DEL1, XI, XJ, XK, XL, THI, THJ, THK, THL, &
 & OI, OJ, OK, OL, RI, RJ, RK, RL, &
 & RIJ, RIK, RLI, RJL, RJK, RKL, THIJ, &
 & THIK, THLI, THJL, THJK, THKL, OIJ, &
 & OIK, OJL, OJK, OLI, OKL, XNIK, XNJL, &
 & XNJK, XNIL, YNIL, YNJK, YNJL, YNIK, &
 & ZNIJ, ZNKL, ZPIJ, ZPKL, THLJ, THIL, THKJ, &
 & THKI, THJI, THLK, VABS, VDIR, &
 & OMEG
 END FUNCTION V2
 REAL(KIND=JWRB) FUNCTION W1(XI, XJ, XK, XL, THI, THJ, THK, THL)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB) :: XI, XJ, XK, XL, THI, THJ, THK, THL, U
 END FUNCTION W1
 REAL(KIND=JWRB) FUNCTION W4(XI, XJ, XK, XL, THI, THJ, THK, THL)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB) :: XI, XJ, XK, XL, THI, THJ, THK, THL, U
 END FUNCTION W4
 REAL(KIND=JWRB) FUNCTION B3(XI, XJ, XK, XL, THI, THJ, THK, THL)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB) :: DEL1, XI, XJ, XK, XL, THI, THJ, THK, THL, &
 & OI, OJ, OK, OL, RI, RJ, RK, RL, &
 & RIJ, RJI, RIK, RKI, RLJ, RJL, RJK, RKJ, &
 & RLI, RIL, RLK, RKL, THIJ, THJI, &
 & THIK, THKI, THLJ, THJL, THJK, THKJ, &
 & THLI, THIL, THLK, THKL, ZIJKL, &
 & VABS, VDIR, A1, A3, W1, OMEG
 END FUNCTION B3
 REAL(KIND=JWRB) FUNCTION B4(XI, XJ, XK, XL, THI, THJ, THK, THL)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB) :: DEL1, XI, XJ, XK, XL, THI, THJ, THK, THL, &
 & OI, OJ, OK, OL, RI, RJ, RK, RL, &
 & RIJ, RIK, RIL, RJL, RJK, RKL, THIJ, THIK, &
 & THIL, THJL, THJK, THLK, THKL, &
 & ZIJKL, VABS, VDIR, &
 & A1, A3, W4, OMEG
 END FUNCTION B4
 REAL(KIND=JWRB) FUNCTION B1(XI, XJ, XK, XL, THI, THJ, THK, THL)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB) :: DEL1, XI, XJ, XK, XL, THI, THJ, THK, THL, &
 & OI, OJ, OK, OL, RI, RJ, RK, RL, &
 & RIJ, RJI, RIK, RKI, RJL, RJK, RLI, &
 & RIL, RKL, THIJ, THJI, THIK, &
 & THKI, THJL, THJK, THLI, THIL, THKL, ZIJKL, &
 & VABS, VDIR, A1, A3, W1, OMEG
 END FUNCTION B1
 REAL(KIND=JWRB) FUNCTION B2(XI, XJ, XK, XL, THI, THJ, THK, THL)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB) :: DEL1, XI, XJ, XK, XL, THI, THJ, THK, THL, &
 & OI, OJ, OK, OL, RI, RJ, RK, RL, &
 & RIJ, RIK, RKI, RJL, RLJ, RJK, RKJ, &
 & RLI, RIL, RKL, THIJ, &
 & THIK, THKI, THJL, THLJ, THJK, THKJ, &
 & THLI, THIL, THKL, ZIJKL, &
 & VABS, VDIR, A1, A3, OMEG
 END FUNCTION B2
 REAL(KIND=JWRB) FUNCTION A1(XI, XJ, XK, THI, THJ, THK)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB) :: DEL1, XI,XJ, XK, THI, THJ, THK, OI, OJ, OK, &
 & OMEG
 END FUNCTION A1
 REAL(KIND=JWRB) FUNCTION A2(XI, XJ, XK, THI, THJ, THK)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB) :: DEL1, XI, XJ, XK, THI, THJ, THK, A1
 END FUNCTION A2
 REAL(KIND=JWRB) FUNCTION A3(XI,XJ,XK,THI,THJ,THK)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB) :: DEL1, OI, OJ, OK, XI, XJ, XK, &
 & THI, THJ, THK, OMEG
 END FUNCTION A3
 REAL(KIND=JWRB) FUNCTION OMEG(X)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB) :: D, XK, X, T
 END FUNCTION OMEG
 REAL(KIND=JWRB) FUNCTION VABS(XI, XJ, THI, THJ)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB) :: XI, XJ, THI, THJ, ARG
 END FUNCTION VABS
 REAL(KIND=JWRB) FUNCTION VDIR(XI, XJ, THI, THJ)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB) :: XI, XJ, THI, THJ, EPS, Y, X
 END FUNCTION VDIR
end interface
