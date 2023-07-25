interface
 SUBROUTINE MUBUF (IU01, BATHY, IU08, NPROPAGS)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 use yowparam , only:&
 & ngx,&
 & ngy
 INTEGER(KIND=JWIM) :: IU01
 INTEGER(KIND=JWIM) :: NPROPAGS
 INTEGER(KIND=JWIM) :: IU08(0:NPROPAGS)
 REAL(KIND=JWRB) :: BATHY(NGX, NGY)
 END SUBROUTINE MUBUF
end interface
