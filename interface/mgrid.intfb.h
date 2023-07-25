interface
 SUBROUTINE MGRID (BATHY)
 use parkind_wave, only:&
 & jwrb
 use yowparam , only:&
 & ngx,&
 & ngy
 REAL(KIND=JWRB), INTENT(INOUT) :: BATHY(NGX,NGY)
 END SUBROUTINE MGRID
end interface
