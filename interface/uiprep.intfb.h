interface
 SUBROUTINE UIPREP (IFORM, LLGRID)
 use parkind_wave, only:&
 & jwim
 INTEGER(KIND=JWIM), INTENT(OUT) :: IFORM
 LOGICAL, INTENT(OUT) :: LLGRID
 END SUBROUTINE UIPREP
end interface