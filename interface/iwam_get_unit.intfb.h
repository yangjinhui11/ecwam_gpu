interface
INTEGER(KIND=JWIM) FUNCTION IWAM_GET_UNIT (KUSO, CDNAME, CDACCESS, CDFORM, KRECL, CDACTION)
 use parkind_wave, only:&
 & jwim
 INTEGER(KIND=JWIM), INTENT(IN) :: KUSO
 CHARACTER(LEN=*), INTENT(IN) :: CDNAME, CDACTION
 CHARACTER(LEN=1), INTENT(IN) :: CDACCESS, CDFORM
 INTEGER(KIND=JWIM), INTENT(IN) :: KRECL
END FUNCTION IWAM_GET_UNIT
end interface
