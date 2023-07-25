      SUBROUTINE IFSTOWAM_NEW (BLK2LOC,NFIELDS, NGPTOTG, NCA, NRA, &
     &                            FIELDS, LWCUR, MASK_IN,&
     &                            NXS,NXE,NYS,NYE,FIELDG ) ! Yangjinhui rewrited in 2023.04.21

        call abor1("IFSTOWAM_NEW should never be called in singe alone model")
      END SUBROUTINE IFSTOWAM_NEW
