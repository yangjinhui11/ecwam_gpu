! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE PROPAGS2_MOD

CONTAINS

SUBROUTINE PROPAGS2_SCC (F1, F3, NINF, NSUP, KIJS, KIJL, NANG, ND3S, ND3E)

! ----------------------------------------------------------------------

!**** *PROPAGS2* -  ADVECTION USING THE CORNER TRANSPORT SCHEME IN SPACE

!*    PURPOSE.
!     --------

!       COMPUTATION OF A PROPAGATION TIME STEP.

!**   INTERFACE.
!     ----------

!       *CALL* *PROPAGS2(F1, F3, NINF, NSUP, KIJS, KIJL, NANG, ND3S, ND3E)*
!          *F1*          - SPECTRUM AT TIME T (with exchange halo).
!          *F3*          - SPECTRUM AT TIME T+DELT
!          *NINF:NSUP+1* - 1st DIMENSION OF F1 and F3
!          *KIJS*        - ACTIVE INDEX OF FIRST POINT
!          *KIJL*        - ACTIVE INDEX OF LAST POINT
!          *NANG*        - NUMBER OF DIRECTIONS
!          *ND3S*        - FREQUENCY INDEX SOLVED BY THIS CALL ND3S:ND3E
!          *ND3E*

!     METHOD.
!     -------

!!  need text !!!!!!!!

!     EXTERNALS.
!     ----------


!     REFERENCE.
!     ----------

!      See IFS Documentation, part VII 

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : COSTH    ,SINTH
      USE YOWSTAT  , ONLY : ICASE    ,IREFRA
      USE YOWTEST  , ONLY : IU06
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KCOR      ,             &
     &            WLATN    ,WLONN    ,WCORN    ,WKPMN    ,WMPMN     ,   &
     &            LLWLATN  ,LLWLONN  ,LLWCORN  ,LLWKPMN  ,LLWMPMN   ,   &
     &            SUMWN    ,                                            &
     &            JXO      ,JYO      ,KCR      ,KPM      ,MPM

! ----------------------------------------------------------------------

      IMPLICIT NONE


      REAL(KIND=JWRB),DIMENSION(NINF:NSUP+1, NANG, ND3S:ND3E), INTENT(IN) :: F1
      REAL(KIND=JWRB),DIMENSION(NINF:NSUP+1, NANG, ND3S:ND3E), INTENT(OUT) :: F3
      INTEGER(KIND=JWIM), INTENT(IN) :: NINF, NSUP
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), INTENT(IN) :: NANG, ND3S, ND3E 


      INTEGER(KIND=JWIM) :: K, M, IJ
      INTEGER(KIND=JWIM) :: IC, ICR, ICL 
      INTEGER(KIND=JWIM) :: KP1, KM1, MM1, MP1, KNS, KEW
      INTEGER(KIND=JWIM) :: JJK, JJY, JJX

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: FJ1, FJ2, FJ3, FJ4, FJ5
      REAL(KIND=JWRB):: RFJ1,RFJ2,RFJ3,RFJ4,RFJ5
! ----------------------------------------------------------------------
!$acc routine vector
!*    SPHERICAL OR CARTESIAN GRID?
!     ----------------------------
      IF (ICASE == 1) THEN

!*      SPHERICAL GRID.
!       ---------------
        IF (IREFRA /= 2 .AND. IREFRA /= 3 ) THEN
!*      WITHOUT DEPTH OR/AND CURRENT REFRACTION.
!       ----------------------------------------
!$acc loop vector collapse(2)  independent
        DO M = ND3S, ND3E
          DO K = 1, NANG
            JJX=JXO(K,1)
            JJY=JYO(K,1)
            JJY=JYO(K,1)
            JJK=KCR(K,1)
              
              DO IJ = KIJS, KIJL
                RFJ1= F1(KLON(IJ,JJX)  ,K  ,M)
                RFJ2= F1(KLAT(IJ,JJY,1),K  ,M)
                RFJ3= F1(KLAT(IJ,JJY,2),K  ,M)
                RFJ4= F1(KCOR(IJ,JJK,1),K  ,M)
                RFJ5= F1(KCOR(IJ,JJK,2),K  ,M)
                F3(IJ,K,M) =                                            &
     &                (1.0_JWRB-SUMWN(IJ,K,M))* F1(IJ           ,K  ,M) &
!    &         + WLONN(IJ,K,M,JXO(K,1)) * F1(KLON(IJ,JXO(K,1))  ,K  ,M) &
!    &         +WLATN(IJ,K,M,JYO(K,1),1)* F1(KLAT(IJ,JYO(K,1),1),K  ,M) &
!    &         +WLATN(IJ,K,M,JYO(K,1),2)* F1(KLAT(IJ,JYO(K,1),2),K  ,M) &
!    &         +       WCORN(IJ,K,M,1,1)* F1(KCOR(IJ,KCR(K,1),1),K  ,M) &
!    &         +       WCORN(IJ,K,M,1,2)* F1(KCOR(IJ,KCR(K,1),2),K  ,M) &
!    &         + WLONN(IJ,K,M,JJX) * F1(KLON(IJ,JJX)  ,K  ,M)           &
!    &         +WLATN(IJ,K,M,JJY,1)* F1(KLAT(IJ,JJY,1),K  ,M)           &
!    &         +WLATN(IJ,K,M,JJY,2)* F1(KLAT(IJ,JJY,2),K  ,M)           &
!    &         +       WCORN(IJ,K,M,1,1)* F1(KCOR(IJ,JJK,1),K  ,M)      &
!    &         +       WCORN(IJ,K,M,1,2)* F1(KCOR(IJ,JJK,2),K  ,M)      &
     &         + WLONN(IJ,K,M,JJX) * RFJ1                            &
     &         +WLATN(IJ,K,M,JJY,1)* RFJ2                            &
     &         +WLATN(IJ,K,M,JJY,2)* RFJ3                            &
     &         +       WCORN(IJ,K,M,1,1)* RFJ4                       &
     &         +       WCORN(IJ,K,M,1,2)* RFJ5

              ENDDO

              DO IC=-1,1,2
                IF (LLWKPMN(K,M,IC)) THEN
                  DO IJ = KIJS, KIJL
                    F3(IJ,K,M) = F3(IJ,K,M)                             &
     &         +      WKPMN(IJ,K,M,IC)* F1(IJ,KPM(K,IC),M)
                  ENDDO!seq
                ENDIF 
              ENDDO !seq

            ENDDO! work
          ENDDO! vector

        ELSE
!*      DEPTH AND CURRENT REFRACTION.
!       -----------------------------

!$acc loop vector
          DO M = ND3S, ND3E
            DO K = 1, NANG
!$acc loop seq
              DO IJ = KIJS, KIJL
                F3(IJ,K,M) = (1.0_JWRB-SUMWN(IJ,K,M))* F1(IJ,K,M)
              ENDDO
!$acc loop seq
              
              DO IC=1,2
                IF (LLWLONN(K,M,IC)) THEN
                  DO IJ = KIJS, KIJL
                    F3(IJ,K,M) = F3(IJ,K,M)                             &
     &           +      WLONN(IJ,K,M,IC)*F1(KLON(IJ,IC) ,K,M)
                  ENDDO
                ENDIF
              ENDDO

!$acc loop seq
              DO ICL=1,2
                DO IC=1,2
                  IF (LLWLATN(K,M,IC,ICL)) THEN
                    DO IJ = KIJS, KIJL
                      F3(IJ,K,M) = F3(IJ,K,M)                           &
     &         +        WLATN(IJ,K,M,IC,ICL)*F1(KLAT(IJ,IC,ICL) ,K,M)
                    ENDDO
                 ENDIF
                ENDDO
              ENDDO

!$acc loop seq
              DO ICL=1,2
                DO ICR=1,4
                  IF (LLWCORN(K,M,ICR,ICL)) THEN
                    DO IJ = KIJS, KIJL
                      F3(IJ,K,M) = F3(IJ,K,M)                           &
     &         +   WCORN(IJ,K,M,ICR,ICL)*F1(KCOR(IJ,KCR(K,ICR),ICL),K,M)
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO

!$acc loop seq
              DO IC=-1,1,2
                IF (LLWKPMN(K,M,IC)) THEN
                  DO IJ = KIJS, KIJL
                    F3(IJ,K,M) = F3(IJ,K,M)                             &
     &         +      WKPMN(IJ,K,M,IC)* F1(IJ,KPM(K,IC),M)
                  ENDDO
                ENDIF
              ENDDO

!$acc loop seq
              DO IC=-1,1,2
                IF (LLWMPMN(K,M,IC)) THEN
                  DO IJ = KIJS, KIJL
                    F3(IJ,K,M) = F3(IJ,K,M)                             &
     &         +      WMPMN(IJ,K,M,IC)* F1(IJ,K  ,MPM(M,IC))
                  ENDDO
                ENDIF
              ENDDO


            ENDDO

          ENDDO

        ENDIF

      ELSE
!*    CARTESIAN GRID.
!     ---------------
        IF (IREFRA == 2 .OR. IREFRA == 3 ) THEN
!*      WITHOUT DEPTH OR/AND CURRENT REFRACTION.
!       ----------------------------------------
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '* PROPAGS2:                              *'
          WRITE (IU06,*) '* CORNER TRANSPORT SCHEME NOT YET READY  *' 
          WRITE (IU06,*) '* FOR  CARTESIAN GRID !                  *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          STOP
        ELSE
!*      DEPTH AND CURRENT REFRACTION.
!       ----------------------------
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '* PROPAGS2:                              *'
          WRITE (IU06,*) '* CORNER TRANSPORT SCHEME NOT YET READY  *' 
          WRITE (IU06,*) '* FOR  CARTESIAN GRID !                  *'
          WRITE (IU06,*) '* FOR DEPTH OR/AND CURRENT REFRACTION !  *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          STOP
        ENDIF

      ENDIF



END SUBROUTINE PROPAGS2_SCC
END MODULE
