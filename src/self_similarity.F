      SUBROUTINE self_similarity(SelfSim,
     &         iAngCen,AlphCen,CCen,RCen,LCen,NCen)
c*****************************************************c
c                                                     c
c     Computes Self Similarity                        c
c                                                     c
c                                                     c
c     INPUT:                                          c
c                                                     c
c       LCen      # total centered functions          c
c       NCen      # non zero centered functions       c
c       iAngCen   Centered ang indices                c
c       AlphCen   Centered exponents                  c
c       CCen      Centered coefficients               c
c       RCen      Centered positions                  c
c                                                     c
c                                                     c
c     OUTPUT:                                         c
c                                                     c
c       SelfSim   Self Similarity                     c
c                                                     c
c                                                     c
c       DCC: 01/20                                    c
c                                                     c
c*****************************************************c
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER  iAngCen(LCen,3),iAng(3),NAngTot
      REAL*8   AlphCen(*),CCen(*),RCen(LCen,3),
     &         SelfSim,Rself(3),
     &         Dself,Sself
      SAVE Zero,Half,One,Two,Four,Eight
      DATA Zero,Half,One,Two,Four,Eight
     &    / 0d0, 0.5d0, 1d0, 2d0, 4d0, 8d0 /

c Starting
      SelfSim = Zero

c  Start loops over functions
      DO i1 = 1,NCen
        DO i2 = i1,NCen
!       DO i2 = 1,NCen
c----------------------------------------------

c  Square distance between centers
        Dself = Zero
        DO j = 1,3
          Dself = Dself + (RCen(i1,j)-RCen(i2,j))**2
        ENDDO

c  Angular part
        DO j = 1,3
          iAng(j) = iAngCen(i1,j)+iAngCen(i2,j)
        ENDDO

c  Compute constants
        Alpha = AlphCen(i1)+AlphCen(i2)
        Cself = exp(-Dself*AlphCen(i1)*AlphCen(i2)/Alpha)
        DO j = 1,3
          Rself(j) = (AlphCen(i1)*RCen(i1,j)+
     &              AlphCen(i2)*RCen(i2,j))/Alpha
        ENDDO

c  Compute overlap integrals
        CALL overlap(iAng,Alpha,Rself,Sself)
        Sself = Sself*CCen(i1)*CCen(i2)*Cself
        IF (i1.ne.i2) THEN
          Sself = Two * Sself
        ENDIF
        SelfSim = SelfSim + Sself
c----------------------------------------------
        ENDDO
      ENDDO
 
      END
