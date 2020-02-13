
      SUBROUTINE make_dens(iAngDens,AlphDens,CDens,RDens,
     &   LDens,NDens,Alpha,CProd,Coef,iAng,ROrb,NOrbX,Tol)
c*****************************************************c
c                                                     c
c       Computes the product of 2 wavefunctions       c
c       to build the density matrix                   c
c                                                     c
c     INPUT:                                          c
c                                                     c
c       NOrbX     wavefunction dimension              c
c       ROrb      # center molec functions            c
c       iAng      angular indices                     c
c       Alpha     exponents                           c
c       CProd     # product coefficients              c
c       LDens     Density leading dimension           c
c                                                     c
c                                                     c
c     OUTPUT:                                         c
c                                                     c
c       iAngDens  Density ang indices                 c
c       AlphDens  Density exponents                   c
c       CDens     Density coefficients                c
c       RDens     Density positions                   c
c       NDens     Density non zero elements           c
c                                                     c
c                                                     c
c       DCC: 05/09                                    c
c                                                     c
c*****************************************************c
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER     ::  iAng(NOrbX,3),iAngFun(3),iAngDens(LDens,3)
      REAL*8      ::  AlphDens(*),CDens(*),Alpha(*),
     &                RFun(3),RDens(LDens,3),CProd(*),Coef(*),
     &                ROrb(NOrbX,3),Tol
      SAVE Zero,Half,One,Two,Four,Eight
      DATA Zero,Half,One,Two,Four,Eight
     &    / 0d0, 0.5d0, 1d0, 2d0, 4d0, 8d0 /

c Initialize some variables
      NDens = 0
      iDens = 1
      iOrb  = 1

      DO i1 = 1,NOrbX
        IJK1 = iAng(i1,1) + iAng(i1,2) + iAng(i1,3)
      DO i2 = i1,NOrbX
        IJK2 = iAng(i2,1) + iAng(i2,2) + iAng(i2,3)

c Squared distance between centers
        Dist = Zero
        DO j = 1,3
          Dist = Dist + (ROrb(i1,j)-ROrb(i2,j))**2
        ENDDO

c Angular part
        DO j = 1,3
          iAngFun(j) = iAng(i1,j) + iAng(i2,j)
        ENDDO

c Product constant
        AlphFun = Alpha(i1) + Alpha(i2)
        CFun = exp(-Dist*Alpha(i1)*Alpha(i2)/AlphFun)
        CFun = Coef(i1) * Coef(i2) * CFun * CProd(iOrb)
        DO j = 1,3
          RFun(j) = (Alpha(i1)*ROrb(i1,j)+Alpha(i2)*ROrb(i2,j))/
     &               AlphFun
        ENDDO
        iOrb = iOrb + 1
        IF (ABS(CFun) .GE. Tol) THEN
c-------------------------------------------
cc
c       S*S product
cc
        IF (IJK1.eq.0.and.IJK2.eq.0) THEN
          CDens(iDens)=Four*CFun
          IF (i1.eq.i2) THEN
            CDens(iDens)=Half*CDens(iDens)
          ENDIF
          AlphDens(iDens)=AlphFun
          DO j=1,3
            iAngDens(iDens,j)=iAngFun(j)
            RDens(iDens,j)=RFun(j)
          ENDDO
          iDens = iDens + 1
          NDens = NDens + 1
        ENDIF
c-------------------------------------------
cc
c       S*P product
cc
        IF ((IJK1+IJK2).eq.1) THEN
c         CFun = CFun
          CDens(iDens) = Zero
          DO  j = 1,3
           CDens(iDens) = CDens(iDens) - Four*CFun*
     .     (dfloat(iAng(i1,j)) * ROrb(i1,j) +
     .     dfloat(iAng(i2,j)) * ROrb(i2,j))
          ENDDO
          CDens(iDens+1) = Four*CFun
          IF (i1.eq.i2) THEN
            CDens(iDens)   = Half * CDens(iDens)
            CDens(iDens+1) = Half * CDens(iDens+1)
          ENDIF
          AlphDens(iDens)   = AlphFun
          AlphDens(iDens+1) = AlphFun
          DO j = 1,3
            iAngDens(iDens,j)   = 0
            iAngDens(iDens+1,j) = iAngFun(j)
            RDens(iDens,j)   = RFun(j)
            RDens(iDens+1,j) = RFun(j)
          ENDDO
          iDens = iDens + 2
          NDens = NDens + 2
        ENDIF
c-------------------------------------------
cc          
c       P*P product
cc        
        IF (IJK1.eq.1 .and. IJK2.eq.1) THEN
c         CFun=CFun
c    ---Funcions S, P i D---
          DO i=1,4
            CDens(iDens+i-1) = Zero
            AlphDens(iDens+i-1) = AlphFun
            DO j=1,3
              RDens(iDens+i-1,j) = RFun(j)
            ENDDO
          ENDDO
          DO j=1,3
            iAngDens(iDens,j) = 0
            iAngDens(iDens+1,j) = iAng(i1,j)
            iAngDens(iDens+2,j) = iAng(i2,j)
            iAngDens(iDens+3,j) = iAngFun(j)
            DO i = 1,3
              CDens(iDens) = CDens(iDens) +
     &        Four*CFun*dfloat(iAng(i1,j)) * ROrb(i1,j)*
     &        dfloat(iAng(i2,i))*ROrb(i2,i)
            ENDDO
            CDens(iDens+1) = CDens(iDens+1)-
     &      Four*CFun*dfloat(iAng(i2,j))*ROrb(i2,j)
            CDens(iDens+2) = CDens(iDens+2)-
     &      Four*CFun*dfloat(iAng(i1,j))*ROrb(i1,j)
          ENDDO
          CDens(iDens+3) = Four * CFun
          IF (i1.eq.i2) THEN
            DO i = 1,4
              CDens(iDens+i-1) = Half*CDens(iDens+i-1)
            ENDDO
          ENDIF
          iDens = iDens + 4
          NDens = NDens + 4
        ENDIF
c-------------------------------------------
cc
c       S*D product
cc
c    ---Funcions S, P i D---
        IF ((IJK1.eq.0 .and. IJK2.eq.2).or.
     &     (IJK1.eq.2 .and. IJK2.eq.0)) THEN
c       ------Case (0 0 2)
          IF ((iAngFun(1).eq.2).or.(iAngFun(2).eq.2).or.
     &       (iAngFun(3).eq.2)) THEN
            CDens(iDens)   = Zero
            CDens(iDens+1) = Zero
            DO i = 1,3
              iAngDens(iDens,i) = 0
              iAngDens(iDens+1,i) = iAngFun(i)/2
              iAngDens(iDens+2,i) = iAngFun(i)
              CDens(iDens) = Four*CFun*ROrb(i1,i)**(iAng(i1,i))
     .                *ROrb(i2,i)**(iAng(i2,i))+CDens(iDens)
              CDens(iDens+1) = -Four*CFun*(ROrb(i1,i)*
     .        dfloat(iAng(i1,i))+ROrb(i2,i)*(dfloat(iAng(i2,i))))+
     .        CDens(iDens+1)
              AlphDens(iDens+i-1) = AlphFun
              DO j = 1,3
                RDens(iDens+i-1,j) = RFun(j)
              ENDDO
            ENDDO
            CDens(iDens+2) = Four*CFun
            iDens = iDens + 3
            NDens = NDens + 3
          ELSE
c       ------Case (1 1 0)
            DO i = 1,4
              DO j = 1,3
                RDens(iDens+i-1,j) = RFun(j)
                iAngDens(iDens+i-1,j) = 0
              ENDDO
              AlphDens(iDens+i-1) = AlphFun
            ENDDO
            DO i = 1,2
              DO j = i+1,3
                IF ((iAngFun(i).eq.1).and.(iAngFun(j).eq.1)) THEN
                  iAngDens(iDens+1,i) = 1
                  iAngDens(iDens+2,j) = 1
                  CDens(iDens) = Four*CFun*(ROrb(i1,i)*
     .       dfloat(iAng(i1,i))+ROrb(i2,i)*dfloat(iAng(i2,i)))*
     .       (ROrb(i1,j)*dfloat(iAng(i1,j))+ROrb(i2,j)*
     .       dfloat(iAng(i2,j)))
                  CDens(iDens+1) = -Four*CFun*(ROrb(i1,j)*
     &       dfloat(iAng(i1,j))+ROrb(i2,j)*dfloat(iAng(i2,j)))
                  CDens(iDens+2) = -Four*CFun*(ROrb(i1,i)*
     &       dfloat(iAng(i1,i))+ROrb(i2,i)*dfloat(iAng(i2,i)))
                ENDIF
              ENDDO
            ENDDO
            CDens(iDens+3) = Four*CFun
            DO i = 1,3
              iAngDens(iDens+3,i) = iAngFun(i)
            ENDDO
            iDens = iDens + 4
            NDens = NDens + 4
          ENDIF
        ENDIF
c-------------------------------------------
cc            
c       P*D product
cc          
c    ---Funcions S, P, D i F---
        IF ((IJK1.eq.1 .and. IJK2.eq.2).or.
     &    (IJK1.eq.2 .and. IJK2.eq.1)) THEN
c       ------Case (1 0 0) (0 1 1)
          IF ((iAngFun(1).eq.1).and.(iAngFun(2).eq.1).and.
     &        (iAngFun(3).eq.1)) THEN
            DO j = 1,3
              DO i = 1,8
                RDens(iDens+i-1,j) = RFun(j)
              ENDDO
              iAngDens(iDens,j) = 0
              iAngDens(iDens+7,j) = 1
              iAngDens(iDens+1,j) = 0
              iAngDens(iDens+2,j) = 0
              iAngDens(iDens+3,j) = 0
              iAngDens(iDens+4,j) = 1
              iAngDens(iDens+5,j) = 1
              iAngDens(iDens+6,j) = 1
              iAngDens(iDens+j,j) = 1
              iAngDens(iDens+3+j,j) = 0
              CDens(iDens+3+j) = -ROrb(i1,j)**(iAng(i1,j))*
     &         ROrb(i2,j)**(iAng(i2,j))
            ENDDO
            CDens(iDens) = -ROrb(i1,1)**(iAng(i1,1))*
     &       ROrb(i2,1)**(iAng(i2,1))*ROrb(i1,2)**(iAng(i1,2))*
     &       ROrb(i1,3)**(iAng(i1,3))*ROrb(i2,2)**(iAng(i2,2))*
     &       ROrb(i2,3)**(iAng(i2,3))
            CDens(iDens+1) = ROrb(i1,2)**(iAng(i1,2))*
     &       ROrb(i1,3)**(iAng(i1,3))*ROrb(i2,2)**(iAng(i2,2))*
     &       ROrb(i2,3)**(iAng(i2,3))
            CDens(iDens+2) = ROrb(i1,1)**(iAng(i1,1))*
     &       ROrb(i1,3)**(iAng(i1,3))*ROrb(i2,1)**(iAng(i2,1))*
     &       ROrb(i2,3)**(iAng(i2,3))
            CDens(iDens+3) = ROrb(i1,1)**(iAng(i1,1))*
     &       ROrb(i1,2)**(iAng(i1,2))*ROrb(i2,1)**(iAng(i2,1))*
     &       ROrb(i2,2)**(iAng(i2,2))
            CDens(iDens+7) = One
            DO i = 1,8
              AlphDens(iDens+i-1) = AlphFun
              CDens(iDens+i-1) = Four*CFun*CDens(iDens+i-1)
            ENDDO
            iDens = iDens + 8
            NDens = NDens + 8
          ENDIF
c
c       ------Case (1 0 0) (1 1 0) and (1 0 0) (0 2 0)
          IF ((iAngFun(1).eq.2).or.(iAngFun(2).eq.2).or.
     .       (iAngFun(3).eq.2)) THEN
            DO j = 1,3
              DO i = 1,6
                iAngDens(iDens+i-1,j) = 0
                RDens(iDens+i-1,j) = RFun(j)
              ENDDO
              iAngDens(iDens+5,j) = iAngFun(j)
              iAngDens(iDens+3,j) = 1
              IF (iAngFun(j).eq.2) THEN
                iAngDens(iDens+4,j) = 2
                iAngDens(iDens+1,j) = 1
                CDens(iDens+3) = -(ROrb(i1,j)*
     .           dfloat(iAng(i1,j))+ROrb(i2,j)*dfloat(iAng(i2,j)))
                CDens(iDens+2) = ROrb(i1,j)**
     .           (iAng(i1,j))*ROrb(i2,j)**(iAng(i2,j))
              ENDIF
              IF (iAngFun(j).eq.1) THEN
                iAngDens(iDens+2,j) = 1
                CDens(iDens+4) = -ROrb(i1,j)**
     .           (iAng(i1,j))*ROrb(i2,j)**(iAng(i2,j))
              ENDIF
              IF (iAngFun(j).eq.0) THEN
                iAngDens(iDens+3,j) = 0
              ENDIF
            ENDDO
            CDens(iDens+1) = CDens(iDens+3)*
     .       CDens(iDens+4)
            CDens(iDens) = CDens(iDens+2)*
     .       CDens(iDens+4)
            CDens(iDens+5) = One
            DO i = 1,6
              AlphDens(iDens+i-1) = AlphFun
              CDens(iDens+i-1) = Four*CFun*CDens(iDens+i-1)
            ENDDO
            iDens=iDens+6
            NDens=NDens+6
          ENDIF
c       ------Case (1 0 0) (2 0 0)
          IF ((iAngFun(1).eq.3).or.(iAngFun(2).eq.3).or.
     .       (iAngFun(3).eq.3)) THEN
            DO i=1,4
              DO j=1,3
                iAngDens(iDens+i-1,j)=(i-1)*iAngFun(j)/3
                RDens(iDens+i-1,j)=RFun(j)
              ENDDO
              AlphDens(iDens+i-1)=AlphFun
            ENDDO
            CDens(iDens)=-ROrb(i1,1)**(iAng(i1,1))*
     .       ROrb(i2,1)**(iAng(i2,1))*ROrb(i1,2)**(iAng(i1,2))*
     .       ROrb(i1,3)**(iAng(i1,3))*ROrb(i2,2)**(iAng(i2,2))*
     .       ROrb(i2,3)**(iAng(i2,3))
            CDens(iDens+1) = Zero
            CDens(iDens+2) = Zero
            DO i = 1,3
              CDens(iDens+1) = ROrb(i1,i)*ROrb(i2,i)*
     .         dfloat(iAng(i1,i)*iAng(i2,i))+CDens(iDens+1)
              CDens(iDens+2) = -(dfloat(iAng(i1,i))*ROrb(i1,i)+
     .         dfloat(iAng(i2,i))*ROrb(i2,i))+CDens(iDens+2)
              IF (iAng(i1,i).eq.2) THEN
                CDens(iDens+1) = CDens(iDens+1)+ROrb(i1,i)**2
              ENDIF
              IF (iAng(i2,i).eq.2) THEN
                CDens(iDens+1) = CDens(iDens+1)+ROrb(i2,i)**2
              ENDIF
            ENDDO
            CDens(iDens+3) = One
          DO i = 1,4
            CDens(iDens+i-1) = Four*CFun*CDens(iDens+i-1)
          ENDDO
          iDens=iDens+4
          NDens=NDens+4
          ENDIF
        ENDIF
c-------------------------------------------
cc
c       D*D product
c    ---S, P, D, F and G functions---
        IF ((IJK1.eq.2) .and. (IJK2.eq.2)) THEN

c       ------Case (2 0 0) (2 0 0)
          IF ((iAngFun(1).eq.4).or.(iAngFun(2).eq.4).or.
     &       (iAngFun(3).eq.4)) THEN
            DO i = 1,5
              DO j = 1,3
                iAngDens(iDens+i-1,j) = (i-1)*iAngFun(j)/4
                RDens(iDens+i-1,j) = RFun(j)
              ENDDO
              AlphDens(iDens+i-1) = AlphFun
            ENDDO
            CDens(iDens) = Four*CFun*ROrb(i1,1)**(iAng(i1,1))*
     &       ROrb(i2,1)**(iAng(i2,1))*ROrb(i1,2)**(iAng(i1,2))*
     &       ROrb(i1,3)**(iAng(i1,3))*ROrb(i2,2)**(iAng(i2,2))*
     &       ROrb(i2,3)**(iAng(i2,3))
            DO i=1,3
              IF (iAngFun(i).eq.4) THEN
                CDens(iDens+1) = -Eight*CFun*(ROrb(i1,i)*
     &           ROrb(i2,i)**2+ROrb(i2,i)*ROrb(i1,i)**2)
                CDens(iDens+2) = Four*CFun*(ROrb(i1,i)**2+
     &           ROrb(i2,i)**2+Four*ROrb(i1,i)*ROrb(i2,i))
                CDens(iDens+3) = -Eight*CFun*(ROrb(i1,i)+
     &           ROrb(i2,i))
              ENDIF
            ENDDO
            CDens(iDens+4) = Four*CFun
            IF (i1.eq.i2) THEN
              DO i=1,5
                CDens(iDens+i-1)=Half*CDens(iDens+i-1)
              ENDDO
            ENDIF
            iDens=iDens+5
            NDens=NDens+5
          ELSE
c
c       ------Case (1 1 0) (2 0 0)
          IF ((iAngFun(1).eq.3).or.(iAngFun(2).eq.3).or.
     &     (iAngFun(3).eq.3)) THEN
            DO j = 1,3
              DO i = 1,8
                RDens(iDens+i-1,j) = RFun(j)
              ENDDO
              iAngDens(iDens,j)   = 0
              iAngDens(iDens+1,j) = 0
              iAngDens(iDens+2,j) = 0
              iAngDens(iDens+4,j) = 0
              iAngDens(iDens+6,j) = 0
              iAngDens(iDens+7,j) = iAngFun(j)
              IF (iAngFun(j).eq.3) THEN
                iAngDens(iDens+1,j) = 1
                iAngDens(iDens+4,j) = 2
                iAngDens(iDens+6,j) = 3
                CDens(iDens+2) = -(ROrb(i1,j)**
     &            iAng(i1,j))*(ROrb(i2,j)**iAng(i2,j))
                CDens(iDens+5) = -(dfloat(iAng(i1,j))*
     &           ROrb(i1,j)+dfloat(iAng(i2,j))*ROrb(i2,j))
              ENDIF
              IF (iAngFun(j).eq.1) THEN
                iAngDens(iDens+2,j) = 1
                CDens(iDens+6) = -(ROrb(i1,j)*
     &           dfloat(iAng(i1,j))+dfloat(iAng(i2,j))*ROrb(i2,j))
              ENDIF
              iAngDens(iDens+3,j) = iAngDens(iDens+1,j)+
     &         iAngDens(iDens+2,j)
              iAngDens(iDens+5,j) = iAngDens(iDens+4,j)+
     &         iAngDens(iDens+2,j)
              IF (iAng(i1,j).eq.2) THEN
                CDens(iDens+3)=(ROrb(i1,j)**2+
     &           Two*ROrb(i1,j)*ROrb(i2,j))
              ENDIF
              IF (iAng(i2,j).eq.2) THEN
                CDens(iDens+3)=(ROrb(i2,j)**2+
     &           Two*ROrb(i1,j)*ROrb(i2,j))
              ENDIF
            ENDDO
            CDens(iDens) = ROrb(i1,1)**(iAng(i1,1))*
     &       ROrb(i2,1)**(iAng(i2,1))*ROrb(i1,2)**(iAng(i1,2))*
     &       ROrb(i1,3)**(iAng(i1,3))*ROrb(i2,2)**(iAng(i2,2))*
     &       ROrb(i2,3)**(iAng(i2,3))
            CDens(iDens+1) = CDens(iDens+3)*CDens(iDens+6)
            CDens(iDens+4) = CDens(iDens+5)*CDens(iDens+6)
            CDens(iDens+7) = One
            DO i = 1,8
              CDens(iDens+i-1) = Four*CFun*CDens(iDens+i-1)
              AlphDens(iDens+i-1)=AlphFun
            ENDDO
            iDens = iDens + 8
            NDens = NDens + 8
         ELSE
c       Aquest else es refereix a que no hi hagi cap 3 a iAngFun
c       ------Case (1 1 0) (0 1 1) i (110)(002)
          IF ((iAngFun(1).eq.1).or.(iAngFun(2).eq.1).or.
     &      (iAngFun(3).eq.1)) THEN
           DO i = 1,12
             AlphDens(iDens+i-1) = AlphFun
             DO j = 1,3
               RDens(iDens+i-1,j) = RFun(j)
             ENDDO
           ENDDO
c       ------Nomes cas (110)(002)
           IF (iAng(i1,1)*iAng(i2,1).eq.0.and.iAng(i1,2)*
     &      iAng(i2,2).eq.0.and.iAng(i1,3)*iAng(i2,3).eq.0) THEN
            ij = 1
            CDens(iDens) = One
            DO j = 1,3
              iAngDens(iDens,j)   = 0
              iAngDens(iDens+1,j) = 0
              iAngDens(iDens+2,j) = 0
              iAngDens(iDens+3,j) = iAngFun(j)-1
              iAngDens(iDens+6,j) = 1
              iAngDens(iDens+7,j) = 0
              iAngDens(iDens+10,j) = 1
              iAngDens(iDens+11,j) = iAngFun(j)
              IF (iAngFun(j).eq.1) THEN
                iAngDens(iDens+ij,j) = 1
                CDens(iDens+7+ij) = -One*
     &          ROrb(i1,j)**(iAng(i1,j))*
     &          ROrb(i2,j)**(iAng(i2,j))
                ij = 2
              ENDIF
              iAngDens(iDens+8,j) = iAngFun(j)-
     &         iAngDens(iDens+1,j)
              iAngDens(iDens+9,j) = iAngFun(j)-
     &         iAngDens(iDens+2,j)
              IF (iAngFun(j).eq.2) THEN
                iAngDens(iDens+6,j) = 0
                iAngDens(iDens+7,j) = 2
                CDens(iDens+6) = ROrb(i1,j)**
     &           (iAng(i1,j))*ROrb(i2,j)**(iAng(i2,j))
                CDens(iDens+10) = -(ROrb(i1,j)*
     &           dfloat(iAng(i1,j))+dfloat(iAng(i2,j))*ROrb(i2,j))
              ENDIF
              iAngDens(iDens+4,j) = 1-iAngDens(iDens+1,j)
              iAngDens(iDens+5,j) = 1-iAngDens(iDens+2,j)
c
              CDens(iDens) = CDens(iDens)*ROrb(i1,j)**
     &         (iAng(i1,j))*ROrb(i2,j)**(iAng(i2,j))
            ENDDO
            CDens(iDens+7) = CDens(iDens+8)*
     &       CDens(iDens+9)
            CDens(iDens+4) = CDens(iDens+8)*
     &       CDens(iDens+10)
            CDens(iDens+5) = CDens(iDens+9)*
     &       CDens(iDens+10)
            CDens(iDens+3) = CDens(iDens+5)*
     &       CDens(iDens+8)
            CDens(iDens+2) = CDens(iDens+6)*
     &       CDens(iDens+8)
            CDens(iDens+1) = CDens(iDens+6)*
     &       CDens(iDens+9)
            CDens(iDens+11) = One
          DO i = 1,12
           CDens(iDens+i-1) = Four*CFun*CDens(iDens+i-1)
          ENDDO
            iDens = iDens + 12
            NDens = NDens + 12
          ELSE
c       ------Nomes cas (110)(011)
          DO j = 1,3
            iAngDens(iDens,j)   = 0
            iAngDens(iDens+1,j) = 0
            iAngDens(iDens+2,j) = 0
            iAngDens(iDens+3,j) = 0
            iAngDens(iDens+4,j) = iAng(i2,j)
            iAngDens(iDens+5,j) = iAng(i1,j)
            iAngDens(iDens+7,j) = 0
            iAngDens(iDens+8,j) = iAngFun(j)
            iAngDens(iDens+9,j) = iAngFun(j)
            iAngDens(iDens+10,j) = 1
            iAngDens(iDens+11,j) = iAngFun(j)
            IF (iAng(i1,j).eq.0) THEN
              iAngDens(iDens+2,j) = 1
              iAngDens(iDens+8,j) = 0
              CDens(iDens+8) = -ROrb(i2,j)
            ENDIF
            IF (iAng(i2,j).eq.0) THEN
              iAngDens(iDens+1,j) = 1
              iAngDens(iDens+9,j) = 0
              CDens(iDens+9) = -ROrb(i1,j)
            ENDIF
            IF (iAngFun(j).eq.2) THEN
              iAngDens(iDens+3,j) = 1
              iAngDens(iDens+7,j) = 2
              CDens(iDens+6) = ROrb(i1,j)*ROrb(i2,j)
              CDens(iDens+10) = -(ROrb(i1,j)+ROrb(i2,j))
            ENDIF
            iAngDens(iDens+6,j) = 1-iAngDens(iDens+3,j)
          ENDDO
          CDens(iDens+4) = CDens(iDens+10)*
     &     CDens(iDens+9)
          CDens(iDens+5) = CDens(iDens+10)*
     &     CDens(iDens+8)
          CDens(iDens+7) = CDens(iDens+8)*
     &     CDens(iDens+9)
          CDens(iDens+1) = CDens(iDens+6)*
     &     CDens(iDens+8)
          CDens(iDens+2) = CDens(iDens+6)*
     &     CDens(iDens+9)
          CDens(iDens+3) = CDens(iDens+5)*
     &     CDens(iDens+9)
          CDens(iDens) = CDens(iDens+6)*
     &     CDens(iDens+7)
          CDens(iDens+11) = One
          DO i = 1,12
            CDens(iDens+i-1) = Four*CFun*CDens(iDens+i-1)
          ENDDO
            iDens = iDens + 12
            NDens = NDens + 12
          ENDIF
c       ------Case (200)(020) and (110)(110)
         ELSE
          DO i = 1,9
            DO j = 1,3
              RDens(iDens+i-1,j) = RFun(j)
            ENDDO
            AlphDens(iDens+i-1) = AlphFun
          ENDDO
c       ------Nomes cas (200)(020)
          IF ((iAng(i1,1).eq.2).or.(iAng(i1,2).eq.2).or.
     &     (iAng(i1,3).eq.2)) THEN
            CDens(iDens+7) = Zero
            CDens(iDens+6) = Zero
            DO j = 1,3
              iAngDens(iDens,j) = 0
              iAngDens(iDens+1,j) = 0
              iAngDens(iDens+2,j) = 0
              iAngDens(iDens+3,j) = iAng(i1,j)
              iAngDens(iDens+4,j) = iAng(i2,j)
              iAngDens(iDens+5,j) = iAngFun(j)/2
              iAngDens(iDens+8,j) = iAngFun(j)
              IF (iAng(i1,j).eq.2) THEN
               iAngDens(iDens+1,j) = 1
               CDens(iDens+4) = ROrb(i1,j)**2
              ENDIF
              IF (iAng(i2,j).eq.2) THEN
               iAngDens(iDens+2,j) = 1
               CDens(iDens+3) = ROrb(i2,j)**2
              ENDIF
              iAngDens(iDens+6,j) = iAngDens(iDens+5,j)+
     &          iAngDens(iDens+2,j)
              iAngDens(iDens+7,j) = iAngDens(iDens+5,j)+
     &          iAngDens(iDens+1,j)
              CDens(iDens+7) = CDens(iDens+7)-ROrb(i2,j)*
     &         dfloat(iAng(i2,j))
              CDens(iDens+6) = CDens(iDens+6)-ROrb(i1,j)*
     &         dfloat(iAng(i1,j))
            ENDDO
            CDens(iDens+8) = One
            CDens(iDens+5) = CDens(iDens+6)*CDens(iDens+7)
            CDens(iDens+1) = CDens(iDens+6)*CDens(iDens+3)
            CDens(iDens+2) = CDens(iDens+4)*CDens(iDens+7)
            CDens(iDens) = CDens(iDens+3)*CDens(iDens+4)
          DO i = 1,9
            CDens(iDens+i-1) = Four*CFun*CDens(iDens+i-1)
          ENDDO
            iDens = iDens + 9
            NDens = NDens + 9
c       ------Nomes cas (110)(110)
          ELSE
          ij = 1
          CDens(iDens+7) = Zero
          CDens(iDens+6) = Zero
          CDens(iDens+4) = Zero
          CDens(iDens+3) = Zero
          DO j = 1,3
            iAngDens(iDens,j)   = 0
            iAngDens(iDens+1,j) = 0
            iAngDens(iDens+2,j) = 0
            iAngDens(iDens+5,j) = iAng(i1,j)
            iAngDens(iDens+8,j) = iAngFun(j)
            IF (iAngFun(j).ne.0) THEN
              iAngDens(iDens+ij,j) = 1
              ij = 2
            ENDIF
            iAngDens(iDens+3,j) = 2*iAngDens(iDens+1,j)
            iAngDens(iDens+4,j) = 2*iAngDens(iDens+2,j)
            iAngDens(iDens+6,j) = iAngFun(j)-iAngDens(iDens+1,j)
            iAngDens(iDens+7,j) = iAngFun(j)-iAngDens(iDens+2,j)
            CDens(iDens+7) = -dfloat(iAngDens(iDens+2,j))*
     &       (ROrb(i1,j)+ROrb(i2,j))+CDens(iDens+7)
            CDens(iDens+6) = -dfloat(iAngDens(iDens+1,j))*
     &       (ROrb(i1,j)+ROrb(i2,j))+CDens(iDens+6)
            CDens(iDens+3) = CDens(iDens+3)+
     &       dfloat(iAngDens(iDens+2,j))*ROrb(i1,j)*ROrb(i2,j)
            CDens(iDens+4) = CDens(iDens+4)+
     &       dfloat(iAngDens(iDens+1,j))*ROrb(i1,j)*ROrb(i2,j)
          ENDDO
          CDens(iDens+8) = One
          CDens(iDens+5) = CDens(iDens+6)*CDens(iDens+7)
          CDens(iDens+2) = CDens(iDens+4)*CDens(iDens+7)
          CDens(iDens+1) = CDens(iDens+3)*CDens(iDens+6)
          CDens(iDens)   = CDens(iDens+3)*CDens(iDens+4)
          DO i = 1,9
            CDens(iDens+i-1) = Four*CFun*CDens(iDens+i-1)
          ENDDO
          IF (i1.eq.i2) THEN
            DO i = 1,9
              CDens(iDens+i-1) = Half*CDens(iDens+i-1)
            ENDDO
          ENDIF
          iDens = iDens + 9
          NDens = NDens + 9
          ENDIF
         ENDIF
        ENDIF
        ENDIF
        ENDIF
c---------------------------------------------------------------------c
        ENDIF  ! Tol test
!10   CONTINUE
c--------------------------------------------------------------
      END DO

      END DO

      END
