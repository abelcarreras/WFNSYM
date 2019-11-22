
      SUBROUTINE aos_product(iAng,Alph,C,R,NAB,NN,
     &   iAngvA,AlphA,CA,RvA,iAngvB,AlphB,CB,RvB)
c*****************************************************c
c                                                     c
c       Computes the product of 2 AOs                 c
c       to build a set of gaussians                   c
c                                                     c
c     INPUT:                                          c
c                                                     c
c       R_        # center functions                  c
c       iAngv_    angular indices                     c
c       Alph_     exponents                           c
c       C_        # product coefficients              c
c       R_        function center                     c
c       NN        product leading dimension           c
c                                                     c
c                                                     c
c     OUTPUT:                                         c
c                                                     c
c       iAng  product angular indices                 c
c       Alph  product exponents                       c
c       C     product coefficients                    c
c       R     product positions                       c
c       NAB   product elements                        c
c                                                     c
c                                                     c
c       DCC: 05/09                                    c
c       DCC: 04/10  WFN-Sym update                    c
c                                                     c
c*****************************************************c
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER  :: iAngvA(3),iAngA(3),iAngvB(3),iAngB(3),
     &  iAng(NN,3),iAngFun(3),iAngX(3)
      REAL*8   :: Alph(NN),C(NN),R(NN,3),RX(3),
     &  RvA(3),RA(3),RvB(3),RB(3),RFun(3)
      SAVE Zero,Half,One,Two,Three,Four,Eight
      DATA Zero,Half,One,Two,Three,Four,Eight
     &    / 0d0, 0.5d0, 1d0, 2d0, 3d0, 4d0, 8d0 /

c Make copies of RvA,RvB and iAngvA,iAngvB 
	DO I = 1,3
	  RA(:)    = RvA(:)
	  RB(:)    = RvB(:)
	  iAngA(:) = iAngvA(:)
	  iAngB(:) = iAngvB(:)
	ENDDO

c Initialize some variables
      NAB = 0
      iAB = 1
	iChange = 0

      ijkA = iAngA(1) + iAngA(2) + iAngA(3)
      ijkB = iAngB(1) + iAngB(2) + iAngB(3)

c Squared distance between centers and angular part
        Dist = Zero
        DO j = 1,3
          Dist = Dist + (RA(j)-RB(j))**2
          iAngFun(j) = iAngA(j) + iAngB(j)
        ENDDO

c Product constant
        AlphFun = AlphA + AlphB
        CFun = DEXP(-Dist*AlphA*AlphB/AlphFun)
        CFun = CA * CB * CFun 
        DO j = 1,3
          RFun(j) = (AlphA*RA(j) + AlphB*RB(j))/AlphFun
        ENDDO

cc
c  S*S product
cc
        IF (ijkA.EQ.0 .AND. ijkB.EQ.0) THEN
          C(iAB)= CFun
          Alph(iAB) = AlphFun
          DO j=1,3
            iAng(iAB,j) = iAngFun(j)
            R(iAB,j) = RFun(j)
          ENDDO
          NAB = NAB + 1

cc
c  S*P product
cc
        ELSEIF ((ijkA+ijkB).EQ.1) THEN
          C(iAB) = Zero
          DO  j = 1,3
           C(iAB) = C(iAB) - CFun*
     &     (dfloat(iAngA(j)) * RA(j) +
     &     dfloat(iAngB(j)) * RB(j))
          ENDDO
          C(iAB+1) = CFun
          Alph(iAB)   = AlphFun
          Alph(iAB+1) = AlphFun
          DO j = 1,3
            iAng(iAB,j)   = 0
            iAng(iAB+1,j) = iAngFun(j)
            R(iAB,j)   = RFun(j)
            R(iAB+1,j) = RFun(j)
          ENDDO
          iAB = iAB + 2
          NAB = NAB + 2

cc          
c  P*P product
cc        
        ELSEIF (ijkA.EQ.1 .AND. ijkB.EQ.1) THEN
cc   ---Funcions S, P i D---
          DO i=1,4
            C(iAB+i-1) = Zero
            Alph(iAB+i-1) = AlphFun
            DO j=1,3
              R(iAB+i-1,j) = RFun(j)
            ENDDO
          ENDDO
          DO j=1,3
            iAng(iAB,j) = 0
            iAng(iAB+1,j) = iAngA(j)
            iAng(iAB+2,j) = iAngB(j)
            iAng(iAB+3,j) = iAngFun(j)
            DO i = 1,3
              C(iAB) = C(iAB) +
     &        CFun*dfloat(iAngA(j)) * RA(j)*
     &        dfloat(iAngB(i))*RB(i)
            ENDDO
            C(iAB+1) = C(iAB+1)-
     &      CFun*dfloat(iAngB(j))*RB(j)
            C(iAB+2) = C(iAB+2)-
     &      CFun*dfloat(iAngA(j))*RA(j)
          ENDDO
          C(iAB+3) = CFun
          iAB = iAB + 4
          NAB = NAB + 4

cc
c  S*D product
cc
cc   ---Funcions S, P i D---
        ELSEIF ((ijkA.EQ.0 .AND. ijkB.EQ.2).OR.
     &     (ijkA.EQ.2 .AND. ijkB.EQ.0)) THEN
cc Case (0 0 2)
          IF ((iAngFun(1).EQ.2).OR.(iAngFun(2).EQ.2).OR.
     &       (iAngFun(3).EQ.2)) THEN
            C(iAB)   = Zero
            C(iAB+1) = Zero
            DO i = 1,3
              iAng(iAB,i) = 0
              iAng(iAB+1,i) = iAngFun(i)/2
              iAng(iAB+2,i) = iAngFun(i)
              IF(iAngFun(i).EQ.2) THEN
	         C(iAB) = CFun*RA(i)**(iAngA(i))
     &                *RB(i)**(iAngB(i))+C(iAB)
               C(iAB+1) = -CFun*(RA(i)*
     &         dfloat(iAngA(i))+RB(i)*(dfloat(iAngB(i))))+
     &         C(iAB+1)
              ENDIF
              Alph(iAB+i-1) = AlphFun
              DO j = 1,3
                R(iAB+i-1,j) = RFun(j)
              ENDDO
            ENDDO
            C(iAB+2) = CFun
            iAB = iAB + 3
            NAB = NAB + 3
          ELSE
cc Case (1 1 0)
            DO i = 1,4
              DO j = 1,3
                R(iAB+i-1,j) = RFun(j)
                iAng(iAB+i-1,j) = 0
              ENDDO
              Alph(iAB+i-1) = AlphFun
            ENDDO
            DO i = 1,2
              DO j = i+1,3
                IF ((iAngFun(i).EQ.1).AND.(iAngFun(j).EQ.1)) THEN
                  iAng(iAB+1,i) = 1
                  iAng(iAB+2,j) = 1
                  C(iAB) = CFun*(RA(i)*
     &       dfloat(iAngA(i))+RB(i)*dfloat(iAngB(i)))*
     &       (RA(j)*dfloat(iAngA(j))+RB(j)*
     &       dfloat(iAngB(j)))
                  C(iAB+1) = -CFun*(RA(j)*
     &       dfloat(iAngA(j))+RB(j)*dfloat(iAngB(j)))
                  C(iAB+2) = -CFun*(RA(i)*
     &       dfloat(iAngA(i))+RB(i)*dfloat(iAngB(i)))
                ENDIF
              ENDDO
            ENDDO
            C(iAB+3) = CFun
            DO i = 1,3
              iAng(iAB+3,i) = iAngFun(i)
            ENDDO
            iAB = iAB + 4
            NAB = NAB + 4
          ENDIF

cc            
c  P*D product
cc          
cc   ---Funcions S, P, D i F---
        ELSEIF ((ijkA.EQ.1 .AND. ijkB.EQ.2).OR.
     &    (ijkA.EQ.2 .AND. ijkB.EQ.1)) THEN
cc Case (1 0 0) (0 1 1)
          IF ((iAngFun(1).EQ.1).AND.(iAngFun(2).EQ.1).
     &        AND.(iAngFun(3).EQ.1)) THEN
            DO j = 1,3
              DO i = 1,8
                R(iAB+i-1,j) = RFun(j)
              ENDDO
              iAng(iAB,j) = 0
              iAng(iAB+7,j) = 1
              iAng(iAB+1,j) = 0
              iAng(iAB+2,j) = 0
              iAng(iAB+3,j) = 0
              iAng(iAB+4,j) = 1
              iAng(iAB+5,j) = 1
              iAng(iAB+6,j) = 1
              iAng(iAB+j,j) = 1
              iAng(iAB+3+j,j) = 0
              C(iAB+3+j) = -RA(j)**(iAngA(j))*
     &         RB(j)**(iAngB(j))
            ENDDO
            C(iAB) = -RA(1)**(iAngA(1))*
     &       RB(1)**(iAngB(1))*RA(2)**(iAngA(2))*
     &       RA(3)**(iAngA(3))*RB(2)**(iAngB(2))*
     &       RB(3)**(iAngB(3))
            C(iAB+1) = RA(2)**(iAngA(2))*
     &       RA(3)**(iAngA(3))*RB(2)**(iAngB(2))*
     &       RB(3)**(iAngB(3))
            C(iAB+2) = RA(1)**(iAngA(1))*
     &       RA(3)**(iAngA(3))*RB(1)**(iAngB(1))*
     &       RB(3)**(iAngB(3))
            C(iAB+3) = RA(1)**(iAngA(1))*
     &       RA(2)**(iAngA(2))*RB(1)**(iAngB(1))*
     &       RB(2)**(iAngB(2))
            C(iAB+7) = One
            DO i = 1,8
              Alph(iAB+i-1) = AlphFun
              C(iAB+i-1) = CFun*C(iAB+i-1)
            ENDDO
            iAB = iAB + 8
            NAB = NAB + 8
          ENDIF
 
cc Case (1 0 0) (1 1 0) and (1 0 0) (0 2 0)
          IF ((iAngFun(1).EQ.2).OR.(iAngFun(2).EQ.2).OR.
     .       (iAngFun(3).EQ.2)) THEN
            DO j = 1,3
              DO i = 1,6
                iAng(iAB+i-1,j) = 0
                R(iAB+i-1,j) = RFun(j)
              ENDDO
              iAng(iAB+5,j) = iAngFun(j)
              iAng(iAB+3,j) = 1
              IF (iAngFun(j).EQ.2) THEN
                iAng(iAB+4,j) = 2
                iAng(iAB+1,j) = 1
                C(iAB+3) = -(RA(j)*
     .           dfloat(iAngA(j))+RB(j)*dfloat(iAngB(j)))
                C(iAB+2) = RA(j)**
     .           (iAngA(j))*RB(j)**(iAngB(j))
              ENDIF
              IF (iAngFun(j).EQ.1) THEN
                iAng(iAB+2,j) = 1
                C(iAB+4) = -RA(j)**
     .           (iAngA(j))*RB(j)**(iAngB(j))
              ENDIF
              IF (iAngFun(j).EQ.0) THEN
                iAng(iAB+3,j) = 0
              ENDIF
            ENDDO
            C(iAB+1) = C(iAB+3)*C(iAB+4)
            C(iAB) = C(iAB+2)*C(iAB+4)
            C(iAB+5) = One
            DO i = 1,6
              Alph(iAB+i-1) = AlphFun
              C(iAB+i-1) = CFun*C(iAB+i-1)
            ENDDO
            iAB=iAB+6
            NAB=NAB+6
          ENDIF
cc Case (1 0 0) (2 0 0)
          IF ((iAngFun(1).EQ.3).OR.(iAngFun(2).EQ.3).OR.
     .       (iAngFun(3).EQ.3)) THEN
            DO i=1,4
              DO j=1,3
                iAng(iAB+i-1,j)=(i-1)*iAngFun(j)/3
                R(iAB+i-1,j)=RFun(j)
              ENDDO
              Alph(iAB+i-1)=AlphFun
            ENDDO
            C(iAB)=-RA(1)**(iAngA(1))*
     .       RB(1)**(iAngB(1))*RA(2)**(iAngA(2))*
     .       RA(3)**(iAngA(3))*RB(2)**(iAngB(2))*
     .       RB(3)**(iAngB(3))
            C(iAB+1) = Zero
            C(iAB+2) = Zero
            DO i = 1,3
              C(iAB+1) = RA(i)*RB(i)*
     .         dfloat(iAngA(i)*iAngB(i))+C(iAB+1)
              C(iAB+2) = -(dfloat(iAngA(i))*RA(i)+
     .         dfloat(iAngB(i))*RB(i))+C(iAB+2)
              IF (iAngA(i).EQ.2) THEN
                C(iAB+1) = C(iAB+1)+RA(i)**2
              ENDIF
              IF (iAngB(i).EQ.2) THEN
                C(iAB+1) = C(iAB+1)+RB(i)**2
              ENDIF
            ENDDO
            C(iAB+3) = One
          DO i = 1,4
            C(iAB+i-1) = CFun*C(iAB+i-1)
          ENDDO
          iAB=iAB+4
          NAB=NAB+4
          ENDIF

cc
c  D*D product
cc
cc   ---S, P, D, F and G functions---
      ELSEIF ((ijkA.EQ.2) .AND. (ijkB.EQ.2)) THEN

cc Case (2 0 0) (2 0 0)
       IF ((iAngFun(1).EQ.4).OR.(iAngFun(2).EQ.4).OR.
     &     (iAngFun(3).EQ.4)) THEN
         DO i = 1,5
           DO j = 1,3
             iAng(iAB+i-1,j) = (i-1)*iAngFun(j)/4
             R(iAB+i-1,j) = RFun(j)
           ENDDO
           Alph(iAB+i-1) = AlphFun
         ENDDO
         C(iAB) = CFun*RA(1)**(iAngA(1))*
     &    RB(1)**(iAngB(1))*RA(2)**(iAngA(2))*
     &    RA(3)**(iAngA(3))*RB(2)**(iAngB(2))*
     &    RB(3)**(iAngB(3))
         DO i=1,3
           IF (iAngFun(i).EQ.4) THEN
             C(iAB+1) = -Two*CFun*(RA(i)*
     &        RB(i)**2+RB(i)*RA(i)**2)
             C(iAB+2) = CFun*(RA(i)**2+
     &        RB(i)**2+Four*RA(i)*RB(i))
             C(iAB+3) = -Two*CFun*(RA(i)+
     &        RB(i))
           ENDIF
         ENDDO
         C(iAB+4) = CFun
         iAB=iAB+5
         NAB=NAB+5

cc Case (1 1 0) (2 0 0)
	 ELSEIF ((iAngFun(1).EQ.3).OR.(iAngFun(2).EQ.3).OR.
     &         (iAngFun(3).EQ.3)) THEN
         DO j = 1,3
           DO i = 1,8
             R(iAB+i-1,j) = RFun(j)
           ENDDO
           iAng(iAB,j)   = 0
           iAng(iAB+1,j) = 0
           iAng(iAB+2,j) = 0
           iAng(iAB+4,j) = 0
           iAng(iAB+6,j) = 0
           iAng(iAB+7,j) = iAngFun(j)
           IF (iAngFun(j).EQ.3) THEN
             iAng(iAB+1,j) = 1
             iAng(iAB+4,j) = 2
             iAng(iAB+6,j) = 3
             C(iAB+2) = -(RA(j)**
     &         iAngA(j))*(RB(j)**iAngB(j))
             C(iAB+5) = -(dfloat(iAngA(j))*
     &        RA(j)+dfloat(iAngB(j))*RB(j))
           ENDIF
           IF (iAngFun(j).EQ.1) THEN
             iAng(iAB+2,j) = 1
             C(iAB+6) = -(RA(j)*
     &        dfloat(iAngA(j))+dfloat(iAngB(j))*RB(j))
           ENDIF
           iAng(iAB+3,j) = iAng(iAB+1,j)+iAng(iAB+2,j)
           iAng(iAB+5,j) = iAng(iAB+4,j)+iAng(iAB+2,j)
           IF (iAngA(j).EQ.2) THEN
             C(iAB+3)=(RA(j)**2+Two*RA(j)*RB(j))
           ENDIF
           IF (iAngB(j).EQ.2) THEN
             C(iAB+3)=(RB(j)**2+Two*RA(j)*RB(j))
           ENDIF
         ENDDO
         C(iAB) = RA(1)**(iAngA(1))*
     &    RB(1)**(iAngB(1))*RA(2)**(iAngA(2))*
     &    RA(3)**(iAngA(3))*RB(2)**(iAngB(2))*
     &    RB(3)**(iAngB(3))
         C(iAB+1) = C(iAB+3)*C(iAB+6)
         C(iAB+4) = C(iAB+5)*C(iAB+6)
         C(iAB+7) = One
         DO i = 1,8
           C(iAB+i-1) = CFun*C(iAB+i-1)
           Alph(iAB+i-1)=AlphFun
         ENDDO
         iAB = iAB + 8
         NAB = NAB + 8

cc Cases (1 1 0) (0 1 1) and (110)(002)
	 ELSEIF ((iAngFun(1).EQ.1).OR.(iAngFun(2).EQ.1).OR.
     &         (iAngFun(3).EQ.1)) THEN
         DO i = 1,12
           Alph(iAB+i-1) = AlphFun
           DO j = 1,3
             R(iAB+i-1,j) = RFun(j)
           ENDDO
         ENDDO

cc Only (110)(002) case
         IF (iAngA(1)*iAngB(1).EQ.0.AND.iAngA(2)*
     &     iAngB(2).EQ.0.AND.iAngA(3)*iAngB(3).EQ.0) THEN
           ij = 1
           C(iAB) = One
           DO j = 1,3
             iAng(iAB,j)   = 0
             iAng(iAB+1,j) = 0
             iAng(iAB+2,j) = 0
             iAng(iAB+3,j) = iAngFun(j)-1
             iAng(iAB+6,j) = 1
             iAng(iAB+7,j) = 0
             iAng(iAB+10,j) = 1
             iAng(iAB+11,j) = iAngFun(j)
             IF (iAngFun(j).EQ.1) THEN
               iAng(iAB+ij,j) = 1
               C(iAB+7+ij) = -One*
     &         RA(j)**(iAngA(j))*RB(j)**(iAngB(j))
               ij = 2
             ENDIF
             iAng(iAB+8,j) = iAngFun(j)-iAng(iAB+1,j)
             iAng(iAB+9,j) = iAngFun(j)-iAng(iAB+2,j)
             IF (iAngFun(j).EQ.2) THEN
               iAng(iAB+6,j) = 0
               iAng(iAB+7,j) = 2
               C(iAB+6) = RA(j)**
     &          (iAngA(j))*RB(j)**(iAngB(j))
               C(iAB+10) = -(RA(j)*
     &          dfloat(iAngA(j))+dfloat(iAngB(j))*RB(j))
             ENDIF
             iAng(iAB+4,j) = 1-iAng(iAB+1,j)
             iAng(iAB+5,j) = 1-iAng(iAB+2,j)
             C(iAB) = C(iAB)*RA(j)**
     &        (iAngA(j))*RB(j)**(iAngB(j))
           ENDDO
           C(iAB+7) = C(iAB+8)*C(iAB+9)
           C(iAB+4) = C(iAB+8)*C(iAB+10)
           C(iAB+5) = C(iAB+9)*C(iAB+10)
           C(iAB+3) = C(iAB+5)*C(iAB+8)
           C(iAB+2) = C(iAB+6)*C(iAB+8)
           C(iAB+1) = C(iAB+6)*C(iAB+9)
           C(iAB+11) = One
           DO i = 1,12
             C(iAB+i-1) = CFun*C(iAB+i-1)
           ENDDO    
           iAB = iAB + 12
           NAB = NAB + 12

         ELSE
cc Only (110)(011) case
          DO j = 1,3
            iAng(iAB,j)   = 0
            iAng(iAB+1,j) = 0
            iAng(iAB+2,j) = 0
            iAng(iAB+3,j) = 0
            iAng(iAB+4,j) = iAngB(j)
            iAng(iAB+5,j) = iAngA(j)
            iAng(iAB+7,j) = 0
            iAng(iAB+8,j) = iAngFun(j)
            iAng(iAB+9,j) = iAngFun(j)
            iAng(iAB+10,j) = 1
            iAng(iAB+11,j) = iAngFun(j)
            IF(iAngA(j).EQ.0) THEN
              iAng(iAB+2,j) = 1
              iAng(iAB+8,j) = 0
              C(iAB+8) = -RB(j)
            ENDIF
            IF(iAngB(j).EQ.0) THEN
              iAng(iAB+1,j) = 1
              iAng(iAB+9,j) = 0
              C(iAB+9) = -RA(j)
            ENDIF
            IF(iAngFun(j).EQ.2) THEN
              iAng(iAB+3,j) = 1
              iAng(iAB+7,j) = 2
              C(iAB+6) = RA(j)*RB(j)
              C(iAB+10) = -(RA(j)+RB(j))
            ENDIF
            iAng(iAB+6,j) = 1-iAng(iAB+3,j)
          ENDDO
          C(iAB+4) = C(iAB+10)*C(iAB+9)
          C(iAB+5) = C(iAB+10)*C(iAB+8)
          C(iAB+7) = C(iAB+8)*C(iAB+9)
          C(iAB+1) = C(iAB+6)*C(iAB+8)
          C(iAB+2) = C(iAB+6)*C(iAB+9)
          C(iAB+3) = C(iAB+5)*C(iAB+9)
          C(iAB) = C(iAB+6)*C(iAB+7)
          C(iAB+11) = One
          DO i = 1,12
            C(iAB+i-1) = CFun*C(iAB+i-1)
          ENDDO
          iAB = iAB + 12
          NAB = NAB + 12
         ENDIF

cc Cases (200)(020) and (110)(110)
        ELSEIF ((iAngFun(1).EQ.2).OR.(iAngFun(2).EQ.2).OR.
     &          (iAngFun(3).EQ.2)) THEN

          DO i = 1,9
            DO j = 1,3
              R(iAB+i-1,j) = RFun(j)
            ENDDO
            Alph(iAB+i-1) = AlphFun
          ENDDO

cc Only (200)(020) case
          IF ((iAngA(1).EQ.2).OR.(iAngA(2).EQ.2).OR.
     &        (iAngA(3).EQ.2)) THEN
            C(iAB+7) = Zero
            C(iAB+6) = Zero
            DO j = 1,3
              iAng(iAB,j) = 0
              iAng(iAB+1,j) = 0
              iAng(iAB+2,j) = 0
              iAng(iAB+3,j) = iAngA(j)
              iAng(iAB+4,j) = iAngB(j)
              iAng(iAB+5,j) = iAngFun(j)/2
              iAng(iAB+8,j) = iAngFun(j)
              IF (iAngA(j).EQ.2) THEN
               iAng(iAB+1,j) = 1
               C(iAB+4) = RA(j)**2
              ENDIF
              IF (iAngB(j).EQ.2) THEN
               iAng(iAB+2,j) = 1
               C(iAB+3) = RB(j)**2
              ENDIF
              iAng(iAB+6,j) = iAng(iAB+5,j)+iAng(iAB+2,j)
              iAng(iAB+7,j) = iAng(iAB+5,j)+iAng(iAB+1,j)
              C(iAB+7) = C(iAB+7)-RB(j)*dfloat(iAngB(j))
              C(iAB+6) = C(iAB+6)-RA(j)*dfloat(iAngA(j))
            ENDDO        
            C(iAB+8) = One
            C(iAB+5) = C(iAB+6)*C(iAB+7)
            C(iAB+1) = C(iAB+6)*C(iAB+3)
            C(iAB+2) = C(iAB+4)*C(iAB+7)
            C(iAB) = C(iAB+3)*C(iAB+4)
            DO i = 1,9
              C(iAB+i-1) = CFun*C(iAB+i-1)
            ENDDO
            iAB = iAB + 9
            NAB = NAB + 9

cc Only (110)(110) case
          ELSE
           ij = 1
           C(iAB+7) = Zero
           C(iAB+6) = Zero
           C(iAB+4) = Zero
           C(iAB+3) = Zero
           DO j = 1,3
            iAng(iAB,j)   = 0
            iAng(iAB+1,j) = 0
            iAng(iAB+2,j) = 0
            iAng(iAB+5,j) = iAngA(j)
            iAng(iAB+8,j) = iAngFun(j)
            IF (iAngFun(j).ne.0) THEN
              iAng(iAB+ij,j) = 1
              ij = 2
            ENDIF
            iAng(iAB+3,j) = 2*iAng(iAB+1,j)
            iAng(iAB+4,j) = 2*iAng(iAB+2,j)
            iAng(iAB+6,j) = iAngFun(j)-iAng(iAB+1,j)
            iAng(iAB+7,j) = iAngFun(j)-iAng(iAB+2,j)
            C(iAB+7) = -dfloat(iAng(iAB+2,j))*
     &       (RA(j)+RB(j))+C(iAB+7)
            C(iAB+6) = -dfloat(iAng(iAB+1,j))*
     &       (RA(j)+RB(j))+C(iAB+6)
            C(iAB+3) = C(iAB+3)+
     &       dfloat(iAng(iAB+2,j))*RA(j)*RB(j)
            C(iAB+4) = C(iAB+4)+
     &       dfloat(iAng(iAB+1,j))*RA(j)*RB(j)
           ENDDO
           C(iAB+8) = One
           C(iAB+5) = C(iAB+6)*C(iAB+7)
           C(iAB+2) = C(iAB+4)*C(iAB+7)
           C(iAB+1) = C(iAB+3)*C(iAB+6)
           C(iAB)   = C(iAB+3)*C(iAB+4)
           DO i = 1,9
             C(iAB+i-1) = CFun*C(iAB+i-1)
           ENDDO
           iAB = iAB + 9
           NAB = NAB + 9
          ENDIF
        ENDIF

c*---------------*c
c*---------------*c
c*  F FUNCTIONS  *c
c*---------------*c
c*---------------*c

cc
c  F*F product
cc
cc   ---S, P, D, F, G, H and I functions---
      ELSEIF ((ijkA.EQ.3) .AND. (ijkB.EQ.3)) THEN

c Order terms: (3,0,0);(2,1,0);(1,1,1)
      IF((iAngB(1).EQ.3 .OR. iAngB(2).EQ.3 .OR. 
     &   iAngB(3).EQ.3) .AND. (iAngA(1).LT.3 .AND.
     &   iAngA(2).LT.3 .AND. iAngA(3).LT.3)) iChange = 1
      IF((iAngB(1).EQ.2 .OR. iAngB(2).EQ.2 .OR.   
     &   iAngB(3).EQ.2) .AND. (iAngA(1).LT.2 .AND.
     &   iAngA(2).LT.2 .AND. iAngA(3).LT.2)) iChange = 1
	IF(iChange.EQ.1) THEN
        iAngX(1:3) = iAngA(1:3)
        iAngA(1:3) = iAngB(1:3)
        iAngB(1:3) = iAngX(1:3)
        RX(1:3) = RA(1:3)
        RA(1:3) = RB(1:3)
        RB(1:3) = RX(1:3)
	ENDIF

c Case [6,0,0]
	IF(iAngFun(1).EQ.6 .OR. iAngFun(2).EQ.6 .OR.
     &   iAngFun(3).EQ.6) THEN

       DO i = 0,6
         DO j = 1,3
           iAng(iAB+i,j) = 0
           R(iAB+i,j) = RFun(j)
         ENDDO
         Alph(iAB+i) = AlphFun
       ENDDO
	 DO i = 1,3
	   IF(iAngFun(i).EQ.6) j1 = i
	 ENDDO
       iAng(iAB  ,j1) = 6
       iAng(iAB+1,j1) = 5
       iAng(iAB+2,j1) = 4
       iAng(iAB+3,j1) = 3
       iAng(iAB+4,j1) = 2
       iAng(iAB+5,j1) = 1
       iAng(iAB+6,j1) = 0

	 C(iAB  ) = CFun
       C(iAB+1) = -CFun*Three*(RA(j1) + RB(j1))
       C(iAB+2) = CFun*Three*
     &  (RA(j1)**2 + Three*RA(j1)*RB(j1) + RB(j1)**2)
       C(iAB+3) = -CFun*( RA(j1)**3 + RB(j1)**3 +
     &  9d0*RA(j1)*RB(j1)*(RA(j1)+RB(j1)) )
       C(iAB+4) = C(iAB+2) * RA(j1) * RB(j1)
       C(iAB+5) = C(iAB+1) * (RA(j1)*RB(j1))**2
       C(iAB+6) = CFun* ( RA(j1)*RB(j1) )**3

       iAB = iAB + 7
       NAB = NAB + 7

c Case [5,1,0]: (3,0,0)(2,1,0)
      ELSEIF(iAngFun(1).EQ.5 .OR. iAngFun(2).EQ.5 .OR.
     &       iAngFun(3).EQ.5) THEN


      DO i = 0,11
        iAng(iAB+i,:) = 0
        R(iAB+i,:)    = RFun(:)
        Alph(iAB+i)   = AlphFun
      ENDDO
	DO i = 1,3
	 IF(iAngFun(i).EQ.5) j1 = i
	 IF(iAngFun(i).EQ.1) j2 = i
	ENDDO
      iAng(iAB  ,j1) = 5
      iAng(iAB  ,j2) = 1
      iAng(iAB+1,j1) = 5
      iAng(iAB+2,j1) = 4
      iAng(iAB+2,j2) = 1
      iAng(iAB+3,j1) = 4
      iAng(iAB+4,j1) = 3
      iAng(iAB+4,j2) = 1
      iAng(iAB+5,j1) = 3
      iAng(iAB+6,j1) = 2
      iAng(iAB+6,j2) = 1
      iAng(iAB+7,j1) = 2
      iAng(iAB+8,j1) = 1
      iAng(iAB+8,j2) = 1
      iAng(iAB+9,j1) = 1
      iAng(iAB+10,j2)= 1

      C(iAB  ) = One
      C(iAB+1) = -RB(j2)
      C(iAB+2) = -(Three*RA(j1)+Two*RB(j1))
      C(iAB+3) = - C(iAB+2) * RB(j2)
      C(iAB+4) = Three*RA(j1)**2 + 6d0*RA(j1)*RB(j1) + RB(j1)**2
      C(iAB+5) = - C(iAB+4) * RB(j2) 
      C(iAB+6) = -(RA(j1)**2+6d0*RA(j1)*RB(j1)+Three*RB(j1)**2)*
     &           RA(j1)
      C(iAB+7) = -C(iAB+6)*RB(j2) 
      C(iAB+8) = (Two*RA(j1)+Three*RB(j1))*RB(j1)*RA(j1)**2
      C(iAB+9) = -C(iAB+8)*RB(j2)
      C(iAB+10) = -RA(j1)*(RA(j1)*RB(j1))**2
      C(iAB+11) = -C(iAB+10)*RB(j2)
	DO i = 0,11
	  C(iAB+i) = CFun * C(iAB+i)
	ENDDO 
	
      iAB = iAB + 12
      NAB = NAB + 12

c Cases [4,n,m]
      ELSEIF(iAngFun(1).EQ.4 .OR. iAngFun(2).EQ.4 .OR.
     &       iAngFun(3).EQ.4) THEN


c cases [4,2,0]: (3,0,0)(1,2,0); (2,1,0)(2,1,0)
	 IF(iAngFun(1).EQ.2 .OR. iAngFun(2).EQ.2 .OR.
     &       iAngFun(3).EQ.2) THEN


       DO i = 0,14
	  DO j = 1,3
         iAng(iAB+i,j) = 0
         R(iAB+i,j)    = RFun(j)
	  ENDDO
        Alph(iAB+i)   = AlphFun
       ENDDO
	 DO i = 1,3
         IF (iAngFun(i).EQ.4) j1 = i
         IF (iAngFun(i).EQ.2) j2 = i
	 ENDDO

       iAng(iAB  ,j1) = 4
       iAng(iAB  ,j2) = 2
       iAng(iAB+1,j1) = 4
       iAng(iAB+1,j2) = 1
       iAng(iAB+2,j1) = 3
       iAng(iAB+2,j2) = 2
       iAng(iAB+3,j1) = 4
       iAng(iAB+4,j1) = 3
       iAng(iAB+4,j2) = 1
       iAng(iAB+5,j1) = 2
       iAng(iAB+5,j2) = 2
       iAng(iAB+6,j1) = 3
       iAng(iAB+7,j1) = 2
       iAng(iAB+7,j2) = 1
       iAng(iAB+8,j1) = 1
       iAng(iAB+8,j2) = 2
       iAng(iAB+9,j1) = 2
       iAng(iAB+10,j1) = 1
       iAng(iAB+10,j2) = 1
       iAng(iAB+11,j2) = 2
       iAng(iAB+12,j1) = 1
       iAng(iAB+13,j2) = 1

       C(iAB  ) = One
       C(iAB+1) = -(iAngA(j2)*RA(j2)+iAngB(j2)*RB(j2)) 
       C(iAB+2) = -(iAngA(j1)*RA(j1)+iAngB(j1)*RB(j1)) 
	 IF(iAngA(j1).EQ.3) THEN
         C(iAB+3) = RB(j2)**2
         C(iAB+5) = Three*RA(j1)*(RA(j1)+RB(j1))
         C(iAB+8) = -(RA(j1)+Three*RB(j1))*RA(j1)**2
         C(iAB+11)= RB(j1)*RA(j1)**3 
	 ELSE
         C(iAB+3) = RA(j2)*RB(j2)
         C(iAB+5) = RA(j1)**2+Four*RA(j1)*RB(j1)+RB(j1)**2
         C(iAB+8) = -Two*(RA(j1)+RB(j1))*RA(j1)*RB(j1)
         C(iAB+11)= (RA(j1)*RB(j1))**2 
	 ENDIF
       C(iAB+4) = C(iAB+1) * C(iAB+2)
       C(iAB+6) = C(iAB+2) * C(iAB+3)
       C(iAB+7) = C(iAB+1) * C(iAB+5)
       C(iAB+9) = C(iAB+3) * C(iAB+5)
       C(iAB+10)= C(iAB+1) * C(iAB+8)
       C(iAB+12)= C(iAB+3) * C(iAB+8)
       C(iAB+13)= C(iAB+1) * C(iAB+11)
       C(iAB+14)= C(iAB+3) * C(iAB+11)

	 DO i =0,14
	   C(iAB+i) = CFun*C(iAB+i)
	 ENDDO
       
      iAB = iAB + 15
      NAB = NAB + 15

c cases [4,1,1]: (3,0,0)(1,1,1); (2,1,0)(2,0,1)
	 ELSE

       DO i = 0,19
         DO j = 1,3
           iAng(iAB+i,j) = 0
           R(iAB+i,j) = RFun(j)
         ENDDO
         Alph(iAB+i) = AlphFun
       ENDDO
       DO i = 1,3
         IF (iAngFun(i).EQ.4) j1 = i
       ENDDO
	 IF(iAngA(j1).EQ.3) THEN
	   j2 = j1 + 1
	   IF(j2.EQ.4) j2 = 1
	 ELSE
	   DO i = 1,3
	     IF(iAngA(i).EQ.1) j2 = i
	   ENDDO
	 ENDIF
	 j3 = 6 - (j1+j2)
	 iAng(iAB  ,j1) = 4
	 iAng(iAB  ,j2) = 1
	 iAng(iAB  ,j3) = 1
	 iAng(iAB+1,j1) = 4
	 iAng(iAB+1,j2) = 1
       iAng(iAB+2,j1) = 4
       iAng(iAB+2,j3) = 1
       iAng(iAB+3,j1) = 3
       iAng(iAB+3,j2) = 1
       iAng(iAB+3,j3) = 1
       iAng(iAB+4,j1) = 4
       iAng(iAB+5,j1) = 3
       iAng(iAB+5,j2) = 1
       iAng(iAB+6,j1) = 3
       iAng(iAB+6,j3) = 1
       iAng(iAB+7,j1) = 2
       iAng(iAB+7,j2) = 1
       iAng(iAB+7,j3) = 1
       iAng(iAB+8,j1) = 3
       iAng(iAB+9,j1) = 2
       iAng(iAB+9,j2) = 1
       iAng(iAB+10,j1) = 2
       iAng(iAB+10,j3) = 1
       iAng(iAB+11,j1) = 1
       iAng(iAB+11,j2) = 1
       iAng(iAB+11,j3) = 1
       iAng(iAB+12,j1) = 2
       iAng(iAB+13,j1) = 1
       iAng(iAB+13,j2) = 1
       iAng(iAB+14,j1) = 1
       iAng(iAB+14,j3) = 1
       iAng(iAB+15,j2) = 1
       iAng(iAB+15,j3) = 1
       iAng(iAB+16,j1) = 1
       iAng(iAB+17,j2) = 1
       iAng(iAB+18,j3) = 1

	 C(iAB  ) = One
	 C(iAB+1) = -RB(j3)
	 C(iAB+2) = -(iAngA(j2)*RA(j2) + iAngB(j2)*RB(j2))
       C(iAB+3) = -(iAngA(j1)*RA(j1) + iAngB(j1)*RB(j1))

       C(iAB+4) = C(iAB+1) * C(iAB+2)
       C(iAB+5) = C(iAB+1) * C(iAB+3)
       C(iAB+6) = C(iAB+2) * C(iAB+3)
	 IF(iAngA(j1).EQ.3) THEN
         C(iAB+7) = Three*RA(j1)*(RA(j1)+RB(j1))
         C(iAB+11)= -(RA(j1)+Three*RB(j1))*RA(j1)**2
         C(iAB+15)= RB(j1)*RA(j1)**3
	 ELSE
         C(iAB+7) = RA(j1)**2+Four*RA(j1)*RB(j1)+RB(j1)**2
         C(iAB+11)= -Two*(RA(j1)+RB(j1))*RA(j1)*RB(j1)
         C(iAB+15)= (RA(j1)*RB(j1))**2
	 ENDIF

       C(iAB+8) = C(iAB+1) * C(iAB+6)
       C(iAB+9) = C(iAB+1) * C(iAB+7)
       C(iAB+10) = C(iAB+2) * C(iAB+7)

       C(iAB+12) = C(iAB+1) * C(iAB+10)
       C(iAB+13) = C(iAB+1) * C(iAB+11)
       C(iAB+14) = C(iAB+2) * C(iAB+11)

       C(iAB+16) = C(iAB+1) * C(iAB+14)
       C(iAB+17) = C(iAB+1) * C(iAB+15)
       C(iAB+18) = C(iAB+2) * C(iAB+15)

       C(iAB+19) = C(iAB+1) * C(iAB+18)

       DO i =0,19
         C(iAB+i) = CFun*C(iAB+i)
       ENDDO

      iAB = iAB + 20
      NAB = NAB + 20

	 ENDIF


c Cases [3,n,m]
      ELSEIF(iAngFun(1).EQ.3 .OR. iAngFun(2).EQ.3 .OR.
     &       iAngFun(3).EQ.3) THEN

c cases [3,2,1]: (3,0,0)(0,2,1);(2,1,0)(1,0,2)=(2,1,0)(0,2,1);(2,1,0)(1,1,1)
       IF(iAngFun(1).EQ.2 .OR. iAngFun(2).EQ.2 .OR.
     &       iAngFun(3).EQ.2) THEN

      DO i = 0,23
        DO j = 1,3
          iAng(iAB+i,j) = 0
          R(iAB+i,j) = RFun(j)
        ENDDO
        Alph(iAB+i) = AlphFun
      ENDDO
	DO i = 1,3
	  IF(iAngFun(i).EQ.3) j1 = i
	  IF(iAngFun(i).EQ.2) j2 = i
	  IF(iAngFun(i).EQ.1) j3 = i
	ENDDO
c change (1,2,0)(2,0,1) by (2,0,1)(1,2,0)
      IF(iAngA(j1).EQ.1 .AND. iAngB(j1).EQ.2) THEN
        iAngX(1:3) = iAngA(1:3)
        iAngA(1:3) = iAngB(1:3)
        iAngB(1:3) = iAngX(1:3)
        RX(1:3) = RA(1:3)
        RA(1:3) = RB(1:3)
        RB(1:3) = RX(1:3)
      ENDIF

	iAng(iAB  ,j1) = 3
	iAng(iAB  ,j2) = 2
	iAng(iAB  ,j3) = 1
      iAng(iAB+1,j1) = 3
      iAng(iAB+1,j2) = 2
      iAng(iAB+2,j1) = 3
      iAng(iAB+2,j2) = 1
      iAng(iAB+2,j3) = 1
      iAng(iAB+3,j1) = 2
      iAng(iAB+3,j2) = 2
      iAng(iAB+3,j3) = 1
      iAng(iAB+4,j1) = 3
      iAng(iAB+4,j2) = 1
      iAng(iAB+5,j1) = 3
      iAng(iAB+5,j3) = 1
      iAng(iAB+6,j1) = 2
      iAng(iAB+6,j2) = 2
      iAng(iAB+7,j1) = 2
      iAng(iAB+7,j2) = 1
      iAng(iAB+7,j3) = 1
      iAng(iAB+8,j1) = 1
      iAng(iAB+8,j2) = 2
      iAng(iAB+8,j3) = 1
      iAng(iAB+9,j1) = 3
      iAng(iAB+10,j1) = 2
      iAng(iAB+10,j2) = 1
      iAng(iAB+11,j1) = 2
      iAng(iAB+11,j3) = 1
      iAng(iAB+12,j1) = 1
      iAng(iAB+12,j2) = 2
      iAng(iAB+13,j2) = 2
      iAng(iAB+13,j3) = 1
      iAng(iAB+14,j1) = 1
      iAng(iAB+14,j2) = 1
      iAng(iAB+14,j3) = 1
      iAng(iAB+15,j1) = 2
      iAng(iAB+16,j2) = 2
      iAng(iAB+17,j1) = 1
      iAng(iAB+17,j2) = 1
      iAng(iAB+18,j1) = 1
      iAng(iAB+18,j3) = 1
      iAng(iAB+19,j2) = 1
      iAng(iAB+19,j3) = 1
      iAng(iAB+20,j1) = 1
      iAng(iAB+21,j2) = 1
      iAng(iAB+22,j3) = 1

      C(iAB  ) = One
      C(iAB+1) = -(iAngA(j3)*RA(j3)+iAngB(j3)*RB(j3))
      C(iAB+2) = -(iAngA(j2)*RA(j2)+iAngB(j2)*RB(j2))
      C(iAB+3) = -(iAngA(j1)*RA(j1)+iAngB(j1)*RB(j1))
      C(iAB+4) = C(iAB+1) * C(iAB+2)
      C(iAB+6) = C(iAB+1) * C(iAB+3)
      C(iAB+7) = C(iAB+2) * C(iAB+3)
      IF(iAngA(j1).EQ.3) THEN
        C(iAB+5) = RB(j2)**2
        C(iAB+8) = Three*RA(j1)**2
        C(iAB+13)= -RA(j1)**3
	ELSEIF(iAngA(j1).EQ.2 .AND. iAngA(j2).EQ.1) THEN
        C(iAB+5) = RA(j2)*RB(j2)
	  C(iAB+8) = RA(j1)*(RA(j1)+Two*RB(j1))
        C(iAB+13)= -RB(j1)*RA(j1)**2
	ELSEIF(iAngA(j1).EQ.2 .AND. iAngA(j2).EQ.0) THEN
        C(iAB+5) = RB(j2)**2
	  C(iAB+8) = RA(j1)*(RA(j1)+Two*RB(j1))
        C(iAB+13)= -RB(j1)*RA(j1)**2
	ELSE
	  write(*,*)"ERROR. Unexpected term in F*F product"
	  STOP
	ENDIF
      C(iAB+ 9) = C(iAB+1) * C(iAB+5)
      C(iAB+10) = C(iAB+1) * C(iAB+7)
      C(iAB+11) = C(iAB+3) * C(iAB+5)
      C(iAB+12) = C(iAB+1) * C(iAB+8)
      C(iAB+14) = C(iAB+2) * C(iAB+8)
      C(iAB+15) = C(iAB+5) * C(iAB+6)
      C(iAB+16) = C(iAB+1) * C(iAB+13)
      C(iAB+17) = C(iAB+1) * C(iAB+14)
      C(iAB+18) = C(iAB+5) * C(iAB+8) 
      C(iAB+19) = C(iAB+2) * C(iAB+13)
      C(iAB+20) = C(iAB+1) * C(iAB+18)
      C(iAB+21) = C(iAB+1) * C(iAB+19)
      C(iAB+22) = C(iAB+5) * C(iAB+13)
      C(iAB+23) = C(iAB+1) * C(iAB+22)
	DO i = 0,23
	  C(iAB+i) = CFun * C(iAB+i)
	ENDDO

      iAB = iAB + 24
      NAB = NAB + 24

c Cases [3,3,0]: (3,0,0)(0,3,0);(2,1,0)(1,2,0)
       ELSE

	  DO i = 1,3
	    IF(iAngA(i).GE.2) j1 = i
	    IF(iAngB(i).GE.2) j2 = i
	  ENDDO    
        DO i = 0,15
          DO j = 1,3
            iAng(iAB+i,j) = 0
            R(iAB+i,j) = RFun(j)
          ENDDO
          Alph(iAB+i) = AlphFun
        ENDDO
        iAng(iAB  ,j1) = 3
        iAng(iAB  ,j2) = 3
        iAng(iAB+1,j1) = 3
        iAng(iAB+1,j2) = 2
        iAng(iAB+2,j1) = 2
        iAng(iAB+2,j2) = 3
        iAng(iAB+3,j1) = 3
        iAng(iAB+3,j2) = 1
        iAng(iAB+4,j1) = 1
        iAng(iAB+4,j2) = 3
        iAng(iAB+5,j1) = 2
        iAng(iAB+5,j2) = 2
        iAng(iAB+6,j1) = 3
        iAng(iAB+7,j2) = 3
        iAng(iAB+8,j1) = 2
        iAng(iAB+8,j2) = 1
        iAng(iAB+9,j1) = 1
        iAng(iAB+9,j2) = 2
        iAng(iAB+10,j1) = 2
        iAng(iAB+11,j2) = 2
        iAng(iAB+12,j1) = 1
        iAng(iAB+12,j2) = 1
        iAng(iAB+13,j1) = 1
        iAng(iAB+14,j2) = 1

        C(iAB  ) = One
        C(iAB+1) = -(iAngA(j2)*RA(j2)+iAngB(j2)*RB(j2))
        C(iAB+2) = -(iAngA(j1)*RA(j1)+iAngB(j1)*RB(j1))
	  IF(iAngA(j1).EQ.3) THEN
         C(iAB+3) = Three * RB(j2)**2
         C(iAB+4) = Three * RA(j1)**2
         C(iAB+6) = -RB(j2)**3
         C(iAB+7) = -RA(j1)**3
	  ELSE
         C(iAB+3) = (Two*RA(j2)+RB(j2))*RB(j2)
         C(iAB+4) = (RA(j1)+Two*RB(j1))*RA(j1)
         C(iAB+6) = -RA(j2)*RB(j2)**2
         C(iAB+7) = -RB(j1)*RA(j1)**2
	  ENDIF
        C(iAB+5) = C(iAB+1) * C(iAB+2)
        C(iAB+8) = C(iAB+2) * C(iAB+3)
        C(iAB+9) = C(iAB+1) * C(iAB+4)
        C(iAB+10) = C(iAB+2) * C(iAB+6)
        C(iAB+11) = C(iAB+1) * C(iAB+7)
        C(iAB+12) = C(iAB+3) * C(iAB+4)
        C(iAB+13) = C(iAB+4) * C(iAB+6)
        C(iAB+14) = C(iAB+3) * C(iAB+7)
        C(iAB+15) = C(iAB+6) * C(iAB+7)

	  DO i = 0,15
	    C(iAB+i) = CFun * C(iAB+i)
	  ENDDO

        iAB = iAB + 16
        NAB = NAB + 16

	ENDIF

c Cases [2,2,2]: (2,1,0)(0,1,2); (1,1,1)(1,1,1)
      ELSEIF(iAngFun(1).EQ.2 .AND. iAngFun(2).EQ.2 .AND.
     &       iAngFun(3).EQ.2) THEN

      DO i = 0,26
        DO j = 1,3
          iAng(iAB+i,j) = 0
          R(iAB+i,j) = RFun(j)
        ENDDO
        Alph(iAB+i) = AlphFun
      ENDDO
	IF(iAngA(1).EQ.1 .AND. iAngA(2).EQ.1) THEN
	  j1 = 1
	  j2 = 2
	ELSE
        DO i = 1,3
          IF(iAngA(i).EQ.2) j1 = i
          IF(iAngB(i).EQ.2) j2 = i
        ENDDO
	ENDIF
	j3 = 6 - (j1+j2)

      iAng(iAB  ,j1) = 2
      iAng(iAB  ,j2) = 2
      iAng(iAB  ,j3) = 2
      iAng(iAB+1,j1) = 2
      iAng(iAB+1,j2) = 2
      iAng(iAB+1,j3) = 1
      iAng(iAB+2,j1) = 2
      iAng(iAB+2,j2) = 1
      iAng(iAB+2,j3) = 2
      iAng(iAB+3,j1) = 1
      iAng(iAB+3,j2) = 2
      iAng(iAB+3,j3) = 2
      iAng(iAB+4,j1) = 2
      iAng(iAB+4,j2) = 2
      iAng(iAB+5,j1) = 2
      iAng(iAB+5,j3) = 2
      iAng(iAB+6,j2) = 2
      iAng(iAB+6,j3) = 2
      iAng(iAB+7,j1) = 2
      iAng(iAB+7,j2) = 1
      iAng(iAB+7,j3) = 1
      iAng(iAB+8,j1) = 1
      iAng(iAB+8,j2) = 2
      iAng(iAB+8,j3) = 1
      iAng(iAB+9,j1) = 1
      iAng(iAB+9,j2) = 1
      iAng(iAB+9,j3) = 2
      iAng(iAB+10,j1) = 2
      iAng(iAB+10,j2) = 1
      iAng(iAB+11,j1) = 2
      iAng(iAB+11,j3) = 1
      iAng(iAB+12,j1) = 1
      iAng(iAB+12,j2) = 2
      iAng(iAB+13,j2) = 2
      iAng(iAB+13,j3) = 1
      iAng(iAB+14,j1) = 1
      iAng(iAB+14,j3) = 2
      iAng(iAB+15,j2) = 1
      iAng(iAB+15,j3) = 2
      iAng(iAB+16,j1) = 1
      iAng(iAB+16,j2) = 1
      iAng(iAB+16,j3) = 1
      iAng(iAB+17,j1) = 2
      iAng(iAB+18,j2) = 2
      iAng(iAB+19,j3) = 2
      iAng(iAB+20,j1) = 1
      iAng(iAB+20,j2) = 1
      iAng(iAB+21,j1) = 1
      iAng(iAB+21,j3) = 1
      iAng(iAB+22,j2) = 1
      iAng(iAB+22,j3) = 1
      iAng(iAB+23,j1) = 1
      iAng(iAB+24,j2) = 1
      iAng(iAB+25,j3) = 1

      C(iAB  ) = One
      C(iAB+1) = -(RA(j3)+RB(j3))
      C(iAB+2) = -(iAngA(j2)*RA(j2)+iAngB(j2)*RB(j2))
      C(iAB+3) = -(iAngA(j1)*RA(j1)+iAngB(j1)*RB(j1))
      C(iAB+4) = RA(j3)*RB(j3)
	IF(iAngA(j1).EQ.2) THEN
	  C(iAB+5) = RB(j2)**2
	  C(iAB+6) = RA(j1)**2
	ELSE
	  C(iAB+5) = RA(j2)*RB(j2)
	  C(iAB+6) = RA(j1)*RB(j1)
	ENDIF
      C(iAB+7) = C(iAB+1) * C(iAB+2)
      C(iAB+8) = C(iAB+1) * C(iAB+3)
      C(iAB+9) = C(iAB+2) * C(iAB+3)
      C(iAB+10) = C(iAB+2) * C(iAB+4)
      C(iAB+11) = C(iAB+1) * C(iAB+5)
      C(iAB+12) = C(iAB+3) * C(iAB+4)
      C(iAB+13) = C(iAB+1) * C(iAB+6)
      C(iAB+14) = C(iAB+3) * C(iAB+5)
      C(iAB+15) = C(iAB+2) * C(iAB+6)
      C(iAB+16) = C(iAB+2) * C(iAB+8)
      C(iAB+17) = C(iAB+4) * C(iAB+5)
      C(iAB+18) = C(iAB+4) * C(iAB+6)
      C(iAB+19) = C(iAB+5) * C(iAB+6)
      C(iAB+20) = C(iAB+4) * C(iAB+9)
      C(iAB+21) = C(iAB+3) * C(iAB+11)
      C(iAB+22) = C(iAB+2) * C(iAB+13)
      C(iAB+23) = C(iAB+3) * C(iAB+17)
      C(iAB+24) = C(iAB+2) * C(iAB+18)
      C(iAB+25) = C(iAB+1) * C(iAB+19)
      C(iAB+26) = C(iAB+6) * C(iAB+17)

      DO i = 0,26
        C(iAB+i) = CFun * C(iAB+i)
      ENDDO

      iAB = iAB + 27
      NAB = NAB + 27

	ELSE
	  WRITE(*,*)'ERROR. Inconsistency in F*F AOs product.'
	  STOP
	ENDIF  ! End F*F cases
cc
c  F*D product
cc
cc   ---S, P, D, F, G and H functions---
      ELSEIF ( ((ijkA.EQ.3).AND.(ijkB.EQ.2)) .OR. 
     &    ((ijkA.EQ.2).AND.(ijkB.EQ.3)) ) THEN

c Order terms: F, D
      IF(ijkB.EQ.3) iChange = 1
      IF(iChange.EQ.1) THEN
        iAngX(1:3) = iAngA(1:3)
        iAngA(1:3) = iAngB(1:3)
        iAngB(1:3) = iAngX(1:3)
        RX(1:3) = RA(1:3)
        RA(1:3) = RB(1:3)
        RB(1:3) = RX(1:3)
      ENDIF
c Case [5,0,0]: (3,0,0)(2,0,0)
      IF(iAngFun(1).EQ.5 .OR. iAngFun(2).EQ.5 .OR.
     &   iAngFun(3).EQ.5) THEN

      DO i = 0,5
        iAng(iAB+i,1:3) = 0
        R(iAB+i,1:3) = RFun(1:3)
        Alph(iAB+i) = AlphFun
      ENDDO
      DO i = 1,3
        IF(iAngFun(i).EQ.5) j1 = i
      ENDDO

	iAng(iAB  ,j1) = 5
	iAng(iAB+1,j1) = 4
	iAng(iAB+2,j1) = 3
	iAng(iAB+3,j1) = 2
	iAng(iAB+4,j1) = 1

      C(iAB  ) = One
      C(iAB+1) = -(Three*RA(j1)+Two*RB(j1))
      C(iAB+2) = Three*RA(j1)**2 + 6d0*RA(j1)*RB(j1) + RB(j1)**2
      C(iAB+3) = -(RA(j1)**2 + 6d0*RA(j1)*RB(j1) 
     &           + Three*RB(j1)**2)*RA(j1)
      C(iAB+4) = (Two*RA(j1) + Three*RB(j1))*RB(j1)*RA(j1)**2
      C(iAB+5) = -(RA(j1)**3) * (RB(j1)**2)
      DO i = 0,5
        C(iAB+i) = CFun * C(iAB+i)
      ENDDO
      iAB = iAB + 6
      NAB = NAB + 6


c Cases [4,1,0]: (3,0,0)(1,1,0); (2,1,0)(2,0,0)
      ELSEIF(iAngFun(1).EQ.4 .OR. iAngFun(2).EQ.4 .OR.
     &       iAngFun(3).EQ.4) THEN

      DO i = 0,9
        iAng(iAB+i,1:3) = 0
        R(iAB+i,1:3) = RFun(1:3)
        Alph(iAB+i) = AlphFun
      ENDDO
      DO i = 1,3
        IF(iAngFun(i).EQ.4) j1 = i
        IF(iAngFun(i).EQ.1) j2 = i
      ENDDO

      iAng(iAB  ,j1) = 4
      iAng(iAB  ,j2) = 1
      iAng(iAB+1,j1) = 4
      iAng(iAB+2,j1) = 3
      iAng(iAB+2,j2) = 1
      iAng(iAB+3,j1) = 3
      iAng(iAB+4,j1) = 2
      iAng(iAB+4,j2) = 1
      iAng(iAB+5,j1) = 2
      iAng(iAB+6,j1) = 1
      iAng(iAB+6,j2) = 1
      iAng(iAB+7,j1) = 1
      iAng(iAB+8,j2) = 1

      C(iAB  ) = One
      C(iAB+1) = -(iAngA(j2)*RA(j2) + iAngB(j2)*RB(j2))
      C(iAB+2) = -(iAngA(j1)*RA(j1) + iAngB(j1)*RB(j1))
      C(iAB+3) = C(iAB+1) * C(iAB+2) 
	IF(iAngA(j1).EQ.3) THEN
        C(iAB+4) = Three*RA(j1)*(RA(j1)+RB(j1))
	  C(iAB+6) = -(RA(j1)+Three*RB(j1))*RA(j1)**2
	ELSE
        C(iAB+4) = RA(j1)**2 + Four*RA(j1)*RB(j1) + RB(j1)**2
	  C(iAB+6) = -Two*RA(j1)*RB(j1)*(RA(j1)+RB(j1))
	ENDIF
      C(iAB+5) = C(iAB+1) * C(iAB+4)
      C(iAB+7) = C(iAB+1) * C(iAB+6)
      C(iAB+8) = (RA(j1)**iAngA(j1)) * (RB(j1)**iAngB(j1))
      C(iAB+9) = C(iAB+1) * C(iAB+8)
      DO i = 0,9
        C(iAB+i) = CFun * C(iAB+i)
      ENDDO
      iAB = iAB + 10
      NAB = NAB + 10

c Cases [3,2,0]: (3,0,0)(0,2,0); (1,2,0)(2,0,0); (2,1,0)(1,1,0)
      ELSEIF(iAngFun(1).EQ.0 .OR. iAngFun(2).EQ.0 .OR.
     &       iAngFun(3).EQ.0) THEN

      DO i = 0,11
        iAng(iAB+i,1:3) = 0
        R(iAB+i,1:3) = RFun(1:3)
        Alph(iAB+i) = AlphFun
      ENDDO
      DO i = 1,3
        IF(iAngFun(i).EQ.3) j1 = i
        IF(iAngFun(i).EQ.2) j2 = i
      ENDDO

      iAng(iAB  ,j1) = 3
      iAng(iAB  ,j2) = 2
      iAng(iAB+1,j1) = 3
      iAng(iAB+1,j2) = 1
      iAng(iAB+2,j1) = 2
      iAng(iAB+2,j2) = 2
      iAng(iAB+3,j1) = 3
      iAng(iAB+4,j1) = 2
      iAng(iAB+4,j2) = 1
      iAng(iAB+5,j1) = 1
      iAng(iAB+5,j2) = 2
      iAng(iAB+6,j1) = 2
      iAng(iAB+7,j2) = 2
      iAng(iAB+8,j1) = 1
      iAng(iAB+8,j2) = 1
      iAng(iAB+9,j1) = 1
      iAng(iAB+10,j2) = 1

      C(iAB  ) = One
      C(iAB+1) = -(iAngA(j2)*RA(j2) + iAngB(j2)*RB(j2))
      C(iAB+2) = -(iAngA(j1)*RA(j1) + iAngB(j1)*RB(j1))
      C(iAB+3) = (RA(j2)**iAngA(j2)) * (RB(j2)**iAngB(j2))
      C(iAB+4) = C(iAB+1) * C(iAB+2) 
	IF(iAngA(j1).EQ.3) THEN
        C(iAB+5) = Three*RA(j1)**2
	ELSEIF(iAngA(j1).EQ.2) THEN
        C(iAB+5) = RA(j1)*(RA(j1)+Two*RB(j1))
	ELSE
        C(iAB+5) = RB(j1)*(RB(j1)+Two*RA(j1))
	ENDIF
      C(iAB+6) = C(iAB+2) * C(iAB+3) 
      C(iAB+7) = -(RA(j1)**iAngA(j1)) * (RB(j1)**iAngB(j1))
      C(iAB+8) = C(iAB+1) * C(iAB+5)
      C(iAB+9) = C(iAB+3) * C(iAB+5)
      C(iAB+10) = C(iAB+1) * C(iAB+7)
      C(iAB+11) = C(iAB+3) * C(iAB+7)
      DO i = 0,11
        C(iAB+i) = CFun * C(iAB+i)
      ENDDO
      iAB = iAB + 12
      NAB = NAB + 12

c Cases [3,1,1]: (3,0,0)(0,1,1); (2,1,0)(1,0,1); (1,1,1)(2,0,0)
      ELSEIF(iAngFun(1).EQ.3 .OR. iAngFun(2).EQ.3 .OR.
     &       iAngFun(3).EQ.3) THEN

      DO i = 0,15
        iAng(iAB+i,1:3) = 0
        R(iAB+i,1:3) = RFun(1:3)
        Alph(iAB+i) = AlphFun
      ENDDO
      DO i = 1,3
        IF(iAngFun(i).EQ.3) j1 = i
      ENDDO
	IF(j1.EQ.3) THEN
	  j2 = 1
	ELSE
	  j2 = j1 + 1
	ENDIF
	j3 = 6 - (j1+j2)

      iAng(iAB  ,j1) = 3
      iAng(iAB  ,j2) = 1
      iAng(iAB  ,j3) = 1
      iAng(iAB+1,j1) = 3
      iAng(iAB+1,j2) = 1
      iAng(iAB+2,j1) = 3
      iAng(iAB+2,j3) = 1
      iAng(iAB+3,j1) = 2
      iAng(iAB+3,j2) = 1
      iAng(iAB+3,j3) = 1
      iAng(iAB+4,j1) = 3
      iAng(iAB+5,j1) = 2
      iAng(iAB+5,j2) = 1
      iAng(iAB+6,j1) = 2
      iAng(iAB+6,j3) = 1
      iAng(iAB+7,j1) = 1
      iAng(iAB+7,j2) = 1
      iAng(iAB+7,j3) = 1
      iAng(iAB+8,j1) = 2
      iAng(iAB+9,j1) = 1
      iAng(iAB+9,j2) = 1
      iAng(iAB+10,j1) = 1
      iAng(iAB+10,j3) = 1
      iAng(iAB+11,j2) = 1
      iAng(iAB+11,j3) = 1
      iAng(iAB+12,j1) = 1
      iAng(iAB+13,j2) = 1
      iAng(iAB+14,j3) = 1

      C(iAB  ) = One
      C(iAB+1) = -(iAngA(j3)*RA(j3) + iAngB(j3)*RB(j3))
      C(iAB+2) = -(iAngA(j2)*RA(j2) + iAngB(j2)*RB(j2))
      C(iAB+3) = -(iAngA(j1)*RA(j1) + iAngB(j1)*RB(j1))
      C(iAB+4) = C(iAB+1) * C(iAB+2) 
      C(iAB+5) = C(iAB+1) * C(iAB+3)
      C(iAB+6) = C(iAB+2) * C(iAB+3)
	IF(iAngA(j1).EQ.3) THEN
        C(iAB+7) = Three*RA(j1)**2
	ELSEIF(iAngA(j1).EQ.2) THEN
        C(iAB+7) = RA(j1)*(RA(j1)+Two*RB(j1))
	ELSE
        C(iAB+7) = RB(j1)*(RB(j1)+Two*RA(j1))
	ENDIF
      C(iAB+8) = C(iAB+1) * C(iAB+6) 
      C(iAB+9) = C(iAB+1) * C(iAB+7) 
      C(iAB+10) = C(iAB+2) * C(iAB+7)
      C(iAB+11) = -(RA(j1)**iAngA(j1)) * (RB(j1)**iAngB(j1)) ! 0^0
      C(iAB+12) = C(iAB+1) * C(iAB+10)
      C(iAB+13) = C(iAB+1) * C(iAB+11)
      C(iAB+14) = C(iAB+2) * C(iAB+11)
      C(iAB+15) = C(iAB+1) * C(iAB+14)
      DO i = 0,15
        C(iAB+i) = CFun * C(iAB+i)
      ENDDO
      iAB = iAB + 16
      NAB = NAB + 16

c Cases [2,2,1]: (2,1,0)(0,0,2); (2,1,0)(0,1,1); (1,1,1)(1,1,0)
	ELSE

      DO i = 0,17
        iAng(iAB+i,1:3) = 0
        R(iAB+i,1:3) = RFun(1:3)
        Alph(iAB+i) = AlphFun
      ENDDO
	j1 = 0
      DO i = 1,3
        IF(iAngA(i).EQ.2) j1 = i
        IF(iAngFun(i).EQ.1) j3 = i
      ENDDO
	IF(j1.EQ.0) j1 = j3 + 1
	IF(j1.EQ.4) j1 = 1
	j2 = 6 - (j1+j3)

      iAng(iAB  ,j1) = 2
      iAng(iAB  ,j2) = 2
      iAng(iAB  ,j3) = 1
      iAng(iAB+1,j1) = 2
      iAng(iAB+1,j2) = 2
      iAng(iAB+2,j1) = 2
      iAng(iAB+2,j2) = 1
      iAng(iAB+2,j3) = 1
      iAng(iAB+3,j1) = 1
      iAng(iAB+3,j2) = 2
      iAng(iAB+3,j3) = 1
      iAng(iAB+4,j1) = 2
      iAng(iAB+4,j2) = 1
      iAng(iAB+5,j1) = 2
      iAng(iAB+5,j3) = 1
      iAng(iAB+6,j1) = 1
      iAng(iAB+6,j2) = 2
      iAng(iAB+7,j2) = 2
      iAng(iAB+7,j3) = 1
      iAng(iAB+8,j1) = 1
      iAng(iAB+8,j2) = 1
      iAng(iAB+8,j3) = 1
      iAng(iAB+9,j1) = 2
      iAng(iAB+10,j2) = 2
      iAng(iAB+11,j1) = 1
      iAng(iAB+11,j2) = 1
      iAng(iAB+12,j1) = 1
      iAng(iAB+12,j3) = 1
      iAng(iAB+13,j2) = 1
      iAng(iAB+13,j3) = 1
      iAng(iAB+14,j1) = 1
      iAng(iAB+15,j2) = 1
      iAng(iAB+16,j3) = 1

      C(iAB  ) = One
      C(iAB+1) = -(iAngA(j3)*RA(j3) + iAngB(j3)*RB(j3))
      C(iAB+2) = -(iAngA(j2)*RA(j2) + iAngB(j2)*RB(j2))
      C(iAB+3) = -(iAngA(j1)*RA(j1) + iAngB(j1)*RB(j1))
      C(iAB+4) = C(iAB+1) * C(iAB+2)
      C(iAB+5) = (RA(j2)**iAngA(j2)) * (RB(j2)**iAngB(j2)) ! 0^0
      C(iAB+6) = C(iAB+1) * C(iAB+3)
      C(iAB+7) = (RA(j1)**iAngA(j1)) * (RB(j1)**iAngB(j1)) ! 0^0
      C(iAB+8) = C(iAB+2) * C(iAB+3)
      C(iAB+9) = C(iAB+1) * C(iAB+5)
      C(iAB+10) = C(iAB+1) * C(iAB+7)
      C(iAB+11) = C(iAB+1) * C(iAB+8)
      C(iAB+12) = C(iAB+3) * C(iAB+5)
      C(iAB+13) = C(iAB+2) * C(iAB+7)
      C(iAB+14) = C(iAB+1) * C(iAB+12)
      C(iAB+15) = C(iAB+1) * C(iAB+13)
      C(iAB+16) = C(iAB+5) * C(iAB+7)
      C(iAB+17) = C(iAB+1) * C(iAB+16)
      DO i = 0,17
        C(iAB+i) = CFun * C(iAB+i)
      ENDDO
      iAB = iAB + 18
      NAB = NAB + 18

	ENDIF

cc
c  F*P product
cc
cc   ---S, P, D, F, and G functions---
      ELSEIF ( ((ijkA.EQ.3).AND.(ijkB.EQ.1)) .OR. 
     &    ((ijkA.EQ.1).AND.(ijkB.EQ.3)) ) THEN

c Order terms: F, P
      IF(ijkB.EQ.3) iChange = 1 
      IF(iChange.EQ.1) THEN
        iAngX(1:3) = iAngA(1:3)
        iAngA(1:3) = iAngB(1:3)
        iAngB(1:3) = iAngX(1:3)
        RX(1:3) = RA(1:3)
        RA(1:3) = RB(1:3)
        RB(1:3) = RX(1:3)
      ENDIF

c Case [4,0,0]: (3,0,0)(1,0,0)
	IF(iAngFun(1).EQ.4 .OR. iAngFun(2).EQ.4 .OR.
     &   iAngFun(3).EQ.4) THEN

      DO i = 0,4
        iAng(iAB+i,1:3) = 0
        R(iAB+i,1:3) = RFun(1:3)
        Alph(iAB+i) = AlphFun
      ENDDO
      DO i = 1,3
        IF(iAngFun(i).EQ.4) j1 = i
      ENDDO

      iAng(iAB  ,j1) = 4
      iAng(iAB+1,j1) = 3
      iAng(iAB+2,j1) = 2
      iAng(iAB+3,j1) = 1

	C(iAB  ) = One
      C(iAB+1) = -(Three*RA(j1)+RB(j1))
      C(iAB+2) = Three*RA(j1)*(RA(j1)+RB(j1))
      C(iAB+3) = -(RA(j1)+Three*RB(j1))*RA(j1)**2
      C(iAB+4) = RB(j1)*RA(j1)**3
      DO i = 0,4
        C(iAB+i) = CFun * C(iAB+i)
      ENDDO
      iAB = iAB + 5
      NAB = NAB + 5

c Cases [3,1,0]: (3,0,0)(0,1,0); (2,1,0)(1,0,0) 
	ELSEIF(iAngFun(1).EQ.3 .OR. iAngFun(2).EQ.3 .OR.
     &       iAngFun(3).EQ.3) THEN

      DO i = 0,7
        iAng(iAB+i,1:3) = 0
        R(iAB+i,1:3) = RFun(1:3)
        Alph(iAB+i) = AlphFun
      ENDDO
      DO i = 1,3
        IF(iAngFun(i).EQ.3) j1 = i
        IF(iAngFun(i).EQ.1) j2 = i
      ENDDO

      iAng(iAB  ,j1) = 3
      iAng(iAB  ,j2) = 1
      iAng(iAB+1,j1) = 3
      iAng(iAB+2,j1) = 2
      iAng(iAB+2,j2) = 1
      iAng(iAB+3,j1) = 2
      iAng(iAB+4,j1) = 1
      iAng(iAB+4,j2) = 1
      iAng(iAB+5,j1) = 1
      iAng(iAB+6,j2) = 1

	C(iAB  ) = One
	C(iAB+1) = -(iAngA(j2)*RA(j2)+iAngB(j2)*RB(j2))
	C(iAB+2) = -(iAngA(j1)*RA(j1)+iAngB(j1)*RB(j1))
	C(iAB+3) = C(iAB+1) * C(iAB+2)
	IF(iAngA(j1).EQ.3) THEN
	  C(iAB+4) = Three*RA(j1)**2
	ELSE
	  C(iAB+4) = RA(j1)*(RA(j1)+Two*RB(j1))
	ENDIF
	C(iAB+5) = C(iAB+1) * C(iAB+4)
	C(iAB+6) = -(RA(j1)**iAngA(j1))*(RB(j1)**iAngB(j1))
	C(iAB+7) = C(iAB+1) * C(iAB+6)
      DO i = 0,7
        C(iAB+i) = CFun * C(iAB+i)
      ENDDO
      iAB = iAB + 8
      NAB = NAB + 8

c Case [2,2,0]: (2,1,0)(0,1,0)
	ELSEIF(iAngFun(1).EQ.0 .OR. iAngFun(2).EQ.0 .OR.
     &       iAngFun(3).EQ.0) THEN

      DO i = 0,8
        iAng(iAB+i,1:3) = 0
        R(iAB+i,1:3) = RFun(1:3)
        Alph(iAB+i) = AlphFun
      ENDDO
      DO i = 1,3
        IF(iAngA(i).EQ.2) j1 = i
        IF(iAngA(i).EQ.1) j2 = i
      ENDDO

      iAng(iAB  ,j1) = 2
      iAng(iAB  ,j2) = 2
      iAng(iAB+1,j1) = 2
      iAng(iAB+1,j2) = 1
      iAng(iAB+2,j1) = 1
      iAng(iAB+2,j2) = 2
      iAng(iAB+3,j1) = 2
      iAng(iAB+4,j2) = 2
      iAng(iAB+5,j1) = 1
      iAng(iAB+5,j2) = 1
      iAng(iAB+6,j1) = 1
      iAng(iAB+7,j2) = 1

	C(iAB  ) = One
	C(iAB+1) = -(RA(j2)+RB(j2))
	C(iAB+2) = -Two*RA(j1)
	C(iAB+3) = RA(j2)*RB(j2)
	C(iAB+4) = RA(j1)**2
	C(iAB+5) = C(iAB+1) * C(iAB+2)
	C(iAB+6) = C(iAB+2) * C(iAB+3)
	C(iAB+7) = C(iAB+1) * C(iAB+4)
	C(iAB+8) = C(iAB+3) * C(iAB+4)
      DO i = 0,8
        C(iAB+i) = CFun * C(iAB+i)
      ENDDO
      iAB = iAB + 9
      NAB = NAB + 9

c Cases [2,1,1]: (2,1,0)(0,0,1); (1,1,1)(1,0,0)
	ELSE

      DO i = 0,11
        iAng(iAB+i,1:3) = 0
        R(iAB+i,1:3) = RFun(1:3)
        Alph(iAB+i) = AlphFun
      ENDDO
      DO i = 1,3
        IF(iAngFun(i).EQ.2) j1 = i
      ENDDO
	IF(j1.EQ.3) THEN
	  j2 = 1
	ELSE
	  j2 = j1+1
	ENDIF
	j3 = 6 - (j1+j2)

      iAng(iAB  ,j1) = 2
      iAng(iAB  ,j2) = 1
      iAng(iAB  ,j3) = 1
      iAng(iAB+1,j1) = 2
      iAng(iAB+1,j2) = 1
      iAng(iAB+2,j1) = 2
      iAng(iAB+2,j3) = 1
      iAng(iAB+3,j1) = 1
      iAng(iAB+3,j2) = 1
      iAng(iAB+3,j3) = 1
      iAng(iAB+4,j1) = 2
      iAng(iAB+5,j1) = 1
      iAng(iAB+5,j2) = 1
      iAng(iAB+6,j1) = 1
      iAng(iAB+6,j3) = 1
      iAng(iAB+7,j2) = 1
      iAng(iAB+7,j3) = 1
      iAng(iAB+8,j1) = 1
      iAng(iAB+9,j2) = 1
      iAng(iAB+10,j3) = 1

	C(iAB  ) = One
	C(iAB+1) = -(iAngA(j3)*RA(j3) + iAngB(j3)*RB(j3))
	C(iAB+2) = -(iAngA(j2)*RA(j2) + iAngB(j2)*RB(j2))
	C(iAB+3) = -(iAngA(j1)*RA(j1) + iAngB(j1)*RB(j1))
	C(iAB+4) = C(iAB+1) * C(iAB+2)
	C(iAB+5) = C(iAB+1) * C(iAB+3)
	C(iAB+6) = C(iAB+2) * C(iAB+3)
	IF(iAngA(j1).EQ.2) THEN
	  C(iAB+7) = RA(j1)**2
	ELSE
	  C(iAB+7) = RA(j1)*RB(j1)
	ENDIF
	C(iAB+8) = C(iAB+1) * C(iAB+6)
	C(iAB+9) = C(iAB+1) * C(iAB+7)
	C(iAB+10) = C(iAB+2) * C(iAB+7)
	C(iAB+11) = C(iAB+1) * C(iAB+10)
      DO i = 0,11
        C(iAB+i) = CFun * C(iAB+i)
      ENDDO
      iAB = iAB + 12
      NAB = NAB + 12

	ENDIF

cc
c  F*S product
cc
cc   ---S, P, D, and F functions---
      ELSEIF( ((ijkA.EQ.3).AND.(ijkB.EQ.0)) .OR.
     &    ((ijkA.EQ.0).AND.(ijkB.EQ.3)) ) THEN

c Order terms: F, S
      IF(ijkB.EQ.3) iChange = 1 
      IF(iChange.EQ.1) THEN
        iAngX(1:3) = iAngA(1:3)
        iAngA(1:3) = iAngB(1:3)
        iAngB(1:3) = iAngX(1:3)
        RX(1:3) = RA(1:3)
        RA(1:3) = RB(1:3)
        RB(1:3) = RX(1:3)
      ENDIF

c case (3,0,0)(0,0,0)
	IF(iAngA(1).EQ.3 .OR. iAngA(2).EQ.3 .OR.
     &   iAngA(3).EQ.3) THEN

      DO i = 0,3
        iAng(iAB+i,1:3) = 0
        R(iAB+i,1:3) = RFun(1:3)
        Alph(iAB+i) = AlphFun
	ENDDO
	DO i = 1,3
	  IF(iAngFun(i).EQ.3) j1 = i
      ENDDO

      iAng(iAB  ,j1) = 3
      iAng(iAB+1,j1) = 2
      iAng(iAB+2,j1) = 1

	C(iAB  ) = CFun
	C(iAB+1) = -CFun * Three * RA(j1)
	C(iAB+2) = -C(iAB+1) * RA(j1)
      C(iAB+3) = -CFun * RA(j1)**3

      iAB = iAB + 4
      NAB = NAB + 4

c case (2,1,0)(0,0,0)
	ELSEIF(iAngA(1).EQ.2 .OR. iAngA(2).EQ.2 .OR.
     &   iAngA(3).EQ.2) THEN

      DO i = 0,5
        iAng(iAB+i,1:3) = 0
        R(iAB+i,1:3) = RFun(1:3)
        Alph(iAB+i) = AlphFun
	ENDDO
	DO i = 1,3
        IF(iAngFun(i).EQ.2) j1 = i
        IF(iAngFun(i).EQ.1) j2 = i
      ENDDO

      iAng(iAB  ,j1) = 2
      iAng(iAB  ,j2) = 1
      iAng(iAB+1,j1) = 2
      iAng(iAB+2,j1) = 1
      iAng(iAB+2,j2) = 1
      iAng(iAB+3,j1) = 1
      iAng(iAB+4,j2) = 1

	C(iAB  ) = One
	C(iAB+1) = -RA(j2)
	C(iAB+2) = -Two * RA(j1)
	C(iAB+3) = C(iAB+1) * C(iAB+2)
	C(iAB+4) = RA(j1)**2
	C(iAB+5) = C(iAB+1) * C(iAB+4)

      DO i = 0,5
        C(iAB+i) = CFun * C(iAB+i)
      ENDDO
      iAB = iAB + 6
      NAB = NAB + 6

c case (1,1,1)(0,0,0)
	ELSE

      DO i = 0,7
        iAng(iAB+i,1:3) = 0
        R(iAB+i,1:3) = RFun(1:3)
        Alph(iAB+i) = AlphFun
      ENDDO

      iAng(iAB  ,1) = 1
      iAng(iAB  ,2) = 1
      iAng(iAB  ,3) = 1
      iAng(iAB+1,1) = 1
      iAng(iAB+1,2) = 1
      iAng(iAB+2,1) = 1
      iAng(iAB+2,3) = 1
      iAng(iAB+3,2) = 1
      iAng(iAB+3,3) = 1
      iAng(iAB+4,1) = 1
      iAng(iAB+5,2) = 1
      iAng(iAB+6,3) = 1

	C(iAB  ) = One
	C(iAB+1) = -RA(3)
	C(iAB+2) = -RA(2)
	C(iAB+3) = -RA(1)
	C(iAB+4) = C(iAB+1) * C(iAB+2)
	C(iAB+5) = C(iAB+1) * C(iAB+3)
	C(iAB+6) = C(iAB+2) * C(iAB+3)
	C(iAB+7) = C(iAB+1) * C(iAB+6)

      DO i = 0,7
        C(iAB+i) = CFun * C(iAB+i)
      ENDDO
      iAB = iAB + 8
      NAB = NAB + 8

      ENDIF
      ENDIF

c Re-order A,B terms if necessary
      IF(iChange.EQ.1) THEN
        iAngX(1:3) = iAngA(1:3)
        iAngA(1:3) = iAngB(1:3)
        iAngB(1:3) = iAngX(1:3)
        RX(1:3) = RA(1:3)
        RA(1:3) = RB(1:3)
        RB(1:3) = RX(1:3)
      ENDIF

      RETURN
      END
