
        SUBROUTINE get_basis(Basis,AtLab,RAt,iAng,Alph,
     & ROrb,COrb,iContr,iLAO,cLAO,NAt,NBas,NOrb,PureCart)
c*****************************************************c
c                                                     c
c       Reads AOs from disk and makes                 c
c       uncontracted basis                            c
c                                                     c
c    INPUT                                            c
c      Basis  atomic basis                            c
c      AtLab  atomic label                            c
c      RAt    atomic position                         c
c                                                     c
c    OUTPUT                                           c
c      ROrb   center  (uncontracted function)         c
c      COrb   coeffs  (uncontracted function)         c
c      iAng   angular indices uncontracted functions  c
c      Alph   exponents (uncontracted functions)      c
c      NC     contracted funcions                     c
c      iLAO   l number for each AO                    c 
c      cLAO   l name   for each AO                    c 
c      iContr basis contraction vector                c
c             contains the # of functions             c
c             in each shell                           c
c                                                     c
c                                                     c
c       DCC: 04/2010                                  c
c       DCC: 06/2010  f-orbitals                      c
c                                                     c
c*****************************************************c
      IMPLICIT REAL*8  (a-h, o-z)
      INTEGER   :: iAng(NOrb,3),iContr(NBas),iAngv(3),iLAO(NBas)
      REAL*8    :: RAt(NAt,3),RAB(3),ROrb(NOrb,3),COrb(NOrb),Alph(NOrb)
      CHARACTER :: Basis*80, OrbType*2,AtLab(*)*2,cLAO(NBas)*4,Aux*255
      LOGICAL   :: PureCart
      SAVE Zero,Half,One,Two,Four
      DATA Zero,Half,One,Two,Four / 0d0,0.5d0,1d0,2d0,4d0  /


c Get the auxiliarity folder
      CALL getenv("WSAUX",Aux)

c Initialize some variables
      ND  = 6   ! cart d-orbitals
      NF  = 10  ! cart f-orbitals
      IF(PureCart) ND = 5 ! pure d-orbitals
      iOrb   = 1
      iCount = 1
      RAB(:) = Zero

      ATOMS: DO k = 1,NAt
        OPEN(Unit=10, Status="old", File=
     &"/"//trim(Aux)//"/basis/"
     &//trim(Basis)//"/"//trim(AtLab(k))//".dat")

        READ(10,*) NShell
        SHELLS: DO iShell = 1,NShell
        READ(10,'(A2)',advance='no') OrbType
        READ(10,*) NC

        iOrb0 = iOrb
c S type functions
        IF (OrbType.eq.'S ') THEN
          xOrbType = One
          iContr(iCount) = NC
          iLAO(iCount) = 0
          cLAO(iCount) = 's   '
          iCount = iCount + 1
          DO j = 1,NC
            READ(10,*) Alph(iOrb),COrb(iOrb)
            DO jj = 1,3
              ROrb(iOrb,jj) = RAt(k,jj)
            ENDDO
            iAng(iOrb,1) = 0
            iAng(iOrb,2) = 0
            iAng(iOrb,3) = 0
            iOrb = iOrb + 1
          ENDDO
          NOrbSh = NC

c P type functions
        ELSEIF (OrbType.eq.'P ') THEN
          xOrbType = DSQRT(3.d0)
          DO j = 0,2
            iContr(iCount+j) = NC
            iLAO(iCount+j) = 1
          ENDDO
          cLAO(iCount)   = 'px  '
          cLAO(iCount+1) = 'py  '
          cLAO(iCount+2) = 'pz  '
          iCount = iCount+3
          DO i = 1,3*NC
            DO j = 1,3
              iAng(iOrb+i-1,j) = 0
              ROrb(iOrb+i-1,j) = RAt(k,j)
            ENDDO
          ENDDO
          DO i = 1,NC
            iAng(iOrb+i-1,1)      = 1
            iAng(iOrb+i-1+NC,2)   = 1
            iAng(iOrb+i-1+2*NC,3) = 1
            READ(10,*) exp1, coef1
            Alph(iOrb+i-1)      = exp1
            Alph(iOrb+i-1+NC)   = exp1
            Alph(iOrb+i-1+2*NC) = exp1
            COrb(iOrb+i-1)      = coef1
            COrb(iOrb+i-1+NC)   = coef1
            COrb(iOrb+i-1+2*NC) = coef1
          ENDDO
          NOrbSh = 3*NC
          iOrb=iOrb+NOrbSh

c SP type functions
        ELSEIF (OrbType.eq.'SP') THEN
          xOrbType = DSQRT(Four)
          iContr(iCount) = NC
          iLAO(iCount)   = 0
          cLAO(iCount)   = 's   '
          DO j = 1,3
            iContr(iCount+j) = NC
            iLAO(iCount+j) = 1
          ENDDO
          cLAO(iCount+1) = 'px  '
          cLAO(iCount+2) = 'py  '
          cLAO(iCount+3) = 'pz  '
          iCount = iCount+4
          DO i = 1,4*NC
            DO j = 1,3
              iAng(iOrb+i-1,j) = 0
              ROrb(iOrb+i-1,j) = RAt(k,j)
            ENDDO
          ENDDO
          DO i=NC,2*NC-1
            iAng(iOrb+i,1)      = 1
            iAng(iOrb+i+NC,2)   = 1
            iAng(iOrb+i+2*NC,3) = 1
          ENDDO
          DO i=1,NC
            READ(10,*) exp1, coef1, coef2
            Alph(iOrb+i-1)      = exp1
            Alph(iOrb+i-1+NC)   = exp1
            Alph(iOrb+i-1+2*NC) = exp1
            Alph(iOrb+i-1+3*NC) = exp1
            COrb(iOrb+i-1)      = coef1
            COrb(iOrb+i-1+NC)   = coef2
            COrb(iOrb+i-1+2*NC) = coef2
            COrb(iOrb+i-1+3*NC) = coef2
          ENDDO
          NOrbSh = 4*NC
          iOrb = iOrb+NOrbSh

c  D type functions
        ELSEIF (OrbType.eq.'D ') THEN
          DO j = 0,ND-1
             iContr(iCount+j) = NC
             iLAO(iCount+j) = 2
          ENDDO
          IF(PureCart) THEN
            cLAO(iCount)   = 'd1  '  
            cLAO(iCount+1) = 'd2  ' 
            cLAO(iCount+2) = 'd3  '
            cLAO(iCount+3) = 'd4  ' 
            cLAO(iCount+4) = 'd5  ' 
          ELSE
            cLAO(iCount)   = 'dx2 '    
            cLAO(iCount+1) = 'dxy '   
            cLAO(iCount+2) = 'dy2 ' 
            cLAO(iCount+3) = 'dxz '   
            cLAO(iCount+4) = 'dyz '   
            cLAO(iCount+5) = 'dz2 '   
          ENDIF
          iCount=iCount+ND
          DO i=1,NC
            READ(10,*) exp1, coef1
            DO j = 0,ND-1
              Alph(iOrb+i-1+j*NC) = exp1
              COrb(iOrb+i-1+j*NC) = coef1
            ENDDO
          ENDDO
          DO i = 1,ND*NC
             DO j = 1,3
               iAng(iOrb+i-1,j) = 0
               ROrb(iOrb+i-1,j) = RAt(k,j)
             ENDDO
          ENDDO
          DO i = 1,NC
            IF(PureCart) THEN
c             Pure d-orbitals order in Q-Chem
              xOrbType = DSQRT(5.d0)
              iAng(iOrb+i-1,1)      = 1  ! xy  (d1)
              iAng(iOrb+i-1,2)      = 1  ! xy  (d1)
              iAng(iOrb+i-1+NC,2)   = 1  ! yz  (d2)
              iAng(iOrb+i-1+NC,3)   = 1  ! yz  (d2)
              iAng(iOrb+i-1+2*NC,3) = 2  ! z2 from 2z2-x2-y2 (d3)
              iAng(iOrb+i-1+3*NC,1) = 1  ! xz  (d4)
              iAng(iOrb+i-1+3*NC,3) = 1  ! xz  (d4)
              iAng(iOrb+i-1+4*NC,1) = 2  ! x2 from x2-y2 (d5)
            ELSE
c             Cartesian d-orbitals order in Q-Chem
              xOrbType = DSQRT(8.d0)
              iAng(iOrb+i-1,1)      = 2  ! xx
              iAng(iOrb+i-1+NC,1)   = 1  ! xy
              iAng(iOrb+i-1+NC,2)   = 1  ! xy
              iAng(iOrb+i-1+2*NC,2) = 2  ! yy
              iAng(iOrb+i-1+3*NC,1) = 1  ! xz 
              iAng(iOrb+i-1+3*NC,3) = 1  ! xz
              iAng(iOrb+i-1+4*NC,2) = 1  ! yz
              iAng(iOrb+i-1+4*NC,3) = 1  ! yz
              iAng(iOrb+i-1+5*NC,3) = 2  ! zz
            ENDIF
          ENDDO
          NOrbSh = ND*NC
          iOrb = iOrb + NOrbSh

c  F type functions
        ELSEIF (OrbType.eq.'F ') THEN
           DO j = 0,NF-1
             iContr(iCount+j) = NC
             iLAO(iCount+j) = 3
           ENDDO
           IF(PureCart) THEN
           ELSE
            cLAO(iCount)   = 'fxxx'
            cLAO(iCount+1) = 'fxxy'
            cLAO(iCount+2) = 'fxyy'
            cLAO(iCount+3) = 'fyyy'
            cLAO(iCount+4) = 'fxxz'
            cLAO(iCount+5) = 'fxyz'
            cLAO(iCount+6) = 'fyyz'
            cLAO(iCount+7) = 'fxzz'
            cLAO(iCount+8) = 'fyzz'
            cLAO(iCount+9) = 'fzzz'
          ENDIF 
          iCount=iCount+NF
          DO i=1,NC
             READ(10,*) exp1, coef1
             DO j = 0,NF-1
               Alph(iOrb+i-1+j*NC) = exp1
               COrb(iOrb+i-1+j*NC) = coef1
             ENDDO
          ENDDO
          DO i = 1,NF*NC
             DO j = 1,3
               iAng(iOrb+i-1,j) = 0
               ROrb(iOrb+i-1,j) = RAt(k,j)
             ENDDO
          ENDDO
          DO i = 1,NC
            IF(PureCart) THEN
c             Pure f-orbitals order in Q-Chem
            ELSE
c             Cartesian f-orbitals order in Q-Chem
              xOrbType = DSQRT(17.36656146d0)
              iAng(iOrb+i-1     ,1) = 3  ! xxx
              iAng(iOrb+i-1+  NC,1) = 2  ! xxy
              iAng(iOrb+i-1+  NC,2) = 1  ! xxy
              iAng(iOrb+i-1+2*NC,1) = 1  ! xyy
              iAng(iOrb+i-1+2*NC,2) = 2  ! xyy
              iAng(iOrb+i-1+3*NC,2) = 3  ! yyy
              iAng(iOrb+i-1+4*NC,1) = 2  ! xxz
              iAng(iOrb+i-1+4*NC,3) = 1  ! xxz
              iAng(iOrb+i-1+5*NC,1) = 1  ! xyz
              iAng(iOrb+i-1+5*NC,2) = 1  ! xyz
              iAng(iOrb+i-1+5*NC,3) = 1  ! xyz
              iAng(iOrb+i-1+6*NC,2) = 2  ! yyz
              iAng(iOrb+i-1+6*NC,3) = 1  ! yyz
              iAng(iOrb+i-1+7*NC,1) = 1  ! xzz
              iAng(iOrb+i-1+7*NC,3) = 2  ! xzz
              iAng(iOrb+i-1+8*NC,2) = 1  ! yzz
              iAng(iOrb+i-1+8*NC,3) = 2  ! yzz
              iAng(iOrb+i-1+9*NC,3) = 3  ! zzz
            ENDIF
          ENDDO
          NOrbSh = NF*NC
          iOrb = iOrb + NOrbSh

        ENDIF

c Include normalization factors
        DO ii = 0,NOrbSh-1
           iA = iOrb0 + ii
           DO j = 1,3
             iAngv(j) = iAng(iA,j)
           ENDDO
           CALL norm_basis(iAngv,Alph(iA),COrb(iA),PureCart)
        ENDDO

c Normalize Shell
        S   = 0.d0
        DO ii = 0,NOrbSh-1
          iA = iOrb0 + ii
          DO jj = ii,NOrbSh-1
            iB = iOrb0 + jj
            CAB = 2.d0
            IF(iA .EQ. iB) CAB = 1.d0
            DO j = 1,3
              iAngv(j) = iAng(iA,j) + iAng(iB,j)
            ENDDO
            AlphAB = Alph(iA) + Alph(iB)
            CALL overlap(iAngv,AlphAB,RAB,SAB)
            S = S + CAB*COrb(iA)*COrb(iB)*SAB
          ENDDO
        ENDDO
  
        S = xOrbType / DSQRT(S)
        DO ii = 0,NOrbSh-1
          iA = iOrb0 + ii
          COrb(iA) = COrb(iA) * S
        ENDDO

        ENDDO SHELLS

        CLOSE(10)
      ENDDO ATOMS

c Check consistency
      IF((iOrb-1).ne.NOrb) THEN
        WRITE(2,*) 'ERROR. NOrb inconsistency in get_basis'
        stop
      ENDIF
      IF((iCount-1).ne.NBas) THEN
        WRITE(2,*) 'ERROR. NBas inconsistency in get_basis'
        stop
      ENDIF

      RETURN
      END



        SUBROUTINE norm_basis(iAng,Alph,COrb,PureCart)
c*****************************************************c
c                                                     c
c       Updates COrb including the its                c
c       normalication factor                          c
c                                                     c
c       INPUT                                         c
c         iAng  :   cartesian indices                 c
c         Alph  :   gaussian exponent                 c
c         COrb  :   gaussian coefficient              c
c         PureCart: pure D-orbitals?                  c
c                                                     c
c       OUTPUT                                        c
c         COrb  :   gaussian coeff updated            c
c                                                     c
c       DCC: 04/10                                    c
c       DCC: 06/10 f-orbitals                         c
c                                                     c
c*****************************************************c
      IMPLICIT REAL*8  (a-h, o-z)
      INTEGER    iAng(3)
      LOGICAL    PureCart
      SAVE Zero,Half,One,Two,Four
      DATA Zero,Half,One,Two,Four / 0d0,0.5d0,1d0,2d0,4d0  /

      PI = Four * DATAN(One)
      iTot = iAng(1) + iAng(2) + iAng(3)

c S orbital
      IF(iTot .eq. 0) THEN
        C = (Two*Alph/PI)**(0.75d0)

c P orbital
      ELSEIF(iTot .eq. 1) THEN
        C = Two*DSQRT(Alph)*(Two*Alph/PI)**0.75d0

c D orbital 
      ELSEIF(iTot .eq. 2) THEN
        C = Four*Alph*(Two*Alph/PI)**(0.75d0)
        IF(.NOT.PureCart) THEN
          DO i=1,3
            IF(iAng(i).eq.2) THEN
              C = C / DSQRT(3.d0)
            ENDIF
          ENDDO
        ELSE
          DO i=1,3
            IF(iAng(i).eq.2) THEN
              C = C / DSQRT(3.d0)
            ENDIF
          ENDDO
          IF(iAng(1) .eq. 2) THEN
            C = C * DSQRT(0.75d0)
          ENDIF
          IF(iAng(3) .eq. 2) THEN
            C = C * DSQRT(0.25d0)
          ENDIF
        ENDIF

c F orbital
      ELSEIF(iTot .eq. 3) THEN
        C = Two**(1.5d0)*(Two*Alph)**(2.25d0)*PI**(0.75d0)
        DO i = 1,3
          IF(iAng(i).eq.3) THEN
            C = C / DSQRT(15.0d0) 
          ELSEIF(iAng(i).eq.2) THEN
            C = C / DSQRT(3.0d0) 
          ENDIF
        ENDDO

      ENDIF

      COrb = COrb * C

      RETURN
      END

      SUBROUTINE lc_cartesian_orbitals(iAngv,LCC,NLC,iAngi)
      INTEGER :: iAngv(3,3),iAngi(3),NLC
      REAL*8  :: LCC(3)

        NLC = 1
        iAngv(1,:) = iAngi
        LCC(1) = 1.d0
        sum = iAngi(1)+iAngi(2)+iAngi(3)
        IF ((iAngi(1).EQ.-2).AND.(sum).EQ.0) THEN
          iAngv(1,:) = (/2,0,0/)
          iAngv(2,:) = (/0,2,0/)
          iAngv(3,:) = (/0,0,2/)
          LCC(:) = (/-1.d0,-1.d0,2.d0/)
          NLC = 3
        ELSEIF ((iAngi(1).EQ.2).AND.(sum).EQ.0) THEN
          iAngv(1,:) = (/2,0,0/)
          iAngv(2,:) = (/0,2,0/)
          iAngv(3,:) = (/0,0,0/)
          LCC(:) = (/1.d0,-1.d0,0.d0/)
          NLC = 2
        ELSEIF ((iAngi(1).EQ.2).AND.(iAngi(3).EQ.1)) THEN
          iAngv(1,:) = (/2,0,1/)
          iAngv(2,:) = (/0,2,1/)
          iAngv(3,:) = (/0,0,0/)
          LCC(:) = (/-1.d0,-1.d0,0.d0/)
          NLC = 2
          IF ((iAngi(2).EQ.-2)) LCC(1) = 1.d0
        ELSEIF (iAngi(1).EQ.3) THEN
          iAngv(1,:) = (/3,0,0/)
          iAngv(2,:) = (/1,2,0/)
          iAngv(3,:) = (/0,0,0/)
          LCC(:) = (/-1.d0,-1.d0,0.d0/)
          NLC = 2
          IF ((iAngi(2).EQ.-2)) LCC(1) = 1.d0
        ELSEIF (ABS(iAngi(2)).EQ.3) THEN
          iAngv(1,:) = (/2,1,0/)
          iAngv(2,:) = (/0,3,0/)
          iAngv(3,:) = (/0,0,0/)
          LCC(:) = (/-1.d0,-1.d0,0.d0/)
          NLC = 2
          IF ((iAngi(2).EQ.-3)) LCC(1) = 1.d0
        ENDIF
      END