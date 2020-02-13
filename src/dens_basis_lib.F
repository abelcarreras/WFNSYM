
        SUBROUTINE dens_basis(iAng,ROrb,R,NCvec,NAt,NOrbX,
     &                       NPrim,NTotShell,NShell,shell_type,NOrb)
c*****************************************************c
c                                                     c
c       Reads AOs from disk                           c
c                                                     c
c       INPUT                                         c
c         AtLab  atomic label                         c
c         R      atomic position                      c
c                                                     c
c       OUTPUT                                        c
c         ROrb   center of cartesian function         c
c         iAng   angular indices cartesian functions  c
c         iOrb                                        c
c         NC     contracted funcions                  c
c         C                                           c
c         Alpha                                       c
c         NCvec  Basis contraction vector             c
c                contains the # of functions          c
c                in each shell                        c
c                                                     c
c                                                     c
c       DCC: 05/09                                    c
c                                                     c
c*****************************************************c
      INTEGER   ::  NPrim(NTotShell),NShell(NAt),
     &              shell_type(NTotShell),iAng(NOrbX,3),NCvec(NOrb)
      REAL*8    ::  R(NAt,3),ROrb(NOrbX,3)
      CHARACTER ::  Aux*255

c Initialize some variables
      iOrb   = 1
      iCount = 1
      iShellC = 0

      DO k = 1,NAt
        DO iShell = 1,NShell(k)
         iShellC = iShellC + 1
         NC = NPrim(iShellC)

c S type functions
         IF  (shell_type(iShellC).eq.0) THEN
          NCvec(iCount) = NC
          iCount = iCount + 1
          DO j = 1,NC
            DO jj = 1,3
              ROrb(iOrb,jj) = R(k,jj)
            ENDDO
            iAng(iOrb,1) = 0
            iAng(iOrb,2) = 0
            iAng(iOrb,3) = 0
            iOrb = iOrb + 1
          ENDDO


c P type functions
         ELSEIF  (shell_type(iShellC).eq.1) THEN
          NCvec(iCount)   = NC
          NCvec(iCount+1) = NC
          NCvec(iCount+2) = NC
          iCount = iCount+3
          DO i = 1,3*NC
            DO j = 1,3
              iAng(iOrb+i-1,j) = 0
              ROrb(iOrb+i-1,j) = R(k,j)
            ENDDO
          ENDDO
          DO i = 1,NC
            iAng(iOrb+i-1,1)      = 1
            iAng(iOrb+i-1+NC,2)   = 1
            iAng(iOrb+i-1+2*NC,3) = 1
          ENDDO
          iOrb=iOrb+3*NC

c SP type functions
         ELSEIF  (shell_type(iShellC).eq.-1) THEN
          NCvec(iCount)   = NC
          NCvec(iCount+1) = NC
          NCvec(iCount+2) = NC
          NCvec(iCount+3) = NC
          iCount = iCount+4
          DO i = 1,4*NC
            DO j = 1,3
              iAng(iOrb+i-1,j) = 0
              ROrb(iOrb+i-1,j) = R(k,j)
            ENDDO
          ENDDO
          DO i=NC,2*NC-1
            iAng(iOrb+i,1)      = 1
            iAng(iOrb+i+NC,2)   = 1
            iAng(iOrb+i+2*NC,3) = 1
          ENDDO
          iOrb = iOrb+4*NC

c  D type functions
         ELSEIF  (shell_type(iShellC).eq.2) THEN
!          WRITE (6,*) 'Cartesian d orbitals not implemented'
!          STOP
           ND = 6
           DO j = 0,ND-1
             NCvec(iCount+j) = NC
           ENDDO
           iCount=iCount+ND
           DO i = 1,ND*NC
              DO j = 1,3
                iAng(iOrb+i-1,j) = 0
                ROrb(iOrb+i-1,j) = R(k,j)
              ENDDO
           ENDDO
           DO i = 1,NC
c            Cartesian d-orbitals order in .fchk
             iAng(iOrb+i-1,1)      = 2  ! xx
             iAng(iOrb+i-1+NC,2)   = 2  ! yy
             iAng(iOrb+i-1+2*NC,3) = 2  ! zz
             iAng(iOrb+i-1+3*NC,1) = 1  ! xy
             iAng(iOrb+i-1+3*NC,2) = 1  ! xy
             iAng(iOrb+i-1+4*NC,1) = 1  ! xz
             iAng(iOrb+i-1+4*NC,3) = 1  ! xz
             iAng(iOrb+i-1+5*NC,2) = 1  ! yz
             iAng(iOrb+i-1+5*NC,3) = 1  ! yz
           ENDDO
!          NOrbSh = ND*NC
           iOrb = iOrb + ND*NC

         ELSEIF  (shell_type(iShellC).eq.-2) THEN
           NAOs=NAOs+5

           NCvec(iCount)   = NC
           NCvec(iCount+1) = NC
           NCvec(iCount+2) = NC
           NCvec(iCount+3) = NC
           NCvec(iCount+4) = NC
           iCount=iCount+5
           DO i = 1,5*NC
             DO j = 1,3
               iAng(iOrb+i-1,j) = 0
               ROrb(iOrb+i-1,j) = R(k,j)
             ENDDO
           ENDDO
           DO i = 1,NC
c Pure d-orbitals order in Q-Chem
            iAng(iOrb+i-1,3)      = 4  ! z2 from 2z2-x2-y2 (d0)
            iAng(iOrb+i-1,1)      = -2 ! x2 from 2z2-x2-y2 (d0)
            iAng(iOrb+i-1,2)      = -2 ! y2 from 2z2-x2-y2 (d0)
            iAng(iOrb+i-1+NC,1)   = 1  ! xz  (d1)
            iAng(iOrb+i-1+NC,3)   = 1  ! xz  (d1)
            iAng(iOrb+i-1+2*NC,2) = 1  ! yz  (d-1)
            iAng(iOrb+i-1+2*NC,3) = 1  ! yz  (d-1)
            iAng(iOrb+i-1+3*NC,1) = 2  ! x2 from x2-y2 (d2)
            iAng(iOrb+i-1+3*NC,2) = -2 ! y2 from x2-y2 (d2)
            iAng(iOrb+i-1+4*NC,1) = 1  ! xy  (d-2)
            iAng(iOrb+i-1+4*NC,2) = 1  ! xy  (d-2)
           ENDDO
           iOrb = iOrb + 5*NC ! pure d-orbitals
c          iOrb = iOrb + 6*NC ! cart d-orbitals

c          NOrb = NOrb + 5*NC ! pure d-orb
c          NOrb = NOrb + 6*NC ! cart d-orb
         ENDIF
        ENDDO
c-------------------------------------------
      ENDDO

      END
