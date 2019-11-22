
      SUBROUTINE make_S(SAB,iAngA,AlphA,CA,RA,
     & iAngB,AlphB,CB,RB,NBas,NBas2,NOrb,iContr,SM)
c*****************************************************c
c                                                     c
c       Computes the overlap matrices between two     c
c       sets of AOs                                   c
c                                                     c
c       INPUT                                         c
c        iAng_   : (x,y,z) exponent indices           c
c        Alph_   : function exponentials              c
c        C_      : function coefficients              c
c        R_      : function positions                 c
c        NOrb    : # uncontracted functions           c
c        NBas    : # basis functions                  c
c        iContr  : contraction vector                 c
c        SM      : symmetry transf matrix             c
c        NN      : max dimension of product (F*F)     c
c                                                     c
c       OUTPUT                                        c
c        SAB     : overlap matrix                     c
c                                                     c
c       DCC: 05/09                                    c
c       DCC: 04/10  WFN-Sym update                    c
c       DCC: 08/10  WFN-Sym F-orbs                    c
c                                                     c
c*****************************************************c
	  IMPLICIT REAL*8  (a-h, o-z)
      INTEGER  iAngA(NOrb,3),iAngB(NOrb,3),
     & iAngvA(3,3),iAngvAi(3),iAngvB(3),iAngBi(3),iContr(NBas)
      REAL*8   AlphA(NOrb),CA(NOrb),RA(NOrb,3),
     &         AlphB(NOrb),CB(NOrb),RB(NOrb,3),
     & AlphvA,CvA,RvA(3),AlphvB,RvB(3),SAB(NBas2),SM(3,3),LCCA(3)
      SAVE Zero, One, Two, Three
      DATA Zero,One,Two,Three / 0.d0,1.d0,2.d0,3.d0 /

      CALL VRload(SAB,NBas2,Zero)

      iA = 1
      DO iMOA = 1,NBas
      DO iCA = 1,iContr(iMOA)

      RvA = RA(iA,:)
      AlphvA = AlphA(iA)
      CvA    = CA(iA)
      CALL lc_cartesian_orbitals(iAngvA,LCCA,NLCA,iAngA(iA,:))
      DO na=1,NLCA
        iB = 1
        DO iMOB = 1,NBas
        iOff = iMOA+(iMOB-1)*NBas
          DO iCB = 1,iContr(iMOB)
            RvB = RB(iB,:)
            AlphvB = AlphB(iB)
            iAngBi = iAngB(iB,:)
            iAngvAi = iAngvA(na,:)
            CALL transform_B(SAB,NBas2,iAngvAi,AlphvA,LCCA(na)*CvA,RvA,
     &      AlphvB,RvB,iOff,iAngBi,NOrb,iB,CB,SM)
            iB = iB + 1
          ENDDO
        ENDDO
      ENDDO

      iA = iA + 1
      ENDDO
      ENDDO

      RETURN
      END

c     Symmetry transform B-orbital
      SUBROUTINE transform_B(SAB,NBas2,iAngvA,AlphvA,CvA,RvA,
     & AlphvB,RvB,iOff,iAngBi,NOrb,iB,CB,SM)
      IMPLICIT REAL*8  (a-h, o-z)
      INTEGER  iAngBi(3), iAngvA(3),iAngvB(3,3),iAng(10,3),iAngvBi(3),
     & iOff
      REAL*8   CB(NOrb),AlphvA,CvA,RvA(3),AlphvB,CvB,RvB(3),
     & SAB(NBas2),stCvB(10),SM(3,3),aux_stCvB(10), LCCB(3)
      SAVE Zero, One, Two, Three
      DATA Zero,One,Two,Three / 0.d0,1.d0,2.d0,3.d0 /

!     aux_stCvB(:) = 0.d0
      CALL VRload(aux_stCvB,10,Zero)
      CALL lc_cartesian_orbitals(iAngvB,LCCB,NLCB,iAngBi)
      DO nb=1,NLCB
        iAngvBi = iAngvB(nb,:)
        CALL transf_orbital(stCvB,iAng,NvB,iAngvBi,
     &  LCCB(nb)*CB(iB),SM)
        aux_stCvB = aux_stCvB + stCvB
      ENDDO

      DO ivB = 1,NvB
        sum = iAng(ivB,1)+iAng(ivB,2)+iAng(ivB,3)
        CvB = aux_stCvB(ivB)
        iAngvBi = iAng(ivB,:)
        CALL om_elements(SABi,27,iAngvA,
     &  AlphvA,CvA,RvA,iAngvBi,AlphvB,CvB,RvB)
        SAB(iOff) = SAB(iOff) + SABi
      ENDDO
      END

      SUBROUTINE om_elements(SABi,NN,iAngvA,
     & AlphvA,CvA,RvA,iAngvB,AlphvB,CvB,RvB)
      INTEGER  iAngAB(NN,3),iAngvA(3),iAngvB(3),iAngvAB(3)
      REAL*8   AlphvA,CvA,RvA(3),AlphvB,CvB,RvB(3),
     & AlphAB(NN),CAB(NN),RAB(NN,3),SABi,RvAB(3), AlphvAB,S
      SAVE Zero, One, Two, Three
      DATA Zero,One,Two,Three / 0.d0,1.d0,2.d0,3.d0 /

c Make the AO * AO product
      SABi = 0
      CALL aos_product(iAngAB,AlphAB,CAB,RAB,NAB,NN,
     &  iAngvA,AlphvA,CvA,RvA,iAngvB,AlphvB,CvB,RvB)
c Compute Ovelap Matrix elements
      DO iAB = 1, NAB
         RvAB(:)    = RAB(iAB,:)
         iAngvAB(:) = iAngAB(iAB,:)
         AlphvAB    = AlphAB(iAB)
         CALL overlap(iAngvAB,AlphvAB,RvAB,S)
         SABi = SABi + CAB(iAB)*S
      ENDDO
      END

      SUBROUTINE transf_orbital(stCvB,iAng,NvB,iAngBi,CBi,SM)
      IMPLICIT REAL*8  (a-h, o-z)
      INTEGER  iAngBi(3),iAng(10,3)
      REAL*8   CBi, stCvB(10),stCvB2(10),SAux(3,3),XAux(3,3),SM(3,3),
     & SFAux(27),SF(27)
      SAVE Zero, One, Two, Three
      DATA Zero,One,Two,Three / 0.d0,1.d0,2.d0,3.d0 /

      CoefD = 0.5d0
      iTotB = iAngBi(1)+iAngBi(2)+iAngBi(3)
      IF(iTotB .EQ. 0) THEN
      NvB = 1
      stCvB(1) = CBi
      iAng(1,:) = 0
      ELSEIF(iTotB .EQ. 1) THEN
      NvB = 3
      DO J = 1,3
          SAux(J,1) = CBi*dfloat(iAngBi(J))
      ENDDO
      CALL AtimsB(stCvB,SM,SAux, 3,1,3, 10,3,3, 1)
      iAng(1:3,:) = 0
      DO J = 1,3
        iAng(J,J) = 1
      ENDDO
      ELSEIF(iTotB .EQ. 2) THEN
      DO I = 1,3
         DO J = 1,I
          SAux(J,I) = CBi*dfloat(iAngBi(I)*iAngBi(J))
          IF(I.NE.J) SAux(J,I) = CoefD*SAux(J,I)
          SAux(I,J) = SAux(J,I)
         ENDDO
         IF(iAngBi(I) .EQ. 2) THEN
           SAux(I,I) = 0.25d0*SAux(I,I)
         ELSE
           SAux(I,I) = 0.d0
         ENDIF
      ENDDO
      CALL AtimsB(XAux,SM,SAux, 3,3,3, 3,3,3, 1)
      CALL AtimsB(SAux,XAux,SM, 3,3,3, 3,3,3, 3)
      ivB = 0
      DO I = 1,3
        DO J = 1,I
          ivB = ivB + 1
          stCvB(ivB) = SAux(J,I)
          IF(I.NE.J) stCvB(ivB) = stCvB(ivB) + SAux(I,J)
        ENDDO
   	  ENDDO
   	    NvB = 6
   	    iAng(:,:) = 0
        iAng(1,1) = 2  ! xx
        iAng(2,1) = 1  ! xy
        iAng(2,2) = 1  ! xy
        iAng(3,2) = 2  ! yy
        iAng(4,1) = 1  ! xz
        iAng(4,3) = 1  ! xz
        iAng(5,2) = 1  ! yz
        iAng(5,3) = 1  ! yz
        iAng(6,3) = 2  ! zz
      ELSEIF(iTotB .EQ. 3) THEN
        NvB = 10
        IJK = 6
        j1 = 1
        j2 = 2
        j3 = 3
      DO I = 1,3
      IF(iAngBi(I).EQ.3) THEN
 	    IJK = 1
 	    j1 = I
      ELSEIF(iAngBi(I).EQ.2) THEN
 	    IJK = 3
 	    j1 = I
 	    DO J = 1,3
 	    IF(iAngBi(J).EQ.1) j2 = J
 	    ENDDO
      ENDIF
      ENDDO
      CALL VRload(SFAux,27,Zero)
      IF(IJK.EQ.1) THEN
          iOff = j1 + 3*(j1-1+3*(j1-1))
          SFAux(iOff) = One
      ELSEIF(IJK.EQ.3) THEN
       CoefF = One / Three
       iOff = j1 + 3*(j1-1+3*(j2-1))
       SFAux(iOff) = CoefF
       iOff = j1 + 3*(j2-1+3*(j1-1))
       SFAux(iOff) = CoefF
       iOff = j2 + 3*(j1-1+3*(j1-1))
       SFAux(iOff) = CoefF
      ELSEIF(IJK.EQ.6) THEN
       CoefF = One / 6.d0
       iOff = j1 + 3*(j2-1+3*(j3-1))
       SFAux(iOff) = CoefF
       iOff = j1 + 3*(j3-1+3*(j2-1))
       SFAux(iOff) = CoefF
       iOff = j2 + 3*(j1-1+3*(j3-1))
       SFAux(iOff) = CoefF
       iOff = j3 + 3*(j1-1+3*(j2-1))
       SFAux(iOff) = CoefF
       iOff = j2 + 3*(j3-1+3*(j1-1))
       SFAux(iOff) = CoefF
       iOff = j3 + 3*(j2-1+3*(j1-1))
       SFAux(iOff) = CoefF
      ENDIF
      iOff = 1
      DO K = 1,3
 	    CALL AtimsB(SAux,SFAux(iOff),SM, 3,3,3, 3,3,3, 3)
        CALL AtimsB(SFAux(iOff),SM,SAux, 3,3,3, 3,3,3, 1)
        iOff = iOff + 9
      ENDDO
      CALL AtimsB(SF,SFAux,SM, 9,3,3, 9,9,3, 3)
      stCvB(1) = SF(1)
      stCvB(2) = SF(2) + SF(4)  + SF(10)
      stCvB(3) = SF(5) + SF(11) + SF(13)
      stCvB(4) = SF(14)
      stCvB(5) = SF(3)  + SF(7)  + SF(19)
      stCvB(6) = SF(6)  + SF(8)  + SF(12) +
     &           SF(16) + SF(20) + SF(22)
      stCvB(7) = SF(15) + SF(17) + SF(23)
      stCvB(8) = SF(9)  + SF(21) + SF(25)
      stCvB(9) = SF(18) + SF(24) + SF(26)
      stCvB(10)= SF(27)
      stCvB(:) = CBi * stCvB(:)

      iAng(:,:) = 0
      iAng(1,1) = 3  ! xxx
      iAng(2,1) = 2  ! xxy
      iAng(2,2) = 1  ! xxy
      iAng(3,1) = 1  ! xyy
      iAng(3,2) = 2  ! xyy
      iAng(4,2) = 3  ! yyy
      iAng(5,1) = 2  ! xxz
      iAng(5,3) = 1  ! xxz
      iAng(6,1) = 1  ! xyz
      iAng(6,2) = 1  ! xyz
      iAng(6,3) = 1  ! xyz
      iAng(7,2) = 2  ! yyz
      iAng(7,3) = 1  ! yyz
      iAng(8,1) = 1  ! xzz
      iAng(8,3) = 2  ! xzz
      iAng(9,2) = 1  ! yzz
      iAng(9,3) = 2  ! yzz
      iAng(10,3)= 3  ! zzz
      ENDIF
      END


      SUBROUTINE overlap_mat(S_AO, AtLab, RAt, NBas, NAt_2, NOrb,
     & NTotShell, NShell, shell_type, NPrim, COrb, Alph)
c*****************************************************c
c                                                     c
c       Computes the overlap matrix                   c
c                                                     c
c       INPUT                                         c
c******************************************************
      IMPLICIT REAL*8  (a-h, o-z)
      INTEGER :: iAng(NOrb,3),iContr(NBas),iAngv(3),iLAO(NBas)
      INTEGER :: NPrim(NTotShell), NShell(NAt_2), shell_type(NTotShell)
      REAL*8    :: RAt(NAt_2,3),ROrb(NOrb,3),COrb(NOrb),Alph(NOrb)
      REAL*8 :: MatId(3,3), S_AO(NBas**2)
      CHARACTER :: cLAO(NBas)*4
      CHARACTER :: AtLab(NAt_2)*2
      SAVE Zero,Half,One,Two,Four
      DATA Zero,Half,One,Two,Four / 0d0,0.5d0,1d0,2d0,4d0  /

      MatId(:,:) = Zero
      DO I = 1,3
        MatId(I,I) = One
      ENDDO

      CALL get_basis(AtLab,RAt,iAng,Alph,
     & ROrb,COrb,iContr,iLAO,cLAO,NAt_2,NBas,NOrb,
     & NPrim, NTotShell, NShell, shell_type)

      CALL make_S(S_AO,iAng,Alph,COrb,ROrb,iAng,Alph,COrb,ROrb,
     & NBas,NBas2,NOrb,iContr,MatId)

      END
