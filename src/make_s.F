
      SUBROUTINE make_S(SAB,iAngA,AlphA,CA,RA,
     & iAngB,AlphB,CB,RB,NBas,NBas2,NOrb,iContr,SM,PureCart)
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
c        PureCart: pure or cartesian d-orbs           c
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
	PARAMETER (NN = 27, NND = 6, NNF = 10)
      INTEGER  iAngA(NOrb,3),iAngB(NOrb,3),iAngAB(NN,3),
     & iAngvA(3),iAngvB(3),iContr(NBas),iAng(NNF,3),iAngvAB(3)
      REAL*8   AlphA(NOrb),CA(NOrb),RA(NOrb,3),
     &         AlphB(NOrb),CB(NOrb),RB(NOrb,3),
     & AlphvA,CvA,RvA(3),AlphvB,CvB,RvB(3),
     & AlphAB(NN),CAB(NN),RAB(NN,3),SAB(NBas2),
     & stCvB(NNF),SAux(3,3),XAux(3,3),SM(3,3),RvAB(3),
     & SFAux(NN),SF(NN)
	LOGICAL  PureCart
      SAVE Zero, One, Two, Three
      DATA Zero,One,Two,Three / 0.d0,1.d0,2.d0,3.d0 /

c Control
	IF(PureCart) THEN
	  WRITE(2,*)"ERROR. Pure F-orbitals not ready"
	  STOP
	ENDIF
	  
c	CoefD = 1.d0/DSQRT(2.d0)
 	CoefD = 0.5d0
      CALL VRload(SAB,NBas2,Zero)

	iA = 1
      DO 10 iMOA = 1,NBas
      DO 10 iCA = 1,iContr(iMOA)

      DO j = 1,3
        RvA(j)    = RA(iA,j)
        iAngvA(j) = iAngA(iA,j)
      ENDDO
      AlphvA = AlphA(iA)
      CvA    = CA(iA)
	
      iB = 1
      DO 20 iMOB = 1,NBas
      DO 20 iCB = 1,iContr(iMOB)

      RvB(:) = RB(iB,:)
      AlphvB = AlphB(iB)

c Symmetry transform B-orbital
      iTotB = iAngB(iB,1)+iAngB(iB,2)+iAngB(iB,3)
      IF(iTotB .EQ. 0) THEN
        NvB = 1
        stCvB(1) = CB(iB)
	  iAng(1,:) = 0
      ELSEIF(iTotB .EQ. 1) THEN
        NvB = 3
	  DO J = 1,3
          SAux(J,1) = CB(iB)*dfloat(iAngB(iB,J))
	  ENDDO
        CALL AtimsB(stCvB,SM,SAux, 3,1,3, NNF,3,3, 1)
	  iAng(1:3,:) = 0
	  DO J = 1,3
	    iAng(J,J) = 1
	  ENDDO
      ELSEIF(iTotB .EQ. 2) THEN
	  NvB = NND
	  DO I = 1,3
         DO J = 1,I
	    SAux(J,I) = CB(iB)*dfloat(iAngB(iB,I)*iAngB(iB,J))
          IF(I.NE.J) SAux(J,I) = CoefD*SAux(J,I)
          SAux(I,J) = SAux(J,I)
         ENDDO
	   IF(iAngB(iB,I) .EQ. 2) THEN
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
        NvB = NNF
	  IJK = 6
	  j1 = 1
	  j2 = 2
	  j3 = 3
	  DO I = 1,3
	    IF(iAngB(iB,I).EQ.3) THEN
		IJK = 1
		j1 = I
	    ELSEIF(iAngB(iB,I).EQ.2) THEN
		IJK = 3
	      j1 = I
		DO J = 1,3
		  IF(iAngB(iB,J).EQ.1) j2 = J
		ENDDO
	    ENDIF
	  ENDDO
	  CALL VRload(SFAux,NN,Zero)
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
     &             SF(16) + SF(20) + SF(22)
        stCvB(7) = SF(15) + SF(17) + SF(23)
        stCvB(8) = SF(9)  + SF(21) + SF(25)
        stCvB(9) = SF(18) + SF(24) + SF(26)
        stCvB(10)= SF(27)

        stCvB(:) = CB(iB) * stCvB(:)

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

c Repeat for each transformed B-orbital
	DO ivB = 1,NvB

	  iAngvB(:) = iAng(ivB,:) 
	  CvB = stCvB(ivB) 
c Make the AO * AO product
	  CALL aos_product(iAngAB,AlphAB,CAB,RAB,NAB,NN,
     & iAngvA,AlphvA,CvA,RvA,
     & iAngvB,AlphvB,CvB,RvB,
     & PureCart)
c Compute Ovelap Matrix elements
        DO iAB = 1,NAB
          RvAB(:)    = RAB(iAB,:)
          iAngvAB(:) = iAngAB(iAB,:)
          AlphvAB    = AlphAB(iAB)
          CALL overlap(iAngvAB,AlphvAB,RvAB,S)
  	    iOff = iMOA+(iMOB-1)*NBas
          SAB(iOff) = SAB(iOff) + CAB(iAB)*S
	  ENDDO
	ENDDO  ! End of transf B-orbitals (ivB loop)

      iB = iB + 1
 20   CONTINUE

      iA = iA + 1
 10   CONTINUE

      RETURN
      END


       SUBROUTINE overlap_mat(S_AO, AtLab, RAt, NBas, NAt_2, NOrb,
     & PureCart, NTotShell, NShell, shell_type, NPrim, COrb, Alph)
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
      LOGICAL   :: PureCart
      SAVE Zero,Half,One,Two,Four
      DATA Zero,Half,One,Two,Four / 0d0,0.5d0,1d0,2d0,4d0  /

      MatId(:,:) = Zero
      DO I = 1,3
        MatId(I,I) = One
      ENDDO

      CALL get_basis(AtLab,RAt,iAng,Alph,
     & ROrb,COrb,iContr,iLAO,cLAO,NAt_2,NBas,NOrb,PureCart,
     & NPrim, NTotShell, NShell, shell_type)

      CALL make_S(S_AO,iAng,Alph,COrb,ROrb,iAng,Alph,COrb,ROrb,
     & NBas,NBas2,NOrb,iContr,MatId,PureCart)

      END
