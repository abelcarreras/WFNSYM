#include "wfnsym.h"

c      PROGRAM Prog_Main
      SUBROUTINE mainlib(NEtot, NEval, NBas, NMo, Norb,
     & NAt_2, NTotShell, iZAt, AtLab, Alph, COrb, NShell,
     & RAt, NPrim, shell_type, igroup_2, ngroup_2, Ca, Cb,
     & RCread_2, VAxis_2, VAxis2_2, iMult, DoOp,
     & OutDim, OutGrim, OutCSM, OutSymLab,
     & OutSDiagA, OutSDiagB, OutWDiagA, OutWDiagB, OutTbl,
     & OutIRLab, OutgIRA, OutgIRB,OutgIRwfA,OutgIRwfB,
     & OutgIRwf, OutSymMat)

c--------------------------------------------------------------
c       David Casanova          06/APR/10	Initializing program
c       David Casanova          28/JUN/10	f-orbitals (starting)
c       David Casanova          22/JUL/10	Print level
c       Abel Carreras            8/SEP/18   Create library
c--------------------------------------------------------------
c	Variables
c
c	unit1		Input. Funcio d'ona obtinguda amb Gaussian
c	unit2		Output. Fitxer general de l'output amb la integracio de la densitat.
c	unit3		Output. Funcio densitat que ens servira pel calcul del solapament
c
c	NAt		# atoms
c	NEtot 	# total electrons
c	NEval 	# valence electrons
c	NBas        # Contracted basis functions
c	NOrb	 	# Uncontracted basis functions
c	COrb	 	Uncontracted AO coefficients            :: (NOrb)
c	ROrb	 	Uncontracted AO centers                 :: (NOrb,3)
c	stROrb	Sym Trans Uncontracted AO centers       :: (NOrb,3)
c	Alph	 	Uncontracted exponents                  :: (NOrb)
c	CA	 	MO coefficients                         :: (NBas2)
c     iContr      AO contraction vector                   :: (NBas)
c     iLAO        AO l vector (value)                     :: (NBas)
c     cLAO        AO l vector (name)                      :: (NBas)
c     RAt         Atomic coordinates                      :: (NAt,3)
c     stRAt       Symmetry Transformed Atomic coordinates :: (NAt,3)
c     S_AO        AO Overlap matrix                       :: (NBas2)
c     S_MO        MO Overlap matrix                       :: (NBas2)
c     iGroup      Group index (identifier)
c     hGroup      # Total symmetry operations
c     dGroup      # Non-redundant operations
c     nGroup      Order of the group
c     nIR	      # Irreducible representations
c     NPrint      Print level  0: Default
c                              1: Print MOs
c                              2: Print S_AO and S_MO matrix
c	HighPrec	Results print level
c			             False: 3 decimals (default)
c					 True : 5 decimals
c--------------------------------------------------------------
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER :: iZAt(NAt_2), NPrim(NTotShell)
      CHARACTER :: AtLab(NAt_2)*2
      REAL*8 :: Alph(Norb),  COrb(Norb), RAt(NAt_2, 3)
      REAL*8 :: CA(NMo*NBas),  CB(NMo*NBas)
      REAL*8 :: RCread_2(3), VAxis_2(3), VAxis2_2(3)
      LOGICAL :: DoOp
c     Output
      INTEGER :: OutDim(3)
      REAL*8 :: OutCSM(100), OutGrim(100)
      CHARACTER :: OutSymLab(100)*6, OutIRLab(100)*3
      REAL*8 :: OutSDiagA(NMo, 100), OutSDiagB(NMo, 100)
      REAL*8 :: OutWDiagA(100), OutWDiagB(100)
      REAL*8 :: OutTbl(100,100)
      REAL*8 :: OutgIRA(NMo,100), OutgIRB(NMo,100)
      REAL*8 :: OutgIRwfA(100), OutgIRwfB(100), OutgIRwf(100)
      REAL*8 :: OutSymMat(100, 3, 3)
c Include files
      INCLUDE 'definitionlib_var.ins'
      INCLUDE 'common_var.ins'
      INCLUDE 'init_var.ins'

c ----------- Added manual arguments-----------
      iout = 0
      igroup = igroup_2
      ngroup = ngroup_2
      CALL SetGroupDim(iGroup, nGroup, dGroup, hGroup)

      NAt = NAt_2
      RCread = RCread_2
      VAxis = VAxis_2
      VAxis2 = VAxis2_2
      basis = 'STO-3G'

      basin = .TRUE.
      molecule = .TRUE.
      natom = .TRUE.
      Unrestricted = .TRUE.
      readRC = .TRUE.
      AXIS = .TRUE.
      AXIS2 = .TRUE.
      Group = 'C'
      SubG = '6V'

      IF (DoOp) THEN
        DoOper = .FALSE.
      ELSE
        DoGroup = .TRUE.
      ENDIF

c -----------------------------------------------

c Check angles in VAxis & VAxis
      IF(.not.AngCHK) THEN
        IF (AXIS) THEN
          Dot = VAxis(1)**2 + VAxis(2)**2 + VAxis(3)**2
          Dot = 1.d0 / DSQRT(Dot)
          VAxis(1) = Dot * VAxis(1)
          VAxis(2) = Dot * VAxis(2)
          VAxis(3) = Dot * VAxis(3)
          AngCHK = .true.
        ELSEIF (ANGLES) THEN
          VAxis(1) = VAxis(1) * PI / 180.d0
          VAxis(2) = VAxis(2) * PI / 180.d0
          AngCHK = .true.
        ENDIF
      ENDIF
      IF(.not.AngCHK2) THEN
        IF (AXIS2) THEN
          Dot = VAxis2(1)**2 + VAxis2(2)**2 + VAxis2(3)**2
          Dot = 1.d0 / DSQRT(Dot)
          VAxis2(:) = Dot * VAxis2(:)
          AngCHK2 = .true.
        ELSEIF (ANGLES2) THEN
          VAxis2(1:2) = VAxis2(1:2) * PI / 180.d0
          AngCHK2 = .true.
        ENDIF
      ENDIF

      ALLOCATE(stRAt(NAt,3))

c Define Axes thresh
      IF(.NOT.Ax_Thresh) IAx_Thresh = 9


c Set up some variables
      NOa = (NEtot + iMult - 1)/2
      NOb = NOa - iMult + 1
      Ncore = (NEtot-NEval)/2
      NOva = NOa - Ncore
      NOvb = NOb - Ncore
      NBas2 = NBas*NBas
      IF(NOa+NOb .NE. NEtot) THEN
        WRITE(iout, *) 'ERROR. Charge/Multiplicity inconsistency'
        STOP
      ENDIF


c Where do we center the molecule?
      IF (NO_CENTER) THEN
c        WRITE(iout,*)'Symmetry analysis at the geometric center'
          RC(:) = Zero
          DO k = 1,NAt
            RC(:) = RC(:) + RAt(k,:)
          ENDDO
          RC(:) = RC(:)/dfloat(NAt)
      ELSEIF (readRC) THEN
c        WRITE(iout,*) 'Symmetry analysis at the input center'
        DO J = 1,3
          RC(J) = ANGSBOHR * RCread(J)
        ENDDO
      ELSEIF (readRN) THEN
c        WRITE(iout,"('Symmetry analysis at the geometric center of:',
c     &        1X,I2,1X,'atoms')") Ncenter
        RC(:) = Zero
        DO k = 1,Ncenter
          RC(:) = RC(:) + RAt(Icenter(k),:)
        ENDDO
        RC(:) = RC(:)/dfloat(Ncenter)
      ELSEIF (readRD) THEN
c        WRITE(iout,*)'Symmetry analysis at the electron density center'
          RC(:) = -RDipole(:)
          DO I = 1,NAt
            RC(:) = RC(:) + dfloat(iZAt(I))*RAt(I,:)
          ENDDO
      ENDIF
      DO J = 1,3
        DO I = 1,NAt
          RAt(I,J) = RAt(I,J) - RC(J)
        ENDDO
      ENDDO

c Define VAxis for Axis_2Atoms and Axis_3atoms cases
      IF (AxVec) THEN
        VAxis(:) = RAt(IAxis(2),:) - RAt(IAxis(1),:)
        Dot = VAxis(1)**2 + VAxis(2)**2 + VAxis(3)**2
        Dot = 1.d0 / DSQRT(Dot)
        VAxis(:) = Dot * VAxis(:)
        AXIS = .TRUE.
      ELSEIF (AxPlane) THEN
        R1(:) = RAt(IAxis(2),:) - RAt(IAxis(1),:)
        R2(:) = RAt(IAxis(3),:) - RAt(IAxis(1),:)
        VAxis(1) = R1(2)*R2(3) - R1(3)*R2(2)
        VAxis(2) = R1(3)*R2(1) - R1(1)*R2(3)
        VAxis(3) = R1(1)*R2(2) - R1(2)*R2(1)
        Dot = VAxis(1)**2 + VAxis(2)**2 + VAxis(3)**2
        Dot = 1.d0 / DSQRT(Dot)
        VAxis(:) = Dot * VAxis(:)
        AXIS = .TRUE.
      ENDIF
c Do the same for secondary axis
      IF (AxVec2) THEN
        VAxis2(:) = RAt(IAxis2(2),:) - RAt(IAxis2(1),:)
        Dot = VAxis2(1)**2 + VAxis2(2)**2 + VAxis2(3)**2
        Dot = 1.d0 / DSQRT(Dot)
        VAxis2(:) = Dot * VAxis2(:)
        AXIS2 = .TRUE.
      ELSEIF (AxPlane2) THEN
        R1(:) = RAt(IAxis2(2),:) - RAt(IAxis2(1),:)
        R2(:) = RAt(IAxis2(3),:) - RAt(IAxis2(1),:)
        VAxis2(1) = R1(2)*R2(3) - R1(3)*R2(2)
        VAxis2(2) = R1(3)*R2(1) - R1(1)*R2(3)
        VAxis2(3) = R1(1)*R2(2) - R1(2)*R2(1)
        Dot = VAxis2(1)**2 + VAxis2(2)**2 + VAxis2(3)**2
        Dot = 1.d0 / DSQRT(Dot)
        VAxis2(:) = Dot * VAxis2(:)
        AXIS2 = .TRUE.
      ENDIF

c Define Symmetry-Operation or Symmetry-Group Matrix
      IF(DoOper) THEN
        ALLOCATE(SymGrp(0,0,0))
        CALL sym_oper(SymMat,iOper,iRot,AXIS,ANGLES,VAxis,iout)
        dGroup = 1
        hGroup = 1
      ELSEIF(DoGroup) THEN
        ALLOCATE(SymGrp(3,3,dGroup))
        CALL sym_group(SymGrp,iGroup,nGroup,dGroup,VAxis,VAxis2,iout)
      ENDIF

c Make operation labels of the group
      ALLOCATE(SymLab(dGroup))
      ALLOCATE(iSymMu(dGroup))
      IF(DoOper) THEN
        IF(iOper.EQ.1) THEN
          SymLab(1) = ' Inv. '
        ELSEIF(iOper.EQ.2) THEN
          SymLab(1) = ' Ref. '
        ELSEIF(iOper.EQ.3) THEN
          line(1:6) = repeat(' ',6)
          WRITE(line(1:2),'(I2)') iRot
          line = adjustl(line)
          SymLab(1) = ' C'//line(1:2)//'  '
        ELSEIF(iOper.EQ.4) THEN
          line(1:6) = repeat(' ',6)
          WRITE(line(1:2),'(I2)') iRot
          line = adjustl(line)
          SymLab(1) = ' S'//line(1:2)//'  '
        ENDIF
        iSymMu(1) = 1
      ELSEIF(DoGroup) THEN
        CALL SymLabels(SymLab)
        CALL SymMultip(iSymMu)
        iTest = 0
        DO iG = 1,dGroup
          ITest = ITest + iSymMu(iG)
        ENDDO
        IF(ITest .NE. hGroup) THEN
          WRITE(iout,*)'ERROR. Operation multiplicity inconsistency'
          STOP
        ENDIF
      ENDIF

c Allocate some space
c     write(iout,*)'allocating space'
      NN = NBas*hGroup
      ALLOCATE(SDiagA(NN),SDiagB(NN))
      ALLOCATE(WDiagA(hGroup),WDiagB(hGroup),WDiag(hGroup))
      ALLOCATE(iAng(NOrb,3),iLAO(NBas),cLAO(NBas))
      ALLOCATE(ROrb(NOrb,3),stROrb(NOrb,3))
      ALLOCATE(iContr(NBas))
c     ALLOCATE(IpivA(NOa),IpivB(NOb))
      ALLOCATE(IpivA(NBas),IpivB(NBas))

c Reading AO Basis Set
c     write(iout,*)'reading AOs'
      CALL get_basis(RAt,iAng,Alph,
     & ROrb,COrb,iContr,iLAO,cLAO,NAt,NBas,NOrb,
     & NPrim, NTotShell, NShell, shell_type)

c Compute AO Overlap Matrix
c     write(iout,*)'computing AO overlap'
      ALLOCATE(S_AO(NBas2),S_MOa(NBas2),S_MOb(NBas2),Xaux(NBas2))
      MatId(:,:) = Zero
      DO I = 1,3
        MatId(I,I) = One
      ENDDO
      CALL make_S(S_AO,iAng,Alph,COrb,ROrb,iAng,Alph,COrb,ROrb,
     & NBas,NBas2,NOrb,iContr,MatId)

#ifdef DEVELOPMENT
c Print AO Overlap Matrix
      line = 'AO Overlap Matrix'
      CALL PrintMatrix(S_AO,cLAO,.TRUE.,.TRUE.,NBas,line,iout)
#endif


c Allocate space for MOs
c      ALLOCATE(CA(NBas2))
      IF(Unrestricted .OR. UHFNO) THEN
c        ALLOCATE(CB(NBas2))
        IF(Unrestricted) THEN
c          WRITE(iout,*)'Unrestricted Wave Function'
        ELSE
          ALLOCATE(WOcc(NBas))
c          WRITE(iout,*)'UHF Natural Orbitals'
        ENDIF
c      ELSE
c        WRITE(iout,*)'Restricted Wave Function'
      ENDIF


c Set group matrices (if necessary)
      IF(DoSymOrb) THEN
        ALLOCATE(SGrpA(NOa*NOa))
        ALLOCATE(SGrpB(0))
        CALL VRload(SGrpA,NOa*NOa,Zero)
c       ALLOCATE(SGrpA(NBAs*NBas))
c       CALL VRload(SGrpA,NBas*NBas,Zero)
        IF(NOb.NE.0 .AND.(NOa.NE.NOb .OR. Unrestricted))THEN
          ALLOCATE(SGrpB(NOb*NOb))
          CALL VRload(SGrpB,NOb*NOb,Zero)
          DO I = 1,NOb
            SGrpB(I+(I-1)*NOb) = One
          ENDDO
        ENDIF
      ENDIF

cgroups
      iOp = 1
      OPERATIONS: DO iG = 1,dGroup

c      WRITE(iout,"('@@@ Operation ',I2,':',1X,A6)") iG,SymLab(iG)
      IF(DoGroup) SymMat(:,:) = SymGrp(:,:,iG)
c      WRITE(iout,*)' '
c      WRITE(iout,*)'Symmetry Transformation matrix'
c      DO i = 1,3
c        WRITE(iout,902) (SymMat(i,j),j=1,3)
c      ENDDO
      DO i = 1,3
        DO j = 1,3
          OutSymMat(iG, i, j) = SymMat(i,j)
        ENDDO
      ENDDO

      CALL AtimsB(stRAt,RAt,SymMat, NAt,3,3, NAt,NAt,3, 3)

c Symmetry Transform ROrb to stROrb
      CALL AtimsB(stROrb,ROrb,SymMat, NOrb,3,3, NOrb,NOrb,3, 3)

c Compute Symmetry Overlap Matrix
c	write(iout,*) 'calling make_S'
      CALL make_S(S_AO,iAng,Alph,COrb,ROrb,iAng,Alph,COrb,stROrb,
     & NBas,NBas2,NOrb,iContr,SymMat)

c Print Symmetry Overlap Matrix
      IF (NPrint .GE. 2) THEN
        line = 'Sym AO Overlap Matrix'
        CALL PrintMatrix(S_AO,cLAO,.TRUE.,.TRUE.,NBas,line,iout)
c Print AOs Symmetry Measures
        WRITE(iout,*)' '
        WRITE(iout,*)'AOs: SOEVs '
        WRITE(iout,*)'---------- '
        DO i = 1,NBas
          WRITE(iout,'(I3,1X,A4,2X,F8.4)') i,cLAO(i),
     &          S_AO(i+(i-1)*NBas)
        ENDDO
        WRITE(iout,*)' '
      ENDIF

c Transform Overlap Matrix to MOs
      CALL AtimsB(Xaux,S_AO,CA, NBas,NBas,NBas, NBas,NBas,NBas, 1)
      CALL AtimsB(S_MOa,CA,Xaux, NBas,NBas,NBas, NBas,NBas,NBas, 2)
      IF(Unrestricted) THEN
        CALL AtimsB(Xaux,S_AO,CB, NBas,NBas,NBas, NBas,NBas,NBas, 1)
        CALL AtimsB(S_MOb,CB,Xaux, NBas,NBas,NBas, NBas,NBas,NBas, 2)
      ELSE
        CALL VRcopy(S_MOb,S_MOa,NBas2)
      ENDIF

c Store the NOx,NOx block of S_MOx in SGrpX (if required)
      IF(DoSymOrb) THEN
        DO I = 0,NOa-1
        DO J = 1,NOa
          SGrpA(J+I*NOa) = SGrpA(J+I*NOa)+S_MOa(J+I*NBas)
        ENDDO
        ENDDO
        CALL writeMat(S_MOa,NBas,iout)
c       CALL VRadd2(SGrpA,S_MOa,NBas2)
        IF(NOb.NE.0 .AND.(NOa.NE.NOb .OR. Unrestricted))THEN
          DO I = 0,NOb-1
          DO J = 1,NOb
            SGrpB(J+I*NOb) = SGrpB(J+I*NOb)+S_MOb(J+I*NBas)
          ENDDO
          ENDDO
        ENDIF
      ENDIF


c Print Sym trans MO Overlap Matrix
      IF (NPrint .GE. 2) THEN
        WRITE(iout,"('sym-trans oper ',I2,':',1X,A6)") iG,SymLab(iG)
        line = 'Sym alpha-MO Overlap Matrix'
        CALL PrintMatrix(S_MOa,cLAO,.FALSE.,.FALSE.,NBas,line,iout)
        IF(Unrestricted) THEN
          line = 'Sym beta-MO Overlap Matrix'
          CALL PrintMatrix(S_MOb,cLAO,.FALSE.,.FALSE.,NBas,line,iout)
        ENDIF
      ENDIF

c Copy the diagonal of S_MO_ to SDiag_
      DO IJ = iOp,iOp+iSymMu(iG)-1
        DO i = 1,NBas
          SDiagA(i+(IJ-1)*NBas) = S_MOa(i+(i-1)*NBas)
          SDiagB(i+(IJ-1)*NBas) = S_MOb(i+(i-1)*NBas)
        ENDDO
      ENDDO

c Print Alpha MOs S value
      IF (NPrint .GE. 2) THEN
        WRITE(iout,*)' '
        WRITE(iout,*)'MOs: SOEVs  alpha  /  beta'
        WRITE(iout,*)'--------------------------'
        DO i = 1,NBas
          WRITE(iout,'(I3,7X,F8.4,2X,F8.4)') i,
     &    S_MOa(i+(i-1)*NBas),S_MOb(i+(i-1)*NBas)
        ENDDO
        WRITE(iout,*)' '
      ENDIF

c Compute the wave function symmetry
      CALL DGETRF(NOa,NOa,S_MOa,NBas,IpivA,Info)
      IF((Unrestricted.OR.(NOa.NE.NOb)) .AND. NOb.GT.0) THEN
        CALL DGETRF(NOb,NOb,S_MOb,NBas,IpivB,Info)
      ELSE
        CALL VRcopy(S_MOb,S_MOa,NBas2)
        IpivB(1:NOa) = IpivA(1:NOa)
      ENDIF
      S_WFa = One
      S_WFb = One
      DO I = 1,NOa
        IF(IpivA(I).EQ.I) THEN
          S_WFa = S_WFa * S_MOa(I+(I-1)*NBas)
        ELSE
          S_WFa = - S_WFa * S_MOa(I+(I-1)*NBas)
        ENDIF
      ENDDO
      IF(NOb.GT.0) THEN
        DO I = 1,NOb
          IF(IpivB(I).EQ.I) THEN
            S_WFb = S_WFb * S_MOb(I+(I-1)*NBas)
          ELSE
            S_WFb = - S_WFb * S_MOb(I+(I-1)*NBas)
          ENDIF
        ENDDO
      ELSE
        S_WFb = Zero
      ENDIF
      DO IJ = iOp,iOp+iSymMu(iG)-1
        WDiagA(IJ) = S_WFa
        WDiagB(IJ) = S_WFb
        IF(NOb.GT.0) THEN
          WDiag(IJ) = S_WFa * S_WFb
        ELSE
          WDiag(IJ) = S_WFa
        ENDIF
      ENDDO

cgroups
        iOp = iOp + iSymMu(iG)
      ENDDO OPERATIONS

      DEALLOCATE(IpivA,IpivB)

c Define level of printing results
c      FMT1 = "('--------')"
c      FMT2 = "(1X,F7.3)"
c      FMT3 = "(2X,A6)"

      IF(DoGroup) THEN
c Build Group Table
      CALL table_dimensions(nIR,iGroup,nGroup,dGroup)
      ALLOCATE(GrpTbl(hGroup,hGroup))
      ALLOCATE(IRLab(nIR),IRMul(nIR))
      CALL group_labels(IRLab,nIR,iGroup,nGroup,iout)
      CALL group_mult(IRMul,nIR,iGroup,nGroup,iout)
      CALL group_table(GrpTbl,nIR,iGroup,nGroup,hGroup,iout)

c Print Group Table
c      WRITE(iout,"(1X,A1,A3,': Ideal Group Table')") Group,SubG
c      WRITE(iout,'(3X)',advance='no')
c      DO i = 1,dGroup
c        WRITE(iout,FMT1,advance='no')
c      ENDDO
c      WRITE(iout,*)
c      WRITE(iout,'(4X)',advance='no')
c      DO i = 1,dGroup
c        WRITE(iout,FMT3,advance='no') SymLab(i)
c      ENDDO
c      WRITE(iout,*)
c      WRITE(iout,'(3X)',advance='no')
c      DO i = 1,dGroup
c        WRITE(iout,FMT1,advance='no')
c      ENDDO
c      WRITE(iout,*)

      jIR = 1
      DO iIR = 1,nIR
c        WRITE(iout,'(A3)',advance='no') IRLab(iIR)
        iOp = 1
        DO iG = 1,dGroup
          IF(IRMul(iIR) .LE. 2) THEN
            Tbl = GrpTbl(jIR,iOp)
          ELSEIF(IRMul(iIR) .EQ. 4) THEN
            Tbl = GrpTbl(jIR,iOp)+GrpTbl(jIR+3,iOp)
          ELSEIF(IRMul(iIR) .EQ. 9) THEN
            Tbl = GrpTbl(jIR,iOp)+GrpTbl(jIR+4,iOp)+GrpTbl(jIR+8,iOp)
          ENDIF
c          WRITE(iout,FMT2,advance='no') Tbl
          OutTbl(iIR, iG) = Tbl
          iOp = iOp + iSymMu(iG)
        ENDDO
c        WRITE(iout,*)
        jIR = jIR + IRMul(iIR)
      ENDDO
c      WRITE(iout,'(3X)',advance='no')
c      DO i = 1,dGroup
c        WRITE(iout,FMT1,advance='no')
c      ENDDO
c      WRITE(iout,*)

      ENDIF ! DoGroup

c*********************c
c Grimme's CSM values c
c*********************c
c Compute Grimme's CSM
      ALLOCATE(csm(hGroup))
      DO iG = 1,hGroup
        csm(iG) = Zero
        DO i = 1,NOa
          csm(iG) = csm(iG) + One-DABS(SDiagA(i+(iG-1)*NBas))
        ENDDO
        DO i = 1,NOb
          csm(iG) = csm(iG) + One-DABS(SDiagB(i+(iG-1)*NBas))
        ENDDO
        csm(iG) = 100.d0 * csm(iG)/dfloat(NOa+NOb)
      ENDDO

c     Set output data
      iOp = 1
      do i = 1, dGroup
        OutGrim(i) = csm(iOp)
        iOp = iOp + iSymMu(i)
      enddo

      iOp = 1
      do i = 1, dGroup
        OutCSM(i) = 100.*(One-WDiag(iOp))
        iOp = iOp + iSymMu(i)
      enddo

      do i = 1, dGroup
         OutSymLab(i) = SymLab(i)
      enddo

      DO i = 1,NMo
        iOp = 1
        DO iG = 1,dGroup
          OutSDiagA(i,iG) = SDiagA(i+(iOp-1)*NBas)
          OutSDiagB(i,iG) = SDiagB(i+(iOp-1)*NBas)
          iOp = iOp + iSymMu(iG)
        ENDDO
      ENDDO

      iOp = 1
      DO iG = 1,dGroup
        OutWDiagA(iG) = WDiagA(iOp)
        OutWDiagB(iG) = WDiagB(iOp)
        iOp = iOp + iSymMu(iG)
      ENDDO

      DO i = 1,nIR
        OutIRLab(i) = IRLab(i)
      ENDDO


      DEALLOCATE(csm)

      IF(DoGroup) THEN

c******************************************c
c Compute Irred. Repr. linear combinations c
c******************************************c
continue
c NEW: MODIFY PROJECTIONS
c Project over Irred. Repr. of the group
      ALLOCATE(gIRA(NBas,nIR),gIRwfA(nIR),gIRwf(nIR))
      ALLOCATE(gIRB(NBas,nIR),gIRwfB(nIR))

c All projections
c     DO iOrb = 1,NBas
c write(iout,"('Projections Orbital: ',I3)") iOrb
c DO iG  = 1,hGroup
c       Xaux(iG) = Zero
c       DO iOp = 0,hGroup-1
c        Xaux(iG) = Xaux(iG)+SDiagA(iOrb+iOp*NBas)*GrpTbl(iG,iOp+1)
c  ENDDO
c  write(iout,"(1x,f9.5)") Xaux(iG)
c ENDDO
c     ENDDO

      DO iOrb = 1,NBas
        ij = 1
        DO iIR = 1,nIR
          ProjA = Zero
          ProjB = Zero
          IF(IRMul(iIR).LE.2) THEN
            DO iG = 0,hGroup-1
              ProjA = ProjA+SDiagA(iOrb+iG*NBas)*GrpTbl(ij,iG+1)
              ProjB = ProjB+SDiagB(iOrb+iG*NBas)*GrpTbl(ij,iG+1)
            ENDDO
            Scale = DFloat(IRMul(iIR)) / DFloat(hGroup)
          ELSEIF(IRMul(iIR).EQ.4) THEN
            DO iG = 0,hGroup-1
              ProjA = ProjA+SDiagA(iOrb+iG*NBas)*GrpTbl(ij,iG+1)
     &           +SDiagA(iOrb+iG*NBas)*GrpTbl(ij+3,iG+1)
              ProjB = ProjB+SDiagB(iOrb+iG*NBas)*GrpTbl(ij,iG+1)
     &                  +SDiagB(iOrb+iG*NBas)*GrpTbl(ij+3,iG+1)
            ENDDO
            Scale = Two / DFloat(hGroup)
          ELSEIF(IRMul(iIR).EQ.9) THEN
            DO iG = 0,hGroup-1
              ProjA = ProjA+SDiagA(iOrb+iG*NBas)*GrpTbl(ij,iG+1)
     &                  +SDiagA(iOrb+iG*NBas)*GrpTbl(ij+4,iG+1)
     &                  +SDiagA(iOrb+iG*NBas)*GrpTbl(ij+8,iG+1)
              ProjB = ProjB+SDiagB(iOrb+iG*NBas)*GrpTbl(ij,iG+1)
     &                  +SDiagB(iOrb+iG*NBas)*GrpTbl(ij+4,iG+1)
     &                  +SDiagB(iOrb+iG*NBas)*GrpTbl(ij+8,iG+1)
            ENDDO
            Scale = 3.d0 / DFloat(hGroup)
          ELSE
            WRITE(iout,*)"ERROR. Wrong IR's multiplicity"
            STOP
          ENDIF
          gIRA(iOrb,iIR) = ProjA * Scale
          gIRB(iOrb,iIR) = ProjB * Scale
          ij = ij + IRMul(iIR)
        ENDDO
      ENDDO

      ij = 1
      DO iIR = 1,nIR
        ProjA = Zero
        ProjB = Zero
        Proj  = Zero
        IF(IRMul(iIR).LE.2) THEN
          DO iG = 1,hGroup
           ProjA = ProjA+WDiagA(iG)*GrpTbl(ij,iG)
           ProjB = ProjB+WDiagB(iG)*GrpTbl(ij,iG)
           Proj  = Proj + WDiag(iG)*GrpTbl(ij,iG)
          ENDDO
          Scale = DFloat(IRMul(iIR)) / DFloat(hGroup)
        ELSEIF(IRMul(iIR).EQ.4) THEN
          DO iG = 1,hGroup
           ProjA = ProjA+WDiagA(iG)*GrpTbl(ij,iG)
     &                  +WDiagA(iG)*GrpTbl(ij+3,iG)
           ProjB = ProjB+WDiagB(iG)*GrpTbl(ij,iG)
     &                  +WDiagB(iG)*GrpTbl(ij+3,iG)
           Proj  = Proj +WDiag(iG)*GrpTbl(ij,iG)
     &                  +WDiag(iG)*GrpTbl(ij+3,iG)
          ENDDO
          Scale = Two / DFloat(hGroup)
        ELSEIF(IRMul(iIR).EQ.9) THEN
          DO iG = 1,hGroup
           ProjA = ProjA+WDiagA(iG)*GrpTbl(ij,iG)
     &                  +WDiagA(iG)*GrpTbl(ij+4,iG)
     &                  +WDiagA(iG)*GrpTbl(ij+8,iG)
           ProjB = ProjB+WDiagB(iG)*GrpTbl(ij,iG)
     &                  +WDiagB(iG)*GrpTbl(ij+4,iG)
     &                  +WDiagB(iG)*GrpTbl(ij+8,iG)
           Proj  = Proj +WDiag(iG)*GrpTbl(ij,iG)
     &                  +WDiag(iG)*GrpTbl(ij+4,iG)
     &                  +WDiag(iG)*GrpTbl(ij+8,iG)
          ENDDO
          Scale = 3.d0 / DFloat(hGroup)
        ELSE
          WRITE(iout,*)"ERROR. Wrong IR's multiplicity"
          STOP
        ENDIF
        gIRwfA(iIR) = ProjA * Scale
        gIRwfB(iIR) = ProjB * Scale
        gIRwf(iIR)  = Proj  * Scale
        ij = ij + IRMul(iIR)
      ENDDO

      ENDIF ! DoGroup

      DO i = 1,NMo
        DO iIR = 1,nIR
c          OutIRLab(i) = IRLab(i)
          OutgIRB(i, iIR)  =  gIRB(i,iIR)
          OutgIRA(i, iIR)  =  gIRA(i,iIR)
        ENDDO
      ENDDO

      DO i = 1,nIR
        OutIRLab(i) = IRLab(i)
        OutgIRwfA(i)  =  gIRwfA(i)
        OutgIRwfB(i)  =  gIRwfB(i)
        OutgIRwf(i)  =  gIRwf(i)
      ENDDO


c ################################
c #  Diagonilze Symmetry Matrix  #
c ################################
      IF(DoSymOrb) THEN

c Allocate some space
      lwork = 3*NOa
c     lwork = 3*NBas
      ALLOCATE(work(lwork))
      ALLOCATE(WOrbA(NOa))
c     ALLOCATE(WOrbA(NBas))

c Scale SGrpA
      Scale = One / DFLOAT(dGroup)
cScale = (EXP(One)-One) / DFLOAT(dGroup)
      CALL VRscale(SGrpA,NOa*NOa,Scale)
c     DO I = 1,NOa
c       SGrpA(I+(I-1)*NOa) = One + SGrpA(I+(I-1)*NOa)
c     ENDDO
c	CALL VRscale(SGrpA,NBas2,Scale)

c Compute the wavefunction symmetry (only alpha)
      CALL VRtrace(Sgroup,SGrpA,NOa)
c      write(iout,*) 'Trace = ',Sgroup
      Sgroup = (exp(Sgroup) - One) / (exp(One) - One)
c      write(iout,*) '< G > = ',Sgroup
      Sgroup = 100.*(One-Sgroup)
c      WRITE(iout,"('**************')")
c      WRITE(iout,"('* S = ',F6.2,' *')") Sgroup
c      WRITE(iout,"('**************')")

c Print SGrpA
      write(iout,*)'SGrpA Matrix'
      CALL writeMat(SGrpA,NOa,iout)
c	CALL writeMat(SGrpA,NBas,iout)

c Diagonalize it
      CALL DSYEV('V','U',NOa,SGrpA,NOa,WOrbA,work,lwork,Info)
c     CALL DSYEV('V','U',NBas,SGrpA,NBas,WOrbA,work,lwork,Info)
      IF(NOb.NE.0 .AND.(NOa.NE.NOb .OR. Unrestricted))THEN
        ALLOCATE(WOrbB(NOb))
        CALL DSYEV('V','U',NOb,SGrpB,NOb,WOrbB,work,lwork,Info)
      ENDIF
      DEALLOCATE(work)

c Transform SGrpX(NOx,NOx) to SOrbX(NBas,NOx)
      ALLOCATE(SOrbA(NBas*NOa))
      CALL AtimsB(SOrbA,CA,SGrpA, NBas,NOa,NOa, NBas,NBas,NOa, 1)

c	ALLOCATE(SOrbA(NBas*NBas))
c	CALL AtimsB(SOrbA,CA,SGrpA, NBas,NBas,NBas, NBas,NBas,NBas, 1)

      DEALLOCATE(SGrpA)
      IF(NOb.NE.0 .AND.(NOa.NE.NOb .OR. Unrestricted))THEN
        ALLOCATE(SOrbB(NBas*NOb))
        CALL AtimsB(SOrbB,CB,SGrpB, NBas,NOb,NOb, NBas,NBas,NOb, 1)
        DEALLOCATE(SGrpB)
      ENDIF

c Print SymOrbs as a Molden file (*.mol)
c      WRITE(iout,*) 'Writing Molden format'
c      CALL make_molden(NAt,iZAt,AtLab,RAt,Basis,input,
c     &     NOa,NBas,SOrbA,WOrbA,PureCart)
      DEALLOCATE(SOrbA,WOrbA)

      IF(NOb.NE.0 .AND.(NOa.NE.NOb .OR. Unrestricted))
     &   DEALLOCATE(SOrbB,WOrbB)

      ENDIF
Cc###################################


 500  FORMAT('     Timmings of WFN-Sym (in seconds)')
 501  FORMAT(' Normalization             :',F9.2)
 502  FORMAT(' Build density             :',F9.2)
 503  FORMAT(' Integrate density         :',F9.2)
 504  FORMAT(' Symmetry transformation   :',F9.2)
 505  FORMAT(' Overlap densities         :',F9.2)

c--------------------------------------------------------------
 100  FORMAT('--------------------------------------------')

 901  FORMAT(13x,i4,16x)
 902  FORMAT(10X,3(1X,F12.8))
 997  FORMAT('Symmetry Overlap  : ',F9.5,3X,F9.5)
 998  FORMAT('Symmetry Measure  : ',F9.5,3X,F9.5)
 999  FORMAT('Density Charge Center (Angs): ',3(f12.8,1x))

 801  FORMAT(2x,i6,3x,3(i2,1x),2x,f15.8,1x,f15.8,3(2x,f9.6))
 802  FORMAT(3x,'N    Tipus',5x,'Exponent',6x,'Coeficient',
     .         17x,'Coordenades')
 803  FORMAT(1x,f9.6,3x,3(1x,f9.6))
 804  FORMAT(i7,3x,3(i1,1x),2(f18.8,3x),3(f12.8))

 790  FORMAT('  *** Wavefunction dimensions ***')
 791  FORMAT('# Total   Electrons            ',I4)
 792  FORMAT('# Core    Electrons            ',I4)
 793  FORMAT('# Valence Electrons            ',I4)
 794  FORMAT('# Occ total   alpha/beta MOs   ',I4,I4)
 795  FORMAT('# Occ valence alpha/beta MOs   ',I4,I4)
 796  FORMAT('# Basis Functions              ',I4)
 797  FORMAT('# Uncontracted Atomic Orbitals ',I4)

 694  FORMAT('*** Density function dimensions ***')
 699  FORMAT('Total gaussians in the density    ',I10)
 399  FORMAT('Non zero gaussians in the density ',I10)
 398  FORMAT('Non zero gauss symm trans density ',I10)
 698  FORMAT('Integrated electron density       ',F15.9)
 697  FORMAT(2X,'      SELF-SIMILARITY:',3X,F15.9)
 696  FORMAT(2X,' INVERSION-SIMILARITY:',3X,F15.9)
 393  FORMAT(2X,' INVERSION C-INDEX   :',3X,F15.9)
 397  FORMAT(2X,' INVERSION MEASURE   :',3X,F15.9)
 695  FORMAT(2X,'REFLECTION-SIMILARITY:',3X,F15.9)
 392  FORMAT(2X,'REFLECTION C-INDEX   :',3X,F15.9)
 396  FORMAT(2X,'REFLECTION MEASURE   :',3X,F15.9)
 395  FORMAT(2X,'  ROTATION MEASURE   :',3X,F15.9)
 391  FORMAT(2X,'  ROTATION C-INDEX   :',3X,F15.9)
 394  FORMAT('Density building CutOff :     1.d-',I2)

c--------------------------------------------------------------

c Free some space
      DEALLOCATE(iAng,iContr)
      DEALLOCATE(S_AO,S_MOa,S_MOb,Xaux)
c      DEALLOCATE(CA)
c      IF(Unrestricted) DEALLOCATE(CB)
      DEALLOCATE(iLAO,cLAO)
      DEALLOCATE(ROrb,stROrb)
      DEALLOCATE(stRAt,Icenter)
      DEALLOCATE(SymGrp,SymLab,iSymMu)
      DEALLOCATE(SDiagA,SDiagB,WDiagA,WDiagB,WDiag)
      IF(DoGroup) THEN
        DEALLOCATE(GrpTbl,IRLab,IRMul)
        DEALLOCATE(gIRA,gIRwfA,gIRwf)
        DEALLOCATE(gIRB,gIRwfB)
      ENDIF

c Saying goodbye
c      WRITE(iout,*)' '
c      WRITE(iout,*)' '
c      WRITE(iout,499)
c      WRITE(iout,*)' '
c      WRITE(iout,498)
c      WRITE(iout,*)' '

 499  FORMAT('    Calculation has finished normally')
 498  FORMAT('  *** Thank you for using Wave-Sym ***')

c Close files
c      CLOSE(1)
c      CLOSE(2)

c     Generate outputs
      OutDim(1) = dGroup
      OutDim(2) = hGroup
      OutDim(3) = nIR

      END


      REAL*8 FUNCTION TimeSec(H,Min,Sec,mSec)

      INTEGER  H, Min, Sec, mSec

      TimeSec = float(3600*H + 60*Min + Sec) +
     &          float(mSec)/1000.

      RETURN
      END
