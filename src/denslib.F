#include "qsym.h"

      SUBROUTINE denslib(AtLab,R,N_elec,LDens,alfa2,coef2,NPrim,NShell,
     &                   shell_type,coef_mol,RCread_2,VAxis_2,VAxis2_2,
     &                   NMOs,NAOs,NOrbX,igroup_2,ngroup_2,DoOper_2,
     &                   NOrb,NAt_2,NTotShell,OutDim,OutCSM_i,OutCSM,
     &                   OutSymLab,OutSassembly)

c--------------------------------------------------------------
c
c   Possible Basis Sets: STO-3G, 6-31G, 4-31G, 6-311G, 3-21G, 6-31Gd
c--------------------------------------------------------------
c       David Casanova          18/MAY/09	Reading Q-Chem MOs (53.0)
c       David Casanova          13/OCT/11	Small modifications
c--------------------------------------------------------------
c	Variables
c
c	NAt		num d'atoms
c	Num_elec	num d'electrons
c	NAOs		num d'orbitals atomics segons la base
c	NOrb	 	num pure      gaussian functions
c	NOrbX	 	num cartesian gaussian functions
c	N_1		orbital molecular que no agafem pel calcul
c   N_2             orbital molecular que no agafem pel calcul
c   vecontb(i)      vector de contraccio de la base. Compte el num de termes de cada capa
c	coef_MO(i,j)    coeficients moleculars. La seva dimensio es (NOrb,Num_elec/2)
c	coef_mol(i)	coeficients moleculars en forma de vector

c   R           Atomic coordinates
c--------------------------------------------------------------
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER :: NOrb,NOrbX,N_elec,NAt_2,igroup_2,ngroup_2
      CHARACTER :: AtLab(NAt_2)*2,str*2
      REAL*8 :: R(NAt_2,3),RT(NAt_2,3),alfa(NOrbX),coef(NOrbX),
     &          alfa2(NOrb),coef2(NOrb),
     &          R_orb(NOrbX,3),CDens(LDens),AlphDens(LDens),
     &          RDens(LDens,3),coef_mol(NMOs*NAOs),MatSym(3,3),
     &          Aux1(3,3),Aux2(3,3),R_fun(3),RCread_2(3),VAxis_2(3),
     &          VAxis2_2(3)
      INTEGER :: i_angl(NOrbX,3),vecontb(NOrb),iAngDens(LDens,3),
     &           NPrim(NTotShell),in(3),shell_type(NTotShell),
     &           NShell(NAt_2)
      INTEGER, ALLOCATABLE :: iAngSym(:,:)
      REAL*8, ALLOCATABLE :: AlphSym(:),CSym(:),RSym(:,:),
     &                       coef_MO(:,:),coeficient(:),S_sym(:)
      LOGICAL :: DoOper_2
c     Output
      INTEGER :: OutDim(3)
      REAL*8  ::  OutCSM_i(100),OutCSM,OutSassembly
      CHARACTER :: OutSymLab(100)*6

c Include files
      INCLUDE 'definitionlib_var.ins'
      INCLUDE 'common_var.ins'
      INCLUDE 'init_var.ins'

c ----------- Added manual arguments-----------
      iout = 0
      igroup = igroup_2
      ngroup = ngroup_2
      NAt = NAt_2
      RCread = RCread_2
      VAxis = VAxis_2
      VAxis2 = VAxis2_2
      TolDens = 1.d-8
      alfa = alfa2
      DoOper = DoOper_2


c Writing Coordinates and transform them to AU unities
      DO k = 1,NAt
        DO j = 1,3
          R(k,j) = R(k,j)*ANGSBOHR
        ENDDO
      ENDDO

c Reading Basis Set
 10   CONTINUE
      CALL dens_basis(i_angl,R_orb,R,vecontb,NAt,
     &                NOrbX,NPrim,NTotShell,NShell,shell_type,NOrb)
      coef = coef2


c Place MOs in Matrix form
      ALLOCATE(coef_MO(NOrbX,N_elec/2))

      ij = 1
      DO i = 1,NMOs
	    iOrb = 0
	    DO j = 1,NAOs
	      DO jj = 1,vecontb(j)
	        coef_MO(iOrb+jj,i) = coef_mol(ij)
	      ENDDO
	      iOrb = iOrb+vecontb(j)
	      ij = ij+1
	    ENDDO
	  ENDDO

c
c Expand x^2-y^2 and z2 only for pure d-orbitals
c
c Introduction of y^2 <- x^2-y^2
c     write(2,*)'expand y2'
      ij = NOrb
      DO i = 1,NOrb
        IF (i_angl(i,1).eq.2 .and. i_angl(i,2).eq.-2) THEN
          alfa(ij+1) = alfa(i)
          coef(ij+1) = -coef(i)*0.8660254038!*DSQRT(0.75D0)
          coef(i) = coef(i)*0.8660254038!*DSQRT(0.75D0)
 	      DO j = 1,NMOs
 	        coef_MO(ij+1,j) = coef_MO(i,j)
	      ENDDO
	      DO j = 1,3
	        R_orb(ij+1,j) = R_orb(i,j)
	      ENDDO
          i_angl(i,1) = 2
          i_angl(i,2) = 0
          i_angl(i,3) = 0
          i_angl(ij+1,1) = 0
          i_angl(ij+1,2) = 2
          i_angl(ij+1,3) = 0
          ij = ij+1
	    END IF


c Introduction of  z^2  =  2z^2-x^2-y^2
c     write(2,*)'expand z2'
        IF (i_angl(i,3).eq.4) THEN
          alfa(ij+1) = alfa(i)
          coef(ij+1) = -coef(i)*0.5D0
          alfa(ij+2) = alfa(i)
          coef(ij+2) = -coef(i)*0.5D0
          coef(i)    = 2.0d0*coef(i)*0.5D0
          DO j = 1,NMOs
            coef_MO(ij+1,j) = coef_MO(i,j)
            coef_MO(ij+2,j) = coef_MO(i,j)
          ENDDO
          DO j = 1,3
            R_orb(ij+1,j) = R_orb(i,j)
            R_orb(ij+2,j) = R_orb(i,j)
          ENDDO
          i_angl(i,1) = 0
          i_angl(i,2) = 0
          i_angl(i,3) = 2
          i_angl(ij+1,1) = 2
          i_angl(ij+1,2) = 0
          i_angl(ij+1,3) = 0
          i_angl(ij+2,1) = 0
          i_angl(ij+2,2) = 2
          i_angl(ij+2,3) = 0
          ij = ij+2

        END IF
      ENDDO

c Compute Normalization Factors
	  DO iOrb = 1,NOrbX
	    DO j = 1,3
            in(j) = i_angl(iOrb,j)
        ENDDO
	    CALL norma(in,alfa(iOrb),Dot)
	    coef(iOrb) = coef(iOrb)*Dot
	  ENDDO

c Product between MO coefficients. Only take needed orbitals
      Nx = (NOrbX)*(NOrbX+1)/2
      ALLOCATE(coeficient(Nx))
      CALL VRload(coeficient,Nx,Zero)

      ij = 1
      DO i_MO = 1,NMOs
       IF(.NOT.OrbCHK) THEN
         k = 1
         DO i1 = 1,NOrbX
           DO i2 = i1,NOrbX
             coeficient(k) = coeficient(k)+
     &        coef_MO(i1,i_MO)*coef_MO(i2,i_MO)
           k = k + 1
           ENDDO
         ENDDO
         ij = ij + 1
       ENDIF
      ENDDO


c Free some space
      DEALLOCATE(coef_MO)
c-------------------------------------------
c  Build density function
      CALL VRload(iAngDens,LDens,Zero)
      CALL VRload(CDens,LDens,Zero)
      CALL VRload(AlphDens,LDens,Zero)
      CALL VRload(RDens,LDens,Zero)
      CALL make_dens(iAngDens,AlphDens,CDens,RDens,LDens,NDens,
     &               alfa,coeficient,coef,i_angl,R_orb,NOrbX,TolDens)
      DEALLOCATE(coeficient)

c Integrate the electron density
      elec_total = Zero
      CALL electronic_integral(elec_total,RC,AlphDens,iAngDens,RDens,
     &                         CDens,LDens,NDens,N_elec)

c Do we wish to center the electron density?
      IF (NO_CENTER) THEN
        CALL VRload(RC,3,Zero)
      ELSEIF (readRC .AND. .NOT. CenCHK) THEN
        DO J = 1,3
          RC(J) = ANGSBOHR * RCread(J)
        ENDDO 
      ENDIF

      Tol = 1.d-8
      IF ((ABS(RC(1)).GT.Tol .OR. ABS(RC(2)).GT.Tol .OR.
     &     ABS(RC(3)).GT.Tol) .AND. .NOT.CenCHK) THEN
        DO J = 1,3
          DO I = 1,NAt
            R(I,J) = R(I,J) - RC(J)
          ENDDO
        ENDDO
 497    FORMAT('Remaking the density at the center...')
        CenCHK = .TRUE.
        GOTO 10
      ENDIF

c**********************************************************c
c  In this version we do NOT move the electron density    *c
c  to its charge center. Instead, we recompute the        *c
c  electron density at RC.                                *c
c**********************************************************c
c  Traslate the electron density to its charge center     *c
c  (The dimension of the density matrices can increase)   *c
c                                                         *c
c     CALL center_dens(ang_cen,alfa_cen,coef_cen,R_cen,   *c
c    & N_cen,RC,iAngDens,AlphDens,CDens,RDens,NDens)      *c
c                                                         *c
c Compute Self and Inversion similarities                 *c
c     CALL self_inversion(S_total,S_inv,                  *c
c    & ang_cen,alfa_cen,coef_cen,R_cen,N_cen)             *c
c**********************************************************c

c compute Self_similarity
      CALL self_similarity(S_total,
     &       iAngDens,AlphDens,CDens,RDens,LDens,NDens)

      CALL SetGroupDim(iGroup, nGroup, dGroup, hGroup)
c Define Symmetry-Operation or Symmetry-Group Matrix
      IF(DoOper) THEN
        dGroup = 1
        hGroup = 1
        ALLOCATE(SymGrp(3,3,dGroup))
!       CALL sym_oper(SymMat,iOper,iRot,AXIS,ANGLES,VAxis,iout)
        CALL sym_oper(SymGrp,iGroup,nGroup,AXIS,ANGLES,VAxis,iout)

      ELSE
        ALLOCATE(SymGrp(3,3,dGroup))
        CALL sym_group(SymGrp,iGroup,nGroup,dGroup,VAxis,VAxis2,iout)
      ENDIF
      ALLOCATE(S_sym(dGroup))

      IF(DoOper) THEN
        IF(iGroup.EQ.1) THEN
          OutSymLab(1) = ' Inv. '
        ELSEIF(iGroup.EQ.2) THEN
          OutSymLab(1) = ' Ref. '
        ELSEIF(iGroup.EQ.3) THEN
          write(str,'(I1)') nGroup
          OutSymLab(1) = ' C'//str//'  '
        ELSEIF(iGroup.EQ.4) THEN
          write(str,'(I1)') nGroup
          OutSymLab(1) = ' S'//str//'  '
        ENDIF
      ELSE
        CALL SymLabels(OutSymLab,iGroup,nGroup,dGroup)
      ENDIF



      LSym = 0
      DO isym = 1,NDens
        isum = iAngDens(isym,1) + iAngDens(isym,2) + iAngDens(isym,3)
        IF (isum .EQ. 0) LSym = LSym + 1
        IF (isum .EQ. 1) LSym = LSym + 3
        IF (isum .EQ. 2) LSym = LSym + 6
        IF (isum .EQ. 3) LSym = LSym + 10
        IF (isum .EQ. 4) LSym = LSym + 15
      ENDDO

      S_sym_total = Zero
      DO iG=1,dGroup
         MatSym(:,:) = SymGrp(:,:,iG)
c Write transformation matrix and new atomic coordinates
         CALL AtimsB(RT,R,MatSym, NAt,3,3, NAt,NAt,3, 1)
         ALLOCATE(iAngSym(LSym,3))
         ALLOCATE(AlphSym(LSym),CSym(LSym),RSym(LSym,3))
         CALL sym_transform(iAngSym,AlphSym,CSym,RSym,LSym,NSym,
     &         MatSym,iAngDens,AlphDens,CDens,RDens,LDens,NDens)
         CALL sym_overlap(iAngDens,AlphDens,CDens,RDens,NDens,LDens,
     &         iAngSym,AlphSym,CSym,RSym,NSym,LSym,S_sym(iG))
         S_sym_total = S_sym_total + S_sym(iG)
!        OutCSM_i(iG) = 100.d0*(One-S_sym(iG)/(dGroup*S_total))
         OutCSM_i(iG) = S_sym(iG)/S_total

         DEALLOCATE(iAngSym)
         DEALLOCATE(AlphSym,CSym,RSym)
      END DO

      OutCSM = 100.d0*(One-S_sym_total/(dGroup*S_total))
!     OutCSM = 100.d0*(One-S_sym_total/(hGroup*S_total))
      OutDim(1) = dGroup
      OutDim(2) = hGroup
      OutDim(3) = nIR
      OutSassembly = S_total

      END

      SUBROUTINE electronic_integral(elec_total,RC,AlphDens,iAngDens,
     &                               RDens,CDens,LDens,NDens,N_elec)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER :: N_elec,ang2(3),i_angl_fun(3)
      REAL*8 :: CDens(LDens),AlphDens(LDens),RDens(LDens,3),
     &          R_fun(3),RC(3)
      INTEGER :: iAngDens(LDens,3)
      SAVE Zero
      DATA Zero / 0d0 /

      DO i = 1,3
        RC(i) = Zero
      ENDDO
      DO i_fun = 1,NDens
        alfa_fun = AlphDens(i_fun)
        i_angl_fun = iAngDens(i_fun,:)
        R_fun = RDens(i_fun,:)
        CALL overlap(i_angl_fun,alfa_fun,R_fun,elec)
        elec_total = elec_total+elec*CDens(i_fun)

        DO k = 1,3
           DO j = 1,3
             ang2(j) = i_angl_fun(j)
           ENDDO
           ang2(k) = ang2(k) + 1
           CALL overlap(ang2,alfa_fun,R_fun,elec)
           RC(k) = RC(k)+elec*CDens(i_fun)/dfloat(N_elec)
        ENDDO
      ENDDO

      END SUBROUTINE


      SUBROUTINE center_charge(R,N_elec,LDens,alfa2,coef2,
     &                   NPrim,NShell,shell_type,coef_mol,
     &                   NMOs,NAOs,NOrbX,NOrb,
     &                   NAt_2,NTotShell,OutRC)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER :: NOrb,NOrbX,N_elec,NAt_2
      REAL*8 :: R(NAt_2,3),alfa(NOrbX),coef(NOrbX),
     &          R_orb(NOrbX,3),CDens(LDens),AlphDens(LDens),
     &          RDens(LDens,3),coef_mol(NMOs*NAOs),
     &          alfa2(NOrb),coef2(NOrb)
      INTEGER :: i_angl(NOrbX,3),vecontb(NOrb),iAngDens(LDens,3),
     &           NPrim(NTotShell),in(3),shell_type(NTotShell),
     &           NShell(NAt_2)
      REAL*8, ALLOCATABLE :: coef_MO(:,:),coeficient(:)
c     Output
      REAL*8 :: OutRC(3)

      INCLUDE 'definitionlib_var.ins'
      INCLUDE 'common_var.ins'
      INCLUDE 'init_var.ins'

c ----------- Added manual arguments-----------
      iout = 0
      TolDens = 1.d-8
      NAt = NAt_2
      alfa = alfa2
      coef = coef2

c Writing Coordinates and transform them to AU unities
      DO k = 1,NAt
        DO j = 1,3
          R(k,j) = R(k,j)*ANGSBOHR
        ENDDO
      ENDDO

c Reading Basis Set
      CALL dens_basis(i_angl,R_orb,R,vecontb,NAt,
     &                NOrbX,NPrim,NTotShell,NShell,shell_type,NOrb)


c Place MOs in Matrix form
      ALLOCATE(coef_MO(NOrbX,N_elec/2))

      ij = 1
      DO i = 1,NMOs
	    iOrb = 0
	    DO j = 1,NAOs
	      DO jj = 1,vecontb(j)
	        coef_MO(iOrb+jj,i) = coef_mol(ij)
	      ENDDO
	      iOrb = iOrb+vecontb(j)
	      ij = ij+1
	    ENDDO
	  ENDDO

c
c Expand x^2-y^2 and z2 only for pure d-orbitals
c
c Introduction of y^2 <- x^2-y^2
c     write(2,*)'expand y2'
      ij = NOrb
      DO i = 1,NOrb
        IF (i_angl(i,1).eq.2 .and. i_angl(i,2).eq.-2) THEN
          alfa(ij+1) = alfa(i)
          coef(ij+1) = -coef(i)*0.8660254038
          coef(i) = coef(i)*0.8660254038
 	      DO j = 1,NMOs
 	        coef_MO(ij+1,j) = coef_MO(i,j)
	      ENDDO
	      DO j = 1,3
	        R_orb(ij+1,j) = R_orb(i,j)
	      ENDDO
          i_angl(i,1) = 2
          i_angl(i,2) = 0
          i_angl(i,3) = 0
          i_angl(ij+1,1) = 0
          i_angl(ij+1,2) = 2
          i_angl(ij+1,3) = 0
          ij = ij+1
	    END IF


c Introduction of  z^2  =  2z^2-x^2-y^2
c     write(2,*)'expand z2'
        IF (i_angl(i,3).eq.4) THEN
          alfa(ij+1) = alfa(i)
          coef(ij+1) = -coef(i)*0.5D0
          alfa(ij+2) = alfa(i)
          coef(ij+2) = -coef(i)*0.5D0
          coef(i)    = 2.0d0*coef(i)*0.5D0
          DO j = 1,NMOs
            coef_MO(ij+1,j) = coef_MO(i,j)
            coef_MO(ij+2,j) = coef_MO(i,j)
          ENDDO
          DO j = 1,3
            R_orb(ij+1,j) = R_orb(i,j)
            R_orb(ij+2,j) = R_orb(i,j)
          ENDDO
          i_angl(i,1) = 0
          i_angl(i,2) = 0
          i_angl(i,3) = 2
          i_angl(ij+1,1) = 2
          i_angl(ij+1,2) = 0
          i_angl(ij+1,3) = 0
          i_angl(ij+2,1) = 0
          i_angl(ij+2,2) = 2
          i_angl(ij+2,3) = 0
          ij = ij+2

        END IF
      ENDDO

c Compute Normalization Factors
	  DO iOrb = 1,NOrbX
	    DO j = 1,3
            in(j) = i_angl(iOrb,j)
        ENDDO
	    CALL norma(in,alfa(iOrb),Dot)
	    coef(iOrb) = coef(iOrb)*Dot
	  ENDDO

c Product between MO coefficients. Only take needed orbitals
      Nx = (NOrbX)*(NOrbX+1)/2
      ALLOCATE(coeficient(Nx))
      CALL VRload(coeficient,Nx,Zero)

      ij = 1
      DO i_MO = 1,NMOs
         k = 1
         DO i1 = 1,NOrbX
           DO i2 = i1,NOrbX
             coeficient(k) = coeficient(k)+
     &        coef_MO(i1,i_MO)*coef_MO(i2,i_MO)
           k = k + 1
           ENDDO
         ENDDO
         ij = ij + 1
      ENDDO


c Free some space
      DEALLOCATE(coef_MO)
c-------------------------------------------
c  Build density function
      CALL VRload(iAngDens,LDens,Zero)
      CALL VRload(CDens,LDens,Zero)
      CALL VRload(AlphDens,LDens,Zero)
      CALL VRload(RDens,LDens,Zero)
      CALL make_dens(iAngDens,AlphDens,CDens,RDens,LDens,NDens,
     &               alfa,coeficient,coef,i_angl,R_orb,NOrbX,TolDens)

      DEALLOCATE(coeficient)

c Integrate the electron density
      elec_total = Zero
      CALL electronic_integral(elec_total,OutRC,AlphDens,iAngDens,
     &                         RDens,CDens,LDens,NDens,N_elec)
      DO i=1,3
      OutRC(i) = OutRC(i)/ANGSBOHR
      END DO

      END SUBROUTINE
