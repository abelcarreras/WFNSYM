
	SUBROUTINE PrintHead(iout)
c
c Prints program header to iout unit
c
      WRITE(iout,*)' ' 
      WRITE(iout,*)'  ******************************************'
      WRITE(iout,*)'  *                                        *'
      WRITE(iout,*)'  *  Welcome to the Wave-Sym Program v1.0  *'
      WRITE(iout,*)'  *                                        *'
      WRITE(iout,*)'  ******************************************'
      WRITE(iout,*)' '
      WRITE(iout,*)'  Symmetry analysis of the electronic wave function'
      WRITE(iout,*)'  by David Casanova, University of Barcelona'
      WRITE(iout,*)' '

	RETURN
	END

      SUBROUTINE PrintOperGrp()
c
c Print available symmetry operations and groups
c

c Header
	CALL PrintHead(6)

c OPERATIONS
      WRITE(6,*)' '
      WRITE(6,*)'---------------------------------'
      WRITE(6,*)' Symmetry Operations (%operation)'
      WRITE(6,*)'---------------------------------'
      WRITE(6,*)'* Identity '
      WRITE(6,*)'  E'
      WRITE(6,*)' '
      WRITE(6,*)'* Inversion '
      WRITE(6,*)'  I'
      WRITE(6,*)' '
      WRITE(6,*)'* Reflection'
      WRITE(6,*)'  R'
      WRITE(6,*)' '
      WRITE(6,*)'* Rotation'
      WRITE(6,*)'  C n       n: rotation order'
      WRITE(6,*)' '
      WRITE(6,*)'* Rotation-Reflection'
      WRITE(6,*)'  S n       n: rotation-reflection order'
      WRITE(6,*)' '

c GROUPS
      WRITE(6,*)'---------------------------------'
      WRITE(6,*)'     Symmetry Groups (%group)'
      WRITE(6,*)'---------------------------------'
      WRITE(6,*)'* Non-Axial '
      WRITE(6,*)'  Cs, Ci'
      WRITE(6,*)' ' 
      WRITE(6,*)'* Single Axis Groups'
      WRITE(6,*)'  Cn, Cnv, Cnh, S2n'
      WRITE(6,*)' ' 
      WRITE(6,*)'* Dihedral Groups'
      WRITE(6,*)'  Dn, Dnd, Dnh'
      WRITE(6,*)' ' 
      WRITE(6,*)'* Cubic Groups'
      WRITE(6,*)'  Td, Oh'
      WRITE(6,*)' ' 

	RETURN
	END


      SUBROUTINE PrintKeys()
c
c Print available keywords
c

c Header
      CALL PrintHead(6)

c Keywords
      WRITE(6,*)'------------------'
      WRITE(6,*)'   KEYWORD LIST   '
      WRITE(6,*)'------------------'
      WRITE(6,*)' '
      WRITE(6,"(' %NAtom               Number of atoms')")
      WRITE(6,"(' %Basis               Atomic basis')")
      WRITE(6,"(' %Operation           Symmetry opertaion')")
      WRITE(6,"(' %Group               Symmetry point group')")
      WRITE(6,"(' %Unrestricted        Unrestricted orbitals')")
      WRITE(6,"(' %Molecule            Molecular coordinate input')")
      WRITE(6,"(' %Center_Oper         Center of symmetry operation')")
      WRITE(6,"(' %Center_Oper_NAtoms  # atoms to define the center')")
      WRITE(6,"(' %No_Center           Center at the origin ')")
      WRITE(6,"(' %Dipole              Electron density center')")
      WRITE(6,"(' %Axis                Principal symmetry axis')")
      WRITE(6,"(' %Axis2               Secondary symmetry axis')")
      WRITE(6,"(' %Angles              Principal symmetry axis')")
      WRITE(6,"(' %Angles2             Secondary symmetry axis')")
      WRITE(6,"(' %Axis_2Atoms         1st axis defined by 2 atoms')")
      WRITE(6,"(' %Axis_3Atoms         1st axis defined by 3 atoms')")
      WRITE(6,"(' %Axis2_2Atoms        2nd axis defined by 2 atoms')")
      WRITE(6,"(' %Axis2_3Atoms        2nd axis defined by 3 atoms')")
      WRITE(6,"(' %AxesThresh          Axis/Axis2 control thresh')")
      WRITE(6,"(' %Print_Level         Print level ')")
      WRITE(6,"(' %HighPrec            Results print level')")
      write(6,*)' '

	RETURN
	END


      SUBROUTINE PrintMatrix(M,cLab,LabC,LabR,N,Title,iUnit)
c Prints Square matrices with labels

      REAL*8      M(*)
      CHARACTER   cLab(N)*4,Title*80
      LOGICAL     LabC,LabR

      inc = 4
      WRITE(iUnit,*) ' '
      WRITE(iUnit,'(A)') Title

      N0 = 1
 100  NF = N0 + inc - 1
      IF(N .LT. NF) NF = N
      WRITE(iUnit,'(15X)',advance='no')
      DO i = N0,NF
        WRITE(iUnit,'(I3,11X)',advance='no') i
      ENDDO
      WRITE(iUnit,'(1X)')
      IF(LabR) THEN
        WRITE(iUnit,'(17X)',advance='no')
        DO i = N0,NF
          WRITE(iUnit,'(A4,10X)',advance='no') cLab(i)
        ENDDO
        WRITE(iUnit,'(7X)')
      ENDIF
      DO i = 1,N
        IF(LabC) THEN
          WRITE(iUnit,'(I3,1X,A4)',advance='no') i,cLab(i)
        ELSE
          WRITE(iUnit,'(I3,1X,A4)',advance='no') i
        ENDIF
        DO j = N0-1,NF-1
          WRITE(iUnit,'(2X,F12.6)',advance='no') M(i+j*N)
        ENDDO
        WRITE(iUnit,'(1X)')
      ENDDO

      N0 = N0 + inc
      IF(N0 .LE. N) GOTO 100

      WRITE(iUnit,*)' '

      RETURN
      END

