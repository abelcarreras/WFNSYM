
          
      SUBROUTINE options(opt,iunit)
c        
c Actives different options depending the keywords
c        
      CHARACTER*80  :: opt0, opt, line, Basis
      CHARACTER     :: Operation*1, Group*1, SubG*3
      REAL*8        :: RCread(3),VAxis(3),VAxis2(3),RDipole(3)
      INTEGER       :: IAxis(3),IAxis2(3),hGroup,dGroup,nGroup
      LOGICAL       :: DO_INV,DO_REF,DO_ROT

c Include files
      INCLUDE 'common_var.ins'

      opt0=opt
      k=index(opt0,' ')
      opt0=opt0(1:k)//repeat(' ',80-k)

c Transform lowercase to uppercase
      do i=2,len_trim(opt0)
        k=ichar(opt0(i:i))
        if (k.ge.97.and.k.le.122) then
          opt0(i:i)=char(k-32)
        endif
      enddo


      select case (opt0)

        case('%NATOM')
          natom = .true.
          k = len_trim(opt0)
          opt0 = repeat(' ',k)//opt(k+1:80)
          read(opt0,*) NAt 

        case('%BASIS')
          basin = .true.
          k = len_trim(opt0)
          line = repeat(' ',k)//opt(k+1:80)
          Basis = adjustl(line)
         
        case('%OPERATION')
          DoOper = .TRUE.
          k = len_trim(opt0)
          line = repeat(' ',k)//opt(k+1:80)
          line = adjustl(line)
          READ(line,'(A1)') Operation
          k = ichar(Operation)
          IF (k.GE.97 .and. k.LE.122) THEN
            Operation = char(k-32)
          ENDIF
          IF(Operation.EQ.'I') THEN
            iOper = 1
          ELSEIF(Operation.EQ.'R') THEN
            iOper = 2
          ELSEIF(Operation.EQ.'C') THEN
            iOper = 3
            line = ' '//line(2:80)
            READ(line,*) iRot
          ELSEIF(Operation.EQ.'S') THEN
            iOper = 4
            line = ' '//line(2:80)
            READ(line,*) iRot
          ELSE
            WRITE(iunit,*) 'ERROR. INCORRECT OPERATION'
            STOP
          ENDIF

c******************c
c* Reading Groups *c
c******************c
        case('%GROUP')
          DoGroup = .TRUE.
          iGroup = -100
          k = len_trim(opt0)
          line = repeat(' ',k)//opt(k+1:80)
          line = adjustl(line)
          READ(line,'(A1)') Group
          k = ichar(Group)
          IF (k.GE.97 .and. k.LE.122) THEN
            Group = char(k-32)
          ENDIF
          line = line(1:79)//'!'     ! line-control
          SubG(1:3) = '   '
          DO i = 1,3
            line = ' '//line(2:80)
            line = adjustl(line)
            READ(line,'(A1)') SubG(i:i)
            k = ichar(SubG(i:i))
            IF (k.GE.97 .AND. k.LE.122) THEN
              SubG(i:i) = char(k-32)
            ELSEIF(k.EQ.33) THEN
              SubG(i:i) = ' '
              GOTO 10    ! Stop when "!" is found
            ENDIF
          ENDDO

 10       IF(Group.EQ.'C' .OR. Group.EQ.'D') THEN
            i = 1
            k1 = ichar(SubG(1:1))
 	     IF (k1.GE.48 .and. k1.LE.57) THEN
	       i = 2
	       k2 = ichar(SubG(2:2))
 	       IF (k2.GE.48 .and. k2.LE.57) THEN
		   i = 3
	         k3 = ichar(SubG(3:3))
		   IF (k3.GE.48 .and. k3.LE.57) THEN
		     WRITE(iunit,*) 'ERROR. GROUP ORDER MUST BE < 100'
	           STOP
	         ENDIF
	       ENDIF
	     ENDIF
c C Groups
	   IF(Group.EQ.'C') THEN
           IF(SubG(i:i).EQ.'H') THEN
		 iGroup = 2  ! Cnh
           ELSEIF(SubG(i:i).EQ.'V') THEN
		 iGroup = 3  ! Cnv
           ELSEIF(SubG(1:3).EQ.'1  ') THEN
	       iGroup = 0   
             nGroup = 1   ! C1
           ELSEIF(SubG(1:3).EQ.'I  ') THEN
             iGroup = 0
             nGroup = 2  ! Ci
           ELSEIF(SubG(1:3).EQ.'S  ') THEN
             iGroup = 0
		 nGroup = 3  ! Cs
           ELSEIF(SubG(1:3).EQ.'INF') THEN
             iGroup = 9
             nGroup = 1  ! Cinf
           ELSEIF(SubG(i:i).EQ.' ') THEN
             i = 4
             iGroup = 1  ! Cn
	     ENDIF
c D Groups
	   ELSE
           IF(SubG(i:i).EQ.'H') THEN
		 iGroup = 5  ! Dnh
           ELSEIF(SubG(i:i).EQ.'D') THEN
             iGroup = 6  ! Dnd
           ELSEIF(SubG(1:3).EQ.'INF') THEN
             iGroup = 9
             nGroup = 2  ! Dinf
           ELSEIF(SubG(i:i).EQ.' ') THEN
             i = 4
             iGroup = 4  ! Dn
	     ENDIF
	   ENDIF

         IF(iGroup.GE.1 .AND. iGroup.LE.7) THEN
           READ(SubG(1:i-1),'(I4)') nGroup
         ENDIF
         IF(Group.EQ.'D' .AND. nGroup.EQ.1) THEN
           WRITE(iunit,997) Group, SubG
           STOP
         ENDIF
c Sn group
         ELSEIF(Group.EQ.'S') THEN
           iGroup = 7
           READ(SubG,*) nGroup
	     IF(nGroup .EQ. 2) THEN
	       Group = 'C'
             SubG(1:3) = 'I  '
             iGroup = 0  ! Sn -> Ci
	     ELSEIF(2*(nGroup/2).NE.nGroup) THEN
	       Group = 'C'
		 iGroup = 1  ! Sn -> Cn
	     ENDIF
           IF(nGroup.GE.100) THEN
	       WRITE(iunit,*) 'ERROR. GROUP ORDER MUST BE < 100'
	       STOP
	     ENDIF
c T,Th,Td group
         ELSEIF(Group.EQ.'T') THEN
           iGroup = 8
           IF(SubG(1:1).EQ.' ') THEN
             nGroup = 1
           ELSEIF(SubG(1:1).EQ.'H') THEN
             nGroup = 2
           ELSEIF(SubG(1:1).EQ.'D') THEN
             nGroup = 3
	     ENDIF
c O,Oh group
         ELSEIF(Group.EQ.'O') THEN
           iGroup = 8
	     IF(SubG(1:1).EQ.' ') THEN
             nGroup = 4
	     ELSEIF(SubG(1:1).EQ.'H') THEN
             nGroup = 5
	     ENDIF
c I,Ih group
         ELSEIF(Group.EQ.'I') THEN
           iGroup = 8
           IF(SubG(1:1).EQ.' ') THEN
             nGroup = 6
           ELSEIF(SubG(1:1).EQ.'H') THEN
             nGroup = 7
           ENDIF
         ELSE
           WRITE(iunit,997) Group, SubG
           STOP
         ENDIF

c Consistency Control
      IF(iGroup.GE.0 .AND. iGroup.LE.9) THEN
        WRITE(iunit,996) Group, SubG
      ELSE
        WRITE(iunit,997) Group, SubG
        STOP
      ENDIF
      IF(iGroup .NE. 7) THEN
       IF(nGroup.LT.2 .OR. nGroup.GT.6) THEN
        WRITE(iunit,*)'ERROR. Irreducible rep. only for n < 7'
        STOP
       ENDIF
      ENDIF


      CALL SetGroupDim(iGroup, nGroup, dGroup, hGroup)

        case('%DCART')
         PureCart = .FALSE.

        case('%UNRESTRICTED')
         Unrestricted = .TRUE.

        case('%HIGHPREC')
         HighPrec = .TRUE.

        case('%MOLECULE')
         molecule = .true.
        
        case('%CENTER_OPER')
          k = len_trim(opt0)
          opt0 = repeat(' ',k)//opt(k+1:80)
          read(opt0,*) RCread(1),RCread(2),RCread(3)
          readRC = .TRUE.

        case('%CENTER_OPER_NATOMS')
         readRN = .true.
         k = len_trim(opt0)
         opt0 = repeat(' ',k)//opt(k+1:80)
         read(opt0,*) Ncenter

        case('%DIPOLE')
          k = len_trim(opt0)
          opt0 = repeat(' ',k)//opt(k+1:80)
          read(opt0,*) RDipole(1),RDipole(2),RDipole(3)
          dipCHK = .TRUE.

        case('%ANGLES')
          k = len_trim(opt0)
          opt0 = repeat(' ',k)//opt(k+1:80)
          read(opt0,*) VAxis(1),VAxis(2)
          ANGLES = .TRUE.

        case('%ANGLES2')
          k = len_trim(opt0)
          opt0 = repeat(' ',k)//opt(k+1:80)
          read(opt0,*) VAxis2(1),VAxis2(2)
          ANGLES2 = .TRUE.

        case('%AXIS')
          k = len_trim(opt0)
          opt0 = repeat(' ',k)//opt(k+1:80)
          read(opt0,*) VAxis(1),VAxis(2),VAxis(3)
          AXIS = .TRUE.

        case('%AXIS2')
          k = len_trim(opt0)
          opt0 = repeat(' ',k)//opt(k+1:80)
          read(opt0,*) VAxis2(1),VAxis2(2),VAxis2(3)
          AXIS2 = .TRUE.

        case('%AXES_THRESH')
         Ax_Thresh = .true.
         k = len_trim(opt0)
         opt0 = repeat(' ',k)//opt(k+1:80)
         read(opt0,*) IAx_Thresh


        case('%AXIS_2ATOMS')
          k = len_trim(opt0)
          opt0 = repeat(' ',k)//opt(k+1:80)
          read(opt0,*) IAxis(1),IAxis(2)
	    AxVec = .TRUE.

        case('%AXIS2_2ATOMS')
          k = len_trim(opt0)
          opt0 = repeat(' ',k)//opt(k+1:80)
          read(opt0,*) IAxis2(1),IAxis2(2)
          AxVec2 = .TRUE.

        case('%AXIS_3ATOMS')
          k = len_trim(opt0)
          opt0 = repeat(' ',k)//opt(k+1:80)
          read(opt0,*) IAxis(1),IAxis(2),IAxis(3)
	    AxPlane = .TRUE.

        case('%AXIS2_3ATOMS')
          k = len_trim(opt0)
          opt0 = repeat(' ',k)//opt(k+1:80)
          read(opt0,*) IAxis2(1),IAxis2(2),IAxis2(3)
          AxPlane2 = .TRUE.

        case('%DO_INVERSION')
          DO_INV = .TRUE.

        case('%DO_REFLECTION')
          DO_INV = .FALSE.
          DO_REF = .TRUE.
        case('%DO_REFLEXION')
          DO_INV = .FALSE.
          DO_REF = .TRUE.

        case('%DO_ROTATION')
          k = len_trim(opt0)
          opt0 = repeat(' ',k)//opt(k+1:80)
          read(opt0,*) iRot
          DO_INV = .FALSE.
          DO_ROT = .TRUE.

        case('%NO_CENTER')
          NO_CENTER = .TRUE.

        case('%PRINT_LEVEL')
         k = len_trim(opt0)
         opt0 = repeat(' ',k)//opt(k+1:80)
         read(opt0,*) NPrint

        case('%PRINT_SYM_ORB')
	    DoSymOrb = .TRUE.

        case('%UHFNO')
          UHFNO = .TRUE.

        case default
          WRITE(iunit,998) opt0
          stop 
      endselect

  998 FORMAT(3X,'* ERROR * Unknown keyword: ',2X,A)
  997 FORMAT(3X,'* ERROR * INCORRECT SYMMETRY GROUP: ',2X,A1,A3)
  996 FORMAT('MEASURES OF THE SYMMETRY GROUP: ',2X,A1,A3)

  999 RETURN
      END

c------------------------------------------------------------------------
        subroutine nextline(iunit)
c
c Aquesta subroutine detecta les linies d'input, evitant linies en blanc i comentaris
c
        integer    :: iunit
        character  :: line*80
        logical    :: ctrl
c
        ctrl=.false.
        do while (.not.ctrl)
         read(iunit,"(a)",end=999) line
         line=adjustl(line)
         if((line(1:1).ne.'!').and.(len_trim(line).ne.0)) then
          backspace(iunit)
          ctrl=.true.
         end if
        end do
c
 999    return
        end subroutine
c------------------------------------------------------------------------

