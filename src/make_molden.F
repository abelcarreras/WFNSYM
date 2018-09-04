
        SUBROUTINE make_molden(NAt,iZAt,AtLab,RAt,Basis,input,
     &             NOcc,NBas,SGrp,WGrp,PureCart)        
c*****************************************************c
c                                                     c
c       Gets some dimensions                          c
c                                                     c
c       INPUT                                         c
c         Basis    atomic basis                       c
c         AtLab    atomic label                       c
c         NAt      number of atoms                    c
c	    RAt      atomic coordinates		      c
c	    SOrb_    sym diagonalized orbitals          c
c	    WOrb_    sym eigenvalues                    c
c         PureCart pure or cartesian D-orbitals       c
c                                                     c
c       OUTPUT                                        c
c         file.mol                                    c
c                                                     c
c       DCC: 09/2010                                  c
c                                                     c
c*****************************************************c
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8,  PARAMETER ::  ANGSBOHR  =  1.889726124993590
      INTEGER   ::  NBas,NAt,iZAt(NAt)
	REAL*8    ::  RAt(NAt,3),SGrp(NBas*NOcc),WGrp(NOcc)
      CHARACTER ::  Basis*80, input*80, line*80, 
     &              OrbType*2,AtLab(NAt)*2,Aux*255
      LOGICAL   ::  PureCart

c Get the auxiliarity folder
      CALL getenv("WSAUX",Aux)

c Control
	IF(PureCart) THEN
	  WRITE(2,*) 'ERROR. Not ready for pure orbitals'
	  STOP
	ENDIF

	Open(Unit=25,file=trim(input)//".mol",status='unknown')

	WRITE(25,*) '======= MOLDEN-FORMATTED INPUT FILE FOLLOWS ======='
	WRITE(25,"('[Molden Format]')")
	WRITE(25,"('[Atoms] (Angs)')")
	DO I = 1,NAt
	  WRITE(25,999) AtLab(I),I,iZAt(I),(RAt(I,J)/ANGSBOHR,J=1,3)
	ENDDO
	WRITE(25,"('[GTO]')")

c Print AO Basis
      DO 10 k = 1,NAt

	WRITE(25,"(I3,4X,'0')") k
	OPEN(Unit=20, Status="old", File=
     &"/"//trim(Aux)//"/basis/"
     &//trim(Basis)//"/"//trim(AtLab(k))//".dat")
      READ(20,*) NShell
      DO iShell = 1,NShell
	  READ(20,'(A2)',advance='no') OrbType
	  READ(20,*) NC, xScl
	  WRITE(25,'(A2,3X,I2,4X,F8.6)') OrbType, NC, xScl
	  DO I = 1,NC
	    READ(20,'(A80)') line
	    WRITE(25,'(A80)') line
	  ENDDO
	ENDDO
	WRITE(25,*) ' '
	CLOSE(20)
 10   CONTINUE

c Print MOs
      WRITE(25,"('[MO]')")
	DO 20 k = 1,NOcc
        WRITE(25,"('Sym=X')")
        WRITE(25,"('Ene= ',F15.6)") WGrp(k)
        WRITE(25,"('Spin=Alpha')")
        WRITE(25,"('Occup=1')")
	  DO 20 iAO = 1,NBas
	    WRITE(25,998) iAO, SGrp(iAO+(k-1)*NBas)
 20	CONTINUE
      WRITE(25,*) ' '

	WRITE(25,"('======= END OF MOLDEN-FORMATTED INPUT FILE =======')")

 998  FORMAT(3X,I3,3X,F15.12)
 999	FORMAT(2X,A2,3X,I3,3X,I2,2X,F14.8,2X,F14.8,2X,F14.8)
	CLOSE(25)

      RETURN
      END

