
        SUBROUTINE make_uhfno(Ca,Cb,WOcc,S,NBas,input)
c*****************************************************c
c                                                     c
c       Builds UHF natural orbitals                   c
c                                                     c
c       INPUT                                         c
c         C_         empty matrices                   c
c         S          AOs overlap matrix               c
c         NBas       AO dimension                     c
c         input      input file                       c
c                                                     c
c       OUTPUT                                        c
c         Ca         UHF natural orbitals             c 
c         WOcc       UHF NOs occupancies              c 
c                                                     c
c       DCC: 11/2012                                  c
c                                                     c
c*****************************************************c
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER   ::  NBas, IOff, reclength
      REAL*8    ::  Ca(NBas**2),Cb(NBas**2),WOcc(NBas),S(NBas**2)
      CHARACTER :: input*80,line*80
      CHARACTER*4 :: cLAO(NBas)
      REAL*8, ALLOCATABLE :: work(:)

	NBas2 = NBas**2
c Open 54.0 file for binary IO (input/output)
      inquire(iolength=reclength) Ca,Cb
c  reclength = reclength*8
c  double    == 8bytes,
c  reclength now equals #bytes( value & grbge )
      OPEN(Unit=16,File=trim(input)//".in.save/54.0",
     &     access='direct',action='read',
     &     recl=reclength,form='unformatted',status="unknown")

c Read density matrices from binary file
      IOff = 1
      READ(16,rec=IOff) Ca,Cb
      CLOSE(16)

c Print Ca and Cb
      line = 'Alpha Density Matrix'
      CALL PrintMatrix(Ca,cLAO,.FALSE.,.FALSE.,NBas,line,2)
      line = 'Beta Density Matrix'
      CALL PrintMatrix(Cb,cLAO,.FALSE.,.FALSE.,NBas,line,2)

c Add alpha and beta densities: Ca = Ca + Cb
	CALL VRadd2(Ca,Cb,NBas2)
      line = 'Total Density Matrix'
      CALL PrintMatrix(Ca,cLAO,.FALSE.,.FALSE.,NBas,line,2)

c ##########################
c #  Diagonilze UHF 1-PDM  #
c ##########################

c Allocate some space
      lwork = 3*NBas
      ALLOCATE(work(lwork))

c Diagonalize it
      CALL DSYGV(2,'V','U',NBas,Ca,NBas,S,NBas,WOcc,work,lwork,Info)
c     CALL DSYEV('V','U',NBas,Ca,NBas,WOcc,work,lwork,Info)
      DEALLOCATE(work)

      RETURN
      END



