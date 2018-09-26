
        SUBROUTINE get_dim(NEtot,NEval,NBas,NOrb,Basis,
     &                     AtLab,NAt,iZAt,PureCart)
c*****************************************************c
c                                                     c
c       Gets some dimensions                          c
c                                                     c
c       INPUT                                         c
c         Basis    atomic basis                       c
c         AtLab    atomic label                       c
c         NAt      number of atoms                    c
c	  PureCart pure or cartesian D-orbitals       c
c                                                     c
c       OUTPUT                                        c
c         NEtot  # total electrons                    c
c         NEval  # valence electrons                  c
c         NBas   # of contracted AOs                  c
c         NOrb   # of uncontracted functions          c
c                                                     c
c       DCC: 04/2010                                  c
c       DCC: 06/2010  f-orbitals                      c
c                                                     c
c*****************************************************c
      INTEGER   ::  NBas,NOrb,NAt,iZAt(NAt)
      CHARACTER ::  Basis*80, OrbType*2,Atom*2,AtLab(*)*2,Aux*255
      LOGICAL   ::  PureCart

c Get the auxiliarity folder
      CALL getenv("WSAUX",Aux)

c Initialize some variables
      ND  = 6   ! cart d-orbitals
      NF  = 10  ! cart f-orbitals
      IF(PureCart) THEN
        ND = 5 ! pure d-orbitals
        NF = 7 ! pure d-orbitals
      END IF
      NEtot = 0
      NEval = 0
      NBas  = 0
      NOrb  = 0

      OPEN(Unit=10, Status="old", File=
     &"/"//trim(Aux)//"/atoms/atomicZ.dat")
      DO k = 1,NAt
         REWIND(10)
         DO iZ = 1,120
          READ(10,'(A2)',advance='no') Atom
          READ(10,*) Ntot, Nval
          IF(trim(Atom) .EQ. trim(AtLab(k))) THEN
            iZAt(k) = Ntot
            NEtot = NEtot + Ntot
            NEval = NEval + Nval
            GOTO 11
          ENDIF
         ENDDO
 11      CONTINUE

        OPEN(Unit=20, Status="old", File=
     &"/"//trim(Aux)//"/basis/"
     &//trim(Basis)//"/"//trim(AtLab(k))//".dat")

        READ(20,*) NShell
        DO iShell = 1,NShell
         READ(20,'(A2)',advance='no') OrbType
         READ(20,*) NC
         IF (OrbType.eq.'S ') THEN
           NBas = NBas + 1
           NOrb = NOrb + NC
         ELSEIF (OrbType.eq.'P ') THEN
           NBas = NBas + 3
           NOrb = NOrb + 3*NC
         ELSEIF (OrbType.eq.'SP') THEN
           NBas = NBas + 4
           NOrb = NOrb + 4*NC
         ELSEIF (OrbType.eq.'D ') THEN
           NBas = NBas + ND
           NOrb = NOrb + ND*NC
         ELSEIF (OrbType.eq.'F ') THEN
           NBas = NBas + NF
           NOrb = NOrb + NF*NC
         ENDIF

         DO j = 1,NC
           READ(20,*)
         ENDDO
        ENDDO
        CLOSE(20)
      ENDDO
      CLOSE(10)
      
      RETURN
      END
