
      SUBROUTINE sym_overlap(iAngA,AlphA,CA,RA,NA,LDA,
     &           iAngB,AlphB,CB,RB,NB,LDB,SAB)
c*****************************************************c
c                                                     c
c       Computes the overlap between two              c
c       electron densities.                           c
c                                                     c
c       INPUT                                         c
c        iAng_   : (x,y,z) exponent indices           c
c        Alph_   : function exponentials              c
c        C_      : function coefficients              c
c        R_      : function positions                 c
c        N_      : number of functions                c
c        LD_     : leading dimension                  c
c                                                     c
c       OUTPUT                                        c
c        SAB     : overlap                            c
c                                                     c
c       DCC: 05/09                                    c
c                                                     c
c*****************************************************c
      parameter   (nn = 1999999)
      INTEGER  iAngA(LDA,3),iAngB(LDB,3),iAngAB(3),
     &         NA,NB,LDA,LDB
      REAL*8   AlphA(*),CA(*),RA(LDA,3),
     &         AlphB(*),CB(*),RB(LDB,3),
     &         Dist,CAB,AlphAB,RAB(3),Test,S,SAB
      SAVE Zero
      DATA Zero / 0d0 /


      Test = 1.d-8
      SAB  = Zero

      DO i1 = 1,NA
      IF (ABS(CA(i1)) .GE. Test) THEN
        DO i2 = 1,NB
        IF (ABS(CB(i2)) .GE. Test) THEN

c Squared distance between centers
      Dist = Zero
      DO j = 1,3
        Dist = Dist + (RA(i1,j)-RB(i2,j))**2
      ENDDO

c Angular part               
      DO j = 1,3
        iAngAB(j) = iAngA(i1,j)+iAngB(i2,j)
      ENDDO

c Product constants    
      AlphAB = AlphA(i1)+AlphB(i2)
      CAB = exp(-Dist*AlphA(i1)*AlphB(i2)/
     &      AlphAB)
      CAB = CA(i1) * CB(i2) * CAB
      DO j = 1,3
        RAB(j) = (AlphA(i1)*RA(i1,j)+
     &        AlphB(i2)*RB(i2,j))/AlphAB
      ENDDO

c Compute overlap
      CALL overlap(iAngAB,AlphAB,RAB,S)
      SAB  =  SAB + CAB*S

        ENDIF
        ENDDO
      ENDIF
      ENDDO

      RETURN
      END
