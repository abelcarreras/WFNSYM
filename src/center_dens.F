      SUBROUTINE center_dens(iAngCen,AlphCen,CCen,RCen,NCen,
     &         RC,iAngDens,AlphDens,CDens,RDens,NDens)
c*****************************************************c
c                                                     c
c       Places the density to its charge center       c
c                                                     c
c     INPUT:                                          c
c                                                     c
c       RC        Density charge center               c
c       NDens     # density functions                 c
c       iAngDens  Density ang indices                 c
c       AlphDens  Density exponents                   c
c       CDens     Density coefficients                c
c       RDens     Density positions                   c
c                                                     c
c                                                     c
c     OUTPUT:                                         c
c                                                     c
c       NCen      # centered functions                c
c       iAngCen   Centered ang indices                c
c       AlphCen   Centered exponents                  c
c       CCen      Centered coefficients               c
c       RCen      Centered positions                  c
c                                                     c
c                                                     c
c       DCC: 05/09                                    c
c                                                     c
c*****************************************************c
      PARAMETER   (NN = 1999999)
      INTEGER   ::  iAngCen(NN,3),iAngDens(NDens,3),NAngTot
      REAL*8    ::  AlphCen(*),CCen(*),RCen(NN,3),
     &          RC(3),AlphDens(*),CDens(*),RDens(NDens,3)
      SAVE Zero,Half,One,Two,Four,Six
      DATA Zero,  Half, One, Two,Four,Six
     &    / 0d0, 0.5d0, 1d0, 2d0, 4d0, 6d0 /

c  Starting      
      iCen = 1
      DO i = 1,NDens
        NAngTot = 0
        DO j = 1,3
          NAngTot = NAngTot+iAngDens(i,j)
        ENDDO
c S type
        IF (NAngTot.eq.0) THEN
         DO j = 1,3
          RCen(iCen,j) = RDens(i,j)-RC(j)
          iAngCen(iCen,j) = iAngDens(i,j)
         ENDDO
         AlphCen(iCen) = AlphDens(i)
         CCen(iCen) = CDens(i)
         iCen = iCen+1
        
c P type
        ELSEIF (NAngTot.eq.1) THEN
         DO j = 1,3
          RCen(iCen,j) = RDens(i,j)-RC(j)
          RCen(iCen+1,j) = RDens(i,j)-RC(j)
          iAngCen(iCen,j) = iAngDens(i,j)
          iAngCen(iCen+1,j) = 0
          IF (iAngDens(i,j).eq.1) THEN
           CCen(iCen+1) = CDens(i)*RC(j)
          ENDIF
         ENDDO
         AlphCen(iCen) = AlphDens(i)
         AlphCen(iCen+1) = AlphDens(i)
         CCen(iCen) = CDens(i)
         iCen = iCen+2
       
c D type
        ELSEIF (NAngTot.eq.2) THEN
         DO j = 1,3
c     (2,0,0)
          IF (iAngDens(i,j).eq.2) THEN
           DO l = 1,3
            AlphCen(iCen+l-1) = AlphDens(i)
            iAngCen(iCen,l) = iAngDens(i,l)
            iAngCen(iCen+1,l) = 0
            iAngCen(iCen+2,l) = 0
            RCen(iCen,l) = RDens(i,l)-RC(l)
            RCen(iCen+1,l) = RDens(i,l)-RC(l)
            RCen(iCen+2,l) = RDens(i,l)-RC(l)
           ENDDO
           iAngCen(iCen+1,j) = 1
           CCen(iCen) = CDens(i)
           CCen(iCen+1) = CDens(i)*Two*RC(j)
           CCen(iCen+2) = CDens(i)*RC(j)**2
           iCen = iCen+3
          ENDIF
         ENDDO
         DO j = 1,2
          IF (iAngDens(i,j).eq.1) THEN
           DO l = j+1,3
            IF (iAngDens(i,l).eq.1) THEN
             DO k = 1,4
              AlphCen(iCen+k-1) = AlphDens(i)
              DO kk = 1,3
               RCen(iCen+k-1,kk) = RDens(i,kk)-RC(kk)
               iAngCen(iCen+k-1,kk) = 0
              ENDDO
             ENDDO
             iAngCen(iCen,j) = 1
             iAngCen(iCen,l) = 1
             iAngCen(iCen+1,j) = 1
             iAngCen(iCen+2,l) = 1
             CCen(iCen) = CDens(i)
             CCen(iCen+1) = CDens(i)*RC(l)
             CCen(iCen+2) = CDens(i)*RC(j)
             CCen(iCen+3) = CDens(i)*RC(l)*RC(j)
             iCen = iCen+4
            ENDIF
           ENDDO
          ENDIF
         ENDDO
       
c F type
        ELSEIF (NAngTot.eq.3) THEN
         DO j = 1,3
c     (3 0 0)
          IF (iAngDens(i,j).eq.3) THEN
           DO k = 1,4
            AlphCen(iCen+k-1) = AlphDens(i)
            DO l = 1,3
             iAngCen(iCen+k-1,l) = 0
             RCen(iCen+k-1,l) = RDens(i,l)-RC(l)
            ENDDO
            iAngCen(iCen+k-1,j) = 4-k
           ENDDO
           CCen(iCen) = CDens(i)
           CCen(iCen+1) = CDens(i)*3.*RC(j)
           CCen(iCen+2) = CDens(i)*3.*(RC(j))**2
           CCen(iCen+3) = CDens(i)*(RC(j))**3
           iCen = iCen+4
          ENDIF
c     (2 1 0)
          IF (iAngDens(i,j).eq.2) THEN
           DO k = 1,6
            DO l = 1,3
             iAngCen(iCen+k-1,l) = 0
             RCen(iCen+k-1,l) = RDens(i,l)-RC(l)
            ENDDO
            AlphCen(iCen+k-1) = AlphDens(i)
           ENDDO
           DO l = 1,3
            iAngCen(iCen,l) = iAngDens(i,l)
            IF (iAngDens(i,l).eq.1) THEN
             iAngCen(iCen+1,j) = 2
             iAngCen(iCen+2,j) = 1
             iAngCen(iCen+2,l) = 1
             iAngCen(iCen+3,j) = 1
             iAngCen(iCen+4,l) = 1
             CCen(iCen+1) = CDens(i)*RC(l)
             CCen(iCen+3) = CDens(i)*Two*RC(j)*RC(l)
             CCen(iCen+5) = CDens(i)*RC(l)*(RC(j))**2
            ENDIF
           ENDDO
           CCen(iCen) = CDens(i)
           CCen(iCen+2) = CDens(i)*Two*RC(j)
           CCen(iCen+4) = CDens(i)*(RC(j))**2
           iCen = iCen+6
          ENDIF
c     (1 1 1)
         ENDDO
          IF ((iAngDens(i,1).eq.1).and.(iAngDens(i,2).eq.1).and.
     &     (iAngDens(i,3).eq.1)) THEN
           DO k = 1,4
            AlphCen(iCen+k-1) = AlphDens(i)
            AlphCen(iCen+k+3) = AlphDens(i)
            DO l = 1,3
             iAngCen(iCen+k-1,l) = 1
             iAngCen(iCen+k+3,l) = 0
             iAngCen(iCen+l,4-l) = 0
             iAngCen(iCen+l+3,l) = 1
             RCen(iCen+k-1,l) = RDens(i,l)-RC(l)
             RCen(iCen+k+3,l) = RDens(i,l)-RC(l)
             CCen(iCen+l) = CDens(i)*RC(4-l)
            ENDDO
           ENDDO
           CCen(iCen) = CDens(i)
         CCen(iCen+4) = CDens(i)*RC(2)*RC(3)
           CCen(iCen+5) = CDens(i)*RC(1)*RC(3)
           CCen(iCen+6) = CDens(i)*RC(1)*RC(2)
           CCen(iCen+7) = CDens(i)*RC(1)*RC(2)*RC(3)
          iCen = iCen+8
          ENDIF
      
c G type
        ELSEIF (NAngTot.eq.4) THEN
         DO j = 1,3
c     (4 0 0)
          IF (iAngDens(i,j).eq.4) THEN
           DO k = 1,5
            AlphCen(iCen+k-1) = AlphDens(i)
            DO l = 1,3
              iAngCen(iCen+k-1,l) = 0
              RCen(iCen+k-1,l) = RDens(i,l)-RC(l)
            ENDDO
            iAngCen(iCen+k-1,j) = 5-k
           ENDDO
           CCen(iCen) = CDens(i)
           CCen(iCen+1) = CDens(i)*Four*RC(j)
           CCen(iCen+2) = CDens(i)*Six*RC(j)**2
           CCen(iCen+3) = CDens(i)*Four*RC(j)**3
           CCen(iCen+4) = CDens(i)*RC(j)**4
           iCen = iCen+5
          ENDIF
c     (3 1 0)
          IF (iAngDens(i,j).eq.3) THEN
           DO k = 1,8
            DO l = 1,3
             iAngCen(iCen+k-1,l) = 0
             RCen(iCen+k-1,l) = RDens(i,l)-RC(l)
            ENDDO
            AlphCen(iCen+k-1) = AlphDens(i)
           ENDDO
           DO l = 1,3
            iAngCen(iCen,l) = iAngDens(i,l)
            IF (iAngDens(i,l).eq.1) THEN
             iAngCen(iCen+1,j) = 3
             iAngCen(iCen+3,j) = 2
             iAngCen(iCen+5,j) = 1
             iAngCen(iCen+2,j) = 2
             iAngCen(iCen+2,l) = 1
             iAngCen(iCen+4,j) = 1
             iAngCen(iCen+4,l) = 1
             iAngCen(iCen+6,l) = 1
             CCen(iCen+1) = CDens(i)*RC(l)
             CCen(iCen+3) = CDens(i)*3.*RC(j)*RC(l)
             CCen(iCen+5) = CDens(i)*3.*(RC(j)**2)*RC(l)
             CCen(iCen+7) = CDens(i)*(RC(j)**3)*RC(l)
            ENDIF
           ENDDO
           CCen(iCen) = CDens(i)
           CCen(iCen+2) = CDens(i)*3.*RC(j)
           CCen(iCen+4) = CDens(i)*3.*RC(j)**2
           CCen(iCen+6) = CDens(i)*RC(j)**3
           iCen = iCen+8
          ENDIF
         ENDDO
         DO j = 1,2
c     (2 2 0)
          IF (iAngDens(i,j).eq.2) THEN
           DO l = j+1,3
            IF (iAngDens(i,l).eq.2) THEN
             DO k = 1,9
              AlphCen(iCen+k-1) = AlphDens(i)
              DO kk = 1,3
               RCen(iCen+k-1,kk) = RDens(i,kk)-RC(kk)
               iAngCen(iCen+k-1,kk) = 0
              ENDDO
             ENDDO
             iAngCen(iCen,j) = 2
             iAngCen(iCen,l) = 2
             iAngCen(iCen+1,j) = 2
             iAngCen(iCen+1,l) = 1
             iAngCen(iCen+2,j) = 1
             iAngCen(iCen+2,l) = 2
             iAngCen(iCen+3,j) = 2
             iAngCen(iCen+4,l) = 2
             iAngCen(iCen+5,j) = 1
             iAngCen(iCen+5,l) = 1
             iAngCen(iCen+6,j) = 1
             iAngCen(iCen+7,l) = 1
             CCen(iCen) = CDens(i)
             CCen(iCen+1) = CDens(i)*Two*RC(l)
             CCen(iCen+2) = CDens(i)*Two*RC(j)
             CCen(iCen+3) = CDens(i)*RC(l)**2
             CCen(iCen+4) = CDens(i)*RC(j)**2
             CCen(iCen+5) = CDens(i)*Four*RC(l)*RC(j)
             CCen(iCen+6) = CDens(i)*Two*RC(j)*RC(l)**2
             CCen(iCen+7) = CDens(i)*Two*RC(l)*RC(j)**2
             CCen(iCen+8) = CDens(i)*(RC(l)*RC(j))**2
             iCen = iCen+9
            ENDIF
           ENDDO
          ENDIF
          IF (iAngDens(i,j).eq.1) THEN
           DO l = j+1,3
            IF (iAngDens(i,l).eq.1) THEN
             DO k = 1,4
              AlphCen(iCen+k-1) = AlphDens(i)
              AlphCen(iCen+k+3) = AlphDens(i)
              AlphCen(iCen+k+7) = AlphDens(i)
              DO kk = 1,3
               RCen(iCen+k-1,kk) = RDens(i,kk)-RC(kk)
               RCen(iCen+k+3,kk) = RDens(i,kk)-RC(kk)
               RCen(iCen+k+7,kk) = RDens(i,kk)-RC(kk)
               iAngCen(iCen+k-1,kk) = iAngDens(i,kk)
               iAngCen(iCen+k+3,kk) = 1
               iAngCen(iCen+k+7,kk) = 0
              ENDDO
             ENDDO
             k = 6-(j+l)
             iAngCen(iCen+1,l) = 0
             iAngCen(iCen+2,j) = 0
             iAngCen(iCen+3,j) = 0
             iAngCen(iCen+3,l) = 0
             iAngCen(iCen+5,l) = 0
             iAngCen(iCen+6,j) = 0
             iAngCen(iCen+7,k) = 0
             iAngCen(iCen+8,k) = 1
             iAngCen(iCen+9,j) = 1
             iAngCen(iCen+10,l) = 1
             CCen(iCen) = CDens(i)
             CCen(iCen+1) = CDens(i)*RC(l)
             CCen(iCen+2) = CDens(i)*RC(j)
             CCen(iCen+3) = CDens(i)*RC(l)*RC(j)
             CCen(iCen+4) = CDens(i)*Two*RC(k)
             CCen(iCen+5) = CDens(i)*Two*RC(k)*RC(l)
             CCen(iCen+6) = CDens(i)*Two*RC(k)*RC(j)
             CCen(iCen+7) = CDens(i)*RC(k)**2
             CCen(iCen+8) = CDens(i)*Two*RC(k)*RC(l)*RC(j)
             CCen(iCen+9) = CDens(i)*RC(l)*RC(k)**2
             CCen(iCen+10) = CDens(i)*RC(j)*RC(k)**2
             CCen(iCen+11) = CDens(i)*RC(j)*RC(l)*RC(k)**2
             iCen = iCen+12
            ENDIF
           ENDDO
          ENDIF
         ENDDO
      
        ENDIF

      ENDDO  

c  Final number of functions in the centered density
      NCen = iCen-1


      RETURN
      END
