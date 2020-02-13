
      SUBROUTINE sym_oper(SymMat,iOper,iRot,Axis,Angles,VAxis,iunit)
c*****************************************************c
c                                                     c
c       Builds the operation matrix to be             c
c       applied to the wave function                  c
c                                                     c
c       INPUT                                         c
c        iOper   : operation index                    c
c        Axis    : axis control                       c
c        VAxis   : orientation axis                   c
c        Angles  : angles control                     c
c        iunit   : unit to write results              c
c                                                     c
c       OUTPUT                                        c
c        SymMat  : symmetry matrix                    c
c                                                     c
c       DCC: 09/10                                    c
c                                                     c
c*****************************************************c
      IMPLICIT REAL*8  (a-h, o-z)
      REAL*8  ::  SymMat(3,3),VAxis(3)
      LOGICAL ::  Axis,Angles

      SAVE Zero,Half,One,Two,Four
      DATA Zero,Half,One,Two,Four / 0d0,.5d0,1d0,2d0,4d0 /

      COMMON /intvar1/ NAt,NPrint,Ncenter,IAxis,IAxis2,IAx_Thresh

      IF (iOper.EQ.1) THEN
        CALL VRload(SymMat,9,Zero)
        DO i = 1,3
          SymMat(i,i) = -One
        ENDDO
      ELSEIF (iOper.EQ.2) THEN
        IF(ANGLES) THEN
          WRITE(iunit,692) VAxis(1),VAxis(2)
        ELSEIF(Axis) THEN
          WRITE(iunit,690) VAxis(1),VAxis(2),VAxis(3)
        ENDIF
        WRITE(iunit,*)'Computing Reflexion Operation'
        CALL reflexion(VAxis,SymMat,ANGLES)
      ELSEIF (iOper.EQ.3) THEN
        IF(ANGLES) THEN
          WRITE(iunit,692) VAxis(1),VAxis(2)
        ELSEIF(Axis) THEN
          WRITE(iunit,690) VAxis(1),VAxis(2),VAxis(3)
        ENDIF
        WRITE(iunit,691) iRot
        CALL rotation(iRot,VAxis,SymMat,ANGLES,.TRUE.)
      ELSEIF (iOper.EQ.4) THEN
        IF(ANGLES) THEN
          WRITE(iunit,692) VAxis(1),VAxis(2)
        ELSEIF(Axis) THEN
          WRITE(iunit,690) VAxis(1),VAxis(2),VAxis(3)
        ENDIF
        WRITE(iunit,693) iRot
        CALL rotation(iRot,VAxis,SymMat,ANGLES,.FALSE.)
      ELSE
        WRITE(iunit,*) 'ERROR. Symmetry Operation not defined'
        STOP
      ENDIF

 690  FORMAT('Symmetry Axis Orientation (x,y,z): (',F7.4,
     &       ',',F7.4,',',F7.4,')')
 691  FORMAT('C',I2,'  Rotation Matrix:')
 692  FORMAT('Symmetry Axis Orientation (theta,phi): (',F6.2,
     &       ',',F6.2,')')
 693  FORMAT('S',I2,'  Improper-Rotation Matrix:')


      RETURN
      END


!     SUBROUTINE sym_group(SymGrp,iGroup,nGroup,dGroup,nIR,
!    &    Axis,Angles,VAxis,Axis2,Angles2,VAxis2,IAx_Thresh,iunit)
      SUBROUTINE sym_group(SymGrp,iGroup,nGroup,dGroup,
     &    VAxis,VAxis2,iunit)
c*****************************************************c
c                                                     c
c       Builds the group matrix to be                 c
c       applied to the wave function                  c
c                                                     c
c       INPUT                                         c
c        iGroup  : group index                        c
c        nGroup  : group order                        c
c        dGroup  : non-redundant operations           c
c        Axis    : axis control                       c
c        VAxis   : orientation axis                   c
c        Angles  : angles control                     c
c        iunit   : unit to write results              c
c                                                     c
c       OUTPUT                                        c
c        SymGrp : matrix with symmetry operations     c
c                                                     c
c       DCC: 09/10                                    c
c                                                     c
c*****************************************************c
      IMPLICIT REAL*8  (a-h, o-z)
      INTEGER ::  dGroup
      REAL*8  ::  SymMat(3,3),SAux0(3,3),SAux1(3,3),
     &            SymGrp(3,3,dGroup),VAxis(3),VAxis2(3),
     &            VAux0(3),VAux1(3),VAux2(3),VAux3(3),VAux4(3)
      LOGICAL ::  Axis,Angles,Axis2,Angles2,PrRot,ImRot

      SAVE Zero,Half,One,Two,Four
      DATA Zero,Half,One,Two,Four / 0d0,.5d0,1d0,2d0,4d0 /

c Set PrRot and ImRot to be used in rotation calls
      PrRot = .TRUE.
      ImRot = .FALSE.

c Set matrices to zero (to be safe)
      SymGrp(:,:,:) = Zero
      Angles = .FALSE.
      Angles2 = .FALSE.
      IAx_Thresh = 9

c Set thresh between axis
      Thresh = 10.d0**(-IAx_Thresh)

c Define identity
      DO I = 1,3
        SymGrp(I,I,1) = One
      ENDDO

c Define the rest
      IF(iGroup.EQ.0) THEN
        IF(nGroup.EQ.2) THEN
          DO I = 1,3
            SymGrp(I,I,2) = -One
          ENDDO
        ELSEIF(nGroup.EQ.3) THEN
        CALL reflexion(VAxis,SymMat,Angles)
        SymGrp(:,:,2) = SymMat(:,:)
        ENDIF

c Cn: Cn,Cnh,Cnv,Dn,Dnh,Dnd
      ELSEIF(iGroup.GE.1 .AND. iGroup.LE.6) THEN
        CALL rotation(nGroup,VAxis,SymMat,ANGLES,.TRUE.)
        SAux0(:,:) = SymMat(:,:) 
        SymGrp(:,:,2) = SAux0(:,:)
        DO iG = 3,1+nGroup/2
          CALL AtimsB(SAux1,SAux0,SymMat, 3,3,3, 3,3,3, 1)
          SAux0(:,:) = SAux1(:,:)
          SymGrp(:,:,iG) = SAux0(:,:)
        ENDDO
c s_h,Sn: Cnh,Dnh
        IF(iGroup.EQ.2 .OR. iGroup.EQ.5) THEN
          iOff = nGroup/2+1
          IF(iGroup.EQ.5) iOff = iOff + nGroup
          CALL reflexion(VAxis,SymMat,Angles)
          SymGrp(:,:,iOff+1) = SymMat(:,:)
          iOff = iOff + 1
          DO iG = 1,nGroup/2
            SAux0(:,:) = SymGrp(:,:,iG+1) 
            CALL AtimsB(SAux1,SymMat,SAux0, 3,3,3, 3,3,3, 1)
            SymGrp(:,:,iG+iOff) = SAux1(:,:)
          ENDDO
        ENDIF
c C2',C2": Dn , Dnh, Dnd &
c s_v,s_d: Cnv, Dnh, Dnd
        IF(iGroup.GE.3 .AND. iGroup.LE.6) THEN
          CALL VRdot(Dot,VAxis,VAxis2,3)
          Dot = DABS(Dot)
          IF(Dot.GT.Thresh) THEN
            WRITE(iunit,"('Axis  = ',3(1x,f9.5))") (VAxis(j),j=1,3)
            WRITE(iunit,"('Axis2 = ',3(1x,f9.5))") (VAxis2(j),j=1,3)
            WRITE(iunit,*)'ERROR. Axes not valid'
            WRITE(iunit,"(' IAx_Thresh = ',I2)") IAx_Thresh
            WRITE(iunit,*)'Modify the group axes or lower AXES_THRESH'
            STOP
          ENDIF
          VAux0(:) = VAxis2(:)
          CALL rotation(2*nGroup,VAxis,SAux0,ANGLES,.TRUE.)
          DO iG = 1,nGroup

c C2',C2": Dn,Dnh,Dnd
            IF(iGroup.GE.4) THEN
              CALL rotation(2,VAux0,SymMat,Angles2,.TRUE.)
              SymGrp(:,:,nGroup/2+1+iG) = SymMat(:,:)
            ENDIF
            CALL AtimsB(VAux1,SAux0,VAux0, 3,1,3, 3,3,3, 1)

c s_v,s_d: Cnv,Dnh,Dnd
            IF(iGroup.EQ.3.OR.iGroup.EQ.5.OR.iGroup.EQ.6) THEN
              IF(iGroup.EQ.3) THEN  ! Cnv
                iOff = nGroup/2+1
                VAux2(:) = VAux0(:)
              ELSEIF(iGroup.EQ.5) THEN  ! Dnh
                iOff = 2*(nGroup+1) - nGroup + 2*(nGroup/2)
                CALL VRcross(VAux2,VAux0,VAxis)
              ELSEIF(iGroup.EQ.6) THEN  ! Dnd 
                iOff = 2*nGroup + 1
                VAux0(:) = VAux1(:) + VAux0(:)
                CALL VRcross(VAux2,VAux0,VAxis)
              ENDIF
              CALL reflexion(VAux2,SymMat,Angles2)
              SymGrp(:,:,iOff+iG) = SymMat(:,:)
            ENDIF

            VAux0(:) = VAux1(:)
          ENDDO
        ENDIF
c Sn: Dnd 
        IF(iGroup.EQ.6) THEN
          iOff = nGroup + nGroup/2 + 1
          iSn = 2*nGroup
          CALL rotation(iSn,VAxis,SymMat,Angles2,.FALSE.)
          SAux0(:,:) = SymGrp(:,:,2)
          SymGrp(:,:,iOff+1) = SymMat(:,:)
          IF(nGroup.GE.3) THEN
            DO iG = 2,(nGroup+1)/2
              CALL AtimsB(SAux1,SAux0,SymMat, 3,3,3, 3,3,3, 1)
              SymGrp(:,:,iOff+iG) = SAux1(:,:)
              SymMat(:,:) = SAux1(:,:)
            ENDDO
          ENDIF
        ENDIF
c Sn
      ELSEIF(iGroup.EQ.7) THEN
        iRot = nGroup
        CALL rotation(nGroup,VAxis,SymMat,ANGLES,.FALSE.)
        SAux0(:,:) = SymMat(:,:)
        SymGrp(:,:,2) = SAux0(:,:)
        DO iG = 3,1+nGroup/2
          CALL AtimsB(SAux1,SAux0,SymMat, 3,3,3, 3,3,3, 1)
          SAux0(:,:) = SAux1(:,:)
          SymGrp(:,:,iG) = SAux0(:,:)
        ENDDO

c T,Th,Td, O,Oh, I,Ih
      ELSEIF(iGroup.EQ.8) THEN
c Td 
        IF(nGroup.EQ.3) THEN
c define: here we define two C3 from z,x axis
c       Dot1 = DATAN(DSQRT(Two))
c       VAux0(:) = VAxis(:)
c       VAxis(:)  =  DSIN(Dot1)*VAxis2(:) + DCOS(Dot1)*VAux0(:)	
c       VAxis2(:) = -DSIN(Dot1)*VAxis2(:) + DCOS(Dot1)*VAux0(:)	

          CALL VRdot(Dot,VAxis,VAxis2,3)
          Dot1 = Two*DATAN(DSQRT(Two))    ! ideal angle between Td axis
          Dot = ABS(Dot-DCOS(Dot1))
c          WRITE(iunit,*)'Diff = ',Dot
          IF(Dot.GT.Thresh) THEN
            WRITE(iunit,"('Axis  = ',3(1x,f15.9))") (VAxis(j),j=1,3)
            WRITE(iunit,"('Axis2 = ',3(1x,f15.9))") (VAxis2(j),j=1,3)
            WRITE(iunit,*)'Input cosine between C3 axes = ',Dot
            WRITE(iunit,*)'ERROR. Axes not valid'
            WRITE(iunit,*)'Ideal cosine between C3 axes = ',DCOS(Dot1)
c            WRITE(iunit,"(' IAx_Thresh = ',I2)") IAx_Thresh
            WRITE(iunit,*)'Modified Td axes: '
            CALL TdAxis2(VAxis,VAxis2)
            WRITE(iunit,"('Axis  = ',3(1x,f15.9))") (VAxis(j),j=1,3)
            WRITE(iunit,"('Axis2 = ',3(1x,f15.9))") (VAxis2(j),j=1,3)
          ENDIF
c          WRITE(iunit,"('Td axis orientation')")
          CALL rotation(3,VAxis,SymMat,.FALSE.,.TRUE.)
c          WRITE(iunit,"('1 C3:  ',3(1x,f15.9))") (VAxis(j),j=1,3)
          SymGrp(:,:,2) = SymMat(:,:)
          CALL rotation(3,VAxis2,SAux0,.FALSE.,.TRUE.)
c          WRITE(iunit,"('2 C3:  ',3(1x,f15.9))") (VAxis2(j),j=1,3)
          SymGrp(:,:,3) = SAux0(:,:)
          CALL AtimsB(VAux0,SymMat,VAxis2, 3,1,3, 3,3,3, 1)
          CALL rotation(3,VAux0,SAux0,.FALSE.,.TRUE.)
c          WRITE(iunit,"('3 C3:  ',3(1x,f15.9))") (VAux0(j),j=1,3)
          SymGrp(:,:,4) = SAux0(:,:)
          CALL AtimsB(VAux1,SymMat,VAux0, 3,1,3, 3,3,3, 1)
          CALL rotation(3,VAux1,SAux0,.FALSE.,.TRUE.)
c          WRITE(iunit,"('4 C3:  ',3(1x,f15.9))") (VAux1(j),j=1,3)
          SymGrp(:,:,5) = SAux0(:,:)
          VAux2(:) = VAxis(:) + VAxis2(:)
          VAux3(:) = VAxis(:) + VAux0(:)
          VAux4(:) = VAxis(:) + VAux1(:)
          CALL rotation(2,VAux2,SAux0,.FALSE.,.TRUE.)
c          WRITE(iunit,"('1 C2:  ',3(1x,f15.9))") (VAux2(j),j=1,3)
          SymGrp(:,:,6) = SAux0(:,:)
          CALL rotation(2,VAux3,SAux0,.FALSE.,.TRUE.)
c          WRITE(iunit,"('2 C2:  ',3(1x,f15.9))") (VAux3(j),j=1,3)
          SymGrp(:,:,7) = SAux0(:,:)
          CALL rotation(2,VAux4,SAux0,.FALSE.,.TRUE.)
c          WRITE(iunit,"('3 C2:  ',3(1x,f15.9))") (VAux4(j),j=1,3)
          SymGrp(:,:,8) = SAux0(:,:)
          CALL rotation(4,VAux2,SAux0,.FALSE.,.FALSE.)
c          WRITE(iunit,"('1 S4:  ',3(1x,f15.9))") (VAux2(j),j=1,3)
          SymGrp(:,:,9) = SAux0(:,:)
          CALL rotation(4,VAux3,SAux0,.FALSE.,.FALSE.)
c          WRITE(iunit,"('2 S4:  ',3(1x,f15.9))") (VAux3(j),j=1,3)
          SymGrp(:,:,10) = SAux0(:,:)
          CALL rotation(4,VAux4,SAux0,.FALSE.,.FALSE.)
c          WRITE(iunit,"('3 S4:  ',3(1x,f15.9))") (VAux4(j),j=1,3)
          SymGrp(:,:,11) = SAux0(:,:)
          CALL VRcross(VAux2,VAxis,VAxis2)
          CALL reflexion(VAux2,SAux0,.FALSE.)
c          WRITE(iunit,"('1 Ref: ',3(1x,f15.9))") (VAux2(j),j=1,3)
          SymGrp(:,:,12) = SAux0(:,:)
          CALL VRcross(VAux2,VAxis,VAux0)
          CALL reflexion(VAux2,SAux0,.FALSE.)
c          WRITE(iunit,"('2 Ref: ',3(1x,f15.9))") (VAux2(j),j=1,3)
          SymGrp(:,:,13) = SAux0(:,:)
          CALL VRcross(VAux2,VAxis,VAux1)
          CALL reflexion(VAux2,SAux0,.FALSE.)
c          WRITE(iunit,"('3 Ref: ',3(1x,f15.9))") (VAux2(j),j=1,3)
          SymGrp(:,:,14) = SAux0(:,:)
          CALL VRcross(VAux2,VAux0,VAux1)
          CALL reflexion(VAux2,SAux0,.FALSE.)
c          WRITE(iunit,"('4 Ref: ',3(1x,f15.9))") (VAux2(j),j=1,3)
          SymGrp(:,:,15) = SAux0(:,:)
          CALL VRcross(VAux2,VAxis2,VAux1)
          CALL reflexion(VAux2,SAux0,.FALSE.)
c          WRITE(iunit,"('5 Ref: ',3(1x,f15.9))") (VAux2(j),j=1,3)
          SymGrp(:,:,16) = SAux0(:,:)
          CALL VRcross(VAux2,VAxis2,VAux0)
          CALL reflexion(VAux2,SAux0,.FALSE.)
c          WRITE(iunit,"('6 Ref: ',3(1x,f15.9))") (VAux2(j),j=1,3)
          SymGrp(:,:,17) = SAux0(:,:)
c Oh  
        ELSEIF(nGroup.EQ.5) THEN
c VAxis,VAxis2 must be perpendicular
          CALL VRdot(Dot,VAxis,VAxis2,3)
          Dot = ABS(Dot)
          IF(Dot.GT.Thresh) THEN
            WRITE(iunit,"('Axis  = ',3(1x,f9.5))") (VAxis(j),j=1,3)
            WRITE(iunit,"('Axis2 = ',3(1x,f9.5))") (VAxis2(j),j=1,3)
            WRITE(iunit,*)'ERROR. Axes not valid'
            WRITE(iunit,"('IAx_Thresh = ',I2)") IAx_Thresh
            WRITE(iunit,*)'Modify the group axes or lower AXES_THRESH'
            STOP
          ENDIF
c Orientation: VAxis = Z, VAxis2 = X, VAux0 = Y.
          CALL VRcross(VAux0,VAxis,VAxis2)
c 2C4,C2 (z,x,y)
          CALL rotation(4,VAxis,SAux0,.FALSE.,PrRot)
          SymGrp(:,:,12) = SAux0(:,:)
          CALL rotation(4,VAxis2,SAux0,.FALSE.,PrRot)
          SymGrp(:,:,13) = SAux0(:,:)
          CALL rotation(4,VAux0,SAux0,.FALSE.,PrRot)
          SymGrp(:,:,14) = SAux0(:,:)
          CALL rotation(2,VAxis,SAux0,.FALSE.,PrRot)
          SymGrp(:,:,15) = SAux0(:,:)
          CALL rotation(2,VAxis2,SAux0,.FALSE.,PrRot)
          SymGrp(:,:,16) = SAux0(:,:)
          CALL rotation(2,VAux0,SAux0,.FALSE.,PrRot)
          SymGrp(:,:,17) = SAux0(:,:)
c inversion
          SymGrp(:,:,18) = Zero  
          DO ij = 1,3
            SymGrp(ij,ij,18) = -One
          ENDDO
c 2S4 (z,x,y)
          CALL rotation(4,VAxis,SAux0,.FALSE.,ImRot)
          SymGrp(:,:,19) = SAux0(:,:)
          CALL rotation(4,VAxis2,SAux0,.FALSE.,ImRot)
          SymGrp(:,:,20) = SAux0(:,:)
          CALL rotation(4,VAux0,SAux0,.FALSE.,ImRot)
          SymGrp(:,:,21) = SAux0(:,:)
c s_h (xy,yz,xz)
          CALL reflexion(VAxis,SAux0,.FALSE.)
          SymGrp(:,:,26) = SAux0(:,:)
          CALL reflexion(VAxis2,SAux0,.FALSE.)
          SymGrp(:,:,27) = SAux0(:,:)
          CALL reflexion(VAux0,SAux0,.FALSE.)
          SymGrp(:,:,28) = SAux0(:,:)
c C2 between 2 C4 axis
          VAux1(:) = VAxis(:) + VAxis2(:) 
          CALL rotation(2,VAux1,SAux0,.FALSE.,PrRot)
          SymGrp(:,:,6)  = SAux0(:,:)
          CALL reflexion(VAux1,SAux0,.FALSE.)
          SymGrp(:,:,29) = SAux0(:,:)
          VAux1(:) = VAxis(:) - VAxis2(:) 
          CALL rotation(2,VAux1,SAux0,.FALSE.,PrRot)
          SymGrp(:,:,7)  = SAux0(:,:)
          CALL reflexion(VAux1,SAux0,.FALSE.)
          SymGrp(:,:,30) = SAux0(:,:)
          VAux1(:) = VAxis(:) + VAux0(:) 
          CALL rotation(2,VAux1,SAux0,.FALSE.,PrRot)
          SymGrp(:,:,8)  = SAux0(:,:)
          CALL reflexion(VAux1,SAux0,.FALSE.)
          SymGrp(:,:,31) = SAux0(:,:)
          VAux1(:) = VAxis(:) - VAux0(:) 
          CALL rotation(2,VAux1,SAux0,.FALSE.,PrRot)
          SymGrp(:,:,9) = SAux0(:,:)
          CALL reflexion(VAux1,SAux0,.FALSE.)
          SymGrp(:,:,32) = SAux0(:,:)
          VAux1(:) = VAxis2(:) + VAux0(:) 
          CALL rotation(2,VAux1,SAux0,.FALSE.,PrRot)
          SymGrp(:,:,10) = SAux0(:,:)
          CALL reflexion(VAux1,SAux0,.FALSE.)
          SymGrp(:,:,33) = SAux0(:,:)
          VAux1(:) = VAxis2(:) - VAux0(:) 
          CALL rotation(2,VAux1,SAux0,.FALSE.,PrRot)
          SymGrp(:,:,11) = SAux0(:,:)
          CALL reflexion(VAux1,SAux0,.FALSE.)
          SymGrp(:,:,34) = SAux0(:,:)
c 2C3,2S6
          VAux1(:) = VAxis(:) + VAxis2(:) + VAux0(:)
          CALL rotation(3,VAux1,SAux0,.FALSE.,PrRot)
          SymGrp(:,:,2) = SAux0(:,:)
          CALL rotation(6,VAux1,SAux0,.FALSE.,ImRot)
          SymGrp(:,:,22) = SAux0(:,:)
          VAux1(:) = VAxis(:) - VAxis2(:) + VAux0(:)
          CALL rotation(3,VAux1,SAux0,.FALSE.,PrRot)
          SymGrp(:,:,3) = SAux0(:,:)
          CALL rotation(6,VAux1,SAux0,.FALSE.,ImRot)
          SymGrp(:,:,23) = SAux0(:,:)
          VAux1(:) = VAxis(:) + VAxis2(:) - VAux0(:)
          CALL rotation(3,VAux1,SAux0,.FALSE.,PrRot)
          SymGrp(:,:,4) = SAux0(:,:)
          CALL rotation(6,VAux1,SAux0,.FALSE.,ImRot)
          SymGrp(:,:,24) = SAux0(:,:)
          VAux1(:) = VAxis(:) - VAxis2(:) - VAux0(:)
          CALL rotation(3,VAux1,SAux0,.FALSE.,PrRot)
          SymGrp(:,:,5) = SAux0(:,:)
          CALL rotation(6,VAux1,SAux0,.FALSE.,ImRot)
          SymGrp(:,:,25) = SAux0(:,:)
    
        ELSE
          WRITE(iunit,*)'ERROR. Only Td and Oh of cubic groups is ready'
          STOP
        ENDIF

c Cinfv, Dinfh
      ELSEIF(iGroup.EQ.9) THEN
        write(iunit,*)'Linear groups not ready'
        stop
      ENDIF 

 690  FORMAT('Symmetry Axis Orientation (x,y,z): (',F7.4,
     &       ',',F7.4,',',F7.4,')')
 691  FORMAT('C',I2,'  Rotation Matrix:')
 692  FORMAT('Symmetry Axis Orientation (theta,phi): (',F6.2,
     &       ',',F6.2,')')
 693  FORMAT('S',I2,'  Improper-Rotation Matrix:')

      RETURN
      END


      SUBROUTINE SymLabels(Lab)
c
c  Builds up a character vector (Lab) containt the 
c  operation list of the point group (iGroup)
c  dGrp defines the lenght of the list

      INTEGER iGroup,nGroup,dGroup,hGroup
      CHARACTER   Lab(dGroup)*6, aux*2

      COMMON /intvar2/  iOper,iGroup,nGroup,dGroup,hGroup

c All groups contain the identity
      Lab(1) = '  E   '

      IF(iGroup.EQ.0) THEN 
        IF(nGroup.EQ.2) THEN
          Lab(2) = '  i   '
        ELSE
          Lab(2) = ' s_h  '
        ENDIF
c Groups with Cn (Cn,Cnh,Cnv,Dn,Dnh,Dnd)
      ELSEIF(iGroup.GE.1 .AND. iGroup.LE.6) THEN 
        IF(nGroup.EQ.2) THEN
          Lab(2) = '  C2  '
        ELSEIF(nGroup.EQ.3) THEN
          Lab(2) = ' 2C3  '
        ELSEIF(nGroup.EQ.4) THEN
          Lab(2) = ' 2C4  '
          Lab(3) = '  C2  '
        ELSEIF(nGroup.EQ.5) THEN
          Lab(2) = '2C5^1 '
          Lab(3) = '2C5^2 '
        ELSEIF(nGroup.EQ.6) THEN
          Lab(2) = ' 2C6  '
          Lab(3) = ' 2C3  '
          Lab(4) = '  C2  '
        ELSEIF (nGroup .EQ. 9) THEN
          Lab(2) = ' 2C9  '
          Lab(3) = ' 2C9^2'
          Lab(4) = ' 2C3  '
          Lab(5) = ' 2C9^4'
        ENDIF
c Groups with s_h and Sn (Cnh,Dnh)
        IF(iGroup.EQ.2 .OR. iGroup.EQ.5) THEN 
        iOff = 0
          IF(iGroup.EQ.5) iOff = nGroup
          IF(nGroup.EQ.2) THEN
            Lab(3+iOff) = ' s_h  '
            Lab(4+iOff) = '  i   '
          ELSEIF(nGroup.EQ.3) THEN
            Lab(3+iOff) = ' s_h  '
            Lab(4+iOff) = ' 2S3  '
          ELSEIF(nGroup.EQ.4) THEN
            Lab(4+iOff) = ' s_h  '
            Lab(5+iOff) = ' 2S4  '
            Lab(6+iOff) = '  i   '
          ELSEIF(nGroup.EQ.5) THEN
            Lab(4+iOff) = ' s_h  '
            Lab(5+iOff) = '2S5^1 '
            Lab(6+iOff) = '2S5^7 '
          ELSEIF(nGroup.EQ.6) THEN
            Lab(5+iOff) = ' s_h  '
            Lab(6+iOff) = ' 2S6  '
            Lab(7+iOff) = ' 2S3  '
            Lab(8+iOff) = '  i   '
          ENDIF
        ENDIF
c Groups with s_v,s_d (Cnv,Dnh)
        IF(iGroup.EQ.3 .OR. iGroup.EQ.5) THEN
          iOff = 0
          IF(iGroup.EQ.5) iOff = nGroup+nGroup/2+1
          IF(nGroup.EQ.2) THEN
            Lab(3+iOff) = ' s_v1 '
            Lab(4+iOff) = ' s_v2 '
          ELSEIF(nGroup.EQ.3) THEN
            Lab(3+iOff) = ' s_v1 '
            Lab(4+iOff) = ' s_v2 '
            Lab(5+iOff) = ' s_v3 '
          ELSEIF(nGroup.EQ.4) THEN
            Lab(4+iOff) = ' s_v1 '
            Lab(5+iOff) = ' s_d1 '
            Lab(6+iOff) = ' s_v2 '
            Lab(7+iOff) = ' s_d2 '
          ELSEIF(nGroup.EQ.5) THEN
            Lab(4+iOff) = ' s_v1 '
            Lab(5+iOff) = ' s_v2 '
            Lab(6+iOff) = ' s_v3 '
            Lab(7+iOff) = ' s_v4 '
            Lab(8+iOff) = ' s_v5 '
          ELSEIF(nGroup.EQ.6) THEN
            Lab(5+iOff)  = ' s_v1 '
            Lab(6+iOff)  = ' s_d1 '
            Lab(7+iOff)  = ' s_v2 '
            Lab(8+iOff)  = ' s_d2 '
            Lab(9+iOff)  = ' s_v3 '
            Lab(10+iOff) = ' s_d3 '
          ENDIF
        ENDIF
c Groups with C2',C2" (Dn,Dnh)
        IF(iGroup.EQ.4 .OR. iGroup.EQ.5) THEN
          IF(nGroup.EQ.2) THEN
            Lab(3:4) = (/" C2'1 "," C2'2 "/)
          ELSEIF(nGroup.EQ.3) THEN
            Lab(3:5) = (/" C2'1 "," C2'2 "," C2'3 "/)
          ELSEIF(nGroup.EQ.4) THEN
            Lab(4:7) = (/" C2'1 ",' C2"1 ',
     &                 " C2'2 ",' C2"2 '/)
          ELSEIF(nGroup.EQ.5) THEN
            Lab(4:8) = (/" C2'1 "," C2'2 ",
     &        " C2'3 "," C2'4 "," C2'5 "/)
          ELSEIF(nGroup.EQ.6) THEN
            Lab(5:10) = (/" C2'1 ",' C2"1 '," C2'2 ",
     &             ' C2"2 '," C2'3 ",' C2"3 '/)
          ENDIF
        ENDIF
c Dnd
        IF(iGroup.EQ.6) THEN 
          iC2 = 2 + nGroup/2
          isd = 2 + 2*nGroup
          DO i = 1,nGroup
            Lab(iC2) = " C2'  "
            Lab(isd) = " s_d  "
            iC2 = iC2 + 1
            isd = isd + 1
          ENDDO
          IF(nGroup.EQ.2) THEN
            Lab(5) = " 2S4  "
          ELSEIF(nGroup.EQ.3) THEN
            Lab(6) = " 2S6  "
            Lab(7) = "  i   "
          ELSEIF(nGroup.EQ.4) THEN
            Lab(8) = " 2S8  "
            Lab(9) = "2S8^3 "
          ELSEIF(nGroup.EQ.5) THEN
            Lab(9)  = "2S10  "
            Lab(10) = "2S10^3"
            Lab(11) = "  i   "
          ELSEIF(nGroup.EQ.6) THEN
            Lab(11) = "2S12  "
            Lab(12) = " 2S4  "
            Lab(13) = "2S12^5"
          ENDIF
        ENDIF
c Sn 
      ELSEIF(iGroup.EQ.7) THEN
        IF(nGroup.EQ.4) THEN
          Lab(2) = " 2S4  "
          Lab(3) = "  C2  "
        ELSEIF(nGroup.EQ.6) THEN
          Lab(2) = " 2S6  "
          Lab(3) = " 2C3  "
          Lab(4) = "  i   "
        ELSEIF(nGroup.EQ.8) THEN
          Lab(2) = " 2S8  "
          Lab(3) = " 2C4  "
          Lab(4) = " 2S8^3"
          Lab(5) = "  C2  "
        ELSEIF(nGroup.EQ.10) THEN
          Lab(2) = " 2S10 "
          Lab(3) = " 2C5  "
          Lab(4) = "2S10^3"
          Lab(5) = " 2C5^3"
          Lab(6) = "  i   "
        ELSEIF(nGroup.EQ.12) THEN
          Lab(2) = " 2S12 "
          Lab(3) = " 2C6  "
          Lab(4) = " 2S4  "
          Lab(5) = " 2C3  "
          Lab(6) = "2S12^5"
          Lab(7) = "  C2  "
        ENDIF
      
c T,Th,Td, O,Oh, I,Ih
      ELSEIF(iGroup.EQ.8) THEN
c Td
        IF(nGroup.EQ.3) THEN
          Lab(2:5)   = ' 2C3  '
          Lab(6:8)   = '  C2  '
          Lab(9:11)  = ' 2S4  '
          Lab(12:17) = ' s_d  '
c Oh
        ELSEIF(nGroup.EQ.5) THEN
          Lab(2:5)   = ' 2C3  '
          Lab(6:11)  = '  C2  '
          Lab(12:14) = ' 2C4  '
          Lab(15:17) = '  C2  '
          Lab(18)    = '  i   '
          Lab(19:21) = ' 2S4  '
          Lab(22:25) = ' 2S6  '
          Lab(26:28) = ' s_h  '
          Lab(29:34) = ' s_d  '
        ELSE
          WRITE(2,*)'ERROR. Only Td in SymLabels'
          STOP
        ENDIF
      ENDIF

      RETURN
      END

      SUBROUTINE SymMultip(Mult)
c
c  Builds up an integer vector (Mult) containt the
c  operation list multiplicity of the point group (iGroup)
c  dGrp defines the lenght of the list

      INTEGER     iGroup,nGroup,dGroup,hGroup,Mult(dGroup)
      LOGICAL     HasSn

      COMMON /intvar2/  iOper,iGroup,nGroup,dGroup,hGroup

      HasSn = .False.
      Mult(1) = 1
      IF(iGroup.EQ.0) THEN
        Mult(:)  = 1
      ELSEIF(iGroup.GE.1 .AND. iGroup.LE.6) THEN
        Mult(2:dGroup) = 1
        IF(nGroup.GE.3) THEN
          IF(iGroup.EQ.2 .OR. iGroup.EQ.5) THEN
            HasSn = .True.
            iOff = 2 + nGroup/2
            IF(iGroup.GE.5) iOff = iOff + nGroup
          ENDIF
          iFinal = (nGroup-1)/2 
          DO i = 1,iFinal
            Mult(i+1) = 2
            IF(HasSn) Mult(i+iOff) = 2
          ENDDO
        ENDIF
        IF(iGroup.EQ.6) THEN
          iOff = nGroup + nGroup/2 + 1
          iFinal = nGroup/2
          DO i = 1,iFinal
            Mult(i+iOff) = 2
          ENDDO
        ENDIF
c Sn
      ELSEIF(iGroup.EQ.7) THEN
        iFinal = dGroup - 1
        Mult(2:iFinal) = 2
        Mult(dGroup) = 1
      ELSEIF(iGroup.EQ.8) THEN
c Td
        IF(nGroup.EQ.3) THEN
          Mult(2:5)   = 2
          Mult(6:8)   = 1
          Mult(9:11)  = 2
          Mult(12:17) = 1
c Oh
        ELSEIF(nGroup.EQ.5) THEN
          Mult(2:5)   = 2
          Mult(6:11)  = 1
          Mult(12:14) = 2
          Mult(15:18) = 1
          Mult(19:25) = 2
          Mult(26:34) = 1
        ELSE
          Mult(:) = 1
        ENDIF
      ENDIF

      RETURN
      END


      SUBROUTINE reflexion(N,MatRef,Ctrl)
c
c  Builds up the reflexion matrix with N normal vector
c  Ctrl controls how the oritentation is given:
c        True  : Theta and Phi angles (in radiants)
c        False : (x,y,z) vector

      IMPLICIT REAL*8 (a-h, o-z)
      REAL*8   MatRef(3,3),N(3)
      LOGICAL  Ctrl
      SAVE     One, Two
      DATA     One, Two / 1.d0, 2.d0  /

      IF (Ctrl) THEN
        Theta = N(1) 
        Phi   = N(2) 
        N(1) = dsin(Theta) * dcos(Phi)
        N(2) = dsin(Theta) * dsin(Phi)
        N(3) = dcos(Theta)
      ELSE
        Dot = 0.0d0
        DO I = 1,3
          Dot = Dot + N(I)**2
        ENDDO
          N(:) = N(:) / DSQRT(Dot)
      ENDIF

      DO 10 I = 1,3
      DO 10 J = 1,3
        MatRef(I,J) = -Two*N(I)*N(J)
        IF (I.EQ.J) THEN
          MatRef(I,J) = MatRef(I,J) + One
        ENDIF
 10   CONTINUE

      RETURN
      END


      SUBROUTINE rotation(iRot,N,MatRot,Ctrl,Pure)
c
c  Builds up the iRot order rotation matrix with N normal vector
c  Ctrl controls how the oritentation is given:
c        True  : Theta and Phi angles (in radiants)
c        False : (x,y,z) vector
c  Pure controls pure and improper rotations
c	   True  : Pure rotation
c	   False : Improper rotation
c
      IMPLICIT REAL*8 (a-h, o-z)
      INTEGER  iRot
      REAL*8   MatRot(3,3), N(3)
      LOGICAL  Ctrl,Pure
      SAVE     One, Two, Four
      DATA     One, Two, Four / 1.d0, 2.d0, 4.d0  /

      PI = Four * DATAN(One)
      Alpha = Two*PI/dfloat(iRot)

      IF (Ctrl) THEN
        Theta = N(1) 
        Phi   = N(2) 
        N(1) = dsin(Theta) * dcos(Phi)
        N(2) = dsin(Theta) * dsin(Phi)
        N(3) = dcos(Theta)
      ELSE
       Dot = 0.0d0
       DO I = 1,3
         Dot = Dot + N(I)**2
       ENDDO
       N(:) = N(:) / DSQRT(Dot)
      ENDIF

c Pure-Rotation or Improper-Rotation?
      IF (Pure) THEN
        TTheta = One-dcos(Alpha)
      ELSE
        TTheta = -One-dcos(Alpha)
      ENDIF

      DO 10 I = 1,3
        DO 10 J = 1,3
          MatRot(I,J) = TTheta * N(I) * N(J)
          IF (I.EQ.J) THEN
            MatRot(I,J) = MatRot(I,J) + dcos(Alpha)
          ENDIF
 10     CONTINUE 
        MatRot(1,2) = MatRot(1,2) + dsin(Alpha)*N(3)
        MatRot(1,3) = MatRot(1,3) - dsin(Alpha)*N(2)
        MatRot(2,1) = MatRot(2,1) - dsin(Alpha)*N(3)
        MatRot(2,3) = MatRot(2,3) + dsin(Alpha)*N(1)
        MatRot(3,1) = MatRot(3,1) + dsin(Alpha)*N(2)
        MatRot(3,2) = MatRot(3,2) - dsin(Alpha)*N(1)
      
        RETURN
      END

      SUBROUTINE TdAxis2(V1,V2)
c
c  Recomputes the 2nd Td axis from VAxis and VAxis2
c
      REAL*8  Aph,TAph,Dot
      REAL*8  Rot(3,3),V(3),V1(3),V2(3)

c Define Td-angle
      Aph = 2.d0*DATAN(DSQRT(2.d0))
      TAph = 1.d0 - DCOS(Aph)

c Define Axis
      CALL VRcross(V,V2,V1)
      Dot = 0.0d0
      DO I = 1,3
        Dot = Dot + V(I)**2
      ENDDO
      V(:) = V(:) / DSQRT(Dot)

c Build the rotation matrix
      DO 10 I = 1,3
        DO 10 J = 1,3
          Rot(I,J) = TAph * V(I) * V(J)
          IF (I.EQ.J) THEN
            Rot(I,J) = Rot(I,J) + DCOS(Aph)
          ENDIF
 10     CONTINUE
        Rot(1,2) = Rot(1,2) + DSIN(Aph)*V(3)
        Rot(1,3) = Rot(1,3) - DSIN(Aph)*V(2)
        Rot(2,1) = Rot(2,1) - DSIN(Aph)*V(3)
        Rot(2,3) = Rot(2,3) + DSIN(Aph)*V(1)
        Rot(3,1) = Rot(3,1) + DSIN(Aph)*V(2)
        Rot(3,2) = Rot(3,2) - DSIN(Aph)*V(1)

c Rotate V1
        CALL AtimsB(V2,Rot,V1, 3,1,3, 3,3,3, 1)

        RETURN
      END

