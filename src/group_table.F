      SUBROUTINE group_table(GrpTbl,nIR,iGroup,nGroup,hGroup,iout)
c*****************************************************c
c                                                     c
c       Returns the table of characters of            c
c       a symmetry point group                        c
c                                                     c
c       COMMON                                        c
c        iGroup  : group index                        c
c        nGroup  : group order                        c
c        hGroup  : group dimension                    c
c                                                     c
c       INPUT                                         c
c        nIR     : # of irreducible rep.              c
c                                                     c
c       OUTPUT                                        c
c        GrpTbl  : group table                        c
c                                                     c
c       DCC: 10/10                                    c
c                                                     c
c*****************************************************c
      IMPLICIT REAL*8  (a-h, o-z)
      INTEGER    iGroup,nGroup,hGroup
      REAL*8     GrpTbl(hGroup,hGroup)

      SAVE Zero,Half,One,Two
      DATA Zero,Half,One,Two / 0d0,0.5d0,1d0,2d0 /

c Set PI value and SQRT(3)/2
      PI   = 4.0d0 * DATAN(One)
      Coef = Half * DSQRT(3.d0) 

c Only up to 6th order (nGroup < 7)
	IF(iGroup .NE. 7) THEN
       IF(nGroup.LT.2 .OR. nGroup.GT.6) THEN
        WRITE(iout,*)'ERROR. Irreducible rep. only for n < 7'
        STOP
       ENDIF
      ENDIF

c***********************c
c Make character tables c
c***********************c
c Nonaxial
      IF(iGroup.EQ.0) THEN
       GrpTbl(:,:) =  One
       GrpTbl(2,2) = -One
       IF(nGroup.EQ.1) THEN
        WRITE(iout,*)'ERROR. Irred. rep. not ready for C1'
        STOP
       ELSEIF(nGroup.EQ.2) THEN
       ELSEIF(nGroup.EQ.3) THEN
       ENDIF
c Cn
      ELSEIF(iGroup.EQ.1) THEN
       IF(nGroup.EQ.2) THEN
        GrpTbl(:,:) =  One
        GrpTbl(2,2) = -One
       ELSEIF(nGroup.EQ.3) THEN
        GrpTbl(1,:) = One
        GrpTbl(2,:) = (/One,-Half,-Half/)
	  GrpTbl(3,:) = (/Zero,Coef,Coef/)
       ELSEIF(nGroup.EQ.4) THEN
        GrpTbl(1,:) = One
        GrpTbl(2,:) = (/One,-One,-One,One/)
        GrpTbl(3,:) = (/One,Zero,Zero,-One/)
        GrpTbl(4,:) = (/Zero,One,-One,Zero/)
       ELSEIF(nGroup.EQ.5) THEN
	  theta  = Two * PI / 5.0d0
	  theta2 = Two * theta
	  Cth  = DCOS(theta)
	  C2th = DCOS(theta2)
	  Sth  = DSIN(theta)
	  S2th = DSIN(theta2)
        GrpTbl(1,:) = One
	  GrpTbl(2,:) = (/One,Cth,Cth,C2th,C2th/)
	  GrpTbl(3,:) = (/Zero,Sth,-Sth,S2th,-S2th/)
	  GrpTbl(4,:) = (/One,C2th,C2th,Cth,Cth/)
	  GrpTbl(5,:) = (/Zero,-S2th,S2th,Sth,-S2th/)
       ELSEIF(nGroup.EQ.6) THEN
        GrpTbl(1,:) = One
	  GrpTbl(2,:) = (/One,-One,-One,One,One,-One/)
	  GrpTbl(3,:) = (/One,Half,Half,-Half,-Half,-One/)
	  GrpTbl(4,:) = (/Zero,Coef,-Coef,Coef,-Coef,Zero/)
        GrpTbl(5,:) = GrpTbl(3,:)*GrpTbl(2,:)
        GrpTbl(6,:) = GrpTbl(4,:)*GrpTbl(2,:)
	 ENDIF
c Cnh
      ELSEIF(iGroup.EQ.2) THEN
       IF(nGroup.EQ.2) THEN
        GrpTbl(1,:) = One
        GrpTbl(2,1) = One
        GrpTbl(2,2) = -One
        GrpTbl(2,3) = -One
        GrpTbl(2,4) = One
        GrpTbl(3,1:2) = One
        GrpTbl(3,3:4) = -One
        GrpTbl(4,:) = GrpTbl(3,:)*GrpTbl(2,:)
       ELSEIF(nGroup.EQ.3) THEN
        GrpTbl(1,:) = One
        GrpTbl(2,:) = (/One,-Half,-Half,One,-Half,-Half/)
	  GrpTbl(3,:) = (/Zero,Coef,-Coef,Zero,Coef,-Coef/)
        GrpTbl(4,:) = (/One,One,One,-One,-One,-One/)
	  GrpTbl(5,:) = GrpTbl(2,:)*GrpTbl(4,:)
	  GrpTbl(6,:) = GrpTbl(3,:)*GrpTbl(4,:)
       ELSEIF(nGroup.EQ.4) THEN
        GrpTbl(1,:) =  One
	  GrpTbl(2,:) = (/One,-One,-One,One,One,-One,-One,One/)
	  GrpTbl(3,:) = (/One,Zero,Zero,-One,-One,Zero,Zero,One/)
	  GrpTbl(4,:) = (/Zero,One,-One,Zero,Zero,-One,One,Zero/)
	  GrpTbl(5,:) = (/One,One,One,One,-One,-One,-One,-One/)
	  GrpTbl(6,:) = GrpTbl(2,:)*GrpTbl(5,:)
	  GrpTbl(7,:) = GrpTbl(3,:)*GrpTbl(5,:)
	  GrpTbl(8,:) = GrpTbl(4,:)*GrpTbl(5,:)
       ELSEIF(nGroup.EQ.5) THEN
        theta  = Two * PI / 5.0d0
        theta2 = Two * theta
        Cth  = DCOS(theta)
        C2th = DCOS(theta2)
        Sth  = DCOS(theta)
        S2th = DCOS(theta2)
	  GrpTbl(1,:) = One
	  GrpTbl(2,:) = (/One,Cth,Cth,C2th,C2th,
     &                  One,Cth,Cth,C2th,C2th/)
	  GrpTbl(3,:) = (/Zero,Sth,-Sth,S2th,-S2th,
     &                  Zero,Sth,-Sth,S2th,-S2th/)
	  GrpTbl(4,:) = (/One,C2th,C2th,Cth,Cth,
     &                  One,C2th,C2th,Cth,Cth/)
        GrpTbl(5,:) = (/Zero,-S2th,S2th,Sth,-Sth,
     &                  Zero,-S2th,S2th,Sth,-Sth/)
        GrpTbl(6,:) = (/(One,i=1,5),(-One,i=6,10)/)
	  GrpTbl(7,:) = GrpTbl(2,:)*GrpTbl(6,:)
	  GrpTbl(8,:) = GrpTbl(3,:)*GrpTbl(6,:)
	  GrpTbl(9,:) = GrpTbl(4,:)*GrpTbl(6,:)
	  GrpTbl(10,:)= GrpTbl(5,:)*GrpTbl(6,:)
       ELSEIF(nGroup.EQ.6) THEN
	  GrpTbl(1,:) = One
	  GrpTbl(2,:) = (/One,-One,-One,One,One,-One,
     &                 -One,One,One,-One,-One,One/)
	  GrpTbl(3,:) = (/One,Half,Half,-Half,-Half,-One,
     &                 -One,-Half,-Half,Half,Half,One/)
	  GrpTbl(4,:) = (/Zero,Coef,-Coef,Coef,-Coef,Zero,
     &                  Zero,-Coef,Coef,-Coef,Coef,Zero/)
	  GrpTbl(5,:) = (/One,(-Half,i=2,5),One,
     &                  One,(-Half,i=2,5),One/)
        GrpTbl(6,:) = (/Zero,-Coef,Coef,Coef,-Coef,Zero,
     &                  Zero,-Coef,Coef,Coef,-Coef,Zero/)
	  GrpTbl(7,:) = (/(One,i=1,6),(-One,i=7,12)/)
        GrpTbl(8,:) = GrpTbl(2,:)*GrpTbl(7,:)
        GrpTbl(9,:) = GrpTbl(3,:)*GrpTbl(7,:)
        GrpTbl(10,:)= GrpTbl(4,:)*GrpTbl(7,:)
        GrpTbl(11,:)= GrpTbl(5,:)*GrpTbl(7,:)
        GrpTbl(12,:)= GrpTbl(6,:)*GrpTbl(7,:)
       ENDIF
c Cnv
      ELSEIF(iGroup.EQ.3) THEN
       IF(nGroup.EQ.2) THEN
        GrpTbl(1,:) = One
	  DO i = 1,2
	    GrpTbl(2,  i) =  One
	    GrpTbl(2,2+i) = -One
	    GrpTbl(3,2*i-1) = One
	    GrpTbl(3,2*i) = -One
	  ENDDO
	  GrpTbl(4,:) = GrpTbl(2,:)*GrpTbl(3,:)
       ELSEIF(nGroup.EQ.3) THEN
        GrpTbl(1,:) = One
        GrpTbl(2,1:3) = One
        GrpTbl(2,4:6) = -One
        GrpTbl(3,1:6) = (/One,-Half,-Half,One,-Half,-Half/)
        GrpTbl(4,1:6) = (/Zero,Coef,-Coef,Zero,-Coef,Coef/)
        GrpTbl(5,1:6) = (/Zero,-Coef,Coef,Zero,-Coef,Coef/)
        GrpTbl(6,1:6) = (/One,-Half,-Half,-One,Half,Half/)
       ELSEIF(nGroup.EQ.4) THEN
        GrpTbl(1,:) = One
	  GrpTbl(2,1:4) = One
	  GrpTbl(2,5:8) = -One
        GrpTbl(3,:) = (/One,-One,-One,One,One,-One,One,-One/)
        GrpTbl(4,:) = (/One,-One,-One,One,-One,One,-One,One/)
        GrpTbl(5,:) = (/One,Zero,Zero,-One,-One,Zero,One,Zero/)
        GrpTbl(6,:) = (/Zero,One,-One,Zero,Zero,One,Zero,-One/)
        GrpTbl(7,:) = (/Zero,-One,One,Zero,Zero,One,Zero,-One/)
        GrpTbl(8,:) = (/One,Zero,Zero,-One,One,Zero,-One,Zero/)
       ELSEIF(nGroup.EQ.5) THEN
        theta  = Two * PI / 5.0d0
        theta2 = Two * theta
        Cth  = DCOS(theta)
        C2th = DCOS(theta2)
        Sth  = DSIN(theta)
        S2th = DSIN(theta2)
        GrpTbl(1,:) = One
        GrpTbl(2,1:5)  = One
        GrpTbl(2,6:10) = -One
        GrpTbl(3,:) = (/One,Cth,Cth,C2th,C2th, 
     &                  One,Cth,C2th,C2th,Cth/)
        GrpTbl(4,:) = (/Zero,Sth,-Sth,S2th,-S2th, 
     &                  Zero,-Sth,-S2th,S2th,Sth/)
        GrpTbl(5,:) = (/Zero,-Sth,Sth,-S2th,S2th,
     &                  Zero,-Sth,-S2th,S2th,Sth/)
        GrpTbl(6,:) = GrpTbl(2,:)*GrpTbl(3,:)
        GrpTbl(7,:) = (/One,C2th,C2th,Cth,Cth,
     &                 -One,-C2th,-Cth,-Cth,-C2th/)
        GrpTbl(8,:) = (/Zero,-S2th,S2th,Sth,-Sth,
     &                  Zero,-S2th,Sth,-Sth,S2th/)
        GrpTbl(9,:) = (/Zero,S2th,-S2th,-Sth,Sth,
     &                  Zero,-S2th,Sth,-Sth,S2th/)
        GrpTbl(10,:) = GrpTbl(2,:)*GrpTbl(7,:)
       ELSEIF(nGroup.EQ.6) THEN
        GrpTbl(1,:) = One
	  GrpTbl(2,1:6)  = One
	  GrpTbl(2,7:12) = -One
	  GrpTbl(3,1:12) = (/One,-One,-One,One,One,-One,
     &                    One,-One,One,-One,One,-One/)
	  GrpTbl(4,:) = GrpTbl(2,:)*GrpTbl(3,:)

        GrpTbl(5,:) = (/One,Half,Half,-Half,-Half,-One,
     &                  One,Half,-Half,-One,-Half,Half/)
        GrpTbl(6,:) = (/Zero,Coef,-Coef,Coef,-Coef,Zero,
     &                  Zero,-Coef,-Coef,Zero,Coef,Coef/)
        GrpTbl(7,:) = (/Zero,-Coef,Coef,-Coef,Coef,Zero,
     &                  Zero,-Coef,-Coef,Zero,Coef,Coef/)
        GrpTbl(8,:) = (/One,Half,Half,-Half,-Half,-One,
     &                 -One,-Half,Half,One,Half,-Half/)
        GrpTbl(9,:)  = GrpTbl(5,:)*GrpTbl(3,:)
        GrpTbl(10,:) = GrpTbl(6,:)*GrpTbl(3,:)
        GrpTbl(11,:) = GrpTbl(7,:)*GrpTbl(3,:)
        GrpTbl(12,:) = GrpTbl(8,:)*GrpTbl(3,:)
       ENDIF
c Dn
      ELSEIF(iGroup.EQ.4) THEN
       IF(nGroup.EQ.2) THEN
	  GrpTbl(1,:) = One
	  GrpTbl(2,:) = (/One,One,-One,-One/)
	  GrpTbl(3,:) = (/One,-One,One,-One/)
	  GrpTbl(4,:) = GrpTbl(2,:)*GrpTbl(3,:)
       ELSEIF(nGroup.EQ.3) THEN
        GrpTbl(1,:) = One
        GrpTbl(2,:) = (/One,One,One, -One,-One,-One/)
        GrpTbl(3,:) = (/One,-Half,-Half,-One,Half,Half/)
        GrpTbl(4,:) = (/Zero,Coef,-Coef,Zero,Coef,-Coef/)
        GrpTbl(5,:) = (/Zero,-Coef,Coef,Zero,Coef,-Coef/)
	  GrpTbl(6,:) = GrpTbl(2,:)*GrpTbl(3,:)
       ELSEIF(nGroup.EQ.4) THEN
        GrpTbl(1,:) = One
        GrpTbl(2,:) = (/One,One,One,One,-One,-One,-One,-One/)
	  GrpTbl(3,:) = (/One,-One,-One,One,One,-One,One,-One/)
	  GrpTbl(4,:) = GrpTbl(2,:)*GrpTbl(3,:)
	  GrpTbl(5,:) = (/One,Zero,Zero,-One, One,Zero,-One,Zero/)
	  GrpTbl(6,:) = (/Zero,One,-One,Zero,Zero,-One,Zero,One/)
	  GrpTbl(7,:) = (/Zero,-One,One,Zero,Zero,-One,Zero,One/)
	  GrpTbl(8,:) = GrpTbl(2,:)*GrpTbl(5,:)
       ELSEIF(nGroup.EQ.5) THEN
        theta  = Two * PI / 5.0d0
        theta2 = Two * theta
        Cth  = DCOS(theta)
        C2th = DCOS(theta2)
        Sth  = DSIN(theta)
        S2th = DSIN(theta2)
        GrpTbl(1,:) = One
        GrpTbl(2,1:5)  = One
        GrpTbl(2,6:10) = -One
        GrpTbl(3,:) = (/One,Cth,Cth,C2th,C2th 
     &               -One,-Cth,-C2th,-C2th,-Cth/)
        GrpTbl(4,:) = (/Zero,Sth,-Sth,S2th,-S2th
     &               Zero,Sth,S2th,-S2th,-Sth/)
        GrpTbl(5,:) = (/Zero,-Sth,Sth,-S2th,S2th
     &               Zero,Sth,S2th,-S2th,-Sth/)
        GrpTbl(6,:) = GrpTbl(2,:)*GrpTbl(3,:)

        GrpTbl(7,:) = (/One,C2th,C2th,Cth,Cth
     &               -One,-C2th,-Cth,-Cth,-C2th/)
        GrpTbl(8,:) = (/Zero,-S2th,S2th,Sth,-Sth
     &               Zero,-S2th,Sth,-Sth,S2th/)
        GrpTbl(9,:) = (/Zero,S2th,-S2th,-Sth,Sth
     &               Zero,-S2th,Sth,-Sth,S2th/)
        GrpTbl(10,:) = GrpTbl(2,:)*GrpTbl(7,:)
       ELSEIF(nGroup.EQ.6) THEN
        GrpTbl(1,:) = One
        GrpTbl(2,1:6)  = One
        GrpTbl(2,7:12) = -One
        GrpTbl(3,:) = (/One,-One,-One,One,One,-One,
     &                  One,-One,One,-One,One,-One/)
        GrpTbl(4,:) = GrpTbl(2,:)*GrpTbl(3,:)

        GrpTbl(5,:) = (/One,Half,Half,-Half,-Half,-One
     &                  One,Half,-Half,-One,-Half,Half/)
        GrpTbl(6,:) = (/Zero,Coef,-Coef,Coef,-Coef,Zero,
     &                  Zero,-Coef,-Coef,Zero,Coef,Coef/)
        GrpTbl(7,:) = (/Zero,-Coef,Coef,-Coef,Coef,Zero,
     &                  Zero,-Coef,-Coef,Zero,Coef,Coef/)
        GrpTbl(8,:)  = GrpTbl(2,:)*GrpTbl(5,:)
        GrpTbl(9,:)  = GrpTbl(3,:)*GrpTbl(5,:)
        GrpTbl(10,:) = GrpTbl(3,:)*GrpTbl(6,:)
        GrpTbl(11,:) = GrpTbl(3,:)*GrpTbl(7,:)
        GrpTbl(12,:) = GrpTbl(3,:)*GrpTbl(8,:)
       ENDIF
c Dnh
      ELSEIF(iGroup.EQ.5) THEN
       IF(nGroup.EQ.2) THEN
        GrpTbl(1,:) = One
        GrpTbl(2,1:4) = (/One,One,-One,-One/)
        GrpTbl(2,5:8) = (/One,One,-One,-One/)
        GrpTbl(3,1:4) = (/One,-One,-One,One/)
        GrpTbl(3,5:8) = (/-One,One,One,-One/)
	  GrpTbl(4,:) = GrpTbl(2,:)*GrpTbl(3,:)
	  GrpTbl(5,1:4) = One
	  GrpTbl(5,5:8) = -One
	  GrpTbl(6,:) = GrpTbl(2,:)*GrpTbl(5,:)
	  GrpTbl(7,:) = GrpTbl(3,:)*GrpTbl(5,:)
	  GrpTbl(8,:) = GrpTbl(4,:)*GrpTbl(5,:)
       ELSEIF(nGroup.EQ.3) THEN
        GrpTbl(1,:) = One
        GrpTbl(2,:)  = (/One,One,One,-One,-One,-One,
     &                   One,One,One,-One,-One,-One/)
        GrpTbl(3,:)  = (/One,-Half,-Half,-One,Half,Half,
     &                   One,-Half,-Half,-One,Half,Half/)
        GrpTbl(4,:)  = (/Zero,Coef,-Coef,Zero,Coef,-Coef,
     &                   Zero,Coef,-Coef,Zero,Coef,-Coef/)
        GrpTbl(5,:)  = (/Zero,-Coef,Coef,Zero,Coef,-Coef,
     &                   Zero,-Coef,Coef,Zero,Coef,-Coef/)
        GrpTbl(6,:) = GrpTbl(2,:)*GrpTbl(3,:)
	  GrpTbl(7,1:6)  = One
	  GrpTbl(7,7:12) = -One
        GrpTbl(8,:)  = GrpTbl(2,:)*GrpTbl(7,:)
        GrpTbl(9,:)  = GrpTbl(3,:)*GrpTbl(7,:)
        GrpTbl(10,:) = GrpTbl(4,:)*GrpTbl(7,:)
        GrpTbl(11,:) = GrpTbl(5,:)*GrpTbl(7,:)
        GrpTbl(12,:) = GrpTbl(6,:)*GrpTbl(7,:)
       ELSEIF(nGroup.EQ.4) THEN
        GrpTbl(1,:) = One
        GrpTbl(2,:) = (/One,One,One,One,-One,-One,-One,-One,
     &                  One,One,One,One,-One,-One,-One,-One/)
        GrpTbl(3,:) = (/One,-One,-One,One,One,-One,One,-One,
     &                  One,-One,-One,One,One,-One,One,-One/)
        GrpTbl(4,:) = (/One,-One,-One,One,-One,One,-One,One,
     &                  One,-One,-One,One,-One,One,-One,One/)
        GrpTbl(5,:) = (/One,Zero,Zero,-One,-One,Zero,One,Zero,
     &                 -One,Zero,Zero,One,One,Zero,-One,Zero/)
        GrpTbl(6,:) = (/Zero,One,-One,Zero,Zero,One,Zero,-One,
     &                  Zero,-One,One,Zero,Zero,-One,Zero,One/)
        GrpTbl(7,:) = GrpTbl(6,:)*GrpTbl(4,:)
        GrpTbl(8,:) = GrpTbl(5,:)*GrpTbl(4,:)
        GrpTbl(9,1:8)  = One
        GrpTbl(9,9:16) = -One
        GrpTbl(10,:) = GrpTbl(2,:)*GrpTbl(9,:)
        GrpTbl(11,:) = GrpTbl(3,:)*GrpTbl(9,:)
        GrpTbl(12,:) = GrpTbl(4,:)*GrpTbl(9,:)
        GrpTbl(13,:) = GrpTbl(5,:)*GrpTbl(9,:)
        GrpTbl(14,:) = GrpTbl(6,:)*GrpTbl(9,:)
        GrpTbl(15,:) = GrpTbl(7,:)*GrpTbl(9,:)
        GrpTbl(16,:) = GrpTbl(8,:)*GrpTbl(9,:)
       ELSEIF(nGroup.EQ.5) THEN
        theta  = Two * PI / 5.0d0
        theta2 = Two * theta
        Cth  = DCOS(theta)
        C2th = DCOS(theta2)
        Sth  = DSIN(theta)
        S2th = DSIN(theta2)
        GrpTbl(1,:) = One
        GrpTbl(2,1:5)   = One
        GrpTbl(2,6:10)  = -One
        GrpTbl(2,11:20) = GrpTbl(2,1:10)
	  GrpTbl(3,1:10) = (/One,Cth,Cth,C2th,C2th,
     &                    -One,-Cth,-Cth,-C2th,-Cth/) 
        GrpTbl(3,11:20) = GrpTbl(3,1:10) 
	  GrpTbl(4,1:10) = (/Zero,Sth,-Sth,S2th,-S2th,
     &                     Zero,Sth,S2th,-S2th,-Sth/)
        GrpTbl(4,11:20) = GrpTbl(4,1:10) 
        GrpTbl(5,1:10) = (/Zero,-Sth,Sth,-S2th,S2th,
     &                     Zero,Sth,S2th,-S2th,-Sth/)
        GrpTbl(5,11:20) = GrpTbl(5,1:10) 
        GrpTbl(6,:) = GrpTbl(2,:)*GrpTbl(3,:)
        GrpTbl(7,1:10) = (/One,C2th,C2th,Cth,Cth,
     &                    -One,-C2th,-Cth,-Cth,-C2th/)
        GrpTbl(7,11:20) = GrpTbl(7,1:10)
        GrpTbl(8,1:10) = (/Zero,-S2th,S2th,Sth,-Sth,
     &                     Zero,-S2th,Sth,-Sth,S2th/)
        GrpTbl(8,11:20) = GrpTbl(8,1:10)
        GrpTbl(9,1:10) = (/Zero,S2th,-S2th,-Sth,Sth,
     &                     Zero,-S2th,Sth,-Sth,S2th/)
        GrpTbl(9,11:20) = GrpTbl(9,1:10)
        GrpTbl(10,:) = GrpTbl(2,:)*GrpTbl(7,:)
	  GrpTbl(11,1:10)  = One
	  GrpTbl(11,11:20) = -One
        GrpTbl(12,:) = GrpTbl(2,:)*GrpTbl(11,:)
        GrpTbl(13,:) = GrpTbl(3,:)*GrpTbl(11,:)
        GrpTbl(14,:) = GrpTbl(4,:)*GrpTbl(11,:)
        GrpTbl(15,:) = GrpTbl(5,:)*GrpTbl(11,:)
        GrpTbl(16,:) = GrpTbl(6,:)*GrpTbl(11,:)
        GrpTbl(17,:) = GrpTbl(7,:)*GrpTbl(11,:)
        GrpTbl(18,:) = GrpTbl(8,:)*GrpTbl(11,:)
        GrpTbl(19,:) = GrpTbl(9,:)*GrpTbl(11,:)
        GrpTbl(20,:) = GrpTbl(10,:)*GrpTbl(11,:)
       ELSEIF(nGroup.EQ.6) THEN
        GrpTbl(1,:) = One
	  GrpTbl(2,1:6)  = One
	  GrpTbl(2,7:12) = -One
	  GrpTbl(2,13:24) = GrpTbl(2,1:12)
        GrpTbl(3,1:12) = (/One,-One,-One,One,One,-One,
     &                     One,-One,One,-One,One,-One/) 
        GrpTbl(3,13:24) = -GrpTbl(3,1:12)
	  GrpTbl(4,:) = GrpTbl(2,:)*GrpTbl(3,:)

	  GrpTbl(5,1:12) = (/One,Half,Half,-Half,-Half,-One,
     &                    -One,-Half,Half,One,Half,-Half/)
	  GrpTbl(5,13:24) = -GrpTbl(5,1:12)
        GrpTbl(6,1:12) = (/Zero,Coef,-Coef,Coef,-Coef,Zero,
     &                     Zero,Coef,Coef,Zero,-Coef,-Coef/)
	  GrpTbl(6,13:24) = -GrpTbl(6,1:12)
        GrpTbl(7,1:12) = (/Zero,-Coef,Coef,-Coef,Coef,Zero,
     &                     Zero,Coef,Coef,Zero,-Coef,-Coef/)
        GrpTbl(7,13:24) = -GrpTbl(7,1:12)
	  GrpTbl(8,:)  = GrpTbl(2,:)*GrpTbl(5,:)
	  GrpTbl(9,:)  = GrpTbl(3,:)*GrpTbl(5,:)
	  GrpTbl(10,:) = GrpTbl(3,:)*GrpTbl(6,:)
	  GrpTbl(11,:) = GrpTbl(3,:)*GrpTbl(7,:)
	  GrpTbl(12,:) = GrpTbl(3,:)*GrpTbl(8,:)
	  GrpTbl(13,1:12)  = One
	  GrpTbl(13,13:24) = -One
        GrpTbl(14,:) = GrpTbl(13,:)*GrpTbl(2,:)
        GrpTbl(15,:) = GrpTbl(13,:)*GrpTbl(3,:)
        GrpTbl(16,:) = GrpTbl(13,:)*GrpTbl(4,:)
        GrpTbl(17,:) = GrpTbl(13,:)*GrpTbl(5,:)
        GrpTbl(18,:) = GrpTbl(13,:)*GrpTbl(6,:)
        GrpTbl(19,:) = GrpTbl(13,:)*GrpTbl(7,:)
        GrpTbl(20,:) = GrpTbl(13,:)*GrpTbl(8,:)
        GrpTbl(21,:) = GrpTbl(13,:)*GrpTbl(9,:)
        GrpTbl(22,:) = GrpTbl(13,:)*GrpTbl(10,:)
        GrpTbl(23,:) = GrpTbl(13,:)*GrpTbl(11,:)
        GrpTbl(24,:) = GrpTbl(13,:)*GrpTbl(12,:)
       ENDIF
c Dnd
      ELSEIF(iGroup.EQ.6) THEN
       IF(nGroup.EQ.2) THEN
        GrpTbl(1,:) = One
	  GrpTbl(2,:) = (/One,One,-One,-One,One,One,-One,-One/)
	  GrpTbl(3,:) = (/One,One,One,One,-One,-One,-One,-One/)
        GrpTbl(4,:) = GrpTbl(2,:)*GrpTbl(3,:)
	  GrpTbl(5,:) = (/One,-One,-One,One,Zero,Zero,Zero/)
        GrpTbl(6,:) = (/Zero,Zero,Zero,Zero,-One,One,-One,One/)
        GrpTbl(7,:) = (/Zero,Zero,Zero,Zero,One,-One,-One,One/)
        GrpTbl(8,:) = GrpTbl(2,:)*GrpTbl(5,:)
       ELSEIF(nGroup.EQ.3) THEN
        GrpTbl(1,:) = One
        GrpTbl(2,1:6) = (/One,One,One,-One,-One,-One/)
	  GrpTbl(2,7:12) = GrpTbl(2,1:6)
        GrpTbl(3,:) = (/One,-Half,-Half,One,-Half,-Half,
     &                  -Half,-Half,One,-Half,One,-Half/)
        GrpTbl(4,:) = (/Zero,Coef,-Coef,Zero,-Coef,Coef,
     &                  -Coef,Coef,Zero,Coef,Zero,-Coef/)
        GrpTbl(5,:) = (/Zero,-Coef,Coef,Zero,-Coef,Coef,
     &                  Coef,-Coef,Zero,Coef,Zero,-Coef/)
        GrpTbl(6,:) = GrpTbl(2,:)*GrpTbl(3,:)
        GrpTbl(7,1:6)  = One
        GrpTbl(7,7:12) = -One
        GrpTbl(8,:)  = GrpTbl(2,:)*GrpTbl(7,:)
        GrpTbl(9,:)  = GrpTbl(3,:)*GrpTbl(7,:)
        GrpTbl(10,:) = GrpTbl(4,:)*GrpTbl(7,:)
        GrpTbl(11,:) = GrpTbl(5,:)*GrpTbl(7,:)
        GrpTbl(12,:) = GrpTbl(6,:)*GrpTbl(7,:)
       ELSEIF(nGroup.EQ.4) THEN
	  Cf = One/DSQRT(Two)
        GrpTbl(1,:) = One
	  GrpTbl(2,1:8) = (/One,One,One,One,-One,-One,-One,-One/)
	  GrpTbl(2,9:16) = GrpTbl(2,1:8)
	  GrpTbl(3,1:8)  = One
	  GrpTbl(3,9:16) = -One
	  GrpTbl(4,:) = GrpTbl(2,:)*GrpTbl(3,:)
	  GrpTbl(5,:) = (/One,Zero,Zero,-One,Zero,One,Zero,-One,
     &                  Cf,Cf,-Cf,-Cf,Cf,Cf,-Cf,-Cf/)
        GrpTbl(6,:) = (/Zero,One,-One,Zero,One,Zero,-One,Zero,
     &                  Cf,-Cf,Cf,-Cf,Cf,-Cf,-Cf,Cf/)
        GrpTbl(7,:) = (/Zero,-One,One,Zero,One,Zero,-One,Zero,
     &                  -Cf,Cf,-Cf,Cf,Cf,-Cf,-Cf,Cf/)
        GrpTbl(8,:)  = GrpTbl(2,:)*GrpTbl(5,:)
        GrpTbl(9,1:8)  = (/One,-One,-One,One,One,-One,One,-One/)
        GrpTbl(9,9:16) = Zero
        GrpTbl(10,1:8)  = Zero
        GrpTbl(10,9:16) = (/-One,One,One,-One,One,-One,One,-One/)
        GrpTbl(11,1:8)  = Zero
        GrpTbl(11,9:16) = (/One,-One,-One,One,One,-One,One,-One/)
        GrpTbl(12,:) = GrpTbl(2,:)*GrpTbl(9,:)
        GrpTbl(13,:) = GrpTbl(3,:)*GrpTbl(5,:)
        GrpTbl(14,:) = GrpTbl(3,:)*GrpTbl(6,:)
        GrpTbl(15,:) = GrpTbl(3,:)*GrpTbl(7,:)
        GrpTbl(16,:) = GrpTbl(3,:)*GrpTbl(8,:)
       ELSEIF(nGroup.EQ.5) THEN
        theta  = Two * PI / 5.0d0
        theta2 = Two * theta
        Cth  = DCOS(theta)
        C2th = DCOS(theta2)
        Sth  = DSIN(theta)
        S2th = DSIN(theta2)
        GrpTbl(1,:) = One
	GrpTbl(2,1:10) = (/One,One,One,One,One,-One,-One,-One,-One,-One/)
	  GrpTbl(2,11:20) = GrpTbl(2,1:10)
	  GrpTbl(3,:) = (/One,Cth,Cth,C2th,C2th,One,Cth,C2th,C2th,Cth,
     &                  C2th,C2th,Cth,Cth,One,C2th,Cth,One,Cth,C2th/)
        GrpTbl(4,:) = (/Zero,Sth,-Sth,S2th,-S2th,Zero,
     &    -Sth,-S2th,S2th,Sth, -S2th,S2th,-Sth,Sth,Zero,
     &     S2th,Sth,Zero,-Sth,-S2th/)
        GrpTbl(5,:) = (/Zero,-Sth,Sth,-S2th,S2th,Zero,
     &     -Sth,-S2th,S2th,Sth, S2th,-S2th,Sth,-Sth,Zero,
     &      S2th,Sth,Zero,-Sth,-S2th/)
        GrpTbl(6,:) = GrpTbl(2,:)*GrpTbl(3,:)
        GrpTbl(7,:) = (/One,C2th,C2th,Cth,Cth,
     &     -One,-C2th,-Cth,-Cth,-C2th,
     &      Cth,Cth,C2th,C2th,One,-Cth,-C2th,-One,-C2th,-Cth/)
        GrpTbl(8,:) = (/Zero,-S2th,S2th,Sth,-Sth,
     &     Zero,-S2th,Sth,-Sth,S2th,
     &     -Sth,Sth,S2th,-S2th,Zero,-Sth,S2th,Zero,-S2th,Sth/)
        GrpTbl(9,:) = (/Zero,S2th,-S2th,-Sth,Sth,
     &     Zero,-S2th,Sth,-Sth,S2th,
     &     Sth,-Sth,-S2th,S2th,Zero,-Sth,S2th,Zero,-S2th,Sth/)
        GrpTbl(10,:) = GrpTbl(2,:)*GrpTbl(7,:)
        GrpTbl(11,1:10) = One
        GrpTbl(11,11:20) = -One
        GrpTbl(12,:) = GrpTbl(11,:)*GrpTbl(2,:)
        GrpTbl(13,:) = GrpTbl(11,:)*GrpTbl(3,:)
        GrpTbl(14,:) = GrpTbl(11,:)*GrpTbl(4,:)
        GrpTbl(15,:) = GrpTbl(11,:)*GrpTbl(5,:)
        GrpTbl(16,:) = GrpTbl(11,:)*GrpTbl(6,:)
        GrpTbl(17,:) = GrpTbl(11,:)*GrpTbl(7,:)
        GrpTbl(18,:) = GrpTbl(11,:)*GrpTbl(8,:)
        GrpTbl(19,:) = GrpTbl(11,:)*GrpTbl(9,:)
        GrpTbl(20,:) = GrpTbl(11,:)*GrpTbl(10,:)
       ELSEIF(nGroup.EQ.6) THEN
        GrpTbl(1,:) = One
	  GrpTbl(2,1:6)   = One
	  GrpTbl(2,7:12)  = -One
	  GrpTbl(2,13:18) = One
	  GrpTbl(2,19:24) = -One
	  GrpTbl(3,1:12)  = One
	  GrpTbl(3,13:24) = -One
        GrpTbl(4,:) = GrpTbl(2,:)*GrpTbl(3,:)
        GrpTbl(5,:) = (/One,Half,Half,-Half,-Half,-One,
     &                 -One,-Half,Half,One,Half,-Half, 
     &                  Coef,Coef,Zero,Zero,-Coef,-Coef
     &                  -Coef,Zero,Coef,Coef,Zero,-Coef/)
        GrpTbl(6,:) = (/Zero,Coef,-Coef,Coef,-Coef,Zero,
     &                  Zero,Coef,Coef,Zero,-Coef,-Coef,
     &                  Half,-Half,One,-One,Half,-Half,
     &                  Half,One,Half,-Half,-One,-Half/)
        GrpTbl(7,:) = (/Zero,-Coef,Coef,-Coef,Coef,Zero,
     &                  Zero,Coef,Coef,Zero,-Coef,-Coef,
     &                 -Half,Half,-One,One,-Half,Half,
     &                  Half,One,Half,-Half,-One,-Half/)
        GrpTbl(8,:) = GrpTbl(2,:)*GrpTbl(5,:)
        GrpTbl(9,:) = (/One,-Half,-Half,-Half,-Half,One,
     &                 -One,Half,Half,-One,Half,Half,
     &                  Half,Half,-One,-One,Half,Half,
     &                 -Half,One,-Half,-Half,One,-Half/)
        GrpTbl(10,:) = (/Zero,-Coef,Coef,Coef,-Coef,Zero,
     &                  Zero,-Coef,Coef,Zero,-Coef,Coef,
     &                 -Coef,Coef,Zero,Zero,Coef,-Coef,
     &                 -Coef,Zero,Coef,-Coef,Zero,Coef/)
        GrpTbl(11,:) = (/Zero,Coef,-Coef,-Coef,Coef,Zero,
     &                  Zero,-Coef,Coef,Zero,-Coef,Coef,
     &                  Coef,-Coef,Zero,Zero,-Coef,Coef,
     &                 -Coef,Zero,Coef,-Coef,Zero,Coef/)
        GrpTbl(12,:) = GrpTbl(9,:)*GrpTbl(2,:)
        GrpTbl(13,1:6) = (/One,-One,-One,One,One,-One/)
        GrpTbl(13,7:18) = Zero
        GrpTbl(13,19:24) = (/-One,One,-One,One,-One,One/)
        GrpTbl(14,1:6) = Zero
        GrpTbl(14,7:18) = (/-One,One,-One,One,-One,One,
     &                       One,-One,-One,One,One,-One/)
        GrpTbl(14,19:24) = Zero
        GrpTbl(15,1:6) = Zero
        GrpTbl(15,7:18) = (/-One,One,-One,One,-One,One,
     &                      -One,One,One,-One,-One,One/)
        GrpTbl(15,19:24) = Zero
        GrpTbl(16,:) = GrpTbl(13,:)*GrpTbl(2,:)
        GrpTbl(17,:) = GrpTbl(3,:)*GrpTbl(9,:)
        GrpTbl(18,:) = GrpTbl(3,:)*GrpTbl(10,:)
        GrpTbl(19,:) = GrpTbl(3,:)*GrpTbl(11,:)
        GrpTbl(20,:) = GrpTbl(3,:)*GrpTbl(12,:)
        GrpTbl(21,:) = GrpTbl(3,:)*GrpTbl(5,:)
        GrpTbl(22,:) = GrpTbl(3,:)*GrpTbl(6,:)
        GrpTbl(23,:) = GrpTbl(3,:)*GrpTbl(7,:)
        GrpTbl(24,:) = GrpTbl(3,:)*GrpTbl(8,:)
       ENDIF
c Sn
      ELSEIF(iGroup.EQ.7) THEN
       IF(nGroup.EQ.4) THEN
	  GrpTbl(1,:) = One 
	  GrpTbl(2,:) = (/One,-One,-One,One/) 
	  GrpTbl(3,:) = (/One,Zero,Zero,-One/) 
	  GrpTbl(4,:) = (/Zero,-One,One,Zero/)
       ELSEIF(nGroup.EQ.6) THEN
	  GrpTbl(1,:) = One 
	  GrpTbl(2,:) = (/One,(-Half,i=2,5),One/) 
	  GrpTbl(3,:) = (/Zero,-Coef,Coef,Coef,-Coef,Zero/)
	  GrpTbl(4,:) = (/One,-One,-One,One,One,-One/)
        GrpTbl(5,:) = GrpTbl(2,:)*GrpTbl(4,:)
        GrpTbl(6,:) = GrpTbl(3,:)*GrpTbl(4,:)
       ELSEIF(nGroup.EQ.8) THEN
	  Cf = One/DSQRT(Two)
	  GrpTbl(1,:) = One 
	  GrpTbl(2,:) = (/One,-One,-One,One,
     &                  One,-One,-One,One/)
	  GrpTbl(3,:) = (/One,Cf,Cf,Zero,
     &                 Zero,-Cf,-Cf,-One/)
        GrpTbl(4,:) = (/Zero,Cf,-Cf,One
     &                 -One,Cf,-Cf,Zero/)
        GrpTbl(5,:) = (/One,Zero,Zero,-One,
     &                 -One,Zero,Zero,One/)
	  GrpTbl(6,:) = (/Zero,-One,One,Zero,
     &                 Zero,One,-One,Zero/)
        GrpTbl(7,:) = GrpTbl(3,:)*GrpTbl(2,:)
        GrpTbl(8,:) = GrpTbl(4,:)*GrpTbl(2,:)
       ELSEIF(nGroup.EQ.10) THEN
        theta  = Two * PI / 5.0d0
        theta2 = Two * theta
        Cth  = DCOS(theta)
        C2th = DCOS(theta2)
        Sth  = DSIN(theta)
        S2th = DSIN(theta2)
	  GrpTbl(1,:) = One 
	  GrpTbl(2,:) = (/One,C2th,C2th,Cth,Cth,
     &                  Cth,Cth,C2th,C2th,One/)
        GrpTbl(3,:) = (/Zero,-S2th,S2th,Sth,-Sth,
     &                 -Sth,Sth,S2th,-S2th,Zero/)
        GrpTbl(4,:) = (/One,Cth,Cth,C2th,C2th,
     &                  C2th,C2th,Cth,Cth,One/)
        GrpTbl(5,:) = (/Zero,-Sth,Sth,-S2th,S2th,
     &                 S2th,-S2th,Sth,-Sth,Zero/)
        GrpTbl(6,:) = (/One,-One,-One,One,One,
     &                 -One,-One,One,One,-One/)
        GrpTbl(7,:) = GrpTbl(2,:)*GrpTbl(6,:)
        GrpTbl(8,:) = GrpTbl(3,:)*GrpTbl(6,:)
        GrpTbl(9,:) = GrpTbl(4,:)*GrpTbl(6,:)
        GrpTbl(10,:)= GrpTbl(5,:)*GrpTbl(6,:)
       ELSEIF(nGroup.EQ.12) THEN
	  GrpTbl(1,:) = One 
	  GrpTbl(2,:) = (/One,-One,-One,One,One,-One,
     &                 -One,One,One,-One,-One,One/) 
        GrpTbl(3,:) = (/One,Coef,Coef,Half,Half,Zero,
     &             Zero,-Half,-Half,-Coef,-Coef,-One/)
	  GrpTbl(4,:) = (/Zero,Half,-Half,Coef,-Coef,One,
     &             -One,Coef,-Coef,Half,-Half,Zero/)
        GrpTbl(5,:) = (/One,Half,Half,-Half,-Half,-One,
     &             -One,-Half,-Half,Half,Half,One/)
        GrpTbl(6,:) = (/Zero,-Coef,Coef,-Coef,Coef,Zero,
     &             Zero,Coef,-Coef,Coef,-Coef,Zero/)
        GrpTbl(7,:) = (/One,Zero,Zero,-One,-One,Zero,
     &             Zero,One,One,Zero,Zero,-One/)
        GrpTbl(8,:) = (/Zero,-One,One,Zero,Zero,One,
     &             -One,Zero,Zero,-One,One,Zero/)
        GrpTbl(9,:) = GrpTbl(5,:)*GrpTbl(2,:)
        GrpTbl(10,:)= GrpTbl(6,:)*GrpTbl(2,:)
        GrpTbl(11,:)= GrpTbl(3,:)*GrpTbl(2,:)
        GrpTbl(12,:)= GrpTbl(4,:)*GrpTbl(2,:)
       ENDIF
c Cubic
      ELSEIF(iGroup.EQ.8) THEN
c Td
       IF(nGroup.EQ.3) THEN
        GrpTbl(1,:) = One
        GrpTbl(2,1:12)  = One
        GrpTbl(2,13:24) = -One
	  GrpTbl(3,:) = (/One,(-Half,i=2,9),(One,i=10,14),
     &   (-Half,i=15,18),One,-Half,-Half,One,-Half,-Half/)
        GrpTbl(4,:) = (/Zero,Coef,-Coef,Coef,-Coef,Coef,-Coef,
     &   Coef,-Coef,(Zero,i=10,14),Coef,Coef,-Coef,-Coef,
     &   Zero,Coef,-Coef,Zero,Coef,-Coef/)
        GrpTbl(5,:) = (/Zero,-Coef,Coef,-Coef,Coef,-Coef,Coef,
     &   Coef,-Coef,(Zero,i=10,14),Coef,Coef,-Coef,-Coef,
     &   Zero,Coef,-Coef,Zero,Coef,-Coef/)
        GrpTbl(6,:) = GrpTbl(2,:)*GrpTbl(3,:)

	  GrpTbl(7,:) = (/One,(Zero,i=2,9),-One,-One,One,
     &   (Zero,i=13,16),One,One,Zero,Zero,-One,Zero,Zero,-One/)
        GrpTbl(8,:) = (/Zero,-One,Zero,-One,Zero,One,Zero,One,Zero,
     &   Zero,Zero,Zero,One,-One,(Zero,i=15,18),
     &   One,Zero,Zero,-One,Zero,Zero/)
        GrpTbl(9,:) = (/Zero,Zero,One,Zero,-One,Zero,One,Zero,-One,
     &   (Zero,i=10,14),One,-One,Zero,Zero,
     &   Zero,-One,Zero,Zero,One,Zero/)
        GrpTbl(10,:) = (/Zero,Zero,-One,Zero,-One,Zero,One,Zero,One,
     &   Zero,Zero,Zero,-One,One,(Zero,i=15,18),
     &   One,Zero,Zero,-One,Zero,Zero/)
        GrpTbl(11,:) = (/One,(Zero,i=2,9),
     &   -One,One,-One,Zero,Zero,One,One,Zero,Zero,
     &   Zero,-One,Zero,Zero,-One,Zero/)
        GrpTbl(12,:) = (/Zero,-One,Zero,One,Zero,One,Zero,-One,Zero,
     &   (Zero,i=10,16),One,-One,Zero,Zero,One,Zero,Zero,-One/)
        GrpTbl(13,:) = (/Zero,One,Zero,-One,Zero,One,Zero,-One,Zero,
     &   (Zero,i=10,14),-One,One,Zero,Zero,
     &    Zero,-One,Zero,Zero,One,Zero/)
        GrpTbl(14,:) = (/Zero,Zero,-One,Zero,One,Zero,One,Zero,-One,
     &   (Zero,i=10,16),-One,One,Zero,Zero,One,Zero,Zero,-One/)
        GrpTbl(15,:) = (/One,(Zero,i=2,9),One,-One,-One,One,One,
     &   (Zero,i=15,18),-One,Zero,Zero,-One,Zero,Zero/)
	  DO i = 7,15
          GrpTbl(i+9,:) = GrpTbl(2,:)*GrpTbl(i,:)
 	  ENDDO
c Oh
       ELSEIF(nGroup.EQ.5) THEN
        GrpTbl(1,:) = One
	  GrpTbl(2,:) = (/(One,i=1,9),(-One,i=10,21),(One,i=22,25),
     &   (-One,i=26,31),(One,i=32,42),(-One,i=43,48)/)
	  GrpTbl(3,:) = (/One,(-Half,i=2,13),(One,i=14,17),
     &   (-Half,i=18,21),(One,i=22,27),(-Half,i=28,39),
     &   (One,i=40,42),(-Half,i=43,46), One,One/)
c                        1      2    3    4    5    6     7    8
        GrpTbl(4,:) = (/Zero,-Coef,Coef,Coef,-Coef,Coef,-Coef,-Coef,
c          9    10   11   12    13   14    15   16
     &    Coef,Coef,Coef,-Coef,-Coef,Zero,Zero,Zero,
c          17    18   19   20    21   22    23   24
     &    Zero,-Coef,-Coef,Coef,Coef,Zero,Zero,Zero,
c          25    26   27   28    29   30    31   32
     &    Zero,Zero,Zero,-Coef,-Coef,Coef,Coef,Coef,
c          33    34   35   36    37   38    39   40
     &    -Coef,-Coef,Coef,-Coef,Coef,Coef,-Coef,Zero,
c          41    42   43   44    45   46    47   48
     &    Zero,Zero,Coef,Coef,-Coef,-Coef,Zero,Zero/)
c                        1      2    3    4    5    6     7    8
        GrpTbl(5,:) = (/Zero,Coef,-Coef,-Coef,Coef,-Coef,Coef,Coef,
c          9    10   11   12    13   14    15   16
     &    -Coef,Coef,Coef,-Coef,-Coef,Zero,Zero,Zero,
c          17    18   19   20    21   22    23   24
     &    Zero,-Coef,-Coef,Coef,Coef,Zero,Zero,Zero,
c          25    26   27   28    29   30    31   32
     &    Zero,Zero,Zero,-Coef,-Coef,Coef,Coef,-Coef,
c          33    34   35   36    37   38    39   40
     &    Coef,Coef,-Coef,Coef,-Coef,-Coef,Coef,Zero,
c          41    42   43   44    45   46    47   48
     &    Zero,Zero,Coef,Coef,-Coef,-Coef,Zero,Zero/)
        GrpTbl(6,:) = GrpTbl(2,:)*GrpTbl(3,:)
        GrpTbl(7,:) = (/One,(Zero,i=2,11),-One,-One,
     &   (Zero,i=14,17),One,One,Zero,Zero,-One,One,-One,
     &   One,Zero,Zero,One,One,(Zero,i=30,39),-One,One,-One,
     &   Zero,Zero,-One,-One,Zero,Zero/)
        GrpTbl(8,:) = (/Zero,One,Zero,Zero,-One,Zero,-One,One,
     &   (Zero,i=9,13),One,-One,One,-One,(Zero,i=18,25),
     &   -One,One,(Zero,i=28,32),One,-One,Zero,-One,
     &    Zero,Zero,One,(Zero,i=40,46),One,-One/)
        GrpTbl(9,:) = (/Zero,Zero,One,-One,Zero,One,Zero,Zero
     &    -One,One,-One,(Zero,i=12,19),-One,One,
     &    (Zero,i=22,29),One,-One,One,Zero,Zero,-One,Zero,
     &     One,-One,(Zero,i=39,42),One,-One,(Zero,i=45,48)/)
        GrpTbl(10,:) = (/Zero,Zero,One,-One,Zero,-One,Zero,Zero
     &    One,(Zero,i=10,13),One,-One,-One,One,(Zero,i=18,25),
     &    One,-One,(Zero,i=28,31),One,Zero,Zero,-One,Zero,
     &   -One,One,(Zero,i=39,46),One,-One/)
        GrpTbl(11,:) = (/One,(Zero,i=2,9),-One,-One,
     &   (Zero,i=12,19),One,One,-One,-One,One,One,
     &   (Zero,i=26,29),One,One,(Zero,i=32,39),-One,-One,One,
     &   -One,-One,(Zero,i=45,48)/)
        GrpTbl(12,:) = (/Zero,One,Zero,Zero,One,Zero,-One,-One,
     &   Zero,Zero,Zero,One,-One,(Zero,i=14,17),One,-One,
     &   (Zero,i=20,27),-One,One,Zero,Zero,Zero,One,One,Zero,-One,
     &   Zero,Zero,-One,(Zero,i=40,44),One,-One,Zero,Zero/)
        GrpTbl(13,:) = (/Zero,One,Zero,Zero,-One,Zero,One,-One,
     &   Zero,One,-One,(Zero,i=12,19),One,-One,(Zero,i=22,29),
     &   -One,One,Zero,One,-One,Zero,One,Zero,Zero,-One,
     &   Zero,Zero,Zero,One,-One,(Zero,i=45,48)/)
        GrpTbl(14,:) = (/Zero,Zero,One,One,Zero,-One,Zero,Zero,
     &   -One,Zero,Zero,One,-One,(Zero,i=14,17),-One,One,
     &   (Zero,i=20,27),One,-One,Zero,Zero,One,Zero,Zero,One,Zero,
     &   -One,-One,(Zero,i=39,44),One,-One,Zero,Zero/)
        GrpTbl(15,:) = (/One,(Zero,i=2,13),-One,-One,One,One,
     &   (Zero,i=18,21),One,-One,-One,One,One,One,(Zero,i=28,39),
     &   One,-One,-One,(Zero,i=43,46),-One,-One/)
        DO i = 7,15
          GrpTbl(i+9,:) = GrpTbl(2,:)*GrpTbl(i,:)
        ENDDO
	  GrpTbl(25,:) = (/(One,i=1,24),(-One,i=25,48)/)
        DO i = 2,24
	    GrpTbl(24+i,:) = GrpTbl(25,:)*GrpTbl(i,:)
	  ENDDO

	 ELSE
         WRITE(*,*)'ERROR. Only Td of Oh (Higher) groups ready'
         STOP
       ENDIF
c Linear
      ELSEIF(iGroup.EQ.9) THEN
        WRITE(iout,*)'ERROR. Irred. rep. not ready for Linear groups'
        STOP
	ENDIF

	RETURN
	END


      SUBROUTINE group_labels(IRLab,nIR,iGroup,nGroup,iout)
c*****************************************************c
c                                                     c
c       Makes Irred. Rep. labels of a                 c
c       symmetry point group                          c
c                                                     c
c       COMMON                                        c
c        iGroup  : group index                        c
c        nGroup  : group order                        c
c                                                     c
c       INPUT                                         c
c        nIR     : # of irreducible rep.              c
c                                                     c
c       OUTPUT                                        c
c        IRLab   : labels of irreducible rep.         c
c                                                     c
c       DCC: 10/10                                    c
c                                                     c
c*****************************************************c
      IMPLICIT REAL*8  (a-h, o-z)
	INTEGER    iGroup,nGroup
      CHARACTER  IRLab(nIR)*3

      SAVE Zero,One,Two
      DATA Zero,One,Two / 0d0,1d0,2d0 /


c Only up to 6th order (nGroup < 7)
      IF(iGroup .NE. 7) THEN
       IF(nGroup.LT.2 .OR. nGroup.GT.6) THEN
        WRITE(iout,*)'ERROR. Irreducible rep. only for n < 7'
        STOP
       ENDIF
      ENDIF

c***********************c
c Make character tables c
c***********************c
c Nonaxial
	IF(iGroup.EQ.0) THEN
	 IF(nGroup.EQ.1) THEN
        WRITE(iout,*)'ERROR. Irred. rep. not ready for C1'
        STOP
	 ELSEIF(nGroup.EQ.2) THEN
	   IRLab(1) = 'Ag '
	   IRLab(2) = 'Au '
	 ELSEIF(nGroup.EQ.3) THEN
	   IRLab(1) = "A' "
	   IRLab(2) = 'A" '
	 ENDIF
c Cn
	ELSEIF(iGroup.EQ.1) THEN
	 IF(nGroup.EQ.2) THEN
         IRLab(1) = 'A  '
         IRLab(2) = 'B  '
	 ELSEIF(nGroup.EQ.3) THEN
         IRLab(1) = 'A  '
         IRLab(2) = 'E  '
	 ELSEIF(nGroup.EQ.4) THEN
         IRLab(1) = 'A  '
         IRLab(2) = 'B  '
         IRLab(3) = 'E  '
	 ELSEIF(nGroup.EQ.5) THEN
         IRLab(1) = 'A  '
         IRLab(2) = 'E1 '
         IRLab(3) = 'E2 '
	 ELSEIF(nGroup.EQ.6) THEN
         IRLab(1) = 'A  '
         IRLab(2) = 'B  '
         IRLab(3) = 'E1 '
         IRLab(4) = 'E2 '
	 ENDIF
c Cnh
	ELSEIF(iGroup.EQ.2) THEN
       IF(nGroup.EQ.2) THEN
         IRLab(1) = 'Ag '
         IRLab(2) = 'Bg '
         IRLab(3) = 'Au ' 
         IRLab(4) = 'Bu '
       ELSEIF(nGroup.EQ.3) THEN
         IRLab(1) = "A' "
         IRLab(2) = "E' "
         IRLab(3) = 'A" '
         IRLab(4) = 'E" '
       ELSEIF(nGroup.EQ.4) THEN
         IRLab(1) = 'Ag '
         IRLab(2) = 'Bg '
         IRLab(3) = 'Eg '
         IRLab(4) = 'Au '
         IRLab(5) = 'Bu '
         IRLab(6) = 'Eu '
       ELSEIF(nGroup.EQ.5) THEN
         IRLab(1) = "A' "
         IRLab(2) = "E1'"
         IRLab(3) = "E2'"
         IRLab(4) = 'A" '
         IRLab(5) = 'E1"'
         IRLab(6) = 'E2"'
       ELSEIF(nGroup.EQ.6) THEN
         IRLab(1) = 'Ag '
         IRLab(2) = 'Bg '
         IRLab(3) = 'E1g'
         IRLab(4) = 'E2g'
         IRLab(5) = 'Au '
         IRLab(6) = 'Bu '
         IRLab(7) = 'E1u'
         IRLab(8) = 'E2u'
       ENDIF
c Cnv
	ELSEIF(iGroup.EQ.3) THEN
       IF(nGroup.EQ.2) THEN
         IRLab(1) = 'A1 '
         IRLab(2) = 'A2 '
         IRLab(3) = 'B1 '
         IRLab(4) = 'B2 '
       ELSEIF(nGroup.EQ.3) THEN
         IRLab(1) = 'A1 '
         IRLab(2) = 'A2 '
         IRLab(3) = 'E  '
       ELSEIF(nGroup.EQ.4) THEN
         IRLab(1) = 'A1 '
         IRLab(2) = 'A2 '
         IRLab(3) = 'B1 '
         IRLab(4) = 'B2 '
         IRLab(5) = 'E  '
       ELSEIF(nGroup.EQ.5) THEN
         IRLab(1) = 'A1 '
         IRLab(2) = 'A2 '
         IRLab(3) = 'E1 '
         IRLab(4) = 'E2 '
       ELSEIF(nGroup.EQ.6) THEN
         IRLab(1) = 'A1 '
         IRLab(2) = 'A2 '
         IRLab(3) = 'B1 '
         IRLab(4) = 'B2 '
         IRLab(5) = 'E1 '
         IRLab(6) = 'E2 '
       ENDIF
c Dn
	ELSEIF(iGroup.EQ.4) THEN
       IF(nGroup.EQ.2) THEN
         IRLab(1) = 'A  '
         IRLab(2) = 'B1 '
         IRLab(3) = 'B2 '
         IRLab(4) = 'B3 '
       ELSEIF(nGroup.EQ.3) THEN
         IRLab(1) = 'A1 '
         IRLab(2) = 'A2 '
         IRLab(3) = 'E  '
       ELSEIF(nGroup.EQ.4) THEN
         IRLab(1) = 'A1 '
         IRLab(2) = 'A2 '
         IRLab(3) = 'B1 '
         IRLab(4) = 'B2 '
         IRLab(5) = 'E  '
       ELSEIF(nGroup.EQ.5) THEN
         IRLab(1) = 'A1 '
         IRLab(2) = 'A2 '
         IRLab(3) = 'E1 '
         IRLab(4) = 'E2 '
       ELSEIF(nGroup.EQ.6) THEN
         IRLab(1) = 'A1 '
         IRLab(2) = 'A2 '
         IRLab(3) = 'B1 '
         IRLab(4) = 'B2 '
         IRLab(5) = 'E1 '
         IRLab(6) = 'E2 '
       ENDIF
c Dnh
	ELSEIF(iGroup.EQ.5) THEN
       IF(nGroup.EQ.2) THEN
         IRLab(1) = 'Ag '
         IRLab(2) = 'B1g'
         IRLab(3) = 'B2g'
         IRLab(4) = 'B3g'
         IRLab(5) = 'Au '
         IRLab(6) = 'B1u'
         IRLab(7) = 'B2u'
         IRLab(8) = 'B3u'
       ELSEIF(nGroup.EQ.3) THEN
         IRLab(1) = "A1'"
         IRLab(2) = "A2'"
         IRLab(3) = "E' "
         IRLab(4) = 'A1"'
         IRLab(5) = 'A2"'
         IRLab(6) = 'E" '
       ELSEIF(nGroup.EQ.4) THEN
         IRLab(1)  = 'A1g'
         IRLab(2)  = 'A2g'
         IRLab(3)  = 'B1g'
         IRLab(4)  = 'B2g'
         IRLab(5)  = 'Eg '
         IRLab(6)  = 'A1u'
         IRLab(7)  = 'A2u'
         IRLab(8)  = 'B1u'
         IRLab(9)  = 'B2u'
         IRLab(10) = 'Eu '
       ELSEIF(nGroup.EQ.5) THEN
         IRLab(1) = "A1'"
         IRLab(2) = "A2'"
         IRLab(3) = "E1'"
         IRLab(4) = "E2'"
         IRLab(5) = 'A1"'
         IRLab(6) = 'A2"'
         IRLab(7) = 'E1"'
         IRLab(8) = 'E2"'
       ELSEIF(nGroup.EQ.6) THEN
         IRLab(1)  = 'A1g'
         IRLab(2)  = 'A2g'
         IRLab(3)  = 'B1g'
         IRLab(4)  = 'B2g'
         IRLab(5)  = 'E1g'
         IRLab(6)  = 'E2g'
         IRLab(7)  = 'A1u'
         IRLab(8) = 'A2u'
         IRLab(9) = 'B1u'
         IRLab(10) = 'B2u'
         IRLab(11) = 'E1u'
         IRLab(12) = 'E2u'
       ENDIF
c Dnd
	ELSEIF(iGroup.EQ.6) THEN
       IF(nGroup.EQ.2) THEN
         IRLab(1) = 'A1 '
         IRLab(2) = 'A2 '
         IRLab(3) = 'B1 '
         IRLab(4) = 'B2 '
         IRLab(5) = 'E  '
       ELSEIF(nGroup.EQ.3) THEN
         IRLab(1) = 'A1g'
         IRLab(2) = 'A2g'
         IRLab(3) = 'Eg '
         IRLab(4) = 'A1u'
         IRLab(5) = 'A2u'
         IRLab(6) = 'Eu '
       ELSEIF(nGroup.EQ.4) THEN
         IRLab(1) = 'A1 '
         IRLab(2) = 'A2 '
         IRLab(3) = 'B1 ' 
         IRLab(4) = 'B2 '
         IRLab(5) = 'E1 '
         IRLab(6) = 'E2 '
         IRLab(7) = 'E3 '
       ELSEIF(nGroup.EQ.5) THEN
         IRLab(1) = 'A1g'
         IRLab(2) = 'A2g'
         IRLab(3) = 'E1g'
         IRLab(4) = 'E2g'
         IRLab(5) = 'A1u'
         IRLab(6) = 'A2u' 
         IRLab(7) = 'E1u'
         IRLab(8) = 'E2u'
       ELSEIF(nGroup.EQ.6) THEN
         IRLab(1) = 'A1 '
         IRLab(2) = 'A2 ' 
         IRLab(3) = 'B1 ' 
         IRLab(4) = 'B2 '
         IRLab(5) = 'E1 '
         IRLab(6) = 'E2 '
         IRLab(7) = 'E3 '
         IRLab(8) = 'E4 '
         IRLab(9) = 'E5 '
       ENDIF
c Sn
	ELSEIF(iGroup.EQ.7) THEN
       IF(nGroup.EQ.4) THEN
         IRLab(1) = 'A  '
         IRLab(2) = 'B  '
         IRLab(3) = 'E  '
       ELSEIF(nGroup.EQ.6) THEN
         IRLab(1) = 'Ag '
         IRLab(2) = 'Eg '
         IRLab(3) = 'Au '
         IRLab(4) = 'Eu '
       ELSEIF(nGroup.EQ.8) THEN
         IRLab(1) = 'A  '
         IRLab(2) = 'B  '
         IRLab(3) = 'E1 '
         IRLab(4) = 'E2 '
         IRLab(5) = 'E3 '
       ELSEIF(nGroup.EQ.10) THEN
         IRLab(1) = 'Ag '
         IRLab(2) = 'E1g'
         IRLab(3) = 'E2g'
         IRLab(4) = 'Au '
         IRLab(5) = 'E1u'
         IRLab(6) = 'E2u'
       ELSEIF(nGroup.EQ.12) THEN
         IRLab(1) = 'A  '
         IRLab(2) = 'B  '
         IRLab(3) = 'E1 '
         IRLab(4) = 'E2 '
         IRLab(5) = 'E3 '
         IRLab(6) = 'E4 '
         IRLab(7) = 'E5 '
	 ELSE
	   WRITE(iout,*)'ERROR. Only considering n = 4,6,8,10,12.'
	   STOP
       ENDIF
c Cubic
	ELSEIF(iGroup.EQ.8) THEN
c Td
       IF(nGroup.EQ.3) THEN
         IRLab(1) = 'A1 '
         IRLab(2) = 'A2 '
         IRLab(3) = 'E  '
         IRLab(4) = 'T1 '
         IRLab(5) = 'T2 '
c Oh
       ELSEIF(nGroup.EQ.5) THEN
         IRLab(1)  = 'A1g'
         IRLab(2)  = 'A2g'
         IRLab(3)  = 'Eg '
         IRLab(4)  = 'T1g'
         IRLab(5)  = 'T2g'
         IRLab(6)  = 'A1u'
         IRLab(7)  = 'A2u'
         IRLab(8)  = 'Eu '
         IRLab(9)  = 'T1u'
         IRLab(10) = 'T2u'
	 ELSE 	 
         WRITE(*,*)'ERROR. Only Td and Oh (Higher) groups ready'
	   STOP
       ENDIF
c Linear
	ELSEIF(iGroup.EQ.9) THEN
        WRITE(iout,*)'ERROR. Irred. rep. not ready for Linear groups'
	  STOP

	ENDIF

	RETURN
	END

      SUBROUTINE group_mult(IRMul,nIR,iGroup,nGroup,iout)
c*****************************************************c
c                                                     c
c       Defines Irred. Rep. multiplicities            c
c                                                     c
c       COMMON                                        c
c        iGroup  : group index                        c
c        nGroup  : group order                        c
c        dGroup  : group dimension                    c
c                                                     c
c       INPUT                                         c
c        nIR     : # of irreducible rep.              c
c                                                     c
c       OUTPUT                                        c
c        IRMul   : multiplicity of irreducible rep.   c
c                                                     c
c       DCC: 12/10                                    c
c                                                     c
c*****************************************************c
      IMPLICIT REAL*8  (a-h, o-z)
      INTEGER    iGroup,nGroup,dGroup,IRMul(nIR)

      SAVE Zero,One,Two
      DATA Zero,One,Two / 0d0,1d0,2d0 /

c Only up to 6th order (nGroup < 7)
      IF(iGroup .NE. 7) THEN
       IF(nGroup.LT.2 .OR. nGroup.GT.6) THEN
        WRITE(iout,*)'ERROR. Irreducible rep. only for n < 7'
        STOP
       ENDIF
      ENDIF

c Define IRMul = 1 by default
      IRMul(:) = 1
c Cn,Cnh
	IF(iGroup.EQ.1 .OR. iGroup.EQ.2) THEN
	 IF(nGroup.EQ.3) THEN
         IRMul(2) = 2
         IF(iGroup.EQ.2) IRMul(4) = 2
       ELSEIF(nGroup.EQ.4) THEN
         IRMul(3) = 2
         IF(iGroup.EQ.2) IRMul(6) = 2
       ELSEIF(nGroup.EQ.5) THEN
         IRMul(2:3) = 2
         IF(iGroup.EQ.2) IRMul(5:6) = 2
       ELSEIF(nGroup.EQ.6) THEN
         IRMul(3:4) = 2
         IF(iGroup.EQ.2) IRMul(7:8) = 2
	 ENDIF
c Cnv, Dn
      ELSEIF(iGroup.EQ.3 .OR. iGroup.EQ.4) THEN
	
       IF(nGroup.EQ.3) THEN
         IRMul(3) = 4
       ELSEIF(nGroup.EQ.4) THEN
         IRMul(5) = 4
       ELSEIF(nGroup.EQ.5) THEN
         IRMul(3:4) = 4
       ELSEIF(nGroup.EQ.6) THEN
         IRMul(5:6) = 4
       ENDIF
c Dnh
      ELSEIF(iGroup.EQ.5) THEN
       IF(nGroup.EQ.3) THEN
         IRMul(3) = 4
         IRMul(6) = 4
       ELSEIF(nGroup.EQ.4) THEN
         IRMul(5)  = 4
         IRMul(10) = 4
       ELSEIF(nGroup.EQ.5) THEN
         IRMul(3:4) = 4
         IRMul(7:8) = 4
       ELSEIF(nGroup.EQ.6) THEN
         IRMul(5:6)   = 4
         IRMul(11:12) = 4
       ENDIF
c Dnd
      ELSEIF(iGroup.EQ.6) THEN
       IF(nGroup.EQ.2) THEN
         IRMul(5) = 4
       ELSEIF(nGroup.EQ.3) THEN
         IRMul(3) = 4
         IRMul(6) = 4
       ELSEIF(nGroup.EQ.4) THEN
         IRMul(5:7) = 4
       ELSEIF(nGroup.EQ.5) THEN
         IRMul(3:4) = 4
         IRMul(7:8) = 4
       ELSEIF(nGroup.EQ.6) THEN
         IRMul(5:9) = 4
       ENDIF
c Sn
      ELSEIF(iGroup.EQ.7) THEN
       IF(nGroup.EQ.4) THEN
         IRMul(3) = 2
       ELSEIF(nGroup.EQ.6) THEN
         IRMul(2) = 2
         IRMul(4) = 2
       ELSEIF(nGroup.EQ.8) THEN
         IRMul(3:5) = 2
       ELSEIF(nGroup.EQ.10) THEN
         IRMul(2:3) = 2
         IRMul(5:6) = 2
       ELSEIF(nGroup.EQ.12) THEN
         IRMul(3:7) = 2
       ENDIF

c Cubic
      ELSEIF(iGroup.EQ.8) THEN
c Td
       IF(nGroup.EQ.3) THEN
         IRMul(3)   = 4
         IRMul(4:5) = 9
c Oh
       ELSEIF(nGroup.EQ.5) THEN
         IRMul(3:5)  = (/4,9,9/)
         IRMul(8:10) = (/4,9,9/)
       ELSE
         WRITE(*,*)'ERROR. Only Td and Oh (Higher) groups ready'
         STOP
       ENDIF
c Linear
      ELSEIF(iGroup.EQ.9) THEN
        WRITE(iout,*)'ERROR. Irred. rep. not ready for Linear groups'
        STOP

      ENDIF

      RETURN
      END


      SUBROUTINE table_dimensions(nIR,iGroup,nGroup,dGroup)
c*****************************************************c
c                                                     c
c       Gets the # of Ireducible Representations      c
c                                                     c
c       INPUT                                         c
c        iGroup  : group index                        c
c        nGroup  : group order                        c
c        dGroup  : group dimension                    c
c                                                     c
c       OUTPUT                                        c
c        nIR     : # of irreducible rep.              c
c                                                     c
c       DCC: 11/10                                    c
c                                                     c
c*****************************************************c
      INTEGER    iGroup,nGroup,dGroup

c Define dimensions 

c Nonaxial
      IF(iGroup.EQ.0) THEN
         nIR = 2
c Cn,Cnh,Sn
      ELSEIF(iGroup.EQ.1 .OR. iGroup.EQ.2 .OR.
     &       iGroup.EQ.7) THEN
        nIR = dGroup 
c Cnv,Dn
      ELSEIF(iGroup.EQ.3 .OR. iGroup.EQ.4) THEN
        nIR = (nGroup-1)/2
        nIR = 2*nGroup - 3*nIR
c Dnh
      ELSEIF(iGroup.EQ.5) THEN
        nIR = (nGroup-1)/2
        nIR = 2*(2*nGroup - 3*nIR)
c Dnd
      ELSEIF(iGroup.EQ.6) THEN
        nIR = nGroup + 3
c Cubic
      ELSEIF(iGroup.EQ.8) THEN
        IF(nGroup.EQ.1) THEN
	    nIR = 3
        ELSEIF(nGroup.EQ.2) THEN
	    nIR = 6
c Td
        ELSEIF(nGroup.EQ.3) THEN
	    nIR = 5
        ELSEIF(nGroup.EQ.4 .OR. nGroup.EQ.6) THEN
	    nIR = 5
        ELSEIF(nGroup.EQ.5 .OR. nGroup.EQ.7) THEN
	    nIR = 10
	  ENDIF
c Linear
      ELSEIF(iGroup.EQ.9) THEN
        WRITE(iout,*)'ERROR. Irred. rep. not ready for Linear groups'
        STOP
      ENDIF

	RETURN
	END
