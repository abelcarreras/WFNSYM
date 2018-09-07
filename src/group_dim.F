c Define Group Dimensions (dGroup and hGroup)

      SUBROUTINE SetGroupDim(iGroup, nGroup, dGroup, hGroup)
      INTEGER       :: iGroup, hGroup,dGroup,nGroup

	   IF(iGroup.EQ.0) THEN
c C1
	     IF(nGroup.EQ.1) THEN
	       dGroup = 1
c Ci,Cs
	     ELSE
	       dGroup = 2
	     ENDIF
	     hGroup = dGroup
c Cn
	   ELSEIF(iGroup.EQ.1) THEN
	     dGroup = 1 + nGroup/2
 	     hGroup = nGroup
c Cnh
	   ELSEIF(iGroup.EQ.2) THEN
	     dGroup = 2*(1 + nGroup/2)
	     hGroup = 2*nGroup
c Cnv,Dn
	   ELSEIF(iGroup.EQ.3 .OR. iGroup.EQ.4) THEN
	     dGroup = nGroup + nGroup/2 + 1
	     hGroup = 2*nGroup
c Dnh
	   ELSEIF(iGroup.EQ.5) THEN
	     dGroup = 2 - (nGroup-2*(nGroup/2))
	     dGroup = 3*nGroup + dGroup
	     hGroup = 4*nGroup
c Dnd
 	   ELSEIF(iGroup.EQ.6) THEN
	     dGroup = 3 * nGroup + 1
	     hGroup = 4*nGroup
c Sn (n = 4,6,8,10,12)
         ELSEIF(iGroup.EQ.7) THEN
           dGroup = 1 + nGroup/2
	     hGroup = nGroup
c Cubic
	   ELSEIF(iGroup.EQ.8) THEN
c Td
	     IF(nGroup.EQ.3) THEN
 	       dGroup = 17
 	       hGroup = 24
c Oh
	     ELSEIF(nGroup.EQ.5) THEN
 	       dGroup = 34
 	       hGroup = 48
	     ELSE
 	       WRITE(iout,*)'ERROR. Only Td and Oh among cubic groups'
 	       STOP
	     ENDIF
	   ENDIF

        END