      SUBROUTINE get_mos(Ca,Cb,NDim,IOff,input,ISpin)
c
      INTEGER   :: IOff, reclength
      REAL*8    :: Ca(NDim),Cb(NDim)
      CHARACTER :: input*80

c Open file for binary IO (input/output)
      inquire(iolength=reclength) Ca,Cb

c     reclength = reclength*8
c double == 8bytes,
c reclength now equals #bytes( value & grbge )

      IOff = IOff + 1
      OPEN(Unit=15,File=trim(input)//".in.save/53.0",
     &     access='direct',action='read',
     &     recl=reclength,form='UNFORMATTED',status="unknown")

c Read variables from binary file
      READ(15,rec=IOff) Ca,Cb
      CLOSE(15)

      RETURN
      END



