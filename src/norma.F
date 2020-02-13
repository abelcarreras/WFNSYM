
        SUBROUTINE norma(i_angl,a,b)

c*****************************************************c
c                                                     c
c       Computes the normalization factor             c
c       of gaussian funcions.                         c
c                                                     c
c       INPUT                                         c
c               i_angl: cartesian indices             c
c               a: gaussian exponent                  c
c                                                     c
c       OUTPUT                                        c
c               b: normalization constant             c
c                                                     c
c       DCC: 05/09                                    c
c                                                     c
c*****************************************************c
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER  ::    i_angl(3),i_total
      REAL *8  ::    a,b,pi

      pi=4.0d0*atan(1.0d0)
      i_total=i_angl(1)+i_angl(2)+i_angl(3)
      if(i_total.eq.0) then
        b=(2.d0*a/pi)**(0.75d0)
      end if
      if(i_total.eq.1) then
        b=(2.d0*a/pi)**(0.75d0)*2.d0*dsqrt(a)
      end if
      if(i_total.eq.2)then
        b=4.d0*a*(2.d0*a/pi)**(0.75d0)
        do i=1,3
          if(i_angl(i).eq.2) then
            b=b*0.5773502692
          end if
        end do
!       if(i_angl(1).eq.2) then
!         b=b*0.8660254038
!       end if
!       if(i_angl(3).eq.2) then
!         b=b*0.5d0
!       end if
      end if

      RETURN
      END

