
        SUBROUTINE overlap(i_angl,a,R,b)
c*****************************************************c
c                                                     c
c       Computes the space integral of gaussian       c
c       funcions.                                     c
c                                                     c
c       INPUT                                         c
c               i_angl: cartesian indices             c
c               a: gaussian exponent                  c
c               R: position of the function           c
c                                                     c
c       OUTPUT                                        c
c               b: integral value (scalar)            c
c                                                     c
c       DCC: 05/09                                    c
c                                                     c
c*****************************************************c
	IMPLICIT REAL*8  (a-h, o-z)
      INTEGER  ::  i_angl(3),i_total
	REAL*8   ::  R(3)
	
      pi = 4.0d0 * DATAN(1.0d0)
      i_total=i_angl(1)+i_angl(2)+i_angl(3)
   
c S function
      IF (i_total.eq.0) THEN
        b=(pi/a)**(1.5d0)

c P function
      ELSEIF (i_total.eq.1) THEN
        b=(dfloat(i_angl(1))*R(1)+dfloat(i_angl(2))*R(2)+
     &    dfloat(i_angl(3))*R(3))*(pi/a)**(1.5d0)

c D function
      ELSEIF (i_total.eq.2) THEN
        if((i_angl(1).eq.2).or.(i_angl(2).
     &    eq.2).or.(i_angl(3).eq.2))then
          b=0.5d0*(dfloat(i_angl(1))*R(1)**2+dfloat(i_angl(2))*
     &    R(2)**2+dfloat(i_angl(3))*R(3)**2)*(pi/a)**(1.5d0)
          b=b+0.5d0*(pi**(1.5d0))/(a**(2.5d0))
        else
          b=0.0d0
          do i=1,2
          do j=i+1,3
            b=b+(dfloat(i_angl(i))*R(i)*
     &           dfloat(i_angl(j))*R(j))
          end do
          end do
          b=b*(pi/a)**(1.5d0)
        end if

c F function
      ELSEIF (i_total.eq.3) THEN
          b=(R(1)**i_angl(1))*(R(2)**i_angl(2))*(R(3)**i_angl(3))
     .      *(pi/a)**(1.5d0)
c       (1,1,1)
          if((i_angl(1).eq.1).and.(i_angl(2).eq.1).
     .       and.(i_angl(3).eq.1)) then
            b=b
          else
            do i=1,3
              if(i_angl(i).eq.1.or.i_angl(i).eq.3)then
                b=b+0.5d0*dfloat(i_angl(i))*R(i)*(pi**(1.5d0))/
     &            (a**(2.5d0))
              end if
            end do
          end if

c G function
      ELSEIF (i_total.eq.4) THEN
          b=(R(1)**i_angl(1))*(R(2)**i_angl(2))*(R(3)**i_angl(3))
     .      *(pi/a)**(1.5d0)
c       (4,0,0)
          if((i_angl(1).eq.4).or.(i_angl(2).eq.4).or.
     .       (i_angl(3).eq.4))then
            b=b+(0.75d0)*(pi**(1.5d0))/(a**(3.5d0))+
     .        3.d0*(R(1)**(i_angl(1)/2))*(R(2)**(i_angl(2)/2))*
     .        (R(3)**(i_angl(3)/2))*(pi**(1.5d0))/(a**(2.5d0))
          end if
c       (2,1,1) 
          if((i_angl(1).ne.0).and.(i_angl(2).ne.0).and.
     .       (i_angl(3).ne.0.))then
            do i=1,2
              if(i_angl(i).eq.1)then
                do j=i+1,3
                 if(i_angl(j).eq.1)then
                  b=b+0.5d0*R(i)*R(j)*(pi**(1.5d0))/(a**(2.5d0))
                 end if
                end do
              end if
            end do
          end if
c       (3,1,0)
          do i=1,3
            if(i_angl(i).eq.3)then
              do j=1,3
                if(i_angl(j).eq.1)then
                  b=b+(1.5d0)*R(i)*R(j)*(pi**(1.5d0))/(a**(2.5d0))
                end if
              end do
            end if
          end do
c       (2,2,0) 
          do i=1,2
            if(i_angl(i).eq.2)then
              do j=i+1,3
                if(i_angl(j).eq.2)then
                  b=b+(0.5d0)*(R(i)**2+R(j)**2)*
     .             (pi**(1.5d0))/(a**(2.5d0))+(0.25d0)*
     .             (pi**(1.5d0))/(a**(3.5d0))
                end if
              end do
            end if
          end do

c H function
      ELSEIF (i_total.eq.5) THEN
         b=(1.d0/a)**(3.5d0)
c       (5,0,0)
         do i=1,3
          if(i_angl(i).eq.5)then
           b=b*((3.75d0)*R(i)+5.d0*a*R(i)**3)
          end if
c       (4,1,0)
          if(i_angl(i).eq.4)then
           do j=1,3
            if(i_angl(j).eq.1)then
             b=b*((0.75d0)*R(j)+3.d0*a*R(j)*R(i)**2)
            end if
           end do
          end if
c       (3,2,0)
          if(i_angl(i).eq.3)then
           do j=1,3
            if(i_angl(j).eq.2)then
             b=b*((0.75d0)*R(i)+(3.d0*R(i)*R(j)**2+R(i)**3)*a*0.5d0)
            end if
           end do
          end if
         end do
         do i=1,2
          do j=i+1,3
c       (3,1,1) (k,i,j)
           if((i_angl(i).eq.1).and.(i_angl(j).eq.1))then
            b=b*a*1.5d0*R(1)*R(2)*R(3)
           end if
c       (2,2,1) (i,j,k)
           if((i_angl(i).eq.2).and.(i_angl(j).eq.2))then
            k=6-(i+j)
            b=b*R(k)*(0.25d0+0.5d0*a*(R(i)**2+R(j)**2))
           end if
          end do
         end do
         b=b+(R(1)**i_angl(1))*(R(2)**i_angl(2))*(R(3)**i_angl(3))/
     .    (a**(1.5d0))
         b=b*pi**(1.5d0)

c I function
      ELSEIF (i_total.eq.6) THEN
         b=(1.d0/a)**(4.5d0)
c       (6,0,0)
         do i=1,3
          if(i_angl(i).eq.6)then
           b=b*7.5d0*(0.25d0+1.5d0*a*R(i)**2+(a*R(i)**2)**2)
          end if
c       (5,1,0) 
          if(i_angl(i).eq.5)then
           do j=1,3
            if(i_angl(j).eq.1)then
             b=b*5.d0*R(j)*(0.75d0*a*R(i)+(a**2)*R(i)**3)
            end if
           end do
          end if
c       (4,2,0) 
          if(i_angl(i).eq.4)then
           do j=1,3
            if(i_angl(j).eq.2)then
             b=b*((0.375d0)+1.5d0*a*(R(i)**2+0.5d0*R(j)**2)+
     .        (0.5d0*R(i)**4+3.d0*(R(i)*R(j))**2)*a**2)
            end if
           end do
          end if
c       (3,2,1) (i,j,k)
          if(i_angl(i).eq.3)then
           do j=1,3
            if(i_angl(j).eq.2)then
             k=6-(i+j)
             b=b*R(k)*(0.75d0*R(i)*a+
     .        0.5d0*(R(i)**3+3.d0*R(i)*R(j)**2)*a**2)
            end if
           end do
          end if
         end do
         do i=1,2
          do j=i+1,3
c       (3,3,0)
           if((i_angl(i).eq.3).and.(i_angl(j).eq.3))then
            b=b*1.5d0*(1.5d0*R(i)*R(j)*a+(R(i)*R(j)**3+
     &        R(j)*R(i)**3)*a**2)
           end if
c       (4,1,1) (k,i,j)
           if((i_angl(i).eq.1).and.(i_angl(j).eq.1))then
            k=6-(i+j)
            b=b*R(i)*R(j)*(0.75d0*a+3.d0*(a*R(k))**2)
           end if
          end do
         end do
c       (2,2,2)
         if((i_angl(1).eq.2).and.(i_angl(2).eq.2).
     &    and.(i_angl(3).eq.2))then
          b=b*(0.125d0+0.25d0*(R(1)**2+R(2)**2+R(3)**2)*a+
     &    0.5d0*((R(1)*R(2))**2+(R(1)*R(3))**2+(R(2)*R(3))**2)*a**2)
         end if
         b=b+(R(1)**i_angl(1))*(R(2)**i_angl(2))*(R(3)**i_angl(3))/
     &    (a**(1.5d0))
         b=b*pi**(1.5d0)

c J function
      ELSEIF (i_total.eq.7) THEN
         b=(1.d0/a)**(4.5d0)
         do i=1,3
c       (7,0,0)
          if(i_angl(i).eq.7)then
           b=b*10.5d0*(1.25d0*R(i)+2.5d0*a*R(i)**3+(a**2)*R(i)**5)
          end if
c       (6,1,0)
          if(i_angl(i).eq.6)then
           do j=1,3
            if(i_angl(j).eq.1)then
             b=b*7.5d0*R(j)*(0.25d0+1.5d0*a*R(i)**2+
     &        (a**2)*R(i)**4)
            end if
           end do
          end if
c       (5,2,0) 
          if(i_angl(i).eq.5)then
           do j=1,3
            if(i_angl(j).eq.2)then
         b=b*((1.875d0)*R(i)+a*(2.5d0*R(i)**3+3.75d0*R(i)*R(j)**2)+
     &        (0.5d0*R(i)**5+5.d0*(R(i)**3)*(R(j)**2))*a**2)
            end if
           end do
          end if
c       (4,3,0) i (4,2,1)
          if(i_angl(i).eq.4)then
           do j=1,3
c       (4,3,0)
            if(i_angl(j).eq.3)then
       b=b*1.5d0*(0.75d0*R(j)+a*(3.d0*R(j)*R(i)**2+0.5d0*R(j)**3)+
     &        (R(j)*R(i)**4+2.d0*(R(i)**2)*(R(j)**3))*(a**2))
            end if
c       (4,2,1) (i,j,k)
            if(i_angl(j).eq.2)then
             k=6-(i+j)
             b=b*1.5d0*R(k)*(0.25d0+a*(R(i)**2+0.5d0*R(j)**2)+
     &        (a**2)*((1.d0/3.d0)*R(i)**4+2.d0*(R(i)*R(j))**2))
            end if
           end do
          end if
         end do
         do i=1,2
          do j=i+1,3
c       (3,3,1) (i,j,k)
           if((i_angl(i).eq.3).and.(i_angl(j).eq.3))then
            k=6-(i+j)
            b=b*1.5d0*R(i)*R(j)*R(k)*a*(1.5d0+a*(R(i)**2+R(j)**2))
           end if
c       (3,2,2) (k,i,j)
           if((i_angl(i).eq.2).and.(i_angl(j).eq.2))then
            k=6-(i+j)
            b=b*R(k)*((0.375d0)+0.25d0*a*(R(k)**2+
     &        3.d0*(R(i)**2+R(j)**2))+0.5d0*((R(i)**2+
     &        R(j)**2)*R(k)**2+3.d0*(R(i)*R(j))**2)*a**2)
           end if
c       (5,1,1) (k,i,j)
           if((i_angl(i).eq.1).and.(i_angl(j).eq.1))then
            k=6-(i+j)
            b=b*5.d0*R(i)*R(j)*R(k)*a*(0.75d0+a*R(k)**2)
           end if
          end do
         end do
         b=b+(R(1)**i_angl(1))*(R(2)**i_angl(2))*(R(3)**i_angl(3))/
     .    (a**(1.5d0))
         b=b*pi**(1.5d0)

c K function
      ELSEIF (i_total.eq.8) THEN
         b=(1.d0/a)**(5.5d0)
         do i=1,3
c       (8,0,0)
          if(i_angl(i).eq.8)then
           b=b*(52.5d0)*(0.125d0+a*R(i)**2+(a**2)*R(i)**4+
     &      (28.d0/105.d0)*(a**3)*R(i)**6)
          end if
c       (7,1,0)
          if(i_angl(i).eq.7)then
           do j=1,3
            if(i_angl(j).eq.1)then
             b=b*a*R(i)*R(j)*10.5d0*(1.25d0+2.5d0*a*R(i)**2+
     &        (a*R(i)**2)**2)
            end if
           end do
          end if
c       (6,2,0)
          if(i_angl(i).eq.6)then
           do j=1,3
            if(i_angl(j).eq.2)then
         b=b*7.5d0*(0.125d0+0.25d0*a*(3.d0*R(i)**2+R(j)**2)+0.5d0*
     &        ((a*R(i))**2)*(3.d0*R(j)**2+R(i)**2)+((1.d0/15.d0)*
     &        R(i)**2+R(j)**2)*(a**3)*R(i)**4)
            end if
           end do
          end if
c       (5,3,0) i (5,2,1)
          if(i_angl(i).eq.5)then
           do j=1,3
c       (5,3,0)
            if(i_angl(j).eq.3)then
       b=b*a*R(i)*R(j)*(5.625d0+3.75d0*(2.d0*R(i)**2+R(j)**2)*a+
     &        (1.5d0*R(i)**4+5.d0*(R(i)*R(j))**2)*a**2)
            end if
c       (5,2,1) (i,j,k)
            if(i_angl(j).eq.2)then
             k=6-(i+j)
       b=b*R(i)*R(k)*a*(1.875d0+a*(2.5d0*R(i)**2+3.75d0*R(j)**2)+
     &        (0.5d0*R(i)**4+5.d0*(R(i)*R(j))**2)*a**2)
            end if
           end do
          end if
c       (4,3,1)
          if(i_angl(i).eq.4)then
           do j=1,3
            if(i_angl(j).eq.3)then
             k=6-(i+j)
       b=b*R(j)*R(k)*a*1.5d0*(0.75d0+(3.d0*R(i)**2+0.5d0*R(j)**2)*a+
     &        (R(i)**4+2.d0*(R(i)*R(j))**2)*a**2)
            end if
           end do 
          end if
         end do 
         do i=1,2
          do j=i+1,3
c       (4,4,0)
           if((i_angl(i).eq.4).and.(i_angl(j).eq.4))then
      b=b*0.75d0*(0.75d0+3.d0*a*(R(i)**2+R(j)**2)+(R(i)**4+R(j)**4+
     &       12.d0*(R(i)*R(j))**2)*a**2+4.d0*((R(i)*R(j))**2)*
     &       (R(i)**2+R(j)**2)*a**3)
           end if
c       (3,3,2) (i,j,k)
           if((i_angl(i).eq.3).and.(i_angl(j).eq.3))then
            k=6-(i+j)
            b=b*R(i)*R(j)*a*(1.125d0+0.75d0*a*(R(i)**2+
     &       R(j)**2+3.d0*R(k)**2)+0.5d0*((R(i)*R(j))**2+
     &       3.d0*(R(i)*R(k))**2+3.d0*(R(j)*R(k))**2)*a**2)
           end if
c       (4,2,2) (k,i,j)
           if((i_angl(i).eq.2).and.(i_angl(j).eq.2))then
            k=6-(i+j)
        b=b*0.25d0*(0.75d0+1.5d0*a*(2.d0*R(k)**2+R(i)**2+R(j)**2)+
     &       (R(k)**4+3.d0*(R(i)*R(j))**2+6.d0*(R(k)*R(i))**2+
     &       6.d0*(R(k)*R(j))**2)*a**2+2.d0*(R(k)**2)*(a**3)*
     &       ((R(i)*R(k))**2+(R(j)*R(k))**2+6.d0*(R(i)*R(j))**2))
           end if
c       (6,1,1) (k,i,j)
           if((i_angl(i).eq.1).and.(i_angl(j).eq.1))then
            k=6-(i+j)
            b=b*a*R(i)*R(j)*7.5d0*(0.25d0+1.5d0*a*R(k)**2+
     &        (R(k)**4)*a**2)
           end if
          end do
         end do
       b=b+(R(1)**i_angl(1))*(R(2)**i_angl(2))*(R(3)**i_angl(3))/
     &    (a**(1.5d0))
         b=b*pi**(1.5d0)
        
      ENDIF

      RETURN
      END

