
      SUBROUTINE sym_transform(iAngSym,AlphSym,CSym,RSym,
     & LSym,NSym,MatSym,
     & iAngdens,AlphDens,CDens,RDens,LDens,NDens)
c*****************************************************c
c                                                     c
c       Transforms the electron density               c
c       upon symmetry operation matrix                c
c       from MatSym,                                  c
c                                                     c
c       INPUT                                         c
c        iAngDens: (x,y,z) exponent indices           c
c        AlphDens: function exponentials              c
c        CDens   : function coefficients              c
c        RDens   : function positions                 c
c        LDens   : total # functions of Dens          c
c        NDens   : non-zero Dens functions            c
c        LSym    : total # functions of Sym           c
c        MatSym  : symmetry transformation            c
c                                                     c
c       OUTPUT                                        c
c        iAngSym : (x,y,z) exponent indices           c
c        AlphSym : function exponentials              c
c        CSym    : function coefficients              c
c        RSym    : function positions                 c
c        NSym    : non-zero Sym functions             c
c                                                     c
c       DCC: 05/09                                    c
c                                                     c
c*****************************************************c
      IMPLICIT REAL*8(a-h,o-z)
      INTEGER  iAngDens(LDens,3),iAngSym(LSym,3)
      REAL*8   AlphDens(LDens),CDens(LDens),RDens(LDens,3),
     &         AlphSym(LSym),CSym(LSym),RSym(LSym,3),
     &         MatSym(3,3),VDens(3),VSym(3)
      SAVE Two,Three,Four
      DATA Two,Three,Four / 2.d0,3.d0,4.d0  /

c Common variables
	INCLUDE 'common_var.ins'

c     Tol = 1.0d-8
      Tol = 0.0d0

      kk = 1
      do  k = 1,NDens
       if (ABS(CDens(k)) .GE. Tol) then
        do j = 1,3            
          VDens(j) = RDens(k,j)
        end do

c Old application of symmetry operation matrix
c       IF (DO_REF) THEN
c        do i=1,3
c          VSym(i)=0.0d0
c          do j=1,3
c            VSym(i)=VSym(i)+MatSym(i,j)*VDens(j)
c          end do
c        end do
c       ELSEIF (DO_ROT) THEN
c        do i=1,3
c          VSym(i)=0.0d0
c          do j=1,3
c            VSym(i)=VSym(i)+MatSym(j,i)*VDens(j)
c          end do
c        end do
c       ENDIF

c New application of symmetry operation matrix
c       write(2,'(3(1x,f8.4))') (VDens(j),j=1,3)
        CALL AtimsB(VSym,MatSym,VDens, 3,1,3, 3,3,3, 2)
c       write(2,'(3(1x,f8.4))') (VSym(j),j=1,3)
c       write(2,*)' '

c S type
      IF (iAngDens(k,1)+iAngDens(k,2)+iAngDens(k,3).eq.0)THEN
         CSym(kk) = CDens(k)
         AlphSym(kk) = AlphDens(k)
         do j = 1,3
          RSym(kk,j) = VSym(j)
          iAngSym(kk,j) = iAngDens(k,j)
         end do
         kk = kk+1

c P type
      ELSEIF(iAngDens(k,1)+iAngDens(k,2)+iAngDens(k,3).eq.1)THEN
         do i = 1,3
          do j = 1,3
           iAngSym(kk-1+i,j) = 0
           RSym(kk-1+i,j) = VSym(j)
          end do
          iAngSym(kk-1+i,i) = 1
          AlphSym(kk-1+i) = AlphDens(k)
          CSym(kk-1+i) = CDens(k)*(MatSym(1,i)*
     . dfloat(iAngDens(k,1))+MatSym(2,i)*dfloat(iAngDens(k,2))+
     . MatSym(3,i)*dfloat(iAngDens(k,3)))
c           kk = kk+1
         end do
         kk = kk+3

c D type
      ELSEIF(iAngDens(k,1)+iAngDens(k,2)+iAngDens(k,3).eq.2)THEN
         do i = 1,3
          do j = 1,3
           iAngSym(kk-1+i,j) = 0
           iAngSym(kk+i+2,j) = 1
           RSym(kk-1+i,j) = VSym(j)
           RSym(kk+i+2,j) = VSym(j)
          end do
          iAngSym(kk-1+i,i) = 2
          iAngSym(kk+i+2,i) = 0
          AlphSym(kk-1+i) = AlphDens(k)
          AlphSym(kk+i+2) = AlphDens(k)
         end do
         do i = 1,3
c (2 0 0) type
          if(iAngDens(k,i).eq.2)then
           do j = 1,3
            CSym(kk-1+j) = CDens(k)*(MatSym(i,j))**2
            CSym(kk+j+2) = Two*(MatSym(i,1)**iAngSym(kk+j+2,1))*
     &      (MatSym(i,2)**iAngSym(kk+j+2,2))*
     &      (MatSym(i,3)**iAngSym(kk+j+2,3))*CDens(k)
           end do
           kk = kk+6
          end if
         end do
         do i = 1,2
          do j = i+1,3
c (1 1 0) type
           if((iAngDens(k,i).eq.1).and.(iAngDens(k,j).eq.1))then
            do ll = 1,3
             CSym(kk) = (MatSym(i,ll)*MatSym(j,ll))*CDens(k)
             kk = kk+1
            end do
            CSym(kk) = (MatSym(i,2)*MatSym(j,3)+MatSym(i,3)*
     &      MatSym(j,2))*CDens(k)
            kk = kk+1
            CSym(kk) = (MatSym(i,1)*MatSym(j,3)+MatSym(i,3)*
     &      MatSym(j,1))*CDens(k)
            kk = kk+1
            CSym(kk) = (MatSym(i,2)*MatSym(j,1)+MatSym(i,1)*
     &      MatSym(j,2))*CDens(k)
            kk = kk+1
           end if
          end do
         end do

c F type
      ELSEIF(iAngDens(k,1)+iAngDens(k,2)+iAngDens(k,3).eq.3)THEN
         do i = 1,3
          do j = 1,3
           iAngSym(kk-1+i,j) = 0
           iAngSym(kk+i+3,j) = 1
           iAngSym(kk+i+6,j) = 1
          end do
          iAngSym(kk-1+i,i) = 3
          iAngSym(kk+3,i) = 1
         end do
         do i = 1,3
          iAngSym(kk+2+2*i,i) = 2
          iAngSym(kk+3+2*i,i) = 2
          iAngSym(kk+3+i,4-i) = 0
          iAngSym(kk+6+i,4-i) = 0
         end do
         do i = 1,10
          AlphSym(kk-1+i) = AlphDens(k)
          do j = 1,3
           RSym(kk-1+i,j) = VSym(j)
          end do
         end do
         do i = 1,3
c (3 0 0) type
          if(iAngDens(k,i).eq.3)then
           do j = 1,3
            CSym(kk-1+j) = CDens(k)*(MatSym(i,j))**3
           end do
           CSym(kk+3) = 6.0d0*MatSym(i,1)*MatSym(i,2)*
     &     MatSym(i,3)*CDens(k)
           do j = 1,6
            CSym(kk+j+3) = Three*CDens(k)*MatSym(i,1)**
     &      iAngSym(kk+j+3,1)*MatSym(i,2)**iAngSym(kk+j+3,2)*
     &      MatSym(i,3)**iAngSym(kk+j+3,3)
           end do
           kk = kk+10
          end if 
c (2 1 0) type
          if(iAngDens(k,i).eq.2)then
           do j = 1,3
            if(iAngDens(k,j).eq.1)then
             do ll = 1,3
                CSym(kk-1+ll) = CDens(k)*MatSym(j,ll)*
     &        MatSym(i,ll)**2
             end do
             CSym(kk+3) = CDens(k)*Two*(MatSym(i,1)*
     &       MatSym(i,2)*MatSym(j,3)+MatSym(i,1)*MatSym(j,2)*
     &       MatSym(i,3)+MatSym(j,1)*MatSym(i,2)*MatSym(i,3))
             CSym(kk+4) = CDens(k)*(MatSym(j,2)*MatSym(i,1)**2+
     &       Two*MatSym(i,1)*MatSym(j,1)*MatSym(i,2))
             CSym(kk+5) = CDens(k)*(MatSym(j,3)*MatSym(i,1)**2+
     &       Two*MatSym(i,1)*MatSym(j,1)*MatSym(i,3))
             CSym(kk+6) = CDens(k)*(MatSym(j,3)*MatSym(i,2)**2+
     &       Two*MatSym(i,2)*MatSym(j,2)*MatSym(i,3))
             CSym(kk+7) = CDens(k)*(MatSym(j,1)*MatSym(i,2)**2+
     &       Two*MatSym(i,2)*MatSym(j,2)*MatSym(i,1))
             CSym(kk+8) = CDens(k)*(MatSym(j,1)*MatSym(i,3)**2+
     &       Two*MatSym(i,3)*MatSym(j,3)*MatSym(i,1))
             CSym(kk+9) = CDens(k)*(MatSym(j,2)*MatSym(i,3)**2+
     &       Two*MatSym(i,3)*MatSym(j,3)*MatSym(i,2))
             kk = kk+10
            end if
           end do
          end if
         end do
c (1 1 1) type
         if((iAngDens(k,1).eq.1).and.(iAngDens(k,2).eq.1).and.
     &     (iAngDens(k,3).eq.1))then
          do i = 1,3
           CSym(kk-1+i) = CDens(k)*MatSym(1,i)*MatSym(2,i)*MatSym(3,i)
          end do
          CSym(kk+3) = CDens(k)*(MatSym(1,1)*MatSym(2,2)*MatSym(3,3)+
     &    MatSym(1,1)*MatSym(2,3)*MatSym(3,2)+MatSym(1,2)*MatSym(2,1)*
     &    MatSym(3,3)+MatSym(1,2)*MatSym(2,3)*MatSym(3,1)+MatSym(1,3)*
     &    MatSym(2,1)*MatSym(3,2)+MatSym(1,3)*MatSym(2,2)*MatSym(3,1))
          do l = 5,10
           do i = 1,3
            do j = 1,3
             if((iAngSym(kk-1+l,i).eq.2).and.
     &          (iAngSym(kk-1+l,j).eq.1))then
              CSym(kk-1+l) = CDens(k)*(MatSym(1,i)*MatSym(2,i)*
     &        MatSym(3,j)+MatSym(1,i)*MatSym(2,j)*MatSym(3,i)+
     &        MatSym(1,j)*MatSym(2,i)*MatSym(3,i))
             end if
            end do
           end do
          end do
          kk = kk+10
         end if

c G type      
      ELSEIF(iAngDens(k,1)+iAngDens(k,2)+iAngDens(k,3).eq.4)THEN
         do i = 1,3
          do j = 1,3
           iAngSym(kk-1+i,j) = 0
           iAngSym(kk+i+2,j) = 2
           iAngSym(kk+i+5,j) = 1
           iAngSym(kk+i+8,j) = 1
           iAngSym(kk+i+11,j) = 1
          end do
         end do
         do i = 1,3
          iAngSym(kk-1+i,i) = 4
          iAngSym(kk+i+5,i) = 2
          iAngSym(kk+i+2,4-i) = 0
          iAngSym(kk+7+2*i,i) = 3
          iAngSym(kk+8+2*i,i) = 3
          iAngSym(kk+8+i,4-i) = 0
          iAngSym(kk+11+i,4-i) = 0
         end do
         do i = 1,15
          AlphSym(kk-1+i) = AlphDens(k)
          do j = 1,3
           RSym(kk-1+i,j) = VSym(j)
          end do
         end do
         do i = 1,3
c (4 0 0) type
          if(iAngDens(k,i).eq.4)then
           do j = 1,3
            CSym(kk-1+j) = CDens(k)*(MatSym(i,j))**4
            CSym(kk+2+j) = 6.0d0*CDens(k)*
     &      MatSym(i,1)**iAngSym(kk+2+j,1)*MatSym(i,2)**
     &      iAngSym(kk+2+j,2)*MatSym(i,3)**iAngSym(kk+2+j,3)
            CSym(kk+5+j) = 12.0d0*CDens(k)*
     &      MatSym(i,1)**iAngSym(kk+5+j,1)*MatSym(i,2)**
     &      iAngSym(kk+5+j,2)*MatSym(i,3)**iAngSym(kk+5+j,3)
            CSym(kk+8+j) = Four*CDens(k)*
     &      MatSym(i,1)**iAngSym(kk+8+j,1)*MatSym(i,2)**
     &      iAngSym(kk+8+j,2)*MatSym(i,3)**iAngSym(kk+8+j,3)
            CSym(kk+11+j) = Four*CDens(k)*
     &      MatSym(i,1)**iAngSym(kk+11+j,1)*MatSym(i,2)**
     &      iAngSym(kk+11+j,2)*MatSym(i,3)**iAngSym(kk+11+j,3)
           end do
           kk = kk+15
          end if
c (3 1 0) type
          if(iAngDens(k,i).eq.3)then
           do j = 1,3
            if(iAngDens(k,j).eq.1)then
             do l = 1,3
              CSym(kk-1+l) = CDens(k)*MatSym(j,l)*
     &        MatSym(i,l)**3
             end do
             do l = 1,2
              do ll = l+1,3
               CSym(kk+l+ll) = Three*CDens(k)*
     &         MatSym(i,l)*MatSym(i,ll)*(MatSym(i,l)*
     &         MatSym(j,ll)+MatSym(i,ll)*MatSym(j,l))
               CSym(kk+11-l-ll) = 6.0d0*CDens(k)*
     &         MatSym(i,l)*MatSym(i,ll)*MatSym(i,6-l-ll)*
     &         MatSym(j,6-l-ll)+Three*CDens(k)*
     &         MatSym(i,6-l-ll)**2*(MatSym(i,l)*
     &         MatSym(j,ll)+MatSym(j,l)*MatSym(i,ll))
               CSym(kk+6+l+ll) = CDens(k)*
     &         MatSym(i,l)**2*(MatSym(i,l)*MatSym(j,ll)+
     &         Three*MatSym(j,l)*MatSym(i,ll))
               CSym(kk+9+l+ll) = CDens(k)*
     &         MatSym(i,ll)**2*(MatSym(i,ll)*MatSym(j,l)+
     &         Three*MatSym(j,ll)*MatSym(i,l))
              end do
             end do
             kk = kk+15
            end if
           end do
          end if
         end do
         do i = 1,2
          do j = i+1,3
c (2 2 0) type
           if((iAngDens(k,i).eq.2).and.(iAngDens(k,j).eq.2))then
            do l = 1,3
             CSym(kk-1+l) = CDens(k)*MatSym(i,l)**2*
     &       MatSym(j,l)**2
            end do
            do l = 1,2
             do ll = l+1,3
              CSym(kk+l+ll) = CDens(k)*(MatSym(i,l)**2*
     &        MatSym(j,ll)**2+MatSym(j,l)**2*MatSym(i,ll)**2+
     &        Four*MatSym(i,l)*MatSym(j,l)*MatSym(i,ll)*
     &        MatSym(j,ll))
              CSym(kk+11-l-ll) = Two*CDens(k)*
     &        (MatSym(i,l)*MatSym(i,ll)*MatSym(j,6-l-ll)**2+
     &        MatSym(j,l)*MatSym(j,ll)*MatSym(i,6-l-ll)**2)+
     &        Four*CDens(k)*MatSym(i,6-l-ll)*
     &        MatSym(j,6-l-ll)*(MatSym(i,l)*MatSym(j,ll)+
     &        MatSym(j,l)*MatSym(i,ll))
              CSym(kk+6+l+ll) = CDens(k)*Two*
     &        MatSym(i,l)*MatSym(j,l)*(MatSym(i,l)*
     &        MatSym(j,ll)+MatSym(i,ll)*MatSym(j,l))
              CSym(kk+9+l+ll) = CDens(k)*Two*
     &        MatSym(i,ll)*MatSym(j,ll)*(MatSym(i,l)*
     &        MatSym(j,ll)+MatSym(i,ll)*MatSym(j,l))
             end do
            end do
            kk = kk+15
           end if
c (2 1 1) type
           if((iAngDens(k,i).eq.1).and.(iAngDens(k,j).eq.1))then
            ij = 6-i-j
            do l = 1,3
             CSym(kk-1+l) = CDens(k)*MatSym(i,l)*
     &       MatSym(j,l)*MatSym(ij,l)**2
            end do
            do l = 1,2
             do ll = l+1,3
              CSym(kk+l+ll) = CDens(k)*(MatSym(ij,l)*
     &        MatSym(j,ll)*(MatSym(ij,l)*MatSym(i,ll)+
     &        Two*MatSym(ij,ll)*MatSym(i,l))+MatSym(ij,ll)*
     &        MatSym(j,l)*(MatSym(ij,ll)*MatSym(i,l)+
     &        Two*MatSym(ij,l)*MatSym(i,ll)))
              CSym(kk+11-l-ll) = CDens(k)*(Two*
     &        (MatSym(i,6-l-ll)*MatSym(j,6-l-ll)*MatSym(ij,l)*
     &        MatSym(ij,ll)+MatSym(ij,6-l-ll)*MatSym(i,6-l-ll)*
     &        (MatSym(ij,l)*MatSym(j,ll)+MatSym(j,l)*
     &        MatSym(ij,ll))+MatSym(ij,6-l-ll)*MatSym(j,6-l-ll)*
     &        (MatSym(ij,l)*MatSym(i,ll)+MatSym(i,l)*
     &        MatSym(ij,ll)))+(MatSym(i,l)*MatSym(j,ll)+
     &        MatSym(j,l)*MatSym(i,ll))*MatSym(ij,6-l-ll)**2)
              CSym(kk+6+l+ll) = CDens(k)*((MatSym(i,l)*
     &        MatSym(j,ll)+MatSym(i,ll)*MatSym(j,l))*
     &        MatSym(ij,l)**2+Two*MatSym(ij,l)*MatSym(i,l)*
     &        MatSym(j,l)*MatSym(ij,ll))
              CSym(kk+9+l+ll) = CDens(k)*((MatSym(i,l)*
     &        MatSym(j,ll)+MatSym(i,ll)*MatSym(j,l))*
     &        MatSym(ij,ll)**2+Two*MatSym(ij,ll)*MatSym(i,ll)*
     &        MatSym(j,ll)*MatSym(ij,l))
             end do
            end do
            kk = kk+15
           end if
          end do
         end do
        end if

       ENDIF
      end do
      NSym = kk-1


      RETURN
      END

