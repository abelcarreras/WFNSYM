
      SUBROUTINE AtimsB(C,A,B,M,N,K,LDC,LDA,LDB,IType)
c     ***************************************************************
c     *                                                             *
c     *  AtimsB constructs the product of two general real          *
c     *         matrices.                                           *
c     *                                                             *
c     *  OUTPUT:                                                    *
c     *     C     - Result of matrix multiply                       *
c     *                                                             *
c     *  INPUT:                                                     *
c     *     A,B   - Matrices to be multiplied                       *
c     *     M     - Number of rows in C and op(A)                   *
c     *     N     - Number of cols in C and op(B)                   *
c     *     K     - Number of cols in op(A) and rows in op(B)       *
c     *               = dimension of "contraction"                  *
c     *     LDC   - Leading dimension of C                          *
c     *     LDA   - Leading dimension of A                          *
c     *     LDB   - Leading dimension of B                          *
c     *     IType - Type of product to construct                    *
c     *                 1 - C = A * B                               *
c     *                 2 - C = A(T) * B                            *
c     *                 3 - C = A * B(T)                            *
c     *                 4 - C = A(T) * B(T)                         *
c     *               +10 - Increment C                             *
c     *              Sign - If < 0, do negative of above            *
c     *                                                             *
c     *     Above, op(X) = X or X(T) as appropriate                 *
c     *                                                             *
c     *  BGJ  (7/93) - Based on routines by CAG and MHG             *
c     *  BGJ (10/94) - BLAS call added                              *
c     *                                                             *
c     ***************************************************************
      IMPLICIT REAL*8(a-h,o-z)
      REAL*8 C(LDC,N),A(LDA,*),B(LDB,*)
      CHARACTER TRANSA,TRANSB

c     C := alpha * op( A ) * op( B ) + beta * C

      IF (IType.gt.0) THEN
         Alpha = 1D0
      ELSE
         Alpha =-1D0
      ENDIF
      JType = ABS(IType)
      IF (JType.gt.10) THEN
         Beta = 1D0
      ELSE
         Beta = 0D0
         CALL VRload(C,LDC*N,0d0)
      ENDIF
      TRANSA = 'N'
      TRANSB = 'N'
      JType = MOD(JType,10)
      IF (JType.eq.2.or.JType.eq.4) TRANSA = 'T'
      IF (JType.eq.3.or.JType.eq.4) TRANSB = 'T'
      CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      RETURN
      END

      SUBROUTINE VRload(A,N,Value)
c
c     Set each element of A to a scalar
c
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 A(N)
      DO i = 1, N-3, 4
         A(i) = Value
         A(i+1) = Value
         A(i+2) = Value
         A(i+3) = Value
      ENDDO
      DO i = i, N
         A(i) = Value
      ENDDO
      RETURN
      END

      SUBROUTINE VRdot(Dot,A,B,N)
c
c     Compute the dot product of A and B
c
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 A(N),B(N)

      Dot = ddot(n,a,1,b,1)

      RETURN
      END

      SUBROUTINE VRcross(C,A,B)
c
c     C = A x B
c
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 C(3),A(3),B(3)
      C(1) = A(2)*B(3) - A(3)*B(2)
      C(2) = A(3)*B(1) - A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)
      RETURN
      END

      SUBROUTINE VRcopy(A,B,N)
c
c     Copy the N first elements of B to A
c
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 A(N),B(N)
      DO i = 1, N-3, 4
         A(i)   = B(i)
         A(i+1) = B(i+1)
         A(i+2) = B(i+2)
         A(i+3) = B(i+3)
      ENDDO
      DO i = i, N
         A(i) = B(i)
      ENDDO
      RETURN
      END

*Deck VRadd
      SUBROUTINE VRadd(C,A,B,N)
c
c     C = A + B
c
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 C(N),A(N),B(N)
      DO i = 1, N-3, 4
         C(i) = A(i) + B(i)
         C(i+1) =  A(i+1) + B(i+1)
         C(i+2) =  A(i+2) + B(i+2)
         C(i+3) =  A(i+3) + B(i+3)
      ENDDO
      DO i = i, N
         C(i) = A(i) + B(i)
      ENDDO
      RETURN
      END

*Deck VRadd2
      SUBROUTINE VRadd2(A, B, N)
c     
c     A = A + B
c     
      implicit none
      INTEGER N
      REAL*8 A(N), B(N)
      INTEGER i

      DO i = 1, N
         A(i) = A(i) + B(i)
      ENDDO

      RETURN
      END

*Deck VRneg
      SUBROUTINE VRneg(A,N)
c
c     A = -A
c
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 A(N)
      DO i = 1, N-3, 4
         A(i) = -A(i)
         A(i+1) = -A(i+1)
         A(i+2) = -A(i+2)
         A(i+3) = -A(i+3)
      ENDDO
      DO i = i, N
         A(i) = -A(i)
      ENDDO
      RETURN
      END

*Deck VRscale
      SUBROUTINE VRscale(A,N,Scale)
c
c     A = scalar * A
c
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 A(N)
      DO i = 1, N-3, 4
         A(i) = Scale * A(i)
         A(i+1) = Scale * A(i+1)
         A(i+2) = Scale * A(i+2)
         A(i+3) = Scale * A(i+3)
      ENDDO
      DO i = i, N
         A(i) = Scale * A(i)
      ENDDO
      RETURN
      END

*Deck VRsub
      SUBROUTINE VRsub(C,A,B,N)
c
c     C = A - B
c
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 C(N),A(N),B(N)
      DO i = 1, N-3, 4
         C(i) = A(i) - B(i)
         C(i+1) = A(i+1) - B(i+1)
         C(i+2) = A(i+2) - B(i+2)
         C(i+3) = A(i+3) - B(i+3)
      ENDDO
      DO i = i, N
         C(i) = A(i) - B(i)
      ENDDO
      RETURN
      END

*Deck VRsub2.  Two arguments for accumulative sub.
      SUBROUTINE VRsub2(A,B,N)
c
c     A = A - B
c
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 A(N),B(N)
      DO i = 1, N-3, 4
         A(i) = A(i) - B(i)
         A(i+1) = A(i+1) - B(i+1)
         A(i+2) = A(i+2) - B(i+2)
         A(i+3) = A(i+3) - B(i+3)
      ENDDO
      DO i = i, N
         A(i) = A(i) - B(i)
      ENDDO
      RETURN
      END

*Deck VRtrace
      SUBROUTINE VRtrace(Trace,A,N)
c
c     Sum the elements of A
c
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 A(N)
      Trace = 0D0
      Do i = 1,N-3,4
         Trace = Trace + A(i)
         Trace = Trace + A(i+1)
         Trace = Trace + A(i+2)
         Trace = Trace + A(i+3)
      ENDDO
      DO i = i,N
         Trace = Trace + A(i)
      ENDDO
      RETURN
      END

*Deck VRnorm
      SUBROUTINE VRnorm(ENorm,A,N)
C
C     Compute the sqrt of the sum of squares
C     For an n x n matrix A (where N = n*n) this gives
C     the Euclidean norm.
C
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 A(N)
      ENorm=0.d0
      DO i=1,N-3,4
         ENorm = ENorm + A(i)*A(i)
         ENorm = ENorm + A(i+1)*A(i+1)
         ENorm = ENorm + A(i+2)*A(i+2)
         ENorm = ENorm + A(i+3)*A(i+3)
      ENDDO
      DO i=i,N
         ENorm = ENorm + A(i)*A(i)
      ENDDO
      RETURN
      END

*Deck
      subroutine writeMat(A,N,iout)
      implicit none
      integer N,i,j,iout
      real*8 A(N,N)

      do i = 1, N
         write(iout, 1000) (A(i,j), j = 1, N)
      enddo
 1000 format(12f9.4)

      return
      end

*Deck TrspSq
      SUBROUTINE TrspSq(A,N)
c     ******************************************************************
c     *                                                                *
c     *           TrspSq transposes an N x N array called A.           *
c     *                                                                *
c     *           PMWG (2/93)                                          *
c     *                                                                *
c     ******************************************************************
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 A(N,N)
      DO 20 I = 2,N
         DO 10 J = 1,I-1
            Temp   = A(J,I)
            A(J,I) = A(I,J)
            A(I,J) = Temp
10       CONTINUE
20    CONTINUE
      RETURN
      END

