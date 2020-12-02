      SUBROUTINE JADAMILU_INTERFACE (Nt,A,JA,IA)
*     simple program to illustrate basic usage of PJD
      IMPLICIT NONE
	  integer          Nt             
*     ..                                                                        
*     .. Array Arguments ..                                                     
      integer            JA(*), IA(*)                                           
      DOUBLE PRECISION   A(*)                                              
*     needed parameters: 
*              N: size of the problem
*         MAXEIG: max. number of wanteg eig (NEIG<=MAXEIG)
*          MAXSP: max. value of MADSPACE
	  integer N
	  PARAMETER (N = Nt)  
      INTEGER  MAXEIG, MAXSP
	  PARAMETER (MAXEIG=3,MAXSP=20)
*     optimal workspace (here MAXSP*MAXSP>MAXEIG)
      INTEGER LX
      PARAMETER(LX=N*(3*MAXSP+MAXEIG+1)+4*MAXSP*MAXSP)
      DOUBLE PRECISION   EIGS(MAXEIG), RES(MAXEIG), X(LX)
*     arguments to pass to the routines
      INTEGER            NEIG, MADSPACE, ISEARCH, NINIT, ICNTL(5)
      INTEGER            ITER, IPRINT, INFO
      DOUBLE PRECISION   SIGMA, TOL, SHIFT, GAP, MEM, DROPTOL
*
*     some local variables
      INTEGER I,K
*     let us define a very simple matrix:
*
*              0     5      0      ...
*              5     1      5       0     ...
*         A =  0     5      2       5      0     ...
*              0     0      5       3      5      0     ...
*             ...   ...   
*
*     we use the storage scheme needed by PJD
*     (only the upper triangular part is referenced)
*      DOUBLE PRECISION A(9)
*      INTEGER IA(6), JA(9)

*	  IA = (/1, 3,5,7,9,10 /)
*	  JA = (/1,2,2,3,3,4,4,5,5 /)
*      A = (/0.0, 5.0, 1.0, 5.0, 2.0, 5.0, 3.0, 5.0, 4.0 /)






*      K=1
*      DO I=1,N
*         IA(I)=K
*         JA(K)=I
*         A(K)=DBLE(I-1)
*         K=K+1
*         IF (I .LT. N) THEN
*            JA(K)=I+1
*            A(K)=5.0D0
*            K=K+1
*         END IF
*      END DO
*      IA(N+1)=K
*
*      PRINT *, 
*     $'================================================================'
*      PRINT *, 
*     $'Computing the smallest eigenvalues with ILU preconditioning'
*      PRINT *, 
*     $'-----------------------------------------------------------'
*
c     set input variables
c     the matrix is already in the required format
c   
c     standard report on standard output
      IPRINT=6
c
c     ISEARCH=0: compute the smallest eigenvalues and no a priori
c                information on the eigenvalues is available
      ISEARCH=0
c     (ISEARCH.eq.0: SIGMA, SHIFT need not to be set)
c     
c     elbow space factor for the fill computed during the ILU
      MEM=20.0
c     tolerence for discarded fill
      DROPTOL=1.d-3
c
c     number of wanted eigenvalues
      NEIG=maxeig
c     no initial approximate eigenvectors
      NINIT=0
c     desired size of the search space
      MADSPACE=maxsp
c     maximum number of iteration steps
      ITER=1000
c     tolerance for the eigenvector residual
      TOL=1.0d-10
c     additional parameters set to default
      ICNTL(1)=0
      ICNTL(2)=0
      ICNTL(3)=0
      ICNTL(4)=0
      ICNTL(5)=0
c     
c     We have the exact matrix: we may call DPJD (double precision version)
c
       CALL DPJD( N,a,ja,ia, EIGS, RES, X, LX, NEIG, 
     $           SIGMA, ISEARCH, NINIT, MADSPACE, ITER, TOL,
     $           SHIFT, DROPTOL, MEM, ICNTL,
     $           IPRINT, INFO, GAP)
c
c
c     release the memory allocated for preconditioning
c     (not needed if the computation really terminates here)
      CALL DPJDCLEANUP
c
      END
