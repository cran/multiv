C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                 C
C  Optimise the variances of a set of groups, by assigning        C
C  the objects in groups such that they are minimally distant     C
C  from group centres.                                            C
C                                                                 C
C  Parameters:                                                    C
C                                                                 C
C  N, M, NG         Numbers of rows, columns, groups,             C
C  A(N,M)           initial data,                                 C
C  MEMGP(N)         group memberships,                            C
C  NGP0             minimum acceptable group cardinality,         C
C  NUMGP(NG)        cardinalities of groups,                      C
C  GPCEN(NG,M)      group centres,                                C
C  COMP(NG)         compactness values for the groups,            C
C  CTOT             sum of these compactnesses,                   C
C  ITMAX            max allowed no. of iterations,                C
C  IERR             error indicator (should be zero).             C
C                                                                 C
C  IERR = 1: invalid group number (<1 or >NG), - is number of     C
C  groups correctly specified?  IERR = 2: a group has < minimum   C
C  allowed number of members, - reduce the number of groups and   C
C  try again.                                                     C
C                                                                 C
C  History                                                        C
C  Author: F. Murtagh, Garching/Munich, Feb. 1986.                C
C  Added parameter ITMAX (formerly hard-wired).  [F.M., Sept. 91] C 
C                                                                 C
C-----------------------------------------------------------------C
        SUBROUTINE MINDST(A,N,M,MEMGP,NGP0,NUMGP,GPCEN,NG,
     X            COMP,CTOT,ITER,ITMAX,IERR)
      implicit double precision (a-h,o-z)
        DIMENSION       A(N,M), MEMGP(N), NUMGP(NG), GPCEN(NG,M),
     X                  COMP(NG)
C
        BIG = 1.0E+30
        ONE = 0.999
        CMAX = BIG
        ITER = 0
  100   ITER = ITER + 1
        IF (ITER.GT.ITMAX) GOTO 500
        CALL GMEANS(A,N,M,MEMGP,NGP0,NUMGP,GPCEN,NG,IERR)
        CALL COMPCT(A,N,M,NG,MEMGP,GPCEN,COMP,CTOT)
        IF (IERR.NE.0) GOTO 500
        IF (NG.LE.1) GOTO 500
        IF (CTOT.GE.CMAX) GOTO 500
        CMAX = CTOT*ONE
        DO 400 I = 1, N
           X = BIG
           DO 300 K = 1, NG
              Y = 0.0
              DO 200 J = 1, M
                 DIFF = GPCEN(K,J) - A(I,J)
                 Y = Y + DIFF*DIFF
  200         CONTINUE
              IF (Y.GE.X) GOTO 300
              X = Y
              ICL = K
  300      CONTINUE
           MEMGP(I) = ICL
  400   CONTINUE
        GOTO 100
C
  500   RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                 C
C  Optimise the variances of a set of groups, by exchanging       C
C  the objects between groups such that they are minimally        C
C  distant from group centres.                                    C
C                                                                 C
C  Parameters:                                                    C
C                                                                 C
C  N, M, NG         Numbers of rows, columns, groups,             C
C  A(N,M)           initial data,                                 C
C  MEMGP(N)         group memberships,                            C
C  NGP0             minimum acceptable group cardinality,         C
C  NUMGP(NG)        cardinalities of groups,                      C
C  GPCEN(NG,M)      group centres,                                C
C  COMP(NG)         compactness values for the groups,            C
C  CTOT             sum of these compactnesses,                   C
C  ITMAX            max allowed no. of iterations,                C
C  IERR             error indicator (should be zero).             C
C                                                                 C
C  IERR = 1: invalid group number (<1 or >NG), - is number of     C
C  groups correctly specified?  IERR = 2: a group has < minimum   C
C  allowed number of members, - reduce the number of groups and   C
C  try again.                                                     C
C                                                                 C
C  History                                                        C
C  Author: F. Murtagh, Garching/Munich, Feb. 1986.                C
C  Added parameter ITMAX (formerly hardwired).  [F.M., Sept. 91]  C
C                                                                 C
C-----------------------------------------------------------------C
        SUBROUTINE EXCH(A,N,M,MEMGP,NGP0,NUMGP,GPCEN,NG,
     X            COMP,CTOT,ITER,ITMAX,IERR)
      implicit double precision (a-h,o-z)
        DIMENSION       A(N,M), MEMGP(N), NUMGP(NG), GPCEN(NG,M),
     X                  COMP(NG)
C
        BIG = 1.0E+30
        ONE = 0.999
        CALL GMEANS(A,N,M,MEMGP,NGP0,NUMGP,GPCEN,NG,IERR)
        CALL COMPCT(A,N,M,NG,MEMGP,GPCEN,COMP,CTOT)
        IF (IERR.NE.0) GOTO 800
        IF (NG.LE.1) GOTO 800
        ITER = 0
        I = 0
        IS = 0
  100   IS = IS + 1
        IF (IS.GT.N) GOTO 800
  200   I = I + 1
        IF (I.LE.N) GOTO 300
        ITER = ITER + 1
        IF (ITER.GT.ITMAX) GOTO 800
        I = 1
  300   ICL = MEMGP(I)
        NUM = NUMGP(ICL)
        IF (NUM.LE.NGP0) GOTO 100
        V = NUM
        EQ = BIG
        DO 600 K = 1, NG
           X = 0.0
           DO 400 J = 1, M
              DIFF = GPCEN(K,J) - A(I,J)
              X = X + DIFF*DIFF
  400      CONTINUE
           IF (K.NE.ICL) GOTO 500
           FRAC1 = V/(V-1.0)
           EP = X*FRAC1
           GOTO 600
  500      FRAC2 = NUMGP(K)
           FRAC = FRAC2/(FRAC2+1.0)
           EK = FRAC*X
           IF (EK.GE.EQ) GOTO 600
           EQ = EK
           IQ = K
           W = FRAC2
  600   CONTINUE
        IF (EQ.GE.EP*ONE) GOTO 100
        IS = 0
        COMP(ICL) = COMP(ICL) - EP
        COMP(IQ) = COMP(IQ) + EQ
        CTOT = CTOT - EP + EQ
        FRAC1 = 1.0/(V-1.0)
        FRAC2 = 1.0/(W+1.0)
        DO 700 J = 1, M
           VAL = A(I,J)
           GPCEN(ICL,J) = (V*GPCEN(ICL,J)-VAL)*FRAC1
           GPCEN(IQ,J) = (W*GPCEN(IQ,J)+VAL)*FRAC2
  700   CONTINUE
        MEMGP(I) = IQ
        NUMGP(ICL) = NUM - 1
        NUMGP(IQ) = NUMGP(IQ) + 1
        GOTO 200
C
  800   CONTINUE
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                 C
C  Standardize to zero mean and unit standard deviation.          C
C                                                                 C
C  Parameters:                                                    C
C                                                                 C
C  N, M, NG         Numbers of rows, columns, groups,             C
C  A(N,M)           initial data, replaced by standardized values.C
C                                                                 C
C  Author: F. Murtagh, Garching/Munich, Feb. 1986.                C
C                                                                 C
C-----------------------------------------------------------------C
        SUBROUTINE STND(A,N,M)
      implicit double precision (a-h,o-z)
        DIMENSION  A(N,M)
C
        DO 500 J = 1, M
           X = 0.0
           DO 100 I = 1, N
              X = X + A(I,J)
  100      CONTINUE
           XBAR = X/FLOAT(N)
           X = 0.0
           DO 200 I = 1, N
              DIFF = A(I,J) - XBAR
              X = X + DIFF*DIFF
  200      CONTINUE
           IF (X.LE.0.0) X = 1.0
           X = 1.0/SQRT(X)
           DO 300 I = 1, N
              A(I,J) = X*(A(I,J)-XBAR)
  300      CONTINUE
  500   CONTINUE
C
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                 C
C  Determine means of the groups.                                 C
C                                                                 C
C  Parameters:                                                    C
C                                                                 C
C  N, M, NG         Numbers of rows, columns, groups,             C
C  A(N,M)           initial data,                                 C
C  MEMGP(N)         group memberships,                            C
C  NGP0             minimum acceptable group cardinality,         C
C  NUMGP(NG)        cardinalities of groups,                      C
C  GPCEN(NG,M)      group centres,                                C
C  IERR             error indicator (should be zero).             C
C                                                                 C
C  Author: F. Murtagh, Garching/Munich, Feb. 1986.                C
C                                                                 C
C-----------------------------------------------------------------C
        SUBROUTINE GMEANS(A,N,M,MEMGP,NGP0,NUMGP,GPCEN,NG,IERR)
      implicit double precision (a-h,o-z)
        DIMENSION       A(N,M), MEMGP(N), NUMGP(NG), GPCEN(NG,M)
C
        DO 200 K = 1, NG
           NUMGP(K) = 0
           DO 100 J = 1, M
              GPCEN(K,J) = 0.0
  100      CONTINUE
  200   CONTINUE           
C
        DO 500 I = 1, N
           ICL = MEMGP(I)
           IF (ICL.GE.1.AND.ICL.LE.NG) GOTO 300
              IERR = 1
              RETURN
  300      CONTINUE
           NUMGP(ICL) = NUMGP(ICL) + 1
           DO 400 J = 1, M
              GPCEN(ICL,J) = GPCEN(ICL,J) + A(I,J)
  400      CONTINUE
  500   CONTINUE
C
        DO 800 K = 1, NG
           NUM = NUMGP(K)
           IF (NUM.GE.NGP0) GOTO 600
              IERR = 2
              RETURN
  600      CONTINUE
           X = 1.0/FLOAT(NUM)
           DO 700 J = 1, M
              GPCEN(K,J) = GPCEN(K,J)*X
  700      CONTINUE
  800   CONTINUE
C
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                 C
C  Determine compactness of the groups (i.e. their variances).    C
C                                                                 C
C  Parameters:                                                    C
C                                                                 C
C  N, M, NG         Numbers of rows, columns, groups,             C
C  A(N,M)           initial data,                                 C
C  MEMGP(N)         group memberships,                            C
C  GPCEN(NG,M)      group centres,                                C
C  COMP(NG)         variances of groups (output),                 C
C  CTOT             sum of these variances.                       C
C                                                                 C
C  Author: F. Murtagh, Garching/Muenchen, Feb. 1986.              C
C                                                                 C
C-----------------------------------------------------------------C
        SUBROUTINE COMPCT(A,N,M,NG,MEMGP,GPCEN,COMP,CTOT)
      implicit double precision (a-h,o-z)
        DIMENSION       A(N,M), MEMGP(N), GPCEN(NG,M), COMP(NG)
C
        CTOT = 0.0
        DO 100 K = 1, NG
           COMP(K) = 0.0
  100   CONTINUE
C
        DO 300 I = 1, N
           ICL = MEMGP(I)
           X = 0.0
           DO 200 J = 1, M
              DIFF = GPCEN(ICL,J) - A(I,J)
              X = X + DIFF*DIFF
  200      CONTINUE
           COMP(ICL) = COMP(ICL) + X
           CTOT = CTOT + X
  300   CONTINUE
C
        RETURN
        END
