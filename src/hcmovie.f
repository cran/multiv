C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  HIERARCHICAL CLUSTERING in M-DIMENSIONAL SPACE. F.Murtagh, Sep. '91.
C  Clusters are determined at each agglomeration, so that they can be graphed.
C
C  Data/variables for all routines:
C
C  DATA(N,M)         input data 
C  DAT2(N*M)         copy of part of DATA (i.e. rows corresp. to a cluster).
C                    In vector form.
C  N,M               nos. of rows and cols. assoc. with input data matrix,
C  IA, IB, CRIT      history of aggloms.; subtree sequence numbers in IA, IB;
C                    the seq. no. of a cluster is the smaller seq. no. of the
C                    agglomerands; CRIT contains criterion value at which
C                    agglomeration takes place; dimsensions of these vectors:
C                    N, of which the first N-1 locations, only, are used.    
C                    With this notation, a test for a non-singleton cluster 
C                    is: a cluster seq. no. must additionally exist at a lower
C                    index location in the IA/IB lists.  This lower index 
C                    location is assocated with a subcluster, or a singleton.
C  MEMBR             vector of length N, used to store cluster cardinalities;
C                    initially all values are 1, indicating that each seq. no.
C                    refers to a singleton cluster.  MEMBR is a real vector,
C                    although storing integer values at all times, for comput-
C                    ational convenience.
C  FLAG              boolean indicator of agglomerable objects/clusters.  Note
C                    that since seq. no. can refer to both singleton-object 
C                    and to non-singleton cluster, FLAG is necessary. 
C  IKLASS            1 (resp. 0) in row i, col. j indicates that object j is
C                    present (resp. absent) in cluster i.
C                    N-1 rows only used. Rows <-> clusters, cols. <-> objects.
C  POTCL             potential cluster: N-length vector of 1s and 0s.
C  POTCL1, POTCL2    agglomerable clusters: each of these is an N-length 
C                    integer-valued boolean indicator vector.
C  ARRAY1(M,M), ARRAY2(M,M), VECT1(M), VECT2(M)   Work arrays and vectors, used
C                    to store the covariance matrix and for other purposes by
C                    the principal components analysis routine, PRCOAN.
C  CENTR(M)          vector of centers of a cluster, used by PRCOAN.
C  V1, V2, VTOT      real scalars used to store variance values.
C
C  Subroutines
C
C  INIT    Initialize
C  GBD     Get best dissimilarity (min. var. criterion, with alpha, version)
C  GBD2    Get best dissimilarity (all other agglomerative criteria)
C  AGG     Agglomerate
C  GNCM    Get new cluster memberships
C  CM      Cluster memberships 
C  AL      Assess linearity (min. var. criterion, with alpha, version)
C  AL2     Assess linearity (all other agglomerative criteria: actually this
C         routine does nothing more than determine memberships of agglomerands)
C  CLLIN   CLuster LINearity: driver for principal components routine
C  PRCOAN  Determines covariance matrix, and coordinates SVD
C  LPCOVCL, LPTRED2, LPTQL2  Used by PRCOAN
C------------------------------------------------------------------------------
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  Initializations
C------------------------------------------------------------------------------
      SUBROUTINE        INIT(MEMBR,FLAG,IKLASS,N,M)
      implicit double precision (a-h,o-z)
      DOUBLE PRECISION              MEMBR(N)
      INTEGER           IKLASS(N,N)
      LOGICAL           FLAG(N)
C
      DO 20 I = 1,N
         MEMBR(I) = 1.
         FLAG(I)  = .TRUE.
         DO 10 J = 1,N-1
            IKLASS(J,I) = 0
 10      CONTINUE
 20   CONTINUE
C
      RETURN
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Determine least diss. using original data.  GBD: Get Best Dissimilarity.
C------------------------------------------------------------------------------
      SUBROUTINE GBD
     .       (MEMBR,FLAG,DLEAST,ILEAST,JLEAST,NCL,IA,IB,IKLASS,POTCL,
     .        array1,array2,vect1,vect2,dat2,N,M,ALPHA,DATA,CENTR)
      implicit double precision (a-h,o-z)
      LOGICAL    FLAG(N)
      DOUBLE PRECISION       MEMBR(N),DATA(N,M),INF,CENTR(M),dat2(1)
      DOUBLE PRECISION       array1(m,m),array2(m,m),vect1(m),vect2(m)
      INTEGER    IA(N), IB(N), IKLASS(N,N), POTCL(N)
C
      INF     = 1.E+20
      DLEAST  = INF
      DO 30 I = 1, N
C-- FLAG controls whether or not seq. no. refers to still active obj./clust.
         IF (.NOT.FLAG(I)) GOTO 700
         DO 20 J = 1,N
            IF (.NOT.FLAG(J)) GOTO 600
            IF (I.EQ.J) GOTO 600
C-- ILAGG and JLAGG - pair of last agglomerated objects/clusters:
            ILAGG = 0
            JLAGG = 0
C-- In the case of clusters of > 2 members, get the member seq. nos.:
C-- Subroutine CM: Cluster Memberships
            IF (MEMBR(I).GT.1) CALL CM(I,ILAGG,NCL,IA,IB,N)
            IF (MEMBR(J).GT.1) CALL CM(J,JLAGG,NCL,IA,IB,N)
C-- Variances of a potential 2-member cluster:
            IF (MEMBR(I).EQ.1.AND.MEMBR(J).EQ.1) THEN
               DISSIM = 0.0
               DO 10 K = 1, M
                  DISSIM = DISSIM + 
     .              (DATA(I,K)-DATA(J,K))*(DATA(I,K)-DATA(J,K))
 10            CONTINUE
               V1   = 0.5*DISSIM
               V2   = 0.0
               VTOT = V1
            ENDIF
C-- Variances of a cluster with potentially > 2 members: 
C-- Subroutine AL: Assess Linearity
            IF (MEMBR(I).GT.1.OR.MEMBR(J).GT.1) THEN
             CALL AL(I,J,ILAGG,JLAGG,IKLASS,POTCL,
     .       array1,array2,vect1,vect2,dat2,V1,V2,VTOT,N,M,DATA,CENTR)
            ENDIF
C-- ALPHA is the linearity vs. compactness coeff. Cf. Murtagh & Raftery, '84.
            XCRIT = ALPHA*V1 + (VTOT-V1)
            IF (XCRIT.GE.DLEAST) GOTO 600
               DLEAST = XCRIT
               ILEAST = I
               JLEAST = J
  600       CONTINUE
   20      CONTINUE
  700    CONTINUE
   30 CONTINUE
C
      RETURN
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Determine least diss. using original data.  GBD: Get Best Dissimilarity.
C     GBD2 -- for all criteria except min. var. one.
C     (In min. var. criterion, allow for linearity-determining; not in other.)
C------------------------------------------------------------------------------
      SUBROUTINE GBD2
     .     (MEMBR,FLAG,DLEAST,ILEAST,JLEAST,NCL,IA,IB,IKLASS,POTCL1,
     .      POTCL2,array1,array2,vect1,vect2,dat2,N,M,DATA,METHOD)
      implicit double precision (a-h,o-z)
      LOGICAL    FLAG(N)
      DOUBLE PRECISION       MEMBR(N),DATA(N,M),INF,dat2(1)
      DOUBLE PRECISION       array1(m,m),array2(m,m),vect1(m),vect2(m)
      INTEGER    IA(N), IB(N), IKLASS(N,N), POTCL1(N), POTCL2(N)
C
      INF     = 1.E+20
      DLEAST  = INF
      DO 30 I = 1, N
C-- FLAG controls whether or not seq. no. refers to still active obj./clust.
         IF (.NOT.FLAG(I)) GOTO 700
         DO 20 J = 1,N
            IF (.NOT.FLAG(J)) GOTO 600
            IF (I.EQ.J) GOTO 600
C-- ILAGG and JLAGG - pair of last agglomerated objects/clusters:
            ILAGG = 0
            JLAGG = 0
C-- In the case of clusters of > 2 members, get the member seq. nos.:
C-- Subroutine CM: Cluster Memberships
            IF (MEMBR(I).GT.1) CALL CM(I,ILAGG,NCL,IA,IB,N)
            IF (MEMBR(J).GT.1) CALL CM(J,JLAGG,NCL,IA,IB,N)
            DMIN = 1000000.0
C-- Variances of a potential 2-member cluster:
            IF (MEMBR(I).EQ.1.AND.MEMBR(J).EQ.1) THEN
               DISSIM = 0.0
               DO 10 K = 1, M
                  DISSIM = DISSIM + 
     .              (DATA(I,K)-DATA(J,K))*(DATA(I,K)-DATA(J,K))
 10            CONTINUE
               DMIN = DISSIM
            ENDIF
            DISMIN = 1000000.0
C-- Variances of a cluster with potentially > 2 members: 
C-- Subroutine AL2: Modeled on `Assess Linearity' routine (AL), but in this
C-- we firstly only det. the clusters' members
            IF (MEMBR(I).GT.1.OR.MEMBR(J).GT.1) THEN
                CALL AL2(I,J,ILAGG,JLAGG,IKLASS,POTCL1,POTCL2,
     .          array1,array2,vect1,vect2,dat2,N,M,DATA)
C-- Now POTCL1 has 1 values to indicate members of cluter 1; ditto for POTCL2.
                IF (METHOD.EQ.2) THEN
C--                Single link
                   DISMIN = 100000.0
                   DO 400 I1 = 1, N
                      IF (POTCL1(I1).EQ.1) THEN
                         DO 380 I2 = 1, N
                            IF (POTCL2(I2).EQ.1) THEN
                               DISSIM = 0.0
                               DO 360 K2 = 1, M
                                  DISSIM=DISSIM+
     .              (DATA(I1,K2)-DATA(I2,K2))*(DATA(I1,K2)-DATA(I2,K2))
 360                           CONTINUE
                               IF (DISSIM.LT.DISMIN) DISMIN = DISSIM
                            ENDIF
 380                     CONTINUE
                      ENDIF
 400               CONTINUE
                ENDIF
C--
                IF (METHOD.EQ.3) THEN
C--                Complete link
                   DISMIN = 0.0
                   DO 500 I1 = 1, N
                      IF (POTCL1(I1).EQ.1) THEN
                         DO 480 I2 = 1, N
                            IF (POTCL2(I2).EQ.1) THEN
                               DISSIM = 0.0
                               DO 460 K2 = 1, M
                                  DISSIM=DISSIM+
     .              (DATA(I1,K2)-DATA(I2,K2))*(DATA(I1,K2)-DATA(I2,K2))
 460                           CONTINUE
                               IF (DISSIM.GT.DISMIN) DISMIN = DISSIM
                            ENDIF
 480                     CONTINUE
                      ENDIF
 500               CONTINUE
                ENDIF
C--
                IF (METHOD.EQ.4) THEN
C--                Average link
                   DISMIN = 0.0
                   DISSIM = 0.0
                   INLINKS = 0
                   DO 590 I1 = 1, N
                      IF (POTCL1(I1).EQ.1) THEN
                         DO 580 I2 = 1, N
                            IF (POTCL2(I2).EQ.1) THEN
                               INLINKS = INLINKS + 1
                               DO 560 K2 = 1, M
                                  DISSIM=DISSIM+
     .              (DATA(I1,K2)-DATA(I2,K2))*(DATA(I1,K2)-DATA(I2,K2))
 560                           CONTINUE
                            ENDIF
 580                     CONTINUE
                      ENDIF
 590                  CONTINUE
                      DISMIN = DISSIM/FLOAT(INLINKS)
                ENDIF
C--
            ENDIF
C
            
            XCRIT = MIN(DMIN,DISMIN)
            IF (XCRIT.GE.DLEAST) GOTO 600
               DLEAST = XCRIT
               ILEAST = I
               JLEAST = J
  600       CONTINUE
   20      CONTINUE
  700    CONTINUE
   30 CONTINUE
C
      RETURN
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Carry out an agglomeration
C------------------------------------------------------------------------------
      SUBROUTINE AGG(ILEAST,JLEAST,DLEAST,NCL,IA,IB,CRIT,MEMBR,FLAG,N)
      implicit double precision (a-h,o-z)
      INTEGER    IA(N), IB(N)
      DOUBLE PRECISION       CRIT(N), MEMBR(N)
      LOGICAL    FLAG(N)
C
      I           = MIN0(ILEAST,JLEAST)
      J           = MAX0(ILEAST,JLEAST)
      ILEAST      = I
      JLEAST      = J
      IA(N-NCL)   = I
      IB(N-NCL)   = J
      CRIT(N-NCL) = DLEAST
C
      RETURN
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  Update cluster memberships array, iklass
C  GNCM: Get New Cluster Memberships
C------------------------------------------------------------------------------
      SUBROUTINE GNCM(I2,J2,NCL,IA,IB,MEMBR,IKLASS,FLAG,N)
      implicit double precision (a-h,o-z)
      INTEGER    IA(N), IB(N), IKLASS(N,N)
      DOUBLE PRECISION       MEMBR(N)
      LOGICAL    FLAG(N)
C
C     We are dealing with clusters I2 and J2.  Are these multiobject clusters,
C     or singletons?  (Remember, the seq. no. values contained in I2 and J2 
C     could be either.)  If I2 (resp. J2) is a singleton, then the value of
C     I2 (resp. J2) is the one member.  If I2 (resp. J2) is a non-singleton
C     cluster, then we must pick up the most recently updated row of array
C     IKLASS.  1 and 0 values in this row indicate memberships in the cluster.
C
C     Use cardinalities vector, MEMBR, to find whether we have a singleton or
C     a non-singleton cluster.
      IAINE = MEMBR(I2)
      IBENJ = MEMBR(J2)
C      
      IF (IAINE.GT.1) GOTO 610
C        Singleton, i.e. one object with seq. no. I2 
         IKLASS(N-NCL,I2)=1
         GOTO 620
  610    CONTINUE
C        Non-singleton, i.e. multi-object cluster
         CALL CM(I2, ILAGG, NCL,IA,IB,N)
         DO 10 ICLST = 1,N
            IF (IKLASS(ILAGG,ICLST).EQ.1) IKLASS(N-NCL,ICLST)=1
   10    CONTINUE
  620 CONTINUE
C
      IF (IBENJ.GT.1) GOTO 630
C        Singleton, i.e. one object with seq. no. J2
         IKLASS(N-NCL,J2)=1
         GOTO 640
  630    CONTINUE
C        Non-singleton cluster
         CALL CM(J2, JLAGG, NCL,IA,IB,N)
         DO 20 ICLST = 1,N
            IF (IKLASS(JLAGG,ICLST).EQ.1) IKLASS(N-NCL,ICLST)=1
   20    CONTINUE
  640 CONTINUE
C
C--  Finally, update MEMBR and FLAG:
      MEMBR(I2)    = MEMBR(I2) + MEMBR(J2)
      FLAG(J2)     = .FALSE.
C
      RETURN
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Given the seq. no. (identifier) of a non-singleton cluster, find where
C     it was most recently met with in the sequence of agglomerations.  The
C     agglomeration number will then provide us with the row of array IKLASS
C     (not used in this subroutine!) which in turn allows us to read off the
C     cluster's members.
C     ICL is the cluster seq. no., ILAGG is the value sought, and 
C     IA and IB relate to the sequence of agglomerations.
C     Subroutine CM: Cluster Memberships
C------------------------------------------------------------------------------
      SUBROUTINE CM(ICL, ILAGG, NCL, IA, IB,N)
      implicit double precision (a-h,o-z)
      INTEGER    IA(N), IB(N)
C
C     What agglom was ICL last met at?
      NAGG = N-NCL-1
      DO 5 I = NAGG,1,-1
         IF (IA(I).EQ.ICL.OR.IB(I).EQ.ICL) GOTO 10
   5  CONTINUE
C--   Should never get to here.
   10 CONTINUE
      ILAGG = I
C
      RETURN
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Subroutine AL: Assess Linearity
C Currently uses principal components analysis
C------------------------------------------------------------------------------
      SUBROUTINE AL(I,J,ILAST,JLAST,IKLASS,POTCL,
     .           array1,array2,vec1,vec2,dat2,V1,V2,VTOT,N,M,DATA,CENTR)
      implicit double precision (a-h,o-z)
      DOUBLE PRECISION       array1(m,m),array2(m,m),vec1(m),vec2(m)
      DOUBLE PRECISION       DATA(N,M),CENTR(M),dat2(1)
      INTEGER    IKLASS(N,N), POTCL(N)
C
      DO 10 IOBJ = 1, N
         POTCL(IOBJ) = 0
 10   CONTINUE
C
C     J is NN of I.  If I is singleton, then ILAST is 0; ditto for J, JLAST.
      IF (ILAST.EQ.0) THEN 
            POTCL(I) = 1
         ELSE 
            DO 20 IOBJ = 1, N
               IF (IKLASS(ILAST,IOBJ).EQ.1) POTCL(IOBJ) = 1
 20         CONTINUE
      ENDIF
      IF (JLAST.EQ.0) THEN 
            POTCL(J) = 1
         ELSE 
            DO 30 IOBJ = 1, N
               IF (IKLASS(JLAST,IOBJ).EQ.1) POTCL(IOBJ) = 1
 30         CONTINUE
      ENDIF
C
C--   Subroutine CLLIN: CLuster LINearity
      CALL CLLIN(POTCL,array1,array2,vec1,vec2,dat2,V1,V2,VTOT,N,M,
     .            DATA,CENTR)
C
      RETURN
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Subroutine AL2: Modeled on AL, Assess Linearity
C Det. clusters' memberships
C------------------------------------------------------------------------------
      SUBROUTINE AL2(I,J,ILAST,JLAST,IKLASS,POTCL1, POTCL2,
     .           array1,array2,vec1,vec2,dat2,N,M,DATA)
      implicit double precision (a-h,o-z)
      DOUBLE PRECISION       array1(m,m),array2(m,m),vec1(m),vec2(m)
      DOUBLE PRECISION       DATA(N,M),dat2(1)
      INTEGER    IKLASS(N,N), POTCL1(N), POTCL2(N)
C
      DO 10 IOBJ = 1, N
         POTCL1(IOBJ) = 0
         POTCL2(IOBJ) = 0
 10   CONTINUE
C
C     J is NN of I.  If I is singleton, then ILAST is 0; ditto for J, JLAST.
      IF (ILAST.EQ.0) THEN 
            POTCL1(I) = 1
         ELSE 
            DO 20 IOBJ = 1, N
               IF (IKLASS(ILAST,IOBJ).EQ.1) POTCL1(IOBJ) = 1
 20         CONTINUE
      ENDIF
      IF (JLAST.EQ.0) THEN 
            POTCL2(J) = 1
         ELSE 
            DO 30 IOBJ = 1, N
               IF (IKLASS(JLAST,IOBJ).EQ.1) POTCL2(IOBJ) = 1
 30         CONTINUE
      ENDIF
C
      RETURN
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Subroutine CLLIN: CLuster LINearity
C Uses principal components analysis
C------------------------------------------------------------------------------
        SUBROUTINE CLLIN(POTCL,array1,array2,vect1,vect2,dat2,V1,V2,
     .                   VTOT,N,M,DATA,CENTR)
      implicit double precision (a-h,o-z)
        INTEGER    POTCL(N)
        DOUBLE PRECISION       DATA(N,M),CENTR(M),vect2(m),dat2(1)
        DOUBLE PRECISION       array1(m,m),array2(m,m),vect1(m)
C
        NROW = 0
        DO 20 I = 1, N
           IF (POTCL(I).EQ.1) THEN
              NROW = NROW + 1
              DO 10 J = 1, M
c                 dat2(NROW,J) = DATA(I,J)
                  dat2(j+(nrow-1)*m) = data(i,j)
  10          CONTINUE
           ENDIF
  20    CONTINUE
C
	METHOD = 2
	IPRINT = 0
C       PCA follows.  Funny name, to distinguish it from other slightly
C       modified replications which might be used in the same S session
	CALL PRCOAN(N,NROW,M,DAT2,METHOD,IPRINT,ARRAY1,VECT1,VECT2,
     .              ARRAY2,IERR,CENTR)
             
C	IF (IERR.NE.0) GOTO 9000
C	GOTO 9900
C 9000	WRITE (6,*) ' ABNORMAL END: IERR =', IERR
 9900	CONTINUE
C
C-- For info: eigenvalues are in VECT1(M-J+1) for J=1,2,...min(NROW,M,7)
C-- 1st e-value is identical to sum of squared values in DAT2(*,1),
C-- 2nd e-value is identical to sum of squared values in DAT2(*,2), etc.
C-- Variance along 1st axis
        V1 = VECT1(M)               
C-- Variance along 2nd axis
        V2 = VECT1(M-1)             
C-- Total variance, - or at least 1st 7 variances
        VTOT = 0.0
        DO 40 J = 1, MIN(NROW,M,7)
           VTOT = VTOT + VECT1(M-J+1)
 40     CONTINUE
C

        RETURN
	END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  Carry out a PRINCIPAL COMPONENTS ANALYSIS (SVD, KARHUNEN-LOEVE EXPANSION).  
C  To call:   CALL PRCOAN(N,NROW,M,DATA,METHOD,IPRINT,A1,W1,W2,A2,IERR) 
C             where
C                                
C  N,NROW,M : integer dims. (NROW is actual # rows in use, out of N) of ...
C  DATA     : input data.              
C             On output, DATA contains in first 7 columns the projections of 
C             the row-points on the first 7 principal components. 
C  A1       : variance-covariance matrix, dimensions M * M.            
C             On output, A1 contains in the first 7 columns the projections of 
C             the column-points on the first 7 principal components. 
C  W1,W2    : real vectors of dimension M (see called routines for use).
C             On output, W1 contains the eigenvalues (read with index M-J+1, 
C             where running index J is J=1,7)
C  A2       : real array of dimensions M * M (see routines for use).   
C  IERR     : error indicator (normally zero).            
C             If IERR > 0, then its value indicates the eigenvalue for which
C             no convergence was obtained.      
C                                                      
C  Inputs here are N, M, DATA, (and IERR).
C                                 
C  HISTORY
C  Adapted from PCA routine written by F. Murtagh in 1986 by F. Murtagh in 
C  October 1990.     
C------------------------------------------------------------------------------
        SUBROUTINE PRCOAN(N,NROW,M,DATA,METHOD,IPRINT,A,
     .                    W,FV1,Z,IERR,CENTR)
      implicit double precision (a-h,o-z)
        DOUBLE PRECISION    DATA(1), A(M,M), W(M), FV1(M), Z(M,M)
        double precision    CENTR(M)
C
C-- Form covariance matrix.
        CALL LPCOVCL(N,M,NROW,DATA,CENTR,A)
C
C-- Carry out eigenreduction.
        M2 = M
        CALL LPTRED2(M,M2,A,W,FV1,Z)
        CALL LPTQL2(M,M2,W,FV1,Z,IERR)
C
        RETURN  
        END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  Det. covariances of cols. First det. means of cols., storing in WORK.
C------------------------------------------------------------------------------
        SUBROUTINE    LPCOVCL(N,M,NROW,DATA,WORK,OUT)
      implicit double precision (a-h,o-z)
        DIMENSION     DATA(1), OUT(M,M), WORK(M)
C
        DO 30 J = 1, M
           WORK(J) = 0.0
           DO 20 I = 1, NROW
               WORK(J) = WORK(J) + DATA(J+(I-1)*M)
   20      CONTINUE
           WORK(J) = WORK(J)/FLOAT(NROW)
   30   CONTINUE
C
C-- Now centre the column points.
        DO 50 I = 1, NROW
           DO 40 J = 1, M
               DATA(J+(I-1)*M) = DATA(J+(I-1)*M)-WORK(J)
   40      CONTINUE
   50   CONTINUE
C
C-- Finally calculate the cross product matrix of the redefined data matrix.
        DO 80 J1 = 1, M
           DO 70 J2 = J1, M
              OUT(J1,J2) = 0.0
              DO 60 I = 1, NROW
                  OUT(J1,J2) = OUT(J1,J2) + 
     .                         DATA(J1+(I-1)*M)*DATA(J2+(I-1)*M)
   60         CONTINUE
              OUT(J2,J1) = OUT(J1,J2)
   70      CONTINUE
   80   CONTINUE
C
        RETURN
        END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Reduce a real, symmetric matrix to a symmetric, tridiagonal matrix. 
C To call:    CALL LPTRED2(NM,N,A,D,E,Z)    where
C                                     
C NM = row dimension of A and Z;      
C N  = order of matrix A (will always be <= NM);
C A  = symmetric matrix of order N to be reduced to tridiag. form;
C D  = vector of dim. N containing, on output, diagonal elts. of
C      tridiagonal matrix.               
C E  = working vector of dim. at least N-1 to contain subdiagonal
C      elements.                          
C Z  = matrix of dims. NM by N containing, on output, orthogonal
C      transformation matrix producing the reduction. 
C                                            
C Normally a call to PTQL2 will follow the call to PTRED2 in order to
C produce all eigenvectors and eigenvalues of matrix A.
C                                        
C Algorithm used: Martin et al., Num. Math. 11, 181-195, 1968. 
C Reference: Smith et al., Matrix Eigensystem Routines - EISPACK
C Guide, Lecture Notes in Computer Science 6, Springer-Verlag, 
C 1976, pp. 489-494.                     
C-----------------------------------------------------------------------------
        SUBROUTINE LPTRED2(NM,N,A,D,E,Z)
      implicit double precision (a-h,o-z)
        DOUBLE PRECISION A(NM,N),D(N),E(N),Z(NM,N)
C
        DO 100 I = 1, N
           DO 100 J = 1, I
              Z(I,J) = A(I,J)
  100   CONTINUE
        IF (N.EQ.1) GOTO 320
        DO 300 II = 2, N
           I = N + 2 - II
           L = I - 1
           H = 0.0
           SCALE = 0.0
           IF (L.LT.2) GOTO 130
           DO 120 K = 1, L
              SCALE = SCALE + ABS(Z(I,K))
  120      CONTINUE
           IF (SCALE.NE.0.0) GOTO 140
  130      E(I) = Z(I,L)
           GOTO 290
  140      DO 150 K = 1, L
              Z(I,K) = Z(I,K)/SCALE
              H = H + Z(I,K)*Z(I,K)
  150      CONTINUE
C
           F = Z(I,L)
           G = -SIGN(SQRT(H),F)
           E(I) = SCALE * G
           H = H - F * G
           Z(I,L) = F - G
           F = 0.0
C
           DO 240 J = 1, L
              Z(J,I) = Z(I,J)/H
              G = 0.0
C             Form element of A*U.
              DO 180 K = 1, J
                 G = G + Z(J,K)*Z(I,K)
  180         CONTINUE
              JP1 = J + 1
              IF (L.LT.JP1) GOTO 220
              DO 200 K = JP1, L
                 G = G + Z(K,J)*Z(I,K)
  200         CONTINUE
C             Form element of P where P = I - U U' / H .
  220         E(J) = G/H
              F = F + E(J) * Z(I,J)
  240      CONTINUE
           HH = F/(H + H)
C          Form reduced A.
           DO 260 J = 1, L
              F = Z(I,J)
              G = E(J) - HH * F
              E(J) = G
              DO 250 K = 1, J
                 Z(J,K) = Z(J,K) - F*E(K) - G*Z(I,K)
  250         CONTINUE
  260      CONTINUE
  290      D(I) = H
  300   CONTINUE
  320   D(1) = 0.0
        E(1) = 0.0
C       Accumulation of transformation matrices.
        DO 500 I = 1, N
           L = I - 1
           IF (D(I).EQ.0.0) GOTO 380
           DO 360 J = 1, L
              G = 0.0
              DO 340 K = 1, L
                 G = G + Z(I,K) * Z(K,J)
  340         CONTINUE
              DO 350 K = 1, L
                 Z(K,J) = Z(K,J) - G * Z(K,I)
  350         CONTINUE
  360      CONTINUE
  380      D(I) = Z(I,I)
           Z(I,I) = 1.0
           IF (L.LT.1) GOTO 500
           DO 400 J = 1, L
              Z(I,J) = 0.0
              Z(J,I) = 0.0
  400      CONTINUE
  500   CONTINUE
C
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Determine eigenvalues and eigenvectors of a symmetric, tridiagonal matrix
C To call:    CALL LPTQL2(NM,N,D,E,Z,IERR)    where
C                             
C NM   = row dimension of Z;    
C N    = order of matrix Z;      
C D    = vector of dim. N containing, on output, eigenvalues;
C E    = working vector of dim. at least N-1;      
C Z    = matrix of dims. NM by N containing, on output, eigenvectors;
C IERR = error, normally 0, but 1 if no convergence.    
C                      
C Normally the call to PTQL2 will be preceded by a call to PTRED2 in 
C order to set up the tridiagonal matrix.   
C                     
C Algorithm used: QL method of Bowdler et al., Num. Math. 11,
C 293-306, 1968.                   
C Reference: Smith et al., Matrix Eigensystem Routines - EISPACK 
C Guide, Lecture Notes in Computer Science 6, Springer-Verlag,
C 1976, pp. 468-474.                 
C-----------------------------------------------------------------------------
        SUBROUTINE LPTQL2(NM,N,D,E,Z,IERR)
      implicit double precision (a-h,o-z)
        DOUBLE PRECISION    D(N), E(N), Z(NM,N)
        DATA    EPS/1.E-12/
C
        IERR = 0
        IF (N.EQ.1) GOTO 1001
        DO 100 I = 2, N
           E(I-1) = E(I)
  100   CONTINUE
        F = 0.0
        B = 0.0
        E(N) = 0.0
C
        DO 240 L = 1, N
           J = 0
           H = EPS * (ABS(D(L)) + ABS(E(L)))
           IF (B.LT.H) B = H
C          Look for small sub-diagonal element.
           DO 110 M = L, N
              IF (ABS(E(M)).LE.B) GOTO 120
C             E(N) is always 0, so there is no exit through
C             the bottom of the loop.
  110      CONTINUE
  120      IF (M.EQ.L) GOTO 220
  130      IF (J.EQ.30) GOTO 1000
           J = J + 1
C          Form shift.
           L1 = L + 1
           G = D(L)
           P = (D(L1)-G)/(2.0*E(L))
           R = SQRT(P*P+1.0)
           D(L) = E(L)/(P+SIGN(R,P))
           H = G-D(L)
C
           DO 140 I = L1, N
              D(I) = D(I) - H
  140      CONTINUE
C
           F = F + H
C          QL transformation.
           P = D(M)
           C = 1.0
           S = 0.0
           MML = M - L
C
           DO 200 II = 1, MML
              I = M - II
              G = C * E(I)
              H = C * P
              IF (ABS(P).LT.ABS(E(I))) GOTO 150
              C = E(I)/P
              R = SQRT(C*C+1.0)
              E(I+1) = S * P * R
              S = C/R
              C = 1.0/R
              GOTO 160
  150         C = P/E(I)
              R = SQRT(C*C+1.0)
              E(I+1) = S * E(I) * R
              S = 1.0/R
              C = C * S
  160         P = C * D(I) - S * G
              D(I+1) = H + S * (C * G + S * D(I))
C             Form vector.
              DO 180 K = 1, N
                 H = Z(K,I+1)
                 Z(K,I+1) = S * Z(K,I) + C * H
                 Z(K,I) = C * Z(K,I) - S * H
  180         CONTINUE
  200      CONTINUE
           E(L) = S * P
           D(L) = C * P
           IF (ABS(E(L)).GT.B) GOTO 130
  220      D(L) = D(L) + F
  240   CONTINUE
C
C       Order eigenvectors and eigenvalues.
        DO 300 II = 2, N
           I = II - 1
           K = I
           P = D(I)
           DO 260 J = II, N
              IF (D(J).GE.P) GOTO 260
              K = J
              P = D(J)
  260      CONTINUE
           IF (K.EQ.I) GOTO 300
           D(K) = D(I)
           D(I) = P
           DO 280 J = 1, N
              P = Z(J,I)
              Z(J,I) = Z(J,K)
              Z(J,K) = P
  280      CONTINUE
  300   CONTINUE
C
        GOTO 1001
C       Set error - no convergence after 30 iterns.
 1000   IERR = L
 1001   RETURN
        END  
C-----------------------------------------------------------------------------









