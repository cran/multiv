C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                            C
C  HIERARCHICAL CLUSTERING using (user-specified) criterion. C
C                                                            C
C  Parameters:                                               C
C                                                            C
C  DATA(N,M)         input data matrix,                      C
C  DISS(LEN)         dissimilarities in lower half diagonal  C
C                    storage; LEN = N.N-1/2,                 C
C  IOPT              clustering criterion to be used,        C
C  IA, IB, CRIT      history of agglomerations; dimensions   C
C                    N, first N-1 locations only used,       C
C  MEMBR, NN, DISNN  vectors of length N, used to store      C 
C                    cluster cardinalities, current nearest  C
C                    neighbour, and the dissimilarity assoc. C
C                    with the latter.                        C
C  FLAG              boolean indicator of agglomerable obj./ C
C                    clusters.                               C
C                                                            C
C  F. Murtagh, ESA/ESO/STECF, Garching, February 1986.       C
C                                                            C
C------------------------------------------------------------C
      SUBROUTINE HC(N,M,LEN,IOPT,DATA,IA,IB,CRIT,MEMBR,NN,DISNN,
     X                FLAG,DISS)
      implicit double precision (a-h,o-z)
      DOUBLE PRECISION DATA(N,M),MEMBR(N),DISS(LEN)
      INTEGER IA(N),IB(N)
      DOUBLE PRECISION CRIT(N)
      DIMENSION NN(N),DISNN(N)
      LOGICAL FLAG(N)
      DOUBLE PRECISION INF
      DATA INF/1.E+20/
C
C  Initializations
C
      DO I=1,N
         MEMBR(I)=1.
         FLAG(I)=.TRUE.
      ENDDO
      NCL=N
C
C  Construct dissimilarity matrix
C
      DO I=1,N-1
         DO J=I+1,N
            IND=IOFFSET(N,I,J)
            DISS(IND)=0.
            DO K=1,M
               DISS(IND)=DISS(IND)+(DATA(I,K)-DATA(J,K))**2
            ENDDO
            IF (IOPT.EQ.1) DISS(IND)=DISS(IND)/2.
C           (Above is done for the case of the min. var. method
C            where merging criteria are defined in terms of variances
C            rather than distances.)
          ENDDO
       ENDDO
C
C  Carry out an agglomeration - first create list of NNs
C
      DO I=1,N-1
         DMIN=INF
         DO J=I+1,N
            IND=IOFFSET(N,I,J)
            IF (DISS(IND).GE.DMIN) GOTO 500
               DMIN=DISS(IND)
               JM=J
  500    CONTINUE
         ENDDO
         NN(I)=JM
         DISNN(I)=DMIN
      ENDDO
C
  400 CONTINUE
C     Next, determine least diss. using list of NNs
      DMIN=INF
      DO I=1,N-1
         IF (.NOT.FLAG(I)) GOTO 600
         IF (DISNN(I).GE.DMIN) GOTO 600
            DMIN=DISNN(I)
            IM=I
            JM=NN(I)
  600    CONTINUE
      ENDDO
      NCL=NCL-1
C
C  This allows an agglomeration to be carried out.
C
      I2=MIN0(IM,JM)
      J2=MAX0(IM,JM)
      IA(N-NCL)=I2
      IB(N-NCL)=J2
      CRIT(N-NCL)=DMIN
C
C  Update dissimilarities from new cluster.
C
      FLAG(J2)=.FALSE.
      DMIN=INF
      DO K=1,N
         IF (.NOT.FLAG(K)) GOTO 800
         IF (K.EQ.I2) GOTO 800
         X=MEMBR(I2)+MEMBR(J2)+MEMBR(K)
         IF (I2.LT.K) THEN
                           IND1=IOFFSET(N,I2,K)
                      ELSE
                           IND1=IOFFSET(N,K,I2)
         ENDIF
         IF (J2.LT.K) THEN
                           IND2=IOFFSET(N,J2,K)
                      ELSE
                           IND2=IOFFSET(N,K,J2)
         ENDIF
         IND3=IOFFSET(N,I2,J2)
         XX=DISS(IND3)
C
C  WARD'S MINIMUM VARIANCE METHOD - IOPT=1.
C
         IF (IOPT.EQ.1) THEN
            DISS(IND1)=(MEMBR(I2)+MEMBR(K))*DISS(IND1)+
     X                 (MEMBR(J2)+MEMBR(K))*DISS(IND2)-
     X                 MEMBR(K)*XX
            DISS(IND1)=DISS(IND1)/X
         ENDIF
C
C  SINGLE LINK METHOD - IOPT=2.
C
         IF (IOPT.EQ.2) THEN
            DISS(IND1)=MIN(DISS(IND1),DISS(IND2))
         ENDIF
C
C  COMPLETE LINK METHOD - IOPT=3.
C
         IF (IOPT.EQ.3) THEN
            DISS(IND1)=MAX(DISS(IND1),DISS(IND2))
         ENDIF
C
C  AVERAGE LINK (OR GROUP AVERAGE) METHOD - IOPT=4.
C
         IF (IOPT.EQ.4) THEN
            DISS(IND1)=(MEMBR(I2)*DISS(IND1)+MEMBR(J2)*DISS(IND2))/
     X                 (MEMBR(I2)+MEMBR(J2))
         ENDIF
C
C  MCQUITTY'S METHOD - IOPT=5.
C
         IF (IOPT.EQ.5) THEN
            DISS(IND1)=0.5*DISS(IND1)+0.5*DISS(IND2)
         ENDIF
C
C  MEDIAN (GOWER'S) METHOD - IOPT=6.
C
         IF (IOPT.EQ.6) THEN
            DISS(IND1)=0.5*DISS(IND1)+0.5*DISS(IND2)-0.25*XX
         ENDIF
C
C  CENTROID METHOD - IOPT=7.
C
         IF (IOPT.EQ.7) THEN
            DISS(IND1)=(MEMBR(I2)*DISS(IND1)+MEMBR(J2)*DISS(IND2)-
     X          MEMBR(I2)*MEMBR(J2)*XX/(MEMBR(I2)+MEMBR(J2)))/
     X          (MEMBR(I2)+MEMBR(J2))
            ENDIF
C
         IF (I2.GT.K) GOTO 800
         IF (DISS(IND1).GE.DMIN) GOTO 800
            DMIN=DISS(IND1)
            JJ=K
  800    CONTINUE
      ENDDO
      MEMBR(I2)=MEMBR(I2)+MEMBR(J2)
      DISNN(I2)=DMIN
      NN(I2)=JJ
C
C  Update list of NNs insofar as this is required.
C
      DO I=1,N-1
         IF (.NOT.FLAG(I)) GOTO 900
         IF (NN(I).EQ.I2) GOTO 850
         IF (NN(I).EQ.J2) GOTO 850
         GOTO 900
  850    CONTINUE
C        (Redetermine NN of I:)
         DMIN=INF
         DO J=I+1,N
            IND=IOFFSET(N,I,J)
            IF (.NOT.FLAG(J)) GOTO 870
            IF (I.EQ.J) GOTO 870
            IF (DISS(IND).GE.DMIN) GOTO 870
               DMIN=DISS(IND)
               JJ=J
  870       CONTINUE
         ENDDO
         NN(I)=JJ
         DISNN(I)=DMIN
  900    CONTINUE
      ENDDO
C
C  Repeat previous steps until N-1 agglomerations carried out.
C
      IF (NCL.GT.1) GOTO 400
C
C
      RETURN
      END
C
C
      FUNCTION IOFFSET(N,I,J)
C  Map row I and column J of upper half diagonal symmetric matrix 
C  onto vector.
      IOFFSET=J+(I-1)*N-(I*(I+1))/2
      RETURN
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                C
C  HIERARCHICAL CLUSTERING using Minimum Variance Criterion,     C
C  using the O(N**2) time Nearest Neighbour Chain algorithm.     C
C                                                                C
C  Parameters:                                                   C
C                                                                C
C  DATA(N,M) :      input data,                                  C
C  IA(N), IB(N), CRIT(N) : sequence of agglomerands and          C
C                   values returned (only locations 1 to N-1     C
C                   are of interest),                            C
C  MEMBR(N), DISS(N), ICHAIN(N) : used in the routines to store  C
C                   cluster cardinalities, nearest neighbour     C    
C                   dissimilarities, and the NN-chain.           C
C  FLAG(N) :        (boolean) used to indicate agglomerable      C
C                   objects and clusters.                        C
C                                                                C
C  Reference: Murtagh, Multidimensional Clustering Algorithms,   C
C             Physica-Verlag, 1985.                              C
C                                                                C
C  F. Murtagh, ESA/ESO/STECF, Garching, February 1986.           C
C                                                                C
C----------------------------------------------------------------C
      SUBROUTINE HCON2(N,M,DATA,IA,IB,CRIT,MEMBR,DISS,ICHAIN,FLAG,istat)
      implicit double precision (a-h,o-z)
      DOUBLE PRECISION    MEMBR(N), DATA(N,M), DISS(N), CRIT(N)
      INTEGER ICHAIN(N), IA(N), IB(N)
      DOUBLE PRECISION INF
      LOGICAL FLAG(N)
      DATA INF/1.E+25/
C     EQUIVALENCE (ICHAIN(1),IA(1)),(DISS(1),CRIT(1))
C
      DO 150 I=1,N
         MEMBR(I)=1
         FLAG(I)=.TRUE.
  150 CONTINUE
      NCL=N
      I1=1
C
C  Start the NN-chain:
C
  200 LEN=N
      ICHAIN(LEN)=I1
      DISS(LEN)=INF
C
C  Determine NN of object I1
C
  300 FLAG(I1)=.FALSE.
C
C  Turn off FLAG so that 0 diss. of I1 with self not obtained.
C
      D=DISS(LEN)
      IF (LEN.LT.N) I2=ICHAIN(LEN+1)
C
C  For identical diss.'s, above ensures that RNN will be found.
C
      CALL DETNN(DATA,FLAG,MEMBR,N,M,I1,I2,D)
      FLAG(I1)=.TRUE.
C
C  If LEN = 1 place obj. I2 as second obj. in NN-chain.
C
      IF (LEN.LT.N) GOTO 350
      LEN=LEN-1
      IF (LEN.LT.N-NCL) GOTO 700
      ICHAIN(LEN)=I2
      DISS(LEN)=D
      GOTO 500
C
C  If LEN < N distinguish between having RNN & continuing NN-chain.
C
  350 CONTINUE
      IF (I2.NE.ICHAIN(LEN+1)) GOTO 400
C
C  Have RNN.
C
      NCL=NCL-1
      CALL AGGLOM(I1,I2,D,DATA,MEMBR,FLAG,IA,IB,CRIT,NCL,N,M)
      LEN=LEN+2
      GOTO 500
  400 CONTINUE
C
C  Grow extra link on NN-chain.
C
      IDUM=ICHAIN(LEN+1)
      FLAG(IDUM)=.FALSE.
      LEN=LEN-1
      IF (LEN.LE.N-NCL) GOTO 700
      ICHAIN(LEN)=I2
      DISS(LEN)=D
      GOTO 500
C
C  Select obj. for continuing to grow (or restarting) NN-chain.
C
  500 CONTINUE
      IF (NCL.EQ.1) GOTO 600
      IF (LEN.EQ.N+1) GOTO 550
      I1=ICHAIN(LEN)
      FLAG(I1)=.TRUE.
      IDUM=ICHAIN(LEN+1)
      IF (LEN.LT.N) FLAG(IDUM)=.TRUE.
C
C  Reestablish agglomerability of objects in NN-chain.
C
      GOTO 300
  550 CALL NEXT(FLAG,I1,N)
      GOTO 200
C
  600 CONTINUE
      RETURN
C
C     SHOULD NEVER REACH HERE...
 700  CONTINUE
C  700 WRITE(6,750)
C  750 FORMAT(' ERROR IN NN-CHAIN ROUTINE - INSUFFICIENT CHAIN SPACE'/)
      istat = 1
      return
      END
C------------------------------------------------------------------------------
      SUBROUTINE DETNN(DATA,FLAG,MEM,N,M,I1,I2,D)
      implicit double precision (a-h,o-z)
C
C  Determine a nearest neighbour.
C
      DOUBLE PRECISION DATA(N,M),MEM(N)
      LOGICAL FLAG(N)
C
      DO 200 I=1,N
         IF (.NOT.FLAG(I)) GOTO 200
         DISS=0.
         DO 100 J=1,M
  100    DISS=DISS+(DATA(I1,J)-DATA(I,J))*(DATA(I1,J)-DATA(I,J))
         DISS=DISS*MEM(I)*MEM(I1)/(MEM(I1)+MEM(I))
         IF (DISS.GE.D) GOTO 200
            D=DISS
            I2=I
  200 CONTINUE
C
      RETURN
      END
C------------------------------------------------------------------------------
      SUBROUTINE AGGLOM(I1,I2,D,DATA,MEM,FLAG,IA,IB,CRIT,NCL,N,M)
      implicit double precision (a-h,o-z)
C
C  Carry out an agglomeration.
C
      DOUBLE PRECISION MEM(N),DATA(N,M),CRIT(N)
      INTEGER IA(N),IB(N)
      LOGICAL FLAG(N)
      INTEGER O1,O2,LB,UB
C
C
      O1=MIN0(I1,I2)
      O2=MAX0(I1,I2)
      DO 100 J=1,M
         DATA(O1,J)=( MEM(O1)*DATA(O1,J)+MEM(O2)*DATA(O2,J) )
     X            / (MEM(O1)+MEM(O2))
         DATA(O2,J)=DATA(O1,J)
  100 CONTINUE
      NAGGL=N-NCL
      MEM(O1)=MEM(O1)+MEM(O2)
      FLAG(O2)=.FALSE.
C
C  Keep sorted list of criterion values: find 1st where new crit. fits.
C
      I=NAGGL-1
  120 IF (D.GE.CRIT(I)) GOTO 140
      I=I-1
      IF (I.GE.1) GOTO 120
C
C  Arriving here must mean that D > all crit. values found so far.
C
      I=0
  140 CONTINUE
C
C  Now, shift rightwards from I+1 to AGGL-1 to make room for new crit.
C
      LB=I+1
      UB=NAGGL-1
      IF (LB.GT.UB) GOTO 180
      J=UB
  160 J1=J+1
      IA(J1)=IA(J)
      IB(J1)=IB(J)
      CRIT(J1)=CRIT(J)
      J=J-1
      IF (J.GE.LB) GOTO 160
  180 CONTINUE
      IA(LB)=O1
      IB(LB)=O2
      CRIT(LB)=D
C
      RETURN
      END
C------------------------------------------------------------------------------
      SUBROUTINE NEXT(FLAG,I1,N)
      implicit double precision (a-h,o-z)
C
C  Determine next agglomerable object/cluster.
C
      LOGICAL FLAG(N)
C
      NXT=I1+1
      IF (NXT.GT.N) GOTO 150
      DO 100 I=NXT,N
         IF (FLAG(I)) GOTO 500
  100 CONTINUE
  150 DO 200 I=1,I1
         IF (FLAG(I)) GOTO 500
  200 CONTINUE
C
C     Should never get to here...
C     STOP
C
  500 I1=I
C
      RETURN
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                               C
C  Given a HIERARCHIC CLUSTERING, described as a sequence of    C
C  agglomerations, prepare the seq. of aggloms. and "horiz."    C
C  order of objects for plotting the dendrogram using S routine C
C  'plclust'.                                                   C
C                                                               C
C  Parameters:                                                  C
C                                                               C
C  IA, IB:       vectors of dimension N defining the agglomer-  C
C                 ations.                                       C
C  IIA, IIB:     used to store IA and IB values differently     C
C                (in form needed for S command `plclust`        C
C  IORDER:       "horiz." order of objects for dendrogram       C
C                                                               C
C  F. Murtagh, ESA/ESO/STECF, Garching, June 1991               C
C                                                               C
C  HISTORY                                                      C
C                                                               C
C  Adapted from routine HCASS, which additionally determines    C
C   cluster assignments at all levels, at extra comput. expense C
C                                                               C
C---------------------------------------------------------------C
      SUBROUTINE HCASS2(N,IA,IB,IORDER,IIA,IIB)
      implicit double precision (a-h,o-z)
      INTEGER IA(N),IB(N),IORDER(N),IIA(N),IIB(N)
C
C     Following bit is to get seq. of merges into format acceptable to plclust
C     I coded clusters as lowest seq. no. of constituents; S's `hclust' codes
C     singletons as -ve numbers, and non-singletons with their seq. nos.
C
      DO 912 I=1,N
         IIA(I)=IA(I)
         IIB(I)=IB(I)
  912 CONTINUE
      DO 915 I=1,N-2
C        In the following, smallest (+ve or -ve) seq. no. wanted
         K=MIN(IA(I),IB(I))
         DO 913 J=I+1, N-1
            IF(IA(J).EQ.K) IIA(J)=-I
            IF(IB(J).EQ.K) IIB(J)=-I
  913    CONTINUE
  915 CONTINUE
      DO 916 I=1,N-1
         IIA(I)=-IIA(I)
         IIB(I)=-IIB(I)
  916 CONTINUE
      DO 917 I=1,N-1
         IF (IIA(I).GT.0.AND.IIB(I).LT.0) THEN
            K = IIA(I)
            IIA(I) = IIB(I)
            IIB(I) = K
         ENDIF
         IF (IIA(I).GT.0.AND.IIB(I).GT.0) THEN
            K1 = MIN(IIA(I),IIB(I))
            K2 = MAX(IIA(I),IIB(I))
            IIA(I) = K1
            IIB(I) = K2
         ENDIF
  917 CONTINUE
C
C
C     NEW PART FOR `ORDER'
C
      IORDER(1) =IIA(N-1)
      IORDER(2) =IIB(N-1)
      LOC=2
      DO 175 I=N-2,1,-1
        DO 169 J=1,LOC
          IF(IORDER(J).EQ.I) THEN
C           REPLACE IORDER(J) WITH IIA(I) AND IIB(I)
            IORDER(J)=IIA(I)
            IF (J.EQ.LOC) THEN
                LOC=LOC+1
                IORDER(LOC)=IIB(I)
                GOTO 171
            ENDIF
            LOC=LOC+1
            DO 95 K=LOC,J+2,-1
               IORDER(K)=IORDER(K-1)
  95        CONTINUE
            IORDER(J+1)=IIB(I)
            GOTO 171
          ENDIF
 169    CONTINUE
C       SHOULD NEVER REACH HERE
 171    CONTINUE
 175  CONTINUE
C 
C
C
      DO 181 I=1,N
         IORDER(I) = -IORDER(I)
 181  CONTINUE
C
C
      RETURN
      END




