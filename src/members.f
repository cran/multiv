C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                               C
C  Given a HIERARCHIC CLUSTERING, described as a sequence of    C
C  agglomerations, derive the assignments into clusters for     C
C  all (except the first and last) levels of the hierarchy.     C
C                                                               C
C  PARAMETERS                                                   C
C                                                               C
C  IA, IB:       vectors of dimension N defining the agglomer-  C
C                 ations.                                       C
C  HVALS:        vector of dim. N, used internally only.        C
C  ICLASS:       array of cluster assignments; dim. N by LEV.   C
C  IIA, IIB:     alt. coding of info. in IA and IB.             C
C                                                               C
C  F. Murtagh, Feb. 1986.                                       C
C                                                               C
C  HISTORY                                                      C
C                                                               C
C  Derived from routine HCASS.F, June 1991.    F.M.             C
C                                                               C
C---------------------------------------------------------------C
      SUBROUTINE ASSGN(N,NPLUS1,IA,IB,ICLASS,HVALS,IIA,IIB)
      implicit double precision (a-h,o-z)
      INTEGER IA(N),IB(N),ICLASS(NPLUS1,N),HVALS(N),IIA(N),IIB(N)
C
C     1 - CONVERT CODING OF SEQ. OF AGGLOMS. FROM S-CONVENTION TO MINE
C     S-CONVENTION: -VE VAL. => OBS. SEQ. NO.; +VE VAL. => CLUSTER REF.
C     MY CONVENTION: ALL VALS. +VE.  IF VAL. ALREADY REFERENCED, THEN 
C     CLUSTER; COL.1 VALS. ALWAYS < COL.2 VALS.
c
c     1st step: for each values, flip in sign
c     Note that n should be # obs. minus 1
      do 100 i = 1, n
         iia(i) = -ia(i)
         iib(i) = -ib(i)
 100  continue
c
c     2nd step: leave references to observations alone (these are now +ve)
c     Clusters: seek 
      do 500 i = 1, n
         if (iia(i).lt.0) then
c           Because it is -ve, iia(i) is a seq. no. of a cluster.
c           Find this cluster, and get the lowest seq. no. of obs. in it.
            iclust = -iia(i)
c           lowest val. among iia(iclust) and iib(iclust)
            iia(i) = iia(iclust)
         endif
c        Same for iib(i)
         if (iib(i).lt.0) then
            iclust = -iib(i)
            iib(i) = iia(iclust)
         endif
c        Ensure that iia(i) < iib(i)
         if (iia(i).gt.iib(i)) then
            itemp = iia(i)
            iia(i) = iib(i)
            iib(i) = itemp
         endif
 500  continue
C
C  Pick out the clusters which the N objects belong to,
C  at levels N-2, N-3, ... 1 of the hierarchy.
C  The clusters are identified by the lowest seq. no. of
C  their members.
C  There are 2, 3, ... N-1 clusters, respectively, for the
C  above levels of the hierarchy.
C
c     NOTE: in following N = # aggloms.; # obs. = N+1
      HVALS(1)=1
      HVALS(2)=IIB(N)
      LOC=3
      DO 59 I=N-1,1,-1
         DO 52 J=1,LOC-1
            IF (IIA(I).EQ.HVALS(J)) GOTO 54
  52     CONTINUE
         HVALS(LOC)=IIA(I)
         LOC=LOC+1
  54     CONTINUE
         DO 56 J=1,LOC-1
            IF (IIB(I).EQ.HVALS(J)) GOTO 58
  56     CONTINUE
         HVALS(LOC)=IIB(I)
         LOC=LOC+1
  58     CONTINUE
  59  CONTINUE
C
      DO 400 LEVEL=1,N-1
         DO 200 I=1,NPLUS1
            ICL=I
            DO 150 ILEV=1,LEVEL
  150       IF (IIB(ILEV).EQ.ICL) ICL=IIA(ILEV)
            NCL = NPLUS1 - LEVEL
            ICLASS(I,NCL-1)=ICL
  200    CONTINUE
  400  CONTINUE
C
      DO 120 I=1,NPLUS1
      DO 120 J=1,N-1
      DO 110 K=2,N
      IF (ICLASS(I,J).NE.HVALS(K)) GOTO 110
         ICLASS(I,J)=K
         GOTO 120
  110 CONTINUE
  120 CONTINUE
C
      RETURN
      END




