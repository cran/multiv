c
      subroutine sammon(n,m,p,x,y,ndis,dstar,d,
     .                  alpha,maxit,diag,iter,tol,err)
c     SAMMON MAPPING (NON-METRIC MULTIDIMENSIONAL SCALING). 
c     Author: F. Murtagh, May 1992.
      implicit double precision (a-h,o-z)
      integer   n,m,p,i,j,k,iter,maxit,diag
      dimension x(n,m), y(n,p), dstar(ndis), d(ndis)
      double precision      alpha,tol,err
c
      iter = 0
c
      c = 0.0
      do i = 1, n-1
         do j = i+1, n
            ind = n*(i-1) - i*(i-1)/2 + j - i
            c = c + dstar(ind)
         enddo
      enddo
c
 100  continue
      iter = iter + 1

      do i = 1, n
         do k = 1, p
c
            dedy   = 0.0
            d2edy2 = 0.0
            do j   = 1, n
c
               if (j.eq.i) goto 200
               if (i.lt.j) then
                           ind = n*(i-1) - i*(i-1)/2 + j - i
               else
                           ind = n*(j-1) - j*(j-1)/2 + i - j
               endif
               if (dstar(ind).le.0.001) goto 200
               if (d(ind).le.0.001) goto 200
c
               xnumer1 = y(i,k)-y(j,k)
               xnumsq  = xnumer1*xnumer1
               xnumer2 = dstar(ind)-d(ind)
               denom1 = dstar(ind)*d(ind)
c
               dedy = dedy + xnumer1*xnumer2/denom1
c
               d2edy2 = d2edy2 + (1.0/denom1)*
     .            (xnumer2-(xnumsq/d(ind))*(1.0+xnumer2/d(ind)))
c
 200           continue
            enddo
            dedy = -dedy
            d2edy2 = abs(-d2edy2)
            deltay = dedy/d2edy2
c
            y(i,k) = y(i,k) - alpha*deltay
         enddo
       enddo
c
       do i = 1, n-1
          do j = i+1, n
            ind = n*(i-1) - i*(i-1)/2 + j - i
             d(ind) = 0.0
             do k = 1, p
                d(ind) = d(ind) + (y(i,k)-y(j,k))*(y(i,k)-y(j,k))
             enddo
             d(ind) = sqrt(d(ind))
          enddo
       enddo
c
       err = 0.0
       do i = 1, n-1
          do j = i+1, n
            ind = n*(i-1) - i*(i-1)/2 + j - i
             if (dstar(ind).gt.0.001)
     .         err = err + 
     .         ((dstar(ind)-d(ind))*(dstar(ind)-d(ind)))/dstar(ind)
          enddo
       enddo
       err = err/c
       if (diag.eq.1) call dblepr(" ",1,err,1)
c
c       write(6,*) ' Error: ',err
       if (err.le.tol) return
       if (iter.ge.maxit) return
       goto 100
c
       end



