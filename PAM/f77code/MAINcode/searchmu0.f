      iteration=0
      niter=100
      do while(iteration.le.niter)
         iteration=iteration+1
         zntot=0.d0
         znd=0.d0
         znp=0.d0
         do i=-L,L
            omega=ome(i)
            zr=omega+xi*0.01d0
            alpha=zr+xmu-ed0
            zeta=zr+xmu-ep0-(V**2/(alpha))
            root = cdsqrt((zeta)**2-d2)
            Gminus=2.d0*one/(zeta+root)
            dosminus=-dimag(Gminus)/pi
            Gd =one/(alpha) + (V**2/alpha**2)*Gminus
            dosd=-dimag(Gd)/pi
            if(omega.le.0)then
               znp=znp+dosminus*dome
               znd=znd+dosd*dome
            endif
         enddo
         
         znp=2.d0*znp
         znd=2.d0*znd
         zntot=znd+znp

*     Store the old value of the parameters for the searching
         nindex1=nindex
         ndelta1=ndelta
*     Check whether the particle number is above, below or in the error respect to the chosen ntot
         if((zntot.ge.nread+nepserr))then
            nindex=-1
         elseif(zntot.le.nread-nepserr)then
            nindex=1
         else
            nindex=0
            niter=1
         endif

*     Start changing the chemical potential value
         if(nindex1+nindex.eq.0)then !avoid loop forth and back
            ndelta=ndelta1/2.d0 !decreasing the step
            xmu=xmu+dfloat(nindex)*ndelta
         else
            ndelta=ndelta1
            xmu=xmu+dfloat(nindex)*ndelta
         endif
         
         write(88,*)iteration,xmu
      enddo
      write(88,*)''
      print*,'----------------------'
      print*,'#iteration=',iteration
      print*,'        mu=',xmu-gzerop,'/',xmu0
      print*,'        nd=',znd
      print*,'        np=',znp
      print*,'      ntot=',zntot,'/',nread
      print*,''


      do i=-L,L
         omega=ome(i)
         zi=xi*omega
         zr=nu(i)+xi*0.01d0
         
         delta=ed0-ep0
         alpha=zr+xmu-ed0
         ialpha=zi+xmu-ed0

         zeta=zr+xmu-ep0-(V**2/(alpha))
         izeta=zi+xmu-ep0-(V**2/(ialpha))
         
         sigmapp=V**2/(alpha)
         isigmapp=V**2/(ialpha)
         
         root = cdsqrt((zeta)**2-d2) !cmplx(gamma**2 -D**2)
         iroot= cdsqrt((izeta)**2-d2) !cmplx(igamma**2 -D**2)

         isq=dimag(iroot)
         w=dimag(zi)
         isig=w*isq/dabs(w*isq)

c     get the p -DOS
         iGminus=2.d0*one/(izeta+isig*iroot)
         Gminus=2.d0*one/(zeta+root) !cmplx((2.d0/(D**2))*(gamma - root))
         dosminus=-dimag(Gminus)/pi
         if(i.eq.0)dosp0=dosminus
c     get the d-DOS
         Gd =one/(alpha) + (V**2/alpha**2)*Gminus
         iGd=one/(ialpha) + (V**2/ialpha**2)*iGminus
         dosd=-dimag(Gd)/pi
         if(i.eq.0)dosd0=dosd
         poles=ome(i)-ep0+xmu-real(sigmapp)
         write(21,*) nu(i),dosminus !  p-DOS
         write(22,*) nu(i),sigmapp
         write(23,*)omega,dimag(iGminus),real(iGminus)
         write(24,*)omega,dimag(isigmapp),real(isigmapp)
         write(25,*) nu(i),dosd !  d-DOS
         write(26,*)omega,dimag(iGd),real(iGd)
         if(abs(poles).le.1.d0)then
            write(27,*)omega,poles
         endif
      enddo
