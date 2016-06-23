      Etot=zero
      Ekin=zero
      Epot=zero
      Ehyb=zero
      Emu=zero
      Eepsi=zero
      free=0.d0
      emax=d
      emin=-d
      depsi=(emax-emin)/dfloat(Lepsi)

      do j=1,Lepsi
         e=emin+dfloat(j)*depsi
         dens=ddens(e,1)
         xfermi=fermi(e,0.d0,beta)
         free=free+depsi*e*dens*xfermi
      enddo
      open(11,file='Efree.analytic')
      write(11,*)temp,free

      do i=1,multi*L
         omega=pi/beta*dfloat(2*i-1)
         zi=xi*omega
         ialpha=zi+xmu-ed0

         izeta0=zi
         izeta=zi+xmu-ep0-(V**2/(ialpha))

!     Green's function, Sigmapp
         greenp0=gfbethe(omega,izeta0,d)
         greenp=gfbethe(omega,izeta,d)
         selfp=V**2/(ialpha)

!     Kinetic energy
         Ekin=Ekin+2.d0*(izeta*greenp-izeta0*greenp0)

!     Hybridization energy
         Ehyb=Ehyb+4.d0*(selfp*greenp)
      enddo

      kag=-real(greenp)*omega**2
      kag0=-real(greenp0)*omega**2
      kag0=abs(kag0)
      kag=abs(kag)
      coda1=(kag0-kag)*pi+
     1     2.d0*(kag*atan(multi*L/kag)-kag0*atan(multi*L/kag0))
      coda1=temp*coda1
      print*,'code =',coda1

      znp=znp-1.d0
      znd=znd-1.d0
      zntot=znp+znd

      Ekin=temp*2.d0*Ekin+2.d0*free+2.d0*coda1
      Ehyb=temp*2.d0*Ehyb
      Eepsi=ed0*znd+ep0*znp
      Emu=-xmu*zntot
      Etot=Ekin+Ehyb+Eepsi+Emu
      
      open(30,file='Etot.analytic',access='append')
      open(31,file='Ekin.analytic',access='append')
      open(32,file='Epot.analytic',access='append')
      open(33,file='Ehyb.analytic',access='append')
      open(34,file='Emu.analytic',access='append')
      open(35,file='Eepsi.analytic',access='append')
      write(30,*)temp,real(Etot)
      write(31,*)temp,real(Ekin)
      write(32,*)temp,real(Epot)
      write(33,*)temp,real(Ehyb)
      write(34,*)temp,real(Emu)
      write(35,*)temp,real(Eepsi)
      do i=30,35
         close(i)
      enddo

