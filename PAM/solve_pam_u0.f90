program get_pam_info
  USE COMMON_VARS
  USE PARSE_CMD
  USE TOOLS
  USE IOTOOLS
  USE FFTGF
  implicit none


  integer                             :: i,L
  real(8)                             :: tpp,tpd,tdd,ep0,ed0,delta,teff,xmu,beta,wmax,eps,dw,de
  real(8)                             :: gzero,gzerop,gzerom,nd,np,ntot,ee,nread,ndelta,nerror,dind,deltaE
  real(8),dimension(:),allocatable    :: wm,wr,en,gdt,gpt,eplus,eminus,epsk
  complex(8)                          :: iw
  complex(8),dimension(:),allocatable :: zeta,alpha,gp,gd,sigp
  logical :: pflag
  namelist/pamvar/L,tpp,tpd,tdd,ep0,ed0,delta,xmu,beta,wmax

  call parse_cmd_variable(L,"L",default=2*8192)
  call parse_cmd_variable(pflag,"PFLAG",default=.false.)
  call parse_cmd_variable(xmu,"XMU",default=0.d0)
  call parse_cmd_variable(beta,"BETA",default=5000.d0)
  call parse_cmd_variable(wmax,"WMAX",default=10.d0)
  call parse_cmd_variable(eps,"EPS",default=1.d-6)
  call parse_cmd_variable(nread,"NREAD",default=0.d0)
  call parse_cmd_variable(nerror,"NERROR",default=1.d-6)
  call parse_cmd_variable(ndelta,"NDELTA",default=0.1d0)
  call parse_cmd_variable(tpp,"TPP",default=0.5d0)
  call parse_cmd_variable(tpd,"TPD",default=1.d0)
  call parse_cmd_variable(tdd,"TDD",default=0.d0)
  call parse_cmd_variable(ep0,"EP0",default=0.d0)
  call parse_cmd_variable(ed0,"ED0",default=0.d0)
  call parse_cmd_variable(xmu,"XMU",default=0.d0)
  delta=ep0-ed0
  write(*,nml=pamvar)

  teff = tpd**2/sqrt(delta**2+4.d0*tpd**2)
  gzerop=0.5d0*(ep0+ed0+sqrt((ed0-ep0)**2 + 4*tpd**2))
  gzerom=0.5d0*(ep0+ed0-sqrt((ed0-ep0)**2 + 4*tpd**2))
  deltaE=gzerop-gzerom
  write(*,"(A,2F14.9)")"Zero_+, Zero_- =",gzerop,gzerom
  write(*,"(A,2F14.9)")"t_eff, W_eff   =",teff/2.d0,teff
  if(delta < 0.d0)gzero=gzerop
  if(delta > 0.d0)gzero=gzerom
  if(delta /= 0.d0)xmu=xmu+gzero

  if(.not.pflag)stop

  allocate(wm(L),wr(L))
  wm = pi/beta*real(2*arange(1,L)-1,8)
  wr = linspace(-wmax,wmax,L,mesh=dw)

  allocate(alpha(L),zeta(L),gp(L),gd(L),sigp(L))
  allocate(gdt(0:L),gpt(0:L))

  if(nread/=0.d0)call search_mu()

  !Get Matsubara axis solution:
  alpha = xi*wm + xmu - ed0
  sigp  = tpd**2/alpha
  zeta  = xi*wm + xmu - ep0 - sigp
  gp    = gfbethe(wm,zeta,2.d0*tpp)
  gd    = one/alpha + (tpd/alpha)**2*gp
  call fftgf_iw2tau(gd,gdt,beta)
  call fftgf_iw2tau(gp,gpt,beta)
  nd=-2.d0*gdt(L)
  np=-2.d0*gpt(L)
  ntot=nd+np

  call splot("Sigmapp_iw.u0",wm,sigp)
  call splot("Gpp_iw.u0",wm,gp)
  call splot("Gdd_iw.u0",wm,gd)



  !Get Real-axis solution:
  alpha=cmplx(wr,eps,8)+xmu-ed0
  sigp =tpd**2/alpha
  zeta =cmplx(wr,eps,8)+xmu-ep0-sigp
  gp   =gfbether(wr,zeta,2.d0*tpp)
  gd   =one/alpha + (tpd/alpha)**2*gp
  call splot("Sigmapp_realw.u0",wr,sigp)
  call splot("Gpp_realw.u0",wr,gp)
  call splot("Gdd_realw.u0",wr,gd)
  call splot("DOSdd.u0",wr,-dimag(gd)/pi)
  call splot("DOSpp.u0",wr,-dimag(gp)/pi)

  allocate(en(L),eplus(L),eminus(L))
  en = linspace(-tpp,tpp,L,mesh=de)
  Eplus=0.5d0*(ed0+ep0+en+sqrt(4.d0*tpd**2+(ep0+en-ed0)**2))
  Eminus=0.5d0*(ed0+ep0+en-sqrt(4.d0*tpd**2+(ep0+en-ed0)**2))
  call splot("HybBands.u0",en,eplus,eminus)
  open(10,file="mu.u0")
  write(10,*)-tpp,xmu
  write(10,*)tpp,xmu
  close(10)
  dind=abs(minval(Eplus)-maxval(Eminus))
  open(10,file='columns.u0')
  write(10,*)"1u, 2nd, 3np, 4ntot, 5delta0, 6delta, 7delta_ind, 8xmu"
  close(10)
  call splot("observables.u0",0.d0,nd,np,ntot,delta,deltaE,dind,xmu)
  open(10,file="Poles.u0")
  do i=1,L
     if(abs(wr(i))<=1.d0)write(10,*)wr(i),wr(i)+xmu-ep0-dreal(sigp(i))
  enddo
  close(10)

contains

  subroutine search_mu()
    real(8) :: ntmp
    integer :: nindex,nindex1,iter
    real(8) :: ndelta1
    logical :: converged
    nindex=0
    iter=0
    converged=.false.
    do while(.not.converged)
       iter=iter+1
       nindex1=nindex
       ndelta1=ndelta
       !
       alpha = xi*wm + xmu - ed0
       sigp  = tpd**2/alpha
       zeta  = xi*wm + xmu - ep0 - sigp
       gp    = gfbethe(wm,zeta,2.d0*tpp)
       gd    = one/alpha + (tpd/alpha)**2*gp
       call fftgf_iw2tau(gd,gdt,beta)
       call fftgf_iw2tau(gp,gpt,beta)
       nd=-2.d0*gdt(L)
       np=-2.d0*gpt(L)
       ntmp=nd+np
       !
       if((ntmp >= nread+nerror))then
          nindex=-1
       elseif(ntmp <= nread-nerror)then
          nindex=1
       else
          nindex=0
       endif
       if(nindex1+nindex==0)then !avoid loop forth and back
          ndelta=ndelta1/exp(1.d0) !decreasing the step       
       else
          ndelta=ndelta1
       endif
       xmu=xmu+real(nindex,8)*ndelta
       write(*,"(A,I3,A,f9.6,A,f9.6,A,f15.12,A,f15.12,L)")"Iter=",iter,"| n=",ntmp," /",nread,&
            "| shift=",nindex*ndelta,"| xmu=",xmu,converged

       converged=(abs(ntmp-nread)<=nerror).OR.(iter>128)
    enddo
  end subroutine search_mu

end program get_pam_info
