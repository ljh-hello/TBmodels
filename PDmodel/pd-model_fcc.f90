! NAME
!   tdd-pam
!
! DESCRIPTION
!   Solve the non-interacting t_dd-PAM model and generate the hamiltoniana matrix H(k), 
!   Model parameters are: t_dd,t_pp,v0/tpd,ep0,ed0.
!   Using the 2d square lattice dispersion. 
!   The matrix is written in the Wannier90 form, as expected by w2CTQMC & ggED code.
!   Output is on file *hkfile.in default.
!
! MODEL Hamiltonian is:
!  |ed0 - 2*tdd*[cos(kx)+cos(ky)],   v0 |
!  |v0 ,   ep0 - 2*tpp*[cos(kx)+cos(ky)]|
!
! OPTIONS
!   nkx=[100]  - Number of points along the x component of the k-grid
!   tpp=[1.0]  - Hopping amplitude between p-electrons
!   tpd=[0.0]  - local hybridization between p- and d-electrons
!   tdd=[0.0]  - Hopping amplitude between d-electrons
!   ep0=[0.0]  - local energy of the p-band
!   ed0=[0.0]  - local energy of the d-band
!   u=[0.0]    - local interaction (required for Double Counting)
!   beta=[100] - Inverse temperature 
!   file=[hkfile.in]- output unit file,default is no write
!   dcflag=[T] - output unit file,default is no write

program tdd_pam
  USE COMMON_VARS
  USE PARSE_CMD
  USE TOOLS
  USE IOTOOLS
  USE FUNCTIONS
  USE TIMER
  USE FFTGF
  USE MATRIX
  implicit none

  integer,parameter    :: L=2000,Norb=2,Ltau=200
  integer              :: i,j,k,iorb,jorb,ik
  real(8)              :: tpp,v0,tpd,tdd,ed0,ep0,u,delta,xmu,beta,eps
  integer              :: dcshift,count
  real(8)              :: epsik,ep,em,fmesh,xmu0,Hdd,Hpp
  real(8)              :: n11,n22
  integer              :: Nkx,Nk
  real(8)              :: ix,iy
  real(8)              :: kx,ky
  real(8),dimension(L) :: wm,wr
  complex(8)           :: Hk(Norb,Norb),fg(L,Norb,Norb),fgr(L,Norb,Norb),w
  complex(8)           :: Uh(2,2),iUh(2,2),Hd(2,2)
  character(len=20)    :: file,nkstring
  logical              :: iexist,ibool,dcflag
  real(8)              :: gzerop,gzerom,gzero
  real(8)              :: Nfk(2,2),Nkk(2,2)
  complex(8),dimension(2,2,L)   :: fgk
  real(8),dimension(2,2,0:Ltau) :: fgkt

  namelist/hkvars/nkx,tpp,tpd,tdd,ep0,ed0,u,xmu,beta,eps,file,dcflag

  nkx=200
  tpp=0.25d0
  tpd=0.4d0
  tdd=0.d0
  ep0=0.d0
  ed0=0.d0
  u=0.d0
  xmu=0.d0
  eps=1.d-3
  beta=100.d0
  file="hkfile.in"
  dcflag=.true.

  inquire(file="inputPDHAM.in",exist=iexist)
  if(iexist)then
     open(10,file="inputPDHAM.in")
     read(10,nml=hkvars)
     close(10)
  endif

  call parse_cmd_variable(nkx,"NKX")
  call parse_cmd_variable(tpp,"TPP")
  call parse_cmd_variable(tpd,"TPD")
  call parse_cmd_variable(tdd,"TDD")
  call parse_cmd_variable(ep0,"EP0")
  call parse_cmd_variable(ed0,"ED0")
  call parse_cmd_variable(u,"U")
  call parse_cmd_variable(xmu,"XMU")
  call parse_cmd_variable(eps,"EPS")
  call parse_cmd_variable(beta,"BETA")
  call parse_cmd_variable(file,"FILE")
  call parse_cmd_variable(dcflag,"DCFLAG")

  dcshift=0 ; if(dcflag)dcshift=1

  write(*,nml=hkvars)

  wm = pi/beta*real(2*arange(1,L)-1,8)
  wr = linspace(-10.d0,10.d0,L,mesh=fmesh)

  Nk=Nkx*Nkx
  print*,"Using Nk=",Nk
  open(50,file=trim(file))
  write(nkstring,*)Nk
  write(50,*)trim(adjustl(trim(Nkstring)))," 1 1 1 1"

  xmu0=xmu
  delta=ep0-ed0
  gzerop=0.5d0*(ep0+ed0+sqrt((Delta)**2 + 4.d0*tpd**2))
  gzerom=0.5d0*(ep0+ed0-sqrt((Delta)**2 + 4.d0*tpd**2))
  gzero=0.d0
  if(delta<0.d0)gzero=gzerop
  if(delta>0.d0)gzero=gzerom
  if(delta /= 0.d0)xmu=xmu+gzero
  print*,xmu,xmu0

  open(101,file="Momentum_distribution.pd")
  fgr=zero ; fg =zero ;ep=0.d0 ; em=0.d0 ;count=0; hdd=0.d0; hpp=0.d0
  call start_progress
  do ix=1,Nkx
     kx = -pi + 2.d0*pi*real(ix-1,8)/real(Nkx,8)
     do iy=1,Nkx
        ky = -pi + 2.d0*pi*real(iy-1,8)/real(Nkx,8)
        epsik   = cos(kx)+cos(ky)
        Hk(1,1) = ed0 - 2.d0*tdd*epsik
        Hk(2,2) = ep0 - 2.d0*tpp*epsik + dble(dcshift)*U/2.d0
        Hk(1,2) = tpd
        Hk(2,1) = tpd
        write(50,"(3(F10.7,1x))")kx,ky,pi
        do iorb=1,Norb
           write(50,"(10(2F10.7,1x))")(Hk(iorb,jorb),jorb=1,Norb)
        enddo

        Hk(2,2) = Hk(2,2)-dble(dcshift)*U/2.d0
        do i=1,L
           w = cmplx(wr(i),eps,8)+xmu
           fgr(i,:,:)=fgr(i,:,:) + inverse_g0k(w,Hk)
           w = xi*wm(i)+xmu
           fg(i,:,:) =fg(i,:,:)  + inverse_g0k(w,Hk)
           fgk(:,:,i) = inverse_g0k(w,Hk)
        enddo


        !The cheapest way to get the momentum distribution, avoiding 
        !weird FFT is to:
        !1) build the diagonal matrix of FD ditribution of the eigen-energies:
        nfk = 0.d0
        nfk(1,1)=fermi(eplus(hk),beta)
        nfk(2,2)=fermi(eminus(hk),beta)
        !2) rotate the matrix back to the initial (non-diagonal) basis:
        Uh = matrix_uk(hk) ; iUh= transpose(Uh)
        nkk = matmul(iUh,matmul(nfk,Uh))        
        write(101,"(4F20.12)")-2.d0*tpp*epsik,nkk(1,1),nkk(2,2),nkk(1,2)

        ! !Get momentum-distribution from G(iw)
        ! do i=1,2
        !    do j=1,2
        !       call fftgf_iw2tau(fgk(i,j,:),fgkt(i,j,0:),beta)
        !       nkk(i,j)=-fgkt(i,j,Ltau)
        !    enddo
        ! enddo
        ! nkk(1,2)=nkk(1,2)-0.5d0
        ! nkk(2,1)=nkk(2,1)-0.5d0
        ! write(100,"(4F20.12)")-2.d0*tpp*epsik,nkk(1,1),nkk(2,2),nkk(1,2)

        ep=ep+eplus(hk)
        em=em+eminus(hk)
        hdd=hdd+Hk(1,1)
        hpp=hpp+Hk(2,2)

        ! !print Hamiltonian values:
        ! write(200,"(4F20.12)")-2.d0*tpp*epsik,dreal(hk(1,1)),dreal(hk(2,2)),dreal(hk(1,2))
        ! !Print Eigenvalues:
        ! hd=0.d0 ; hd(1,1)=eplus(hk) ; hd(2,2)=eminus(hk)
        ! write(201,"(4F20.12)")-2.d0*tpp*epsik,dreal(hd(1,1)),dreal(hd(2,2)),dreal(hd(1,2))
        ! !Diagonalize Hamiltonian and print diagonal terms, must be equal to 201
        ! Uh = matrix_uk(hk) ; iUh= transpose(Uh)
        ! hd = matmul(matmul(Uh,hk),iUh)
        ! write(202,"(4F20.12)")-2.d0*tpp*epsik,dreal(hd(1,1)),dreal(hd(2,2)),dreal(hd(1,2))
        count=count+1
        call progress_bar(count,Nk)
     enddo
  enddo
  close(101)
  call stop_progress
  fgr= fgr/real(Nk,8)
  fg = fg/real(Nk,8)
  ep = ep/real(Nk,8)
  em = em/real(Nk,8)
  hdd= hdd/real(Nk,8)
  hpp= hpp/real(Nk,8)

  open(10,file="pam_shift")
  rewind(10)
  write(10,"(F20.12)")gzero
  close(10)

  ik = 0
  open(10,file="Eigenbands.pd")
  !From \Gamma=(0,0) to X=(pi,0): 100 steps
  do ix=1,100
     ik=ik+1
     kx = 0.d0 + pi*real(ix-1,8)/100.d0
     ky = 0.d0
     epsik   = cos(kx)+cos(ky)
     Hk(1,1) = ed0 - 2.d0*tdd*epsik
     Hk(2,2) = ep0 - 2.d0*tpp*epsik
     Hk(1,2) = tpd
     Hk(2,1) = tpd
     write(10,*)ik,eplus(hk),eminus(hk)
  enddo
  !From X=(pi,0) to M=(pi,pi): 100 steps
  do iy=1,100
     ik=ik+1
     kx = pi
     ky = 0.d0 + pi*real(iy-1,8)/100.d0
     epsik   = cos(kx)+cos(ky)
     Hk(1,1) = ed0 - 2.d0*tdd*epsik
     Hk(2,2) = ep0 - 2.d0*tpp*epsik
     Hk(1,2) = tpd
     Hk(2,1) = tpd
     write(10,*)ik,eplus(hk),eminus(hk)
  enddo
  !From M=(pi,pi) to \Gamma=(0,0): 100 steps
  do ix=1,100
     ik=ik+1
     iy=ix
     kx = pi - pi*real(ix-1,8)/100.d0
     ky = pi - pi*real(iy-1,8)/100.d0
     epsik   = cos(kx)+cos(ky)
     Hk(1,1) = ed0 - 2.d0*tdd*epsik
     Hk(2,2) = ep0 - 2.d0*tpp*epsik
     Hk(1,2) = tpd
     Hk(2,1) = tpd
     write(10,*)ik,eplus(hk),eminus(hk)
  enddo
  close(10)

  call splot("DOSdd.pd",wr,-dimag(fgr(:,1,1))/pi)
  call splot("DOSpp.pd",wr,-dimag(fgr(:,2,2))/pi)
  call splot("Gdd_iw.pd",wm,fg(:,1,1))
  call splot("Gpp_iw.pd",wm,fg(:,2,2))

  n11 = -2.d0*sum(dimag(fgr(:,1,1))*fermi(wr(:),beta))*fmesh/pi
  n22 = -2.d0*sum(dimag(fgr(:,2,2))*fermi(wr(:),beta))*fmesh/pi
  open(10,file="observables.pd")
  write(10,"(14F20.12)")tpp,tpd,tdd,xmu,u,n11,n22,n11+n22,&
       ep,em,abs(ep-em),gzerop,gzerom,abs(gzerop-gzerom)
  close(10)

  write(*,"(A,2F14.9)")"Center of mass of the bands   =",hdd,hpp
  write(*,"(A,2F14.9)")"Center of mass of hyb bands   =",ep,em
  write(*,"(A,2F14.9)")"Zero_+, Zero_-                =",gzerop,gzerom
  write(*,"(A,3F14.9)")"Occupations                   =",n11,n22,n11+n22
  if(dcflag)then
     write(*,"(A,2F14.9)")"U / DC shift of p-electrons   =",U,U/2.d0
     write(*,"(A,2F14.9)")"p-level before/after DC shift =",ep0,ep0+dble(dcshift)*U/2.d0
  end if


  call msg("Remember to open the file:"//trim(file))
  call msg("Parameters for "//trim(file)//" are in +paramaters4_"//trim(file))
  !reset xmu to original value
  xmu=xmu0
  open(10,file="parameters4_"//trim(file))
  write(10,nml=hkvars)
  close(10)
  close(50)




contains

  function inverse_g0k(iw,hk) result(g0k)
    integer                     :: i
    complex(8),dimension(2,2)   :: hk
    complex(8)                  :: iw
    complex(8),dimension(2,2)   :: g0k
    complex(8)                  :: delta,ppi,vmix
    g0k=zero
    delta = iw - hk(1,1)
    ppi   = iw - hk(2,2)
    vmix  = -hk(1,2)
    g0k(1,1) = one/(delta - abs(vmix)**2/ppi)
    g0k(2,2) = one/(ppi - abs(vmix)**2/delta)
    g0k(1,2) = -vmix/(ppi*delta - abs(vmix)**2)
    g0k(2,1) = conjg(g0k(1,2))
  end function inverse_g0k

  function eplus(hk)
    complex(8),dimension(2,2) :: hk
    real(8)                   :: eplus
    eplus = hk(1,1)+hk(2,2) + sqrt(abs(hk(1,1)-hk(2,2))**2 + 4.d0*hk(1,2)*hk(2,1) )
    eplus = eplus/2.d0
  end function eplus

  function eminus(hk)
    complex(8),dimension(2,2) :: hk
    real(8)                   :: eminus
    eminus = hk(1,1)+hk(2,2) - sqrt(abs(hk(1,1)-hk(2,2))**2 + 4.d0*hk(1,2)*hk(2,1) )
    eminus = eminus/2.d0
  end function eminus


  function matrix_uk(hk) result(uk)
    complex(8),dimension(2,2) :: hk
    complex(8),dimension(2,2) :: uk
    complex(8) :: delta0,ek,u,v
    delta0 = -(hk(1,1)-hk(2,2))
    ek     = sqrt( delta0**2 + 4.d0*hk(1,2)*hk(2,1) )
    u = sqrt(0.5d0*(1.d0 - delta0/ek))
    v = sqrt(0.5d0*(1.d0 + delta0/ek))
    uk(1,1) = u
    uk(2,2) = u
    uk(1,2) = v
    uk(2,1) =-v
  end function matrix_uk

end program tdd_pam


