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
!   nky=[100]  - Number of points along the y component of the k-grid
!   nkz=[1]    - Number of points along the z component of the k-grid (not used)
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
  use fftgf
  implicit none

  integer,parameter    :: L=2000,Norb=2,Ltau=200
  integer              :: i,j,k,iorb,jorb,ik
  real(8)              :: tpp,v0,tpd,tdd,ed0,ep0,u,delta,xmu,beta,eps
  integer              :: dcshift,count
  real(8)              :: epsik,ep,em,fmesh,xmu0,Dpp,Ddd,dos,de
  real(8)              :: n11,n22
  integer              :: Nkx,Nky,Nkz,Nk
  real(8)              :: ix,iy,iz
  real(8)              :: kx,ky,kz
  real(8),dimension(L) :: wm,wr
  complex(8)           :: Hk(Norb,Norb),fg(L,Norb,Norb),fgr(L,Norb,Norb),w
  character(len=20)    :: file,nkstring
  logical              :: iexist,ibool,dcflag
  real(8)              :: gzerop,gzerom,gzero,nkk(2,2)

  complex(8),dimension(2,2,L)   :: fgk
  real(8),dimension(2,2,0:Ltau) :: fgkt

  namelist/hkvars/nkx,nky,nkz,tpp,tpd,tdd,ep0,ed0,u,xmu,beta,eps,file,dcflag

  nkx=200
  nky=200
  nkz=1
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
  call parse_cmd_variable(nky,"NKY")
  call parse_cmd_variable(nkz,"NKZ")
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

  Nk=Nkx!*Nky*Nkz
  call msg("Using Nk="//txtfy(Nk))
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


  de = 2.d0/real(Nk,8)
  fgr=zero ; fg =zero ; ep=0.d0 ; em=0.d0 ; count=0
  call start_timer
  do ix=1,Nk
     kx = -1.d0 + real(ix-1,8)*de
     Hk(1,1) = ed0 - 2.d0*tdd*kx
     Hk(2,2) = ep0 - 2.d0*tpp*kx + dble(dcshift)*U/2.d0
     Hk(1,2) = tpd
     Hk(2,1) = Hk(1,2)
     write(50,"(3(F10.7,1x))")kx
     do iorb=1,Norb
        write(50,"(10(2F10.7,1x))")(Hk(iorb,jorb),jorb=1,Norb)
     enddo
     Hk(2,2) = Hk(2,2)-dble(dcshift)*U/2.d0
     dos = dens_bethe(kx,1.d0)*de
     do i=1,L
        w = cmplx(wr(i),eps,8)+xmu
        fgr(i,:,:)=fgr(i,:,:) + inverse_g0k(w,Hk)*dos
        w = xi*wm(i)+xmu
        fg(i,:,:) =fg(i,:,:)  + inverse_g0k(w,Hk)*dos
     enddo
     count=count+1
     call eta(count,Nk)
  enddo
  call stop_timer


  open(10,file="pam_shift")
  rewind(10)
  write(10,"(F20.12)")gzero
  close(10)

  call splot("DOSdd.pd",wr,-dimag(fgr(:,1,1))/pi)
  call splot("DOSpp.pd",wr,-dimag(fgr(:,2,2))/pi)
  call splot("Gdd_iw.pd",wm,fg(:,1,1))
  call splot("Gpp_iw.pd",wm,fg(:,2,2))

  n11 = -2.d0*sum(dimag(fgr(:,1,1))*fermi(wr(:),beta))*fmesh/pi
  n22 = -2.d0*sum(dimag(fgr(:,2,2))*fermi(wr(:),beta))*fmesh/pi
  open(10,file="observables.pd")
  write(10,"(14F20.12)")tpp,tpd,tdd,xmu,u,n11,n22,n11+n22,&
       gzerop,gzerom,abs(gzerop-gzerom)
  close(10)

  write(*,"(A,2F14.9)")"Zero_+, Zero_-                =",gzerop,gzerom
  if(dcflag)then
     write(*,"(A,2F14.9)")"U / DC shift of p-electrons   =",U,U/2.d0
     write(*,"(A,2F14.9)")"p-level before/after DC shift =",ep0,ep0+dble(dcshift)*U/2.d0
  end if


  call msg("Remember to open the file:"//trim(file))
  call msg("Parameters for "//trim(file)//" are in +paramaters4_"//trim(file))  
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
    eminus = hk(1,1)+hk(2,2) -sqrt(abs(hk(1,1)-hk(2,2))**2 + 4.d0*hk(1,2)*hk(2,1) )
    eminus = eminus/2.d0
  end function eminus

end program tdd_pam


