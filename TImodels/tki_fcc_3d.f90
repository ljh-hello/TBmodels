! DESCRIPTION
!   Solve the non-interacting TKI model and generate the hamiltoniana matrix H(k), 
!   Model parameters are: ed,ef,lambda,ed0.
!
! MODEL Hamiltonian is:
!
! |     E_d(k)              &         [lambda*V(k)]^+  |
! |    [lambda*V(k)]        &        E_f(k)            |
!

! E_d(k) = ed0 - ed*[Cos(kx)+Cos(ky)]\sigma_0
! E_f(k) = -ef*[Cos(kx)+Cos(ky)]\sigma_0
! V(k)  = [2Sin(kx),2Sin(ky),0].dot.[\sigma_x,\sigma_y,\sigma_z]
include "TIGHT_BINDING.f90"
program tki
  USE CONSTANTS
  USE ARRAYS
  USE MATRIX
  USE IOTOOLS
  USE TIMER
  USE DMFT_TOOLS
  USE COLORS
  USE TIGHT_BINDING
  implicit none

  integer,parameter                       :: L=1000,Norb=4
  integer                                 :: Nktot,Nkpath,Nkx,Nky,Nkz
  integer                                 :: i,j,k,ik,iorb,jorb
  integer                                 :: ix,iy,iz
  real(8)                                 :: kx,ky,kz
  real(8),dimension(:),allocatable        :: kxgrid,kygrid,kzgrid
  real(8),dimension(:,:),allocatable      :: kpath
  complex(8),dimension(:,:,:),allocatable :: Hk
  !model parameters:
  real(8)                                 :: ed,ef,ed0,lambda
  real(8)                                 :: xmu,beta,eps

  real(8),dimension(L)                    :: wm,wr
  real(8)                                 :: fmesh,Ndens(Norb)
  complex(8)                              :: w
  complex(8)                              :: Hloc(Norb,Norb)
  complex(8)                              :: Gmats(Norb,Norb,L),Greal(Norb,Norb,L)
  character(len=20)                       :: file,nkstring
  logical                                 :: iexist,ibool,dcflag



  call parse_input_variable(Nkx,"NKX","input.conf",default=10)
  call parse_input_variable(Nky,"NKY","input.conf",default=10)
  call parse_input_variable(Nkz,"NKZ","input.conf",default=10)
  call parse_input_variable(Nkpath,"NKPATH","input.conf",default=500)
  call parse_input_variable(ed,"ed","input.conf",default=1.d0)
  call parse_input_variable(ef,"ef","input.conf",default=0.2d0)
  call parse_input_variable(ed0,"ed0","input.conf",default=0.d0)
  call parse_input_variable(lambda,"LAMBDA","input.conf",default=0.4d0)
  call parse_input_variable(xmu,"XMU","input.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","input.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","input.conf",default=1000.d0)
  call parse_input_variable(file,"FILE","input.conf",default="hkfile_tki.in")

  Nktot=Nkx*Nky*Nkz
  allocate(Hk(Nktot,Norb,Norb))
  allocate(kxgrid(Nkx),kygrid(Nky),kzgrid(Nkz))

  write(*,*) "Using Nk_total="//txtfy(Nktot)
  Greal = zero
  Gmats = zero 
  Hloc  = zero
  kxgrid = kgrid(Nkx)
  kygrid = kgrid(Nky)
  kzgrid = kgrid(Nkz)

  call start_timer

  Hk = build_hk_model(Nktot,Norb,hk_model,kxgrid,kygrid,kzgrid)
  call write_hk_w90(trim(file),Norb,2,1,Hk,kxgrid,kygrid,kzgrid)
  Hloc=sum(Hk(:,:,:),dim=1)/Nktot
  call write_hloc(Hloc,"tki_hloc.dat")

  wm = pi/beta*real(2*arange(1,L)-1,8)
  wr = linspace(-10d0,10d0,L,mesh=fmesh)
  do ik=1,Nktot
     do i=1,L
        w = dcmplx(wr(i),eps)+xmu
        Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(w,Hk(ik,:,:))/Nktot
        w = xi*wm(i)+xmu
        Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k(w,Hk(ik,:,:))/Nktot
     enddo
     call eta(ik,Nktot)
  enddo
  call stop_timer


  do iorb=1,Norb
     Ndens(iorb)=get_density_fromFFT(Gmats(iorb,iorb,:),beta)
     call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_iw.nint",wm,Gmats(iorb,iorb,:))
     call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_realw.nint",wr,&
          -dimag(Greal(iorb,iorb,:)),dreal(Greal(iorb,iorb,:)))
  enddo


  Nktot=6*Nkpath
  allocate(kpath(8,3))
  kpath(1,:)=[0,0,0]!G
  kpath(2,:)=[1,0,0]!X
  kpath(3,:)=[1,1,0]!M
  kpath(4,:)=[0,0,0]!G
  kpath(5,:)=[1,1,1]!R
  kpath(6,:)=[1,0,1]!M
  kpath(7,:)=[0,0,1]!X
  kpath(8,:)=[1,1,1]!R
  kpath=kpath*pi
  call solve_along_BZpath(Nkpath,kpath,Norb,Hk_model,&
       colors_name=[character(len=20) :: 'red','red','blue','blue'],&
       points_name=[character(len=20) :: "G","X","M","G","R","M","X","R"],&
       file="Eigenband.nint")


  open(10,file="observables.nint")
  write(10,"(10F20.12)")(Ndens(iorb),iorb=1,Norb),sum(Ndens)
  close(10)
  write(*,"(A,10F14.9)")"Occupations                   =",(Ndens(iorb),iorb=1,Norb),sum(Ndens)

contains


  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: kx,ky,kz
    real(8)                   :: ck
    real(8)                   :: cx,cy,cz,sx,sy,sz
    complex(8),dimension(N,N) :: hk
    complex(8),dimension(2,2) :: tau0,taux,tauy,tauz
    if(size(kpoint)/=3)stop "hk_model: error in kpoint dimensions"
    if(N/=4)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    kz=kpoint(3)
    !Pauli matrices
    tau0=zero;tau0(1,1)=one;tau0(2,2)=one
    taux=zero;taux(1,2)=one;taux(2,1)=one
    tauy=zero;tauy(1,2)=-xi;tauy(2,1)=xi
    tauz=zero;tauz(1,1)=one;tauz(2,2)=-one
    !
    cx=cos(kx);cy=cos(ky);cz=cos(kz)
    sx=sin(kx);sy=sin(ky);sz=sin(kz)
    ck=cx+cy+cz
    !
    Hk          = zero
    Hk(1:2,1:2) = (ed0+ed*ck)*tau0
    Hk(3:4,3:4) = -ef*ck*tau0      !this must have inverted sign to get a TI state
    !Hk(1:2,3:4) = -2d0*lambda*(sx*taux + sy*tauy + sz*tauz)
    Hk(1:2,3:4) =-2d0*lambda*(sx*transpose(conjg(taux)) + sy*transpose(conjg(tauy)) + sz*transpose(conjg(tauz)))
    Hk(3:4,1:2) = -2d0*lambda*(sx*taux + sy*tauy + sz*tauz)
  end function hk_model


  function inverse_g0k(iw,hk) result(g0k)
    complex(8)                  :: iw
    complex(8),dimension(4,4)   :: hk
    complex(8),dimension(4,4)   :: g0k
    g0k = -hk
    forall(i=1:4)g0k(i,i) = iw + xmu + g0k(i,i)
    call matrix_inverse(g0k)
  end function inverse_g0k



  function get_density_fromFFT(giw,beta) result(n)
    complex(8),dimension(:) :: giw
    real(8)                 :: gtau(size(giw))
    real(8)                 :: beta,n
    call fftgf_iw2tau(giw,gtau,beta,[0.d0,1.d0,-xmu,0d0])
    n = -2.d0*gtau(size(giw))
  end function get_density_fromFFT

end program tki


