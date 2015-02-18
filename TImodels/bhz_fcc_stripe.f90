program bhz_fcc_stripe
  USE SF_CONSTANTS
  USE SF_ARRAYS
  USE SF_LINALG
  USE SF_IOTOOLS
  USE SF_TIMER
  USE DMFT_MISC
  USE DMFT_TIGHT_BINDING
  USE DMFT_PARSE_INPUT
  USE DMFT_FFTGF
  USE DMFT_FFTAUX
  implicit none

  integer,parameter                       :: L=2048,Norb=2,Nso=2*Norb
  integer                                 :: Nk,Nktot,Nkpath,Nkx,Ly,Npts,Nslat
  integer                                 :: i,j,k,ik,iorb,jorb
  integer                                 :: ix,iy,iz
  real(8)                                 :: kx,ky,kz
  real(8),dimension(:),allocatable        :: kxgrid
  real(8),dimension(:,:),allocatable      :: kpath
  complex(8),dimension(:,:,:),allocatable :: Hk,Hkr

  real(8)                                 :: mh,rh,lambda,delta
  real(8)                                 :: xmu,beta,eps
  real(8),dimension(L)                    :: wm,wr
  real(8)                                 :: n(Nso)
  complex(8)                              :: w,Hloc(Nso,Nso)
  complex(8)                              :: Gmats(Nso,Nso,L),Greal(Nso,Nso,L)
  character(len=20)                       :: file,nkstring
  logical                                 :: iexist,ibool,dcflag


  call parse_input_variable(nkx,"NKX","inputBHZ.conf",default=25)
  call parse_input_variable(Ly,"Ly","inputBHZ.conf",default=10)
  call parse_input_variable(nkpath,"NKPATH","inputBHZ.conf",default=500)
  call parse_input_variable(mh,"MH","inputBHZ.conf",default=3.d0)
  call parse_input_variable(rh,"RH","inputBHZ.conf",default=0.d0)
  call parse_input_variable(lambda,"LAMBDA","inputBHZ.conf",default=0.3d0)
  call parse_input_variable(delta,"DELTA","inputBHZ.conf",default=0.d0)
  call parse_input_variable(xmu,"XMU","inputBHZ.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputBHZ.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputBHZ.conf",default=100.d0)
  call parse_input_variable(file,"FILE","inputBHZ.conf",default="hkfile_bhz.in")
  call save_input_file("inputBHZ.conf")

  !START SOLVING THE INHOMO PROBLEM: EDGE STATES IN THE STRIPE GEOMETRY.
  !Generate the Hamiltonian for the full BZ (a set of 1D lines)
  Nslat=Ly*Nso
  allocate(kxgrid(Nkx))
  allocate(Hkr(Nslat,Nslat,Nkx))
  Hkr = build_Hkr_model(hkr_model,Ly,Nso,kxgrid,[0d0],[0d0],pbc=.false.)
  call write_hk_w90("Hkrfile_BHZ.data",&
       No=Nslat,&
       Nd=Norb,&
       Np=0,&
       Nineq=Ly,&
       Hk=Hkr,&
       kxgrid=kxgrid,kygrid=[0d0],kzgrid=[0d0])



  Npts=3
  allocate(Kpath(Npts,1))
  kpath(1,:)=[-1]*pi
  kpath(2,:)=[ 0]*pi
  kpath(3,:)=[ 1]*pi
  call solve_HkR_along_BZpath(hkr_model,Ly,Nso,kpath,Nkpath,"hkr_Eigenbands.nint",.false.)



contains



  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: kx,ky
    complex(8),dimension(N,N) :: hk
    if(N/=4)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    Hk          = zero
    Hk(1:2,1:2) = hk_bhz2x2(kx,ky)
    Hk(3:4,3:4) = conjg(hk_bhz2x2(-kx,-ky))
  end function hk_model
  function hk_bhz2x2(kx,ky) result(hk)
    real(8)                   :: kx,ky
    complex(8),dimension(2,2) :: hk
    hk = (mh-cos(kx)-cos(ky))*pauli_tau_z + lambda*sin(kx)*pauli_tau_x + lambda*sin(ky)*pauli_tau_y 
  end function hk_bhz2x2




  function hkr_model(kpoint,Nlat,N,pbc) result(Hrk)
    real(8),dimension(:)                :: kpoint
    real(8)                             :: kx
    integer                             :: Nlat,N
    complex(8),dimension(Nlat*N,Nlat*N) :: Hrk
    integer                             :: i,Idmin,Idmax,Itmin,Itmax
    logical                             :: pbc
    kx=kpoint(1)
    Hrk=zero
    do i=1,Nlat
       Idmin=1+(i-1)*N
       Idmax=i*N
       Hrk(Idmin:Idmax,Idmin:Idmax)=h0_rk_bhz(kx,N)
    enddo
    do i=1,Nlat-1
       Idmin=1+(i-1)*N
       Idmax=i*N
       Itmin=i*N+1
       Itmax=(i+1)*N
       Hrk(Idmin:Idmax,Itmin:Itmax)=t0_rk_bhz(N)
       Hrk(Itmin:Itmax,Idmin:Idmax)=conjg(transpose(t0_rk_bhz(N)))
    enddo
    if(pbc)then
       Itmin=(Nlat-1)*N+1
       Itmax= Nlat*N
       Hrk(1:N,Itmin:Itmax)=conjg(transpose(t0_rk_bhz(N)))!T0
       Hrk(Itmin:Itmax,1:N)=t0_rk_bhz(N)!conjg(transpose(T0))
    endif
  end function hkr_model


  function h0_rk_bhz(kx,N) result(H)
    real(8)                    :: kx
    integer                    :: N
    complex(8),dimension(N,N)  :: H
    complex(8),dimension(4,4) :: gamma1,gamma5
    gamma1=kronecker_product_pauli_matrices(pauli_sigma_x,pauli_sigma_z)
    gamma5=kronecker_product_pauli_matrices(pauli_sigma_z,pauli_sigma_0)
    H = (mh-cos(kx))*gamma5 + lambda*sin(kx)*gamma1
  end function h0_rk_bhz

  function t0_rk_bhz(N) result(H)
    integer                    :: N
    complex(8),dimension(N,N)  :: H
    complex(8),dimension(4,4) :: gamma2,gamma5
    gamma2=kronecker_product_pauli_matrices(-pauli_sigma_y,pauli_sigma_0)
    gamma5=kronecker_product_pauli_matrices(pauli_sigma_z,pauli_sigma_0)
    H = xi*0.5d0*lambda*gamma2 + 0.5d0*gamma5
  end function T0_rk_bhz



  function inverse_g0k(iw,hk) result(g0k)
    complex(8)                  :: iw
    complex(8),dimension(4,4)   :: hk
    complex(8),dimension(4,4)   :: g0k
    g0k=zero
    if(rh==0.AND.delta==0)then
       g0k(1:2,1:2) = inverse_g0k2x2(iw,hk(1:2,1:2))
       g0k(3:4,3:4) = inverse_g0k2x2(iw,hk(3:4,3:4))
    else
       g0k = -hk
       forall(i=1:4)g0k(i,i) = iw + xmu + g0k(i,i)
       call matrix_inverse(g0k)
    endif
  end function inverse_g0k
  !
  function inverse_g0k2x2(iw,hk) result(g0k)
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
  end function inverse_g0k2x2

end program bhz_fcc_stripe


