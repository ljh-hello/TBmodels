! DESCRIPTION
!   Solve the non-interacting BHZ model with AFM 2x2 basis 
!   generate the hamiltoniana matrix H(k), 

program bhz_afm
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                       :: L=2048,Norb=2,Nspin=2,Nlat=2,Nso=Nspin*Norb,Nlso=Nlat*Nso
  integer                                 :: Nk,Nktot,Nkpath,Nkx,Npts
  integer                                 :: i,j,k,ik,iorb,jorb
  integer                                 :: ix,iy,iz
  real(8)                                 :: kx,ky,kz
  real(8),dimension(:),allocatable        :: kxgrid
  real(8),dimension(:,:),allocatable      :: kpath,ktrims,ddk
  complex(8),dimension(:,:,:),allocatable :: Hk
  real(8),dimension(:),allocatable        :: Wtk
  real(8)                                 :: mh,rh,lambda,delta,wmax
  real(8)                                 :: xmu,beta,eps,Ekin,Eloc
  real(8),dimension(L)                    :: wm,wr
  real(8)                                 :: n(Nlso)
  complex(8)                              :: w
  complex(8)                              :: Hloc(Nlso,Nlso)
  complex(8)                              :: Hlocsite(Nlat,Nso,Nso)
  complex(8)                              :: Gmats(Nlso,Nlso,L),Greal(Nlso,Nlso,L),Sfoo(Nlso,Nlso,L)
  character(len=20)                       :: file,nkstring
  logical                                 :: iexist,ibool,dcflag
  integer                                 :: icount,jcount,isite,jsite,ispin,jspin,row,col
  !Dirac matrices:
  complex(8),dimension(4,4)               :: Gamma1,Gamma2,Gamma3,Gamma4,Gamma5

  call parse_input_variable(nkx,"NKX","inputBHZ.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputBHZ.conf",default=500)
  call parse_input_variable(mh,"MH","inputBHZ.conf",default=3.d0)
  call parse_input_variable(rh,"RH","inputBHZ.conf",default=0.d0)
  call parse_input_variable(lambda,"LAMBDA","inputBHZ.conf",default=0.3d0)
  call parse_input_variable(delta,"DELTA","inputBHZ.conf",default=0.d0)
  call parse_input_variable(xmu,"XMU","inputBHZ.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputBHZ.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputBHZ.conf",default=1000.d0)
  call parse_input_variable(wmax,"WMAX","inputBHZ.conf",default=10d0)
  call parse_input_variable(file,"FILE","inputBHZ.conf",default="hkfile_bhz_afm.in")
  call save_input_file("inputBHZ.conf")

  Gamma1 = kron_pauli(pauli_z,pauli_x)
  Gamma2 =-kron_pauli(pauli_0,pauli_y)
  Gamma3 = kron_pauli(pauli_x,pauli_x)
  Gamma4 = kron_pauli(pauli_y,pauli_x)
  Gamma5 = kron_pauli(pauli_0,pauli_z)

  Nktot=Nkx*Nkx
  allocate(Hk(Nlso,Nlso,Nktot))
  allocate(Wtk(Nktot))
  allocate(kxgrid(Nkx))
  write(*,*) "Using Nk_total="//txtfy(Nktot)

  Hloc  = zero
  kxgrid = kgrid(Nkx)
  Hk = TB_build_model(hk_model,Nlso,kxgrid,kxgrid,[0d0])
  Wtk = 1d0/Nktot
  Hloc=sum(Hk(:,:,:),dim=3)/Nktot


  !Build the local GF:
  Gmats=zero
  Greal=zero
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-wmax,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(eps,"eps")
  Sfoo =zero
  call dmft_gloc_matsubara(Hk,Wtk,Gmats,Sfoo,iprint=1)
  Sfoo =zero
  call dmft_gloc_realaxis(Hk,Wtk,Greal,Sfoo,iprint=1)
  do iorb=1,Nlso
     n(iorb) = fft_get_density(Gmats(iorb,iorb,:),beta)
  enddo
  !plot observables
  open(10,file="observables.nint")
  write(10,"(20F20.12)")(n(iorb),iorb=1,Nlso),sum(n)
  close(10)
  write(*,"(A,20F14.9)")"Occupations =",(n(iorb),iorb=1,Nlso),sum(n)



  !solve along the standard path in the 2D BZ.
  Npts=4
  allocate(kpath(Npts,3))
  kpath(1,:)=kpoint_Gamma
  kpath(2,:)=kpoint_X1
  kpath(3,:)=kpoint_M1
  kpath(4,:)=kpoint_Gamma
  call solve_Hk_along_BZpath(Hk_model,Nlso,kpath,Nkpath,&
       colors_name=[red1,blue1,red1,blue1, red1,blue1,red1,blue1],&
       points_name=[character(len=20) :: 'G', 'X', 'M', 'G'],&
       file="Eigenbands_afm.nint")



contains



  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)          :: kpoint
    integer                       :: N
    real(8)                       :: kx,ky
    complex(8),dimension(N,N)     :: hk
    complex(8),dimension(N,N)     :: h0,tk
    complex(8),dimension(Nso,Nso) :: M
    complex(8),dimension(Nso,Nso) :: tx,ty,thx,thy
    !
    if(N/=Nlso)stop "hk_model error: N != Nlso" 
    kx = kpoint(1)
    ky = kpoint(2)
    !
    !
    M  = Mh*Gamma5
    tx = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma1
    thx= -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma1
    !
    ty = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma2
    thy= -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma2
    !
    ! H2 =  | m1                       tx + tx^+.e^i.2.kx + ty^+.e^i.(kx+ky) + ty^+.e^i.(kx-ky) |
    !       | tx^+ + tx.e^-i.2.kx + ty.e^-i.(kx+ky)+ ty^+.e^-i.(kx-ky)          m2              |
    !
    hk(1:4,1:4)    = M
    hk(1:4,5:8)    = tx  + thx*exp(xi*2*kx) + thy*exp(xi*(kx+ky)) + ty*exp(xi*(kx-ky))
    !
    hk(5:8,1:4)    = thx + tx*exp(-xi*2*kx) + ty*exp(-xi*(kx+ky)) + thy*exp(-xi*(kx-ky))
    hk(5:8,5:8)    = M
    !
  end function hk_model



end program bhz_afm







