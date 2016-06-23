! DESCRIPTION
!   Solve the non-interacting BHZ model with AFM 2x2 basis 
!   generate the hamiltoniana matrix H(k), 

program bhz_afm
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                       :: L=2048,Norb=2,Nspin=2,Nlat=8,Nso=Nspin*Norb,Nlso=Nlat*Nso
  integer                                 :: Nk,Nktot,Nkpath,Nkx,Npts
  integer                                 :: i,j,k,ik,iorb,jorb
  real(8),dimension(:),allocatable        :: kxgrid
  real(8),dimension(:,:),allocatable      :: kpath,ktrims,ddk
  complex(8),dimension(:,:,:),allocatable :: Hk
  real(8)                                 :: mh,rh,lambda,delta
  real(8)                                 :: xmu,beta,eps,Ekin,Eloc
  real(8),dimension(L)                    :: wm,wr
  real(8)                                 :: n(Nlso)
  complex(8)                              :: w
  complex(8)                              :: Hloc(Nlso,Nlso)
  complex(8)                              :: Hlocsite(Nlat,Nso,Nso)
  complex(8)                              :: Gmats(Nlso,Nlso,L),Greal(Nlso,Nlso,L)
  character(len=20)                       :: file,nkstring
  logical                                 :: iexist,ibool,dcflag
  integer                                 :: icount,jcount,isite,jsite,ispin,jspin,row,col
  !Dirac matrices:
  complex(8),dimension(4,4)               :: Gamma1,Gamma2,Gamma3,Gamma4,Gamma5

  call parse_input_variable(nkx,"NKX","inputBHZ.conf",default=14)
  call parse_input_variable(nkpath,"NKPATH","inputBHZ.conf",default=500)
  call parse_input_variable(mh,"MH","inputBHZ.conf",default=3.d0)
  call parse_input_variable(rh,"RH","inputBHZ.conf",default=0.d0)
  call parse_input_variable(lambda,"LAMBDA","inputBHZ.conf",default=0.3d0)
  call parse_input_variable(delta,"DELTA","inputBHZ.conf",default=0.d0)
  call parse_input_variable(xmu,"XMU","inputBHZ.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputBHZ.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputBHZ.conf",default=1000.d0)
  call parse_input_variable(file,"FILE","inputBHZ.conf",default="hkfile_bhz_afm.in")
  call save_input_file("inputBHZ.conf")

  if(Norb/=2)stop "Norb != 2"
  if(Nspin/=2)stop "Nspin != 2"
  if(Nso/=4)stop "Nso != 4"
  if(Nlso/=32)stop "Nlso != 32"

  Gamma1 = kron_pauli(pauli_z,pauli_x)
  Gamma2 =-kron_pauli(pauli_0,pauli_y)
  Gamma3 = kron_pauli(pauli_x,pauli_x)
  Gamma4 = kron_pauli(pauli_y,pauli_x)
  Gamma5 = kron_pauli(pauli_0,pauli_z)

  Nktot=Nkx*Nkx*Nkx
  allocate(Hk(Nlso,Nlso,Nktot))
  allocate(kxgrid(Nkx))
  write(*,*) "Using Nk_total="//txtfy(Nktot)

  Greal = zero
  Gmats = zero 
  Hloc  = zero
  kxgrid = kgrid(Nkx)
  Hk = build_hk_model(hk_model,Nlso,kxgrid,kxgrid,kxgrid)
  call write_hk_w90(trim(file),Nlso,&
       Nd=Nso,&
       Np=0,   &
       Nineq=8,&
       hk=Hk,  &
       kxgrid=kxgrid,&
       kygrid=kxgrid,&
       kzgrid=kxgrid)
  Hloc=sum(Hk(:,:,:),dim=3)/Nktot




  !Build the local GF:
  wm = pi/beta*real(2*arange(1,L)-1,8)
  wr = linspace(-10.d0,10.d0,L)
  do ik=1,Nktot
     do i=1,L
        w = dcmplx(wr(i),eps)+xmu
        Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(w,Hk(:,:,ik))/Nktot
        w = xi*wm(i)+xmu
        Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k(w,Hk(:,:,ik))/Nktot
     enddo
  enddo
  do iorb=1,Nlso
     call splot("Gloc_lso"//reg(txtfy(iorb))//"_iw.nint",wm,Gmats(iorb,iorb,:))
     call splot("Gloc_lso"//reg(txtfy(iorb))//"_realw.nint",wr,&
          -dimag(Greal(iorb,iorb,:))/pi,dreal(Greal(iorb,iorb,:)))
     n(iorb) = fft_get_density(Gmats(iorb,iorb,:),beta)
  enddo
  !plot observables
  open(10,file="observables.nint")
  write(10,"(33F20.12)")(n(iorb),iorb=1,Nlso),sum(n)
  close(10)
  write(*,"(A,33F14.9)")"Occupations =",(n(iorb),iorb=1,Nlso),sum(n)



  !solve along the standard path in the 2D BZ.
  Npts=8
  allocate(kpath(Npts,3))
  kpath(1,:)=[0,1,0]!X
  kpath(2,:)=[0,0,0]!G
  kpath(3,:)=[1,1,0]!M
  kpath(4,:)=[1,1,1]!R
  kpath(5,:)=[0,0,1]!Z
  kpath(6,:)=[1,0,1]!A
  kpath(7,:)=[0,0,0]!G
  kpath(8,:)=[0,0,1]!Z
  kpath=kpath*pi
  call solve_Hk_along_BZpath(Hk_model,Nlso,kpath,Nkpath,&
       colors_name=[character(len=10) :: 'red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue'],&
       points_name=[character(len=20) :: "X","G","M","R","Z","A","G","Z"],&
       file="Eigenbands_afm.nint")



contains



  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)              :: kpoint
    integer                           :: N
    real(8)                           :: kx,ky,kz
    complex(8),dimension(N,N)         :: hk
    complex(8),dimension(4*Nso,4*Nso) :: h2d,hz,hzh
    complex(8),dimension(Nso,Nso)     :: M
    complex(8),dimension(Nso,Nso)     :: tx,ty,tz,thx,thy,thz
    !
    if(N/=Nlso)stop "hk_model error: N != Nlso" 
    kx = 2d0*kpoint(1)
    ky = 2d0*kpoint(2)
    kz = 2d0*kpoint(3)
    M  = Mh*Gamma5
    tx = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma1
    ty = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma2
    tz = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma3
    thx= -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma1
    thy= -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma2
    thz= -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma3
    !
    ! Hk =  | Hk^2d   Hk^z  |
    !       | Hk^z^+  Hk^2d |
    ! 
    ! Hk^2d =  | m1                 tx + tx^+.e^-ikx    0                   ty^+ + ty.e^iky |
    !          | tx^+ + tx.e^-ikx   m2                  ty^+ + ty.e^-iky    0               |
    !          | 0                  ty + ty^+e^-iky     m3                  tx^+ + tx.e^ikx |
    !          | ty + ty^+e^-iky    0                   tx + tx^+.e^-ikx    m4              |
    !
    ! Hk^z  =  | tz^+ + tz.e^-ikz   0                   0                   0               |
    !          | 0                  tz^+ + tz.e^-ikz    0                   0               |
    !          | 0                  0                   tz^+ + tz.e^-ikz    0               |
    !          | 0                  0                   0                   tz^+ + tz.e^-ikz|
    !
    h2d(1:4,1:4)    = M
    h2d(1:4,5:8)    = tx  + thx*exp(-xi*kx)
    h2d(1:4,9:12)   = zero
    h2d(1:4,13:16)  = thy + ty*exp(xi*ky)
    !
    h2d(5:8,1:4)    = thx + tx*exp(xi*kx)
    h2d(5:8,5:8)    = M
    h2d(5:8,9:12)   = thy + ty*exp(xi*ky) 
    h2d(5:8,13:16)  = zero
    !
    h2d(9:12,1:4)   = zero
    h2d(9:12,5:8)   = ty  + thy*exp(-xi*ky)
    h2d(9:12,9:12)  = M
    h2d(9:12,13:16) = thx + tx*exp(xi*kx)
    !
    h2d(13:16,1:4)  = ty  + thy*exp(-xi*ky)
    h2d(13:16,5:8)  = zero
    h2d(13:16,9:12) = tx  + thx*exp(-xi*kx)
    h2d(13:16,13:16)= M
    !
    Hz              = zero
    Hz(1:4,1:4)     = thz + tz*exp(-xi*kz)
    Hz(5:8,5:8)     = thz + tz*exp(-xi*kz)
    Hz(9:12,9:12)   = thz + tz*exp(-xi*kz)
    Hz(13:16,13:16) = thz + tz*exp(-xi*kz)
    !
    Hzh             = zero
    Hzh(1:4,1:4)    = tz + thz*exp(xi*kz)
    Hzh(5:8,5:8)    = tz + thz*exp(xi*kz)
    Hzh(9:12,9:12)  = tz + thz*exp(xi*kz)
    Hzh(13:16,13:16)= tz + thz*exp(xi*kz)
    !
    Hk(1:16,1:16)   = H2d
    Hk(1:16,17:32)  = Hz
    Hk(17:32,1:16)  = Hzh
    Hk(17:32,17:32) = H2d
  end function hk_model



  function inverse_g0k(zeta,hk) result(g0k)
    complex(8)                      :: zeta
    complex(8),dimension(Nlso,Nlso) :: hk
    complex(8),dimension(Nlso,Nlso) :: g0k
    g0k = -hk
    forall(i=1:Nlso)g0k(i,i)=zeta - hk(i,i)
    call matrix_inverse(g0k)
  end function inverse_g0k





end program bhz_afm







