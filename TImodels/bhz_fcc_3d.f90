program bhz_3d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                       :: L=2048,Norb=2,Nspin=2,Nso=Nspin*Norb
  integer                                 :: Nk,Nktot,Nkpath,Nkx,Npts,z2(4)
  integer                                 :: i,j,k,ik,iorb,jorb
  integer                                 :: ix,iy,iz
  real(8)                                 :: kx,ky,kz
  real(8),dimension(:),allocatable        :: kxgrid
  real(8),dimension(:,:),allocatable      :: kpath,ktrims
  complex(8),dimension(:,:,:),allocatable :: Hk
  real(8)                                 :: ez,mh,rh,lambda,delta,lz
  real(8)                                 :: xmu,beta,eps,Ekin,Eloc
  real(8),dimension(L)                    :: wm,wr
  real(8)                                 :: n(Nso)
  complex(8)                              :: w,Hloc(Nso,Nso)
  complex(8)                              :: Gmats(Nso,Nso,L),Greal(Nso,Nso,L)
  character(len=20)                       :: file,nkstring
  logical                                 :: iexist,ibool,dcflag,iener


  call parse_input_variable(nkx,"NKX","inputBHZ.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputBHZ.conf",default=500)
  call parse_input_variable(ez,"ez","inputBHZ.conf",default=1d0)
  call parse_input_variable(lz,"lz","inputBHZ.conf",default=1d0)
  call parse_input_variable(mh,"MH","inputBHZ.conf",default=3.d0)
  call parse_input_variable(lambda,"LAMBDA","inputBHZ.conf",default=0.3d0)
  call parse_input_variable(xmu,"XMU","inputBHZ.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputBHZ.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputBHZ.conf",default=1000.d0)
  call parse_input_variable(iener,"IENER","inputBHZ.conf",default=.false.)
  call parse_input_variable(file,"FILE","inputBHZ.conf",default="hkfile_bhz.in")
  call save_input_file("inputBHZ.conf")


  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
  Nktot=Nkx*Nkx*Nkx
  allocate(Hk(Nso,Nso,Nktot))
  allocate(kxgrid(Nkx))
  write(*,*) "Using Nk_total="//txtfy(Nktot)
  Greal = zero
  Gmats = zero 
  Hloc  = zero
  kxgrid = kgrid(Nkx)
  Hk = build_hk_model(hk_model,Nso,kxgrid,kxgrid,kxgrid)
  call write_hk_w90(trim(file),Nso,&
       Nd=Norb,&
       Np=1,   &
       Nineq=1,&
       hk=Hk,  &
       kxgrid=kxgrid,&
       kygrid=kxgrid,&
       kzgrid=kxgrid)
  Hloc=sum(Hk(:,:,:),dim=3)/Nktot
  call write_Hloc(Hloc)



  !Evaluate the Z2 index:
  !STRONG TI
  z2(1) = z2_number(reshape( [ [0,0,0] , [1,0,0] , [1,1,0] , [0,1,0] , [0,1,1] , [0,0,1] , [1,0,1] , [1,1,1] ] , [3,8] )*pi)
  !WEAK TI
  !K=1: n_1=1, n_2,3=0,1
  z2(2) = z2_number(reshape( [ [1,0,0] , [1,1,0] , [1,1,1] , [1,0,1] ] , [3,4])*pi)
  !K=2: n_2=1, n_1,2=0,1
  z2(3) = z2_number(reshape( [ [0,1,0] , [0,1,1] , [1,1,1] , [1,1,0] ] , [3,4])*pi)
  !k=3: n_3=1, n_1,2=0,1
  z2(4) = z2_number(reshape( [ [0,0,1] , [0,1,1] , [1,1,1] , [1,0,1] ] , [3,4])*pi)
  open(100,file="z2_invariant.nint")
  write(100,*)z2
  close(100)


  !solve along the standard path in the 2D BZ.
  Npts=8
  allocate(kpath(Npts,3))
  ! kpath(1,:)=[0,0,0]!G
  ! kpath(2,:)=[1,0,0]!X
  ! kpath(3,:)=[1,1,0]!M
  ! kpath(4,:)=[0,0,0]!G
  ! kpath(5,:)=[1,1,1]!R
  ! kpath(6,:)=[1,0,1]!M
  ! kpath(7,:)=[0,0,1]!Z
  ! kpath(8,:)=[1,1,1]!R
  ! kpath(9,:)=[0,0,0]!G
  ! kpath(10,:)=[0,0,1]!Z
  kpath(1,:)=[0,1,0]!X
  kpath(2,:)=[0,0,0]!G
  kpath(3,:)=[1,1,0]!M
  kpath(4,:)=[1,1,1]!R
  kpath(5,:)=[0,0,1]!Z
  kpath(6,:)=[1,0,1]!A
  kpath(7,:)=[0,0,0]!G
  kpath(8,:)=[0,0,1]!Z
  kpath=kpath*pi
  call solve_Hk_along_BZpath(Hk_model,Nso,kpath,Nkpath,&
       colors_name=[character(len=20) :: 'red','blue','red','blue'],&
                                !points_name=[character(len=20) :: "G","X","M","G","R","M","Z","R","G","Z"],&
       points_name=[character(len=20) :: "X","G","M","R","Z","A","G","Z"],&
       file="Eigenband.nint")








  !Build the local GF:
  wm = pi/beta*real(2*arange(1,L)-1,8)
  wr = linspace(-10.d0,10.d0,L)
  call start_timer()
  do ik=1,Nktot
     do i=1,L
        w = dcmplx(wr(i),eps)+xmu
        Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(w,Hk(:,:,ik))/Nktot
        w = xi*wm(i)+xmu
        Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k(w,Hk(:,:,ik))/Nktot
     enddo
     call eta(ik,Nktot)
  enddo
  call stop_timer()
  do iorb=1,Nso
     call splot("Gloc_l"//reg(txtfy(iorb))//reg(txtfy(iorb))//"_iw.nint",wm,Gmats(iorb,iorb,:))
     call splot("Gloc_l"//reg(txtfy(iorb))//reg(txtfy(iorb))//"_realw.nint",wr,&
          -dimag(Greal(iorb,iorb,:))/pi,dreal(Greal(iorb,iorb,:)))
     n(iorb) = fft_get_density(Gmats(iorb,iorb,:),beta)
  enddo
  open(10,file="observables.nint")
  write(10,"(10F20.12)")(n(iorb),iorb=1,Nso),sum(n)
  close(10)
  write(*,"(A,10F14.9)")"Occupations =",(n(iorb),iorb=1,Nso),sum(n)



  if(iener)then
     Ekin=get_kinetic_energy(Hk,L)
     Eloc=get_local_energy(Hk,L)
     open(10,file="energy.nint")
     write(10,"(3F20.12)")Ekin,Eloc,Ekin-Eloc
     close(10)
     write(*,"(A,F14.9)")"<K>           =",Ekin
     write(*,"(A,F14.9)")"<E0>          =",Eloc
  endif
  deallocate(Kpath,kxgrid,Hk)





contains


  function z2_number(ktrims) result(z2)
    real(8),dimension(:,:),intent(in)       :: ktrims
    complex(8),dimension(:,:,:),allocatable :: Htrims
    complex(8),dimension(:),allocatable     :: Delta
    integer                                 :: z2
    integer                                 :: i,j,Ntrim,itrim,Nocc
    !
    Ntrim=size(Ktrims,2)
    allocate(Htrims(Nso,Nso,Ntrim))
    allocate(Delta(Ntrim))
    !
    do itrim=1,Ntrim
       Htrims(:,:,itrim) = hk_model(Ktrims(:,itrim),Nso)
       Delta(itrim)=-sign(1d0,dreal(Htrims(1,1,itrim)))
    enddo
    !
    z2=product(Delta(:))
    if(z2>0)then
       z2=0
    else
       z2=1
    end if
    !
  end function z2_number



  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: kx,ky,kz
    complex(8),dimension(N,N) :: hk
    if(N/=4)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    kz=kpoint(3)
    !
    Hk          = zero
    Hk(1:2,1:2) = &
         (Mh - cos(kx) - cos(ky) - ez*cos(kz))*pauli_tau_z +&
         lambda*sin(kx)*pauli_tau_x + lambda*sin(ky)*pauli_tau_y
    Hk(3:4,3:4) = conjg( &
         (Mh-cos(-kx) - cos(-ky) - ez*cos(-kz))*pauli_tau_z +&
         lambda*sin(-kx)*pauli_tau_x + lambda*sin(-ky)*pauli_tau_y)
    Hk(1:2,3:4) = lz*lambda*sin(kz)*pauli_tau_x
    Hk(3:4,1:2) = lz*lambda*sin(kz)*pauli_tau_x
  end function hk_model

  function inverse_g0k(iw,hk) result(g0k)
    complex(8)                  :: iw
    complex(8),dimension(4,4)   :: hk
    complex(8),dimension(4,4)   :: g0k
    g0k = (iw + xmu)*zeye(4)-hk
    call inv(g0k)
  end function inverse_g0k


  function get_kinetic_energy(Hk,Liw) result(ed_Ekin)
    integer                                  :: Lk,No,Liw
    integer                                  :: i,ik,iorb
    complex(8),dimension(:,:,:)              :: Hk ![Nso][Nso][Nk]
    real(8),dimension(size(Hk,3))            :: Wtk
    !
    real(8),dimension(:),allocatable         :: wm
    complex(8),dimension(:,:),allocatable    :: Ak,Bk
    complex(8),dimension(:,:),allocatable    :: Ck,Zk
    complex(8),dimension(:,:),allocatable    :: Zeta,Gk,Tk
    real(8)                                  :: Tail0,Tail1,spin_degeneracy
    !
    real(8)                                  :: H0,ed_Ekin
    !
    No = size(Hk,1)
    Lk = size(Hk,3)
    if(No/=size(Hk,2))stop "get_kinetic_energy: size(Hk,1)!=size(Hk,2) [Norb_total]"
    if(Lk/=size(Wtk))stop "get_kinetic_energy: size(Wtk)!=size(Hk,3) [L_k]"
    !
    allocate(wm(Liw))
    allocate(Ak(No,No),Bk(No,No),Ck(No,No),Zk(No,No),Zeta(No,No),Gk(No,No),Tk(No,No))
    !
    wm = pi/beta*dble(2*arange(1,Liw)-1)
    Wtk=1d0/Lk
    !
    H0=0d0
    Zk=0d0 ; forall(i=1:No)Zk(i,i)=1d0
    do ik=1,Lk
       Bk=-Hk(:,:,ik)
       do i=1,Liw
          Gk = (xi*wm(i)+xmu)*Zk(:,:) - Hk(:,:,ik)
          select case(No)
          case default
             call matrix_inverse(Gk)
          case(1)
             Gk = 1d0/Gk
          end select
          Tk = Zk(:,:)/(xi*wm(i)) - Bk(:,:)/(xi*wm(i))**2
          Ck = matmul(Hk(:,:,ik),Gk - Tk)
          H0 = H0 + Wtk(ik)*trace_matrix(Ck,No)
       enddo
    enddo
    spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0=H0/beta*2.d0*spin_degeneracy
    !
    Tail0=0d0
    Tail1=0d0
    do ik=1,Lk
       Ak= Hk(:,:,ik)
       Bk=-Hk(:,:,ik)
       Ck= matmul(Ak,Bk)
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak,No)
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck,No)
    enddo
    Tail0=spin_degeneracy*Tail0
    Tail1=spin_degeneracy*Tail1*beta
    ed_Ekin=H0+Tail0+Tail1
    deallocate(wm,Ak,Bk,Ck,Zk,Zeta,Gk,Tk)
  end function get_kinetic_energy


  function get_local_energy(Hk,Liw) result(ed_Eloc)
    integer                                  :: Lk,No,Liw
    integer                                  :: i,ik,iorb
    complex(8),dimension(:,:,:)              :: Hk ![Nso][Nso][Nk]
    real(8),dimension(size(Hk,3))            :: Wtk
    !
    real(8),dimension(:),allocatable         :: wm
    complex(8),dimension(:,:),allocatable    :: Ak,Bk,Hloc
    complex(8),dimension(:,:),allocatable    :: Ck,Zk
    complex(8),dimension(:,:),allocatable    :: Zeta,Gk,Tk
    real(8)                                  :: Tail0,Tail1,spin_degeneracy
    !
    real(8)                                  :: H0,ed_Eloc
    !
    No = size(Hk,1)
    Lk = size(Hk,3)
    if(No/=size(Hk,2))stop "get_kinetic_energy: size(Hk,1)!=size(Hk,2) [Norb_total]"
    if(Lk/=size(Wtk))stop "get_kinetic_energy: size(Wtk)!=size(Hk,3) [L_k]"
    !
    allocate(wm(Liw))
    allocate(Ak(No,No),Bk(No,No),Ck(No,No),Zk(No,No),Zeta(No,No),Gk(No,No),Tk(No,No),Hloc(No,No))
    !
    wm = pi/beta*dble(2*arange(1,Liw)-1)
    Wtk=1d0/Lk
    !
    Hloc=sum(Hk,3)/Lk
    H0=0d0
    Zk=0d0 ; forall(i=1:No)Zk(i,i)=1d0
    do ik=1,Lk
       Bk=-Hk(:,:,ik)
       do i=1,Liw
          Gk = (xi*wm(i)+xmu)*Zk(:,:) - Hk(:,:,ik)
          select case(No)
          case default
             call matrix_inverse(Gk)
          case(1)
             Gk = 1d0/Gk
          end select
          Tk = Zk(:,:)/(xi*wm(i)) - Bk(:,:)/(xi*wm(i))**2
          Ck = matmul(Hloc,Gk - Tk)
          H0 = H0 + Wtk(ik)*trace_matrix(Ck,No)
       enddo
    enddo
    spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0=H0/beta*2.d0*spin_degeneracy
    !
    Tail0=0d0
    Tail1=0d0
    do ik=1,Lk
       Ak= Hloc
       Bk=-Hk(:,:,ik)
       Ck= matmul(Ak,Bk)
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak,No)
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck,No)
    enddo
    Tail0=spin_degeneracy*Tail0
    Tail1=spin_degeneracy*Tail1*beta
    ed_Eloc=H0+Tail0+Tail1
    deallocate(wm,Ak,Bk,Ck,Zk,Zeta,Gk,Tk,Hloc)
  end function get_local_energy


  function trace_matrix(M,dim) result(tr)
    integer                       :: dim
    complex(8),dimension(dim,dim) :: M
    complex(8) :: tr
    integer                       :: i
    tr=dcmplx(0d0,0d0)
    do i=1,dim
       tr=tr+M(i,i)
    enddo
  end function trace_matrix


end program bhz_3d


