! DESCRIPTION
!   Solve the non-interacting BHZ model and generate the hamiltoniana matrix H(k), 
!   Model parameters are: t_dd,v0/tpd,ep0,ed0.
!   Using the 2d square lattice dispersion. 
!   The matrix is written in the Wannier90 form, as expected by w2CTQMC & ggED code.
!   Output is on file *hkfile.in default.
!
! MODEL Hamiltonian is:
!
! |     h^{2x2}(k)              &         hso^{2x2}(k)        |
! |      [hso^{2x2}]*(k)        &        [h^{2x2}]*(-k)       |
!
!
! h^{2x2}(k):=
!
! | m-(Cos{kx}+Cos{ky})         & \lambda*(Sin{kx}-i*Sin{ky}) |
! | \lambda*(Sin{kx}+i*Sin{ky}) & -m+(Cos{kx}+Cos{ky})        |
!
! hso^{2x2}(k):=
! | xi*rh*(sin(kx)-xi*sin(ky))  &         \delta              |
! |         -\delta             &             0               |
program bhz_fcc
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                       :: L=2048,Norb=2,Nspin=2,Nso=Nspin*Norb
  integer                                 :: Nk,Nktot,Nkpath,Nkx,Npts
  integer                                 :: i,j,k,ik,iorb,jorb
  integer                                 :: ix,iy,iz
  real(8)                                 :: kx,ky,kz
  real(8),dimension(:),allocatable        :: kxgrid
  real(8),dimension(:,:),allocatable      :: kpath
  complex(8),dimension(:,:,:),allocatable :: Hk

  real(8)                                 :: mh,rh,lambda,delta
  real(8)                                 :: xmu,beta,eps,Ekin,Eloc
  real(8),dimension(L)                    :: wm,wr
  real(8)                                 :: n(Nso)
  complex(8)                              :: w,Hloc(Nso,Nso)
  complex(8)                              :: Gmats(Nso,Nso,L),Greal(Nso,Nso,L)
  character(len=20)                       :: file,nkstring
  logical                                 :: iexist,ibool,dcflag,iener


  call parse_input_variable(nkx,"NKX","inputBHZ.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputBHZ.conf",default=500)
  call parse_input_variable(mh,"MH","inputBHZ.conf",default=3.d0)
  call parse_input_variable(rh,"RH","inputBHZ.conf",default=0.d0)
  call parse_input_variable(lambda,"LAMBDA","inputBHZ.conf",default=0.3d0)
  call parse_input_variable(delta,"DELTA","inputBHZ.conf",default=0.d0)
  call parse_input_variable(xmu,"XMU","inputBHZ.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputBHZ.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputBHZ.conf",default=1000.d0)
  call parse_input_variable(iener,"IENER","inputBHZ.conf",default=.false.)
  call parse_input_variable(file,"FILE","inputBHZ.conf",default="hkfile_bhz.in")
  call save_input_file("inputBHZ.conf")


  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
  Nktot=Nkx*Nkx
  allocate(Hk(Nso,Nso,Nktot))
  allocate(kxgrid(Nkx))
  write(*,*) "Using Nk_total="//txtfy(Nktot)
  Greal = zero
  Gmats = zero 
  Hloc  = zero
  kxgrid = kgrid(Nkx)
  Hk = build_hk_model(hk_model,Nso,kxgrid,kxgrid,[0d0])
  call write_hk_w90(trim(file),Nso,&
       Nd=Norb,&
       Np=1,   &
       Nineq=1,&
       hk=Hk,  &
       kxgrid=kxgrid,&
       kygrid=kxgrid,&
       kzgrid=[0d0])
  Hloc=sum(Hk(:,:,:),dim=3)/Nktot
  call write_Hloc(Hloc)

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
  do iorb=1,Nso
     call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_iw.nint",wm,Gmats(iorb,iorb,:))
     call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_realw.nint",wr,&
          -dimag(Greal(iorb,iorb,:))/pi,dreal(Greal(iorb,iorb,:)))
     n(iorb) = fft_get_density(Gmats(iorb,iorb,:),beta)
  enddo


  !solve along the standard path in the 2D BZ.
  Npts=5
  allocate(kpath(Npts,3))
  kpath(1,:)=kpoint_X1
  kpath(2,:)=kpoint_Gamma
  kpath(3,:)=kpoint_M1
  kpath(4,:)=kpoint_X1
  kpath(5,:)=kpoint_Gamma
  call solve_Hk_along_BZpath(Hk_model,Nso,kpath,Nkpath,&
       colors_name=[character(len=20) :: 'red','blue','red','blue'],&
       points_name=[character(len=20) :: 'X', 'G', 'M', 'X', 'G'],&
       file="Eigenband.nint")
  !plot observables

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
  ! function hk_model(kpoint,N) result(hk)
  !   real(8),dimension(:)      :: kpoint
  !   integer                   :: N
  !   real(8)                   :: kx,ky
  !   complex(8),dimension(N,N) :: hk
  !   if(size(kpoint)/=3)stop "hk_model: error in kpoint dimensions"
  !   if(N/=4)stop "hk_model: error in N dimensions"
  !   kx=kpoint(1)
  !   ky=kpoint(2)
  !   Hk          = zero
  !   Hk(1:2,1:2) = hk_bhz2x2(kx,ky)
  !   Hk(3:4,3:4) = conjg(hk_bhz2x2(-kx,-ky))
  !   Hk(1,4) = -delta ; Hk(4,1)=-delta
  !   Hk(2,3) =  delta ; Hk(3,2)= delta
  !   Hk(1,3) = xi*rh*(sin(kx)-xi*sin(ky))
  !   Hk(3,1) =-xi*rh*(sin(kx)+xi*sin(ky))
  ! end function hk_model
  ! function hk_bhz2x2(kx,ky) result(hk)
  !   real(8)                   :: kx,ky
  !   complex(8),dimension(2,2) :: hk
  !   real(8)                   :: epsik
  !   epsik   = cos(kx)+cos(ky)
  !   hk(1,1) = mh - epsik
  !   hk(2,2) =-conjg(hk(1,1))!-mh + epsik
  !   hk(1,2) = lambda*(sin(kx)-xi*sin(ky))
  !   hk(2,1) = lambda*(sin(kx)+xi*sin(ky))
  ! end function hk_bhz2x2


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

  ! function Eigk(hk) result(eig)
  !   complex(8),dimension(4,4) :: hk
  !   real(8),dimension(4)      :: eig
  !   !call matrix_diagonalize(hk,eig)
  !   ! call matrix_diagonalize(hk(1:2,1:2),eig(1:2))
  !   ! call matrix_diagonalize(hk(3:4,3:4),eig(3:4))
  !   eig(1:2) = eigk2x2(hk(1:2,1:2))
  !   eig(3:4) = eigk2x2(hk(3:4,3:4))
  ! end function Eigk

  ! function Eigk2x2(hk) result(eig)
  !   complex(8),dimension(2,2) :: hk
  !   real(8),dimension(2)      :: eig
  !   complex(8)                :: delta0,vk2,sqrt0,u,v
  !   delta0 = (hk(1,1)-hk(2,2))
  !   vk2    = 4.d0*hk(1,2)*hk(2,1)
  !   sqrt0  = sqrt( abs(delta0)**2 + vk2 )
  !   eig(1)=hk(1,1)+hk(2,2) + sqrt0
  !   eig(2)=hk(1,1)+hk(2,2) - sqrt0
  !   eig = eig/2.d0
  !   if(sqrt0==0.d0)then
  !      u=sqrt(0.5d0)
  !      v=u
  !   else
  !      u = sqrt(0.5d0*(1.d0 + delta0/sqrt0))
  !      v = sqrt(0.5d0*(1.d0 - delta0/sqrt0))
  !   endif
  !   hk(1,1) = u
  !   hk(2,1) = -v
  !   hk(1,2) = v
  !   hk(2,2) = u
  ! end function Eigk2x2


  ! function matrix_uk(hk) result(uk)
  !   complex(8),dimension(2,2) :: hk
  !   complex(8),dimension(2,2) :: uk
  !   complex(8) :: delta0,ek,u,v
  !   delta0 = -(hk(1,1)-hk(2,2))
  !   ek     = sqrt( delta0**2 + 4.d0*hk(1,2)*hk(2,1) )
  !   u = sqrt(0.5d0*(1.d0 - delta0/ek))
  !   v = sqrt(0.5d0*(1.d0 + delta0/ek))
  !   uk(1,1) = u
  !   uk(1,2) = v
  !   uk(2,1) =-v
  !   uk(2,2) = u
  ! end function matrix_uk



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


end program bhz_fcc


