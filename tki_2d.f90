! DESCRIPTION
!   Solve the non-interacting TKI model and generate the hamiltoniana matrix H(k), 
!   Model parameters are: ed,ef,lambda,ed0.
!
! MODEL Hamiltonian is:
!
! |     E_d(k)              &         [lambda*V(k)]^+  |
! |    [lambda*V(k)]        &        E_f(k)            |
!

! E_d(k) =    - 2*td*[Cos(kx)+Cos(ky)]\sigma_0
! E_f(k) = Ef + 2*alpha*td*[Cos(kx)+Cos(ky)]\sigma_0
! V(k)  = lambda*[2Sin(kx),2Sin(ky),0].dot.[\sigma_x,\sigma_y,\sigma_z]
program tki_2d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                       :: L=1024,Norb=2,Nspin=2,Nso=Nspin*Norb
  integer                                 :: Nk,Nktot,Nkpath,Nkx,Npts
  integer                                 :: i,j,k,ik,iorb,jorb
  real(8)                                 :: kx,ky
  real(8),dimension(:),allocatable        :: kxgrid
  complex(8),dimension(:,:,:),allocatable :: Hk
  !model parameters:
  real(8)                                 :: ts,alpha,ef,lambda
  real(8)                                 :: xmu,beta,eps


  real(8)                                 :: chern,z2
  real(8),dimension(:,:,:),allocatable    :: nkgrid
  real(8),dimension(L)                    :: wm,wr
  real(8)                                 :: fmesh,Ndens(Nso)
  complex(8)                              :: w
  complex(8)                              :: Hloc(Nso,Nso)
  complex(8)                              :: Gmats(Nso,Nso,L),Greal(Nso,Nso,L)
  character(len=20)                       :: file,nkstring
  logical                                 :: iexist,ibool,dcflag
  real(8),dimension(:,:),allocatable      :: kpath,ktrims


  call parse_input_variable(Nkx,"NKX","input.conf",default=10)
  call parse_input_variable(Nkpath,"NKPATH","input.conf",default=300)
  call parse_input_variable(ts,"ts","input.conf",default=0.5d0)
  call parse_input_variable(alpha,"ALPHA","input.conf",default=0.2d0,comment="BHZ needs: alpha=1")
  call parse_input_variable(ef,"ef","input.conf",default=-4.d0,comment="BHZ needs: ef=-4*m")
  call parse_input_variable(lambda,"LAMBDA","input.conf",default=0.3d0)
  call parse_input_variable(xmu,"XMU","input.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","input.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","input.conf",default=1000.d0)
  call parse_input_variable(file,"FILE","input.conf",default="hkfile_tki.in")
  call save_input_file("input.conf")


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
  call write_Hloc(Hloc,"tki_hloc.dat")


  call start_timer
  wm = pi/beta*real(2*arange(1,L)-1,8)
  wr = linspace(-10d0,10d0,L,mesh=fmesh)
  do ik=1,Nktot
     do i=1,L
        w = dcmplx(wr(i),eps)+xmu
        Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(w,Hk(:,:,ik))/Nktot
        w = xi*wm(i)+xmu
        Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k(w,Hk(:,:,ik))/Nktot
     enddo
     call eta(ik,Nktot)
  enddo
  call stop_timer



  do iorb=1,Nso
     Ndens(iorb)=get_density_fromFFT(Gmats(iorb,iorb,:),beta)
     call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_iw.nint",wm,Gmats(iorb,iorb,:))
     call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_realw.nint",wr,&
          -dimag(Greal(iorb,iorb,:)),dreal(Greal(iorb,iorb,:)))
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
       colors_name=[character(len=20) :: 'red','red','blue','blue'],&
       points_name=[character(len=20) :: 'X', 'G', 'M', 'X', 'G'],&
       file="Eigenband.nint")


  open(10,file="density.nint")
  write(10,"(10F20.12)")(Ndens(iorb),iorb=1,Nso),sum(Ndens)
  close(10)
  write(*,"(A,10F14.9)")"Occupations =",(Ndens(iorb),iorb=1,Nso),sum(Ndens)



  allocate(ktrims(2,4))
  ktrims=reshape( [ [0d0,0d0] , [0d0,pi] , [pi,0d0] , [pi,pi] ] , shape(ktrims))
  z2 = z2_number(ktrims,[3,4])
  open(100,file="z2_invariant.nint")
  write(100,*) z2
  close(100)


contains


  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: kx,ky,kz
    real(8)                   :: tf
    complex(8),dimension(N,N) :: hk
    if(N/=4)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    kz=0d0
    tf=-alpha*ts !this must have inverted sign to get a TI state
    !
    Hk(1:2,1:2) = (0d0- 2*ts*(cos(kx)+cos(ky)))*pauli_0 !- xmu*pauli_0
    Hk(3:4,3:4) = (ef - 2*tf*(cos(kx)+cos(ky)))*pauli_0 !- xmu*pauli_0
    Hk(1:2,3:4) = -2d0*lambda*(sin(kx)*pauli_x + sin(ky)*pauli_y)
    Hk(3:4,1:2) = -2d0*lambda*(sin(kx)*pauli_x + sin(ky)*pauli_y)
    !
  end function hk_model


  function inverse_g0k(iw,hk) result(g0k)
    complex(8)                  :: iw
    complex(8),dimension(4,4)   :: hk
    complex(8),dimension(4,4)   :: g0k
    g0k = (iw+xmu)*eye(4)-hk 
    call inv(g0k)
  end function inverse_g0k


  function get_density_fromFFT(giw,beta) result(n)
    complex(8),dimension(:) :: giw
    real(8)                 :: gtau(size(giw))
    real(8)                 :: beta,n
    call fft_gf_iw2tau(giw,gtau,beta)!,[0.d0,1.d0,-xmu,0d0])
    n = -2.d0*gtau(size(giw))
  end function get_density_fromFFT


  function z2_number(ktrims,band_indices) result(z2)
    real(8),dimension(:,:),intent(in)       :: ktrims
    integer,dimension(:),intent(in)         :: band_indices
    complex(8),dimension(:,:,:),allocatable :: Htrims
    real(8),dimension(:,:),allocatable      :: Etrims
    complex(8),dimension(:),allocatable     :: Delta
    real(8)                                 :: z2
    integer                                 :: i,j,Ntrim,itrim,Nocc
    !
    Ntrim=size(Ktrims,2)
    Nocc = size(band_indices)
    allocate(Htrims(Nso,Nso,Ntrim),Etrims(Nocc,Ntrim))
    allocate(Delta(Ntrim))
    !
    do itrim=1,Ntrim
       Htrims(:,:,itrim) = hk_model(Ktrims(:,itrim),Nso)
       Delta(itrim)=-sign(1d0,dreal(Htrims(1,1,itrim)))
    enddo
    z2=product(Delta(:))
    if(z2>0)then
       z2=0d0
    else
       z2=1d0
    end if
  end function z2_number


end program tki_2d


