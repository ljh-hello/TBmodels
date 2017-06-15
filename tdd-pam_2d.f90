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
!  |ed0 - 2*tdd*[cos(kx)+cos(ky)],   tpd-4*v0*sin(kx)*sin(ky) |
!  |tpd-4*v0*sin(kx)*sin(ky)        ,   ep0 - 2*tpp*[cos(kx)+cos(ky)]|
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
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                       :: Norb=2,Nspin=1,Nso=Nspin*Norb
  integer                                 :: Nk,Nktot,Nkpath,Nkx,Npts,L
  integer                                 :: i,j,k,ik,iorb,jorb
  integer                                 :: ix,iy,iz
  real(8)                                 :: kx,ky,kz
  real(8),dimension(:,:),allocatable      :: kpath,ktrims
  complex(8),dimension(:,:,:),allocatable :: Hk
  real(8),dimension(:),allocatable        :: Wtk

  real(8)                                 :: tpp,alpha,tpd,v0,ep0,ed0
  real(8)                                 :: xmu,beta,eps,Eout(2)
  real(8)                                 :: xmu0,delta,gzerop,gzerom,gzero
  real(8)                                 :: dens(Nso)
  complex(8)                              :: Hloc(Nso,Nso)
  complex(8),dimension(:,:,:),allocatable :: Gmats,Greal,Sfoo
  character(len=20)                       :: file



  call parse_input_variable(nkx,"NKX","inputPAM.conf",default=200)
  call parse_input_variable(nkpath,"NKPATH","inputPAM.conf",default=500)
  call parse_input_variable(L,"L","inputPAM.conf",default=2048)
  call parse_input_variable(tpp,"TPP","inputPAM.conf",default=0.25d0)
  call parse_input_variable(tpd,"TPD","inputPAM.conf",default=0.4d0)
  call parse_input_variable(alpha,"ALPHA","inputPAM.conf",default=0d0)
  call parse_input_variable(v0 ,"V0","inputPAM.conf",default=0d0)
  call parse_input_variable(ep0,"EP0","inputPAM.conf",default=0d0)
  call parse_input_variable(ed0,"ED0","inputPAM.conf",default=0d0)
  call parse_input_variable(xmu,"XMU","inputPAM.conf",default=0d0)
  call parse_input_variable(eps,"EPS","inputPAM.conf",default=1d-3)
  call parse_input_variable(beta,"BETA","inputPAM.conf",default=100d0)
  call parse_input_variable(file,"FILE","inputPAM.conf",default="hkfile_pam.in")
  call save_input_file("inputPAM.conf")
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-5d0,"wini")
  call add_ctrl_var(5d0,"wfin")
  call add_ctrl_var(eps,"eps")

  Nktot=Nkx*Nkx
  write(*,*) "Using Nk_total="//txtfy(Nktot)

  xmu0=xmu
  delta=ep0-ed0
  gzerop=0.5d0*(ep0+ed0+sqrt((Delta)**2 + 4.d0*tpd**2))
  gzerom=0.5d0*(ep0+ed0-sqrt((Delta)**2 + 4.d0*tpd**2))
  gzero=0.d0
  if(delta<0.d0)gzero=gzerop
  if(delta>0.d0)gzero=gzerom
  if(delta /= 0.d0)xmu=xmu+gzero
  write(*,*)'shift mu to (from) = ',xmu,'(',xmu-gzero,')'
  write(*,*)'shift is           = ',gzero

  open(10,file="pam_shift.tb")
  rewind(10)
  write(10,"(F20.12)")gzero
  close(10)


  allocate(Hk(Nso,Nso,Nktot))
  allocate(Wtk(Nktot))
  call TB_set_bk(bkx=[pi2,0d0],bky=[0d0,pi2])

  call TB_build_model(Hk,hk_model,Nso,[Nkx,Nkx])
  Wtk = 1d0/Nktot

  
  call TB_write_hk(Hk,trim(file),Nso,&
       Nd=Norb,Np=1,Nineq=1,&
       Nkvec=[Nkx,Nkx])


  !GET LOCAL PART OF THE HAMILTONIAN
  Hloc=sum(Hk,dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  call TB_write_Hloc(Hloc)


  !SOLVE ALONG A PATH IN THE BZ.
  Npts=4
  allocate(kpath(Npts,3))
  kpath(1,:)=kpoint_Gamma
  kpath(2,:)=kpoint_M1
  kpath(3,:)=kpoint_X1
  kpath(4,:)=kpoint_Gamma
  call TB_Solve_model(Hk_model,Nso,kpath,Nkpath,&
       colors_name=[red1,blue1],&
       points_name=[character(len=20) :: 'G', 'M', 'X', 'G'],&
       file="Eigenband.tb")



  !Build the local GF:
  allocate(Gmats(Nso,Nso,L))
  allocate(Greal(Nso,Nso,L))
  allocate(Sfoo(Nso,Nso,L))
  Gmats=zero
  Greal=zero
  Sfoo =zero
  call dmft_gloc_matsubara(Hk,Wtk,Gmats,Sfoo,iprint=1)
  Sfoo =zero
  call dmft_gloc_realaxis(Hk,Wtk,Greal,Sfoo,iprint=1)
  do iorb=1,Nso
     dens(iorb) = fft_get_density(Gmats(iorb,iorb,:),beta)
  enddo
  !plot observables
  open(10,file="observables.nint")
  write(10,"(20F20.12)")(dens(iorb),iorb=1,Nso),sum(dens)
  close(10)
  write(*,"(A,20F14.9)")"Occupations =",(dens(iorb),iorb=1,Nso),sum(dens)

  Sfoo = zero
  Eout = dmft_kinetic_energy(Hk,Wtk,Sfoo)
  print*,Eout



contains


  function Hk_model(kpoint,N) result(Hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: kx,ky,epsik,vpsik
    complex(8),dimension(N,N) :: Hk
    if(N/=2)stop "hk_model error: wrong N dimensions"
    kx = kpoint(1)
    ky = kpoint(2)
    !
    epsik = cos(kx)+cos(ky)
    vpsik = sin(kx)*sin(ky)
    !
    Hk(1,1) = ed0 - 2.d0*alpha*tpp*epsik
    Hk(2,2) = ep0 - 2.d0*tpp*epsik
    Hk(1,2) = tpd - 4.d0*v0*vpsik
    Hk(2,1) = tpd - 4.d0*v0*vpsik
    !
  end function Hk_model



end program tdd_pam


