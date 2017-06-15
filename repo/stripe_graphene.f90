program graphene
  USE SCIFOR
  USE VECTORS
  implicit none


  integer,parameter                  :: L=1000
  integer                            :: i,j,k,ik,ix,iy,Nlev
  real(8)                            :: ts,tsp,phi,delta,xmu,beta,eps,wmax
  integer                            :: Nkx,Nky,Nk,ipath
  real(8)                            :: kx,ky
  real(8),dimension(2)               :: Kvec
  real(8),dimension(2)               :: a1,a2,a3
  real(8),dimension(2)               :: nn1,nn2,nn3
  complex(8),allocatable             :: Hk(:,:),fg(:,:,:),fgr(:,:,:)
  real(8),dimension(:,:),allocatable :: fk
  real(8),dimension(:),allocatable   :: kxgrid,kygrid
  logical                            :: iexist,gflag
  real(8),dimension(L)               :: wm,wr
  real(8),dimension(:,:),allocatable :: KPath
  complex(8) :: iw
  namelist/hkvars/nkx,nky,ts,tsp,xmu,beta,eps,Nlev,gflag

  Nlev=10
  nkx=0
  nky=50
  ts=1.d0
  tsp=0.d0
  phi=0.d0
  delta=0.d0
  xmu=0.d0
  wmax=10.d0
  eps=3.d-2
  beta=100.d0
  gflag=.true.

  inquire(file="inputGRAPHENE.in",exist=iexist)
  if(iexist)then
     open(10,file="inputGRAPHENE.in")
     read(10,nml=hkvars)
     close(10)
  else
     open(10,file="default.inputGRAPHENE.in")
     write(10,nml=hkvars)
     close(10)
     stop "inputFILE not found: write  default.inputGRAPHENE.in"
  endif

  call parse_cmd_variable(nlev,"NLEV")
  call parse_cmd_variable(nkx,"NKX")
  call parse_cmd_variable(nky,"NKY")
  call parse_cmd_variable(ts,"TS")
  call parse_cmd_variable(tsp,"TSP")
  call parse_cmd_variable(xmu,"XMU")
  call parse_cmd_variable(wmax,"wmax")
  call parse_cmd_variable(eps,"EPS")
  call parse_cmd_variable(beta,"BETA")
  call parse_cmd_variable(gflag,"GFLAG")

  open(20,file="parameters_graphene.hc")
  write(20,nml=hkvars)
  close(20)

  wm = pi/beta*real(2*arange(1,L)-1,8)
  wr = linspace(-wmax,wmax,L)

  Nk=Nky
  call msg("Using Nk="//txtfy(Nk))

  allocate(kxgrid(Nkx),kygrid(Nky),fk(Nk,Nlev))
  allocate(Hk(Nlev,Nlev),fg(L,Nlev,Nlev),fgr(L,Nlev,Nlev))

  call start_progress

  !Honeycomb lattice basis: alat=1
  a1=[3.d0, sqrt(3.d0)]/2.d0
  a2=[3.d0,-sqrt(3.d0)]/2.d0
  a3=a2-a1

  !nearest-neighbor displacements:
  nn1=[ 1.d0/2.d0, sqrt(3.d0)/2.d0]
  nn2=[ 1.d0/2.d0,-sqrt(3.d0)/2.d0]
  nn3=[-1.d0     ,0.d0]

  if(gflag)write(*,*)"Calculating the GF:"
  fg =zero
  fgr=zero
  do iy=1,Nk
     ky = -pi + 2.d0*pi*real(iy-1,8)/real(Nky,8)
     kygrid(iy)=ky
     Kvec=[0.d0,ky]

     if(gflag)then
        Hk = get_hk(Nlev,kvec)
        do i=1,L
           iw = cmplx(wr(i),eps,8)+xmu
           fgr(i,:,:)=fgr(i,:,:) + inverse_g0k(Nlev,iw,Hk)
           iw = xi*wm(i)+xmu
           fg(i,:,:) =fg(i,:,:)  + inverse_g0k(Nlev,iw,Hk)
        enddo
     endif

     Hk = get_hk(Nlev,kvec)
     call matrix_diagonalize(Hk,fk(iy,:))
     call progress(iy,Nky)
  enddo

  call stop_progress

  open(100,file="Eigk.hc")
  do i=1,Nlev
     do iy=1,Nky
        write(100,*)kygrid(iy),fk(iy,i)
     enddo
     write(100,*)""
  enddo
  close(100)

  if(gflag)then
     fgr= fgr/real(Nk,8)
     fg = fg/real(Nk,8)
     open(101,file="DOS.hc")
     open(102,file="G_iw.hc")
     do i=1,Nlev
        do iy=1,L
           write(101,*)wr(iy),-dimag(fgr(iy,i,i))/pi-0.1d0*dble(i-1)
           write(102,*)wm(iy),dimag(fg(iy,i,i)),dreal(fg(iy,i,i))
        enddo
        write(100,*)""
        write(101,*)""
        write(102,*)""
     enddo
     close(101)
     close(102)
  endif

contains

  function get_hk(N,kpnt) result(hk)
    real(8),dimension(2)      :: kpnt
    complex(8),dimension(N,N) :: hk
    real(8)                   :: Rlat(2)
    integer                   :: i,j,N
    real(8)                   :: arg
    Rlat=[0.d0,sqrt(3.d0)]
    arg = dot_product(kpnt,Rlat)
    Hk=zero
    do i=1,Nlev-1,2
       Hk(i,i+1)=-ts*(1.d0+exp(-xi*arg))      
    enddo
    do i=2,Nlev-1,2
       Hk(i,i+1)=-ts
    enddo
    forall(i=1:N,j=1:N,j<i)Hk(i,j)=conjg(Hk(j,i))
  end function get_hk

  function inverse_g0k(N,iw,hk) result(g0k)
    integer                   :: i,N
    complex(8),dimension(N,N) :: hk
    complex(8)                :: iw
    complex(8),dimension(N,N) :: g0k
    g0k=-hk
    forall(i=1:N)g0k(i,i)=iw+g0k(i,i)
    call matrix_inverse(g0k)
  end function inverse_g0k


end program graphene
