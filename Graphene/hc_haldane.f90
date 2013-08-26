program graphene
  USE SCIFOR
  USE VECTORS
  implicit none


  integer,parameter                  :: L=2000,Norb=2
  integer                            :: i,j,k,ik,ix,iy,iorb,jorb
  real(8)                            :: ts,tsp,phi,delta,xmu,beta,eps,wmax,Mh
  integer                            :: Nkx,Nky,Nk,ipath
  real(8)                            :: kx,ky
  real(8),dimension(2) :: Kvec
  real(8),dimension(2) :: a1,a2,a3
  real(8),dimension(2) :: nn1,nn2,nn3
  complex(8)                         :: Hk(Norb,Norb),fg(L,Norb,Norb),fgr(L,Norb,Norb),w
  real(8),dimension(:,:),allocatable :: fk
  real(8),dimension(:),allocatable   :: kxgrid,kygrid
  logical                            :: iexist
  real(8),dimension(L)               :: wm,wr
  real(8),dimension(:,:),allocatable :: KPath
  namelist/hkvars/nkx,nky,ts,tsp,xmu,beta,eps

  nkx=50
  nky=50
  ts=1.d0
  tsp=1.d0/3.d0/sqrt(3.d0)
  Mh=0.d0
  phi=0.d0
  xmu=0.d0
  wmax=10.d0
  eps=3.d-2
  beta=100.d0

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

  call parse_cmd_variable(nkx,"NKX")
  call parse_cmd_variable(nky,"NKY")
  call parse_cmd_variable(ts,"TS")
  call parse_cmd_variable(tsp,"TSP")
  call parse_cmd_variable(phi,"phi")
  call parse_cmd_variable(Mh,"Mh")
  call parse_cmd_variable(xmu,"XMU")
  call parse_cmd_variable(wmax,"wmax")
  call parse_cmd_variable(eps,"EPS")
  call parse_cmd_variable(beta,"BETA")

  open(20,file="parameters_graphene.hc")
  write(20,nml=hkvars)
  close(20)

  wm = pi/beta*real(2*arange(1,L)-1,8)
  wr = linspace(-wmax,wmax,L)

  Nk=Nkx*Nky
  call msg("Using Nk="//txtfy(Nk))

  allocate(kxgrid(Nkx),kygrid(Nky),fk(Nkx,Nky))


  call start_progress

  !Honeycomb lattice basis: alat=1
  a1=[3.d0, sqrt(3.d0)]/2.d0
  a2=[3.d0,-sqrt(3.d0)]/2.d0
  a3=a2-a1

  !nearest-neighbor displacements:
  nn1=[ 1.d0/2.d0, sqrt(3.d0)/2.d0]
  nn2=[ 1.d0/2.d0,-sqrt(3.d0)/2.d0]
  nn3=[-1.d0     ,0.d0]



  ik=0
  do ix=1,Nkx
     kx = -pi + 2.d0*pi*real(ix-1,8)/real(Nkx,8)
     kxgrid(ix)=kx
     do iy=1,Nky
        ky = -pi + 2.d0*pi*real(iy-1,8)/real(Nky,8)
        kygrid(iy)=ky
        ik=ik+1
        Kvec=[kx,ky]
        Hk = get_hk(kvec)

        do i=1,L
           w = cmplx(wr(i),eps,8)+xmu
           fgr(i,:,:)=fgr(i,:,:) + inverse_g0k(w,Hk)
           w = xi*wm(i)+xmu
           fg(i,:,:) =fg(i,:,:)  + inverse_g0k(w,Hk)
        enddo

        fk(ix,iy)=eplus(hk)
        call progress(ik,Nk)
     enddo
  enddo
  close(101)
  call stop_progress

  fgr= fgr/real(Nk,8)
  fg = fg/real(Nk,8)

  call splot3d("Eigk.hc",kxgrid,kygrid,fk(:,:))
  call splot("DOS.hc",wr,-dimag(fgr(:,1,1))/pi)
  call splot("G_iw.hc",wm,fg(:,1,1))


  allocate(Kpath(4,2))
  KPath(1,:)=[0.0,0.0]
  KPath(2,:)=[1.d0/3.d0, 1.d0/3.d0/sqrt(3.d0)]*2.d0*pi
  KPath(3,:)=[1.d0/3.d0,-1.d0/3.d0/sqrt(3.d0)]*2.d0*pi
  KPath(4,:)=[0.0,0.0]
  open(10,file="KPath.hc")
  do ipath=1,3
     do j=1,100
        Kvec = Kpath(ipath,:)+(Kpath(ipath+1,:)-Kpath(ipath,:))*dble(j)/dble(100)
        Hk=get_hk(Kvec)
        write(10,*)(ipath-1)*100+j,eplus(hk),eminus(hk)
     enddo
  enddo
  close(10)

contains

  function get_hk(kpnt) result(hk)
    real(8),dimension(2) :: kpnt
    complex(8),dimension(2,2) :: hk
    real(8)    :: arg1,arg2,arg3
    complex(8) :: fkp,epsk
    arg1 = dot_product(kpnt,nn1)
    arg2 = dot_product(kpnt,nn2)
    arg3 = dot_product(kpnt,nn3)
    fkp = exp(-xi*arg1)+exp(-xi*arg2)+exp(-xi*arg3)
    arg1 = dot_product(kpnt,a1)
    arg2 = dot_product(kpnt,a2)
    arg3 = dot_product(kpnt,a3)
    Hk(1,1) =  Mh-2.d0*tsp*(cos(arg1) + cos(arg2) + cos(arg3))*exp(-xi*pi*phi)
    Hk(2,2) = -Mh-2.d0*tsp*(cos(arg1) + cos(arg2) + cos(arg3))*exp(xi*pi*phi)
    Hk(1,2) = -ts*fkp
    Hk(2,1) = -ts*conjg(fkp)
  end function get_hk

  function inverse_g0k(iw,hk) result(g0k)
    integer                     :: i
    complex(8),dimension(2,2)   :: hk
    complex(8)                  :: iw
    complex(8),dimension(2,2)   :: g0k
    complex(8)                  :: delta,ppi,vmix12,vmix21
    g0k=zero
    delta = iw - hk(1,1)
    ppi   = iw - hk(2,2)
    vmix12  = -hk(1,2)
    vmix21  = -hk(2,1)
    g0k(1,1) = one/(delta - vmix12*vmix21/ppi)
    g0k(2,2) = one/(ppi - vmix12*vmix21/delta)
    g0k(1,2) = -vmix12/(ppi*delta - vmix12*vmix21)
    g0k(2,1) = -vmix21/(ppi*delta - vmix12*vmix21)!conjg(g0k(1,2))
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
    eminus = hk(1,1)+hk(2,2) - sqrt(abs(hk(1,1)-hk(2,2))**2 + 4.d0*hk(1,2)*hk(2,1) )
    eminus = eminus/2.d0
  end function eminus

end program graphene
