! NAME
!   hm_2bands_hyb
!
! DESCRIPTION
!   Solve the non-interacting HM 2bands with hybridizations 
!   and generate the hamiltoniana matrix H(k), 
!   Model parameters are: t_dd,lambda,Mh.
!   Using the 3d cubic lattice dispersion. 
!   Output is on file *hkfile.in default.
!
! MODEL Hamiltonian is:
!  | Mh - e0*[cos(kx)+cos(ky)+cos(kz)]  ,        lambda*(cos(kx)-cos(ky))*cos(kz)  |
!  | lambda*(cos(kx)-cos(ky))*cos(kz)   ,       -Mh - e0*[cos(kx)+cos(ky)+cos(kz)] |
!

program hm2bhyb
  USE SCIFOR
  implicit none

  integer,parameter    :: L=2000,Norb=2,Lkstep=100
  integer              :: i,j,k
  integer              :: iorb,jorb
  integer              :: ik
  real(8)              :: ix,iy,iz
  real(8)              :: kx,ky,kz
  character(len=32)    :: finput
  character(len=20)    :: file
  integer              :: Lk,Nk,Nkx,Nky,Nkz
  real(8)              :: e0,lambda,Mh,u,xmu,beta,eps
  real(8)              :: wmax,fmesh
  complex(8)           :: zeta
  complex(8)           :: Hk(Norb,Norb)
  complex(8)           :: Gmats(Norb,Norb,L),Greal(Norb,Norb,L)
  real(8)              :: Eval(Norb),Ndens(Norb)
  real(8),allocatable  :: kpath(:,:)
  real(8)              :: wm(L),wr(L)
  logical              :: iexist,ibool,dcflag

  call parse_cmd_variable(finput,"FINPUT",default="inputHM2B.conf")
  call parse_input_variable(Nkx,"NKX",finput,default=20)
  call parse_input_variable(Nky,"NKY",finput,default=20)
  call parse_input_variable(Nkz,"NKZ",finput,default=20)
  call parse_input_variable(e0,"E0",finput,default=1.d0)
  call parse_input_variable(lambda ,"LAMBDA",finput,default=0.d0)
  call parse_input_variable(Mh,"MH",finput,default=0.d0)
  call parse_input_variable(U,"U",finput,default=0.d0)
  call parse_input_variable(xmu,"XMU",finput,default=0.d0)
  call parse_input_variable(beta,"BETA",finput,default=1000.d0)
  call parse_input_variable(file,"FILE",finput,default="hkfile.in")
  call parse_input_variable(eps,"EPS",finput,default=6.d-2)
  call parse_input_variable(wmax,"WMAX",finput,default=10.d0)
  call save_input_file(finput)

  wm = pi/beta*real(2*arange(1,L)-1,8)
  wr = linspace(-wmax,wmax,L,mesh=fmesh)

  Nk=Nkx*Nky*Nkz
  print*,"Using Nk=",Nk
  open(50,file=trim(file))
  write(50,*)"#",Nk,Norb,0,0


  Gmats=zero
  Greal=zero
  ik=0
  call start_progress
  do ix=1,Nkx
     kx = -pi + 2.d0*pi*(ix-1)/Nkx
     do iy=1,Nky
        ky = -pi + 2.d0*pi*(iy-1)/Nky
        do iz=1,Nkz
           kz = -pi + 2.d0*pi*(iz-1)/Nkz
           !
           Hk = Hk_model(kx,ky,kz)
           !
           write(50,"(3(F10.7,1x))")kx,ky,kz
           do iorb=1,Norb
              write(50,"(10(2F10.7,1x))")(Hk(iorb,jorb),jorb=1,Norb)
           enddo
           !
           do i=1,L
              zeta = dcmplx(wr(i),eps)+xmu
              Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(zeta,Hk)
              zeta = xi*wm(i)+xmu
              Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k(zeta,Hk)
           enddo
           !
           ik=ik+1
           call progress_bar(ik,Nk)
        enddo
     enddo
  enddo
  call stop_progress
  Gmats= Gmats/dble(Nk)
  Greal= Greal/dble(Nk)

  Lk=6*Lkstep
  allocate(kpath(Lk,3))
  ik = 0
  !1: (0,0,0)-->(0,0,pi)
  do iz=1,Lkstep
     ik = ik+1
     kx = 0.d0
     ky = 0.d0
     kz = 0.d0 + pi*(iz-1)/Lkstep
     kpath(ik,:)=[kx,ky,kz]
  enddo
  !2: (0,0,pi)-->(0,pi,pi)
  do iy=1,Lkstep
     ik = ik+1
     kx = 0.d0
     ky = 0.d0 + pi*(iy-1)/Lkstep
     kz = pi
     kpath(ik,:)=[kx,ky,kz]
  enddo
  !3: (0,pi,pi)-->(0,pi,0)
  do iz=1,Lkstep
     ik=ik+1
     kx = 0.d0
     ky = pi
     kz = pi - pi*(iz-1)/Lkstep
     kpath(ik,:)=[kx,ky,kz]
  enddo
  !4: (0,pi,0)-->(pi,pi,0)
  do ix=1,Lkstep
     ik  =ik+1
     kx = 0.d0 + pi*(ix-1)/Lkstep
     ky = pi
     kz = 0.d0
     kpath(ik,:)=[kx,ky,kz]
  enddo
  !5: (pi,pi,0)-->(pi,pi,pi)
  do iz=1,Lkstep
     ik = ik+1
     kx = pi
     ky = pi
     kz = 0.d0 + pi*(iz-1)/Lkstep
     kpath(ik,:)=[kx,ky,kz]
  enddo
  !6: (pi,pi,pi)-->(0,0,0)
  do ix=1,Lkstep
     ik = ik+1
     kx = pi - pi*(ix-1)/Lkstep
     ky = pi - pi*(ix-1)/Lkstep
     kz = pi - pi*(ix-1)/Lkstep
     kpath(ik,:)=[kx,ky,kz]
  enddo

  open(10,file="Eigenbands.nint")
  do ik=1,Lk
     Hk = Hk_model(kpath(ik,1),kpath(ik,2),kpath(ik,3))
     call Eigensolve(Hk,Eval)
     write(10,'(I,10F18.12)')ik,(Eval(iorb),iorb=1,Norb)
  enddo
  close(10)


  do iorb=1,Norb
     Ndens(iorb)=get_density_fromFFT(Gmats(iorb,iorb,:),beta)
     call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_iw.nint",wm,Gmats(iorb,iorb,:))
     call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_realw.nint",wr,Greal(iorb,iorb,:))
  enddo


  open(10,file="observables.nint")
  write(10,"(14F20.12)")e0,lambda,Mh,xmu,u,(Ndens(iorb),iorb=1,Norb),sum(Ndens)
  close(10)
  write(*,"(A,10F14.9)")"Occupations                   =",(Ndens(iorb),iorb=1,Norb),sum(Ndens)




contains


  function Hk_model(kx,ky,kz) result(Hk)
    real(8)                         :: kx,ky,kz
    real(8)                         :: epsik,vpsik
    complex(8),dimension(Norb,Norb) :: Hk
    epsik = cos(kx)+cos(ky)+cos(kz)
    vpsik = (cos(kx)-cos(ky))*cos(kz)
    Hk(1,1) = Mh - e0*epsik
    Hk(2,2) =-Mh - e0*epsik
    Hk(1,2) = lambda*vpsik
    Hk(2,1) = lambda*vpsik
  end function Hk_model


  function inverse_g0k(zeta,hk) result(g0k)
    complex(8)                      :: zeta
    complex(8),dimension(Norb,Norb) :: hk
    complex(8),dimension(Norb,Norb) :: g0k
    complex(8)                      :: h11,h22,h12
    integer                         :: i
    g0k=zero
    h11 = zeta - hk(1,1)
    h22 = zeta - hk(2,2)
    h12 =      - hk(1,2)
    g0k(1,1) = one/(h11 - abs(h12)**2/h22)
    g0k(2,2) = one/(h22 - abs(h12)**2/h11)
    g0k(1,2) = -h12/(h11*h22 - abs(h12)**2)
    g0k(2,1) = conjg(g0k(1,2))
  end function inverse_g0k


  subroutine Eigensolve(hk,Evals)
    complex(8),dimension(Norb,Norb) :: hk
    real(8),dimension(Norb)         :: evals
    real(8)                         :: eplus,eminus
    complex(8),dimension(2,2)       :: uk
    complex(8)                      :: delta00,ek,u,v
    !Evals
    eplus   = hk(1,1)+hk(2,2) + sqrt(abs(hk(1,1)-hk(2,2))**2 + 4.d0*hk(1,2)*hk(2,1) )
    eplus   = eplus/2.d0
    eminus  = hk(1,1)+hk(2,2) -sqrt(abs(hk(1,1)-hk(2,2))**2 + 4.d0*hk(1,2)*hk(2,1) )
    eminus  = eminus/2.d0
    !Eigenvectors
    delta00  = -(hk(1,1)-hk(2,2))
    ek      = sqrt( delta00**2 + 4.d0*hk(1,2)*hk(2,1) )
    u       = sqrt(0.5d0*(1.d0 - delta00/ek))
    v       = sqrt(0.5d0*(1.d0 + delta00/ek))
    uk(1,1) = u
    uk(2,2) = u
    uk(1,2) = v
    uk(2,1) =-v
    !
    Hk = Uk
    Evals(1)=Eplus;Evals(2)=Eminus
  end subroutine Eigensolve


  !ADDITIONAL ROUTINES:
  function get_density_fromFFT(giw,beta) result(n)
    complex(8),dimension(:) :: giw
    real(8)                 :: gtau(0:size(giw))
    real(8)                 :: beta,n
    call fftgf_iw2tau(giw,gtau,beta)
    n = -2.d0*gtau(size(giw))
  end function get_density_fromFFT

end program hm2bhyb


