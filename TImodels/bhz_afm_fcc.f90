! DESCRIPTION
!   Solve the non-interacting BHZ model with AFM 2x2 basis 
!   generate the hamiltoniana matrix H(k), 
!   The matrix is written in the Wannier90 form, as expected by w2CTQMC & ED code.
!   Output is on file *hkfile.in default.

program bhz_afm
  USE CONSTANTS
  USE PARSE_INPUT
  USE IOFILE
  USE ARRAYS
  USE MATRIX
  USE TIMER
  implicit none

  integer,parameter    :: L=1000,Norb=4*4
  integer              :: i,j,k,ik,iorb,jorb
  real(8)              :: mh,rh,lambda,delta
  real(8)              :: xmu,beta,eps
  integer              :: dcshift,count
  real(8)              :: epsik,fmesh,n(Norb),eig(Norb)
  integer              :: Nkx,Nk
  real(8)              :: ix,iy
  real(8)              :: kx,ky
  real(8),dimension(L) :: wm,wr
  complex(8)           :: w,Hk(Norb,Norb),Hloc(Norb,Norb),Hksite(Norb,Norb),Hlocsite(Norb,Norb)
  complex(8)           :: fg(L,Norb,Norb),fgr(L,Norb,Norb)
  character(len=20)    :: file,nkstring
  logical              :: iexist,ibool,dcflag
  integer :: icount,jcount,isite,jsite,ispin,jspin,row,col


  call parse_input_variable(nkx,"NKX","inputBHZ.in",default=25)
  call parse_input_variable(mh,"MH","inputBHZ.in",default=3.d0)
  call parse_input_variable(rh,"RH","inputBHZ.in",default=0.d0)
  call parse_input_variable(lambda,"LAMBDA","inputBHZ.in",default=0.3d0)
  call parse_input_variable(delta,"DELTA","inputBHZ.in",default=0.d0)
  call parse_input_variable(xmu,"XMU","inputBHZ.in",default=0.d0)
  call parse_input_variable(eps,"EPS","inputBHZ.in",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputBHZ.in",default=1000.d0)
  call parse_input_variable(file,"FILE","inputBHZ.in",default="hkfile_bhz.in")



  wm = pi/beta*real(2*arange(1,L)-1,8)
  wr = linspace(-10.d0,10.d0,L,mesh=fmesh)

  Nk=Nkx*Nkx
  write(*,*)"Using Nk="//txtfy(Nk)
  open(50,file=trim(file))
  open(51,file="jan_order_"//trim(file))
  write(nkstring,*)Nk
  write(50,*)trim(adjustl(trim(Nkstring)))," 2 2" !Nk, Nwannier, Nbands

  call start_progress
  fgr=zero ; fg =zero ;count=0; Hloc=zero ; Hlocsite=zero
  do ix=1,Nkx
     kx = -pi + 2.d0*pi*dble(ix-1)/dble(Nkx)
     do iy=1,Nkx
        ky = -pi + 2.d0*pi*dble(iy-1)/dble(Nkx)
        Hk     = hk_bhz_afm(mh,lambda,kx,ky)
        Hksite = iso2site(Hk)
        write(50,"(3(F10.7,1x))")kx,ky,pi
        write(51,"(3(F10.7,1x))")kx,ky,pi
        do iorb=1,Norb
           write(50,"(100(2F10.7,1x))")(Hksite(iorb,jorb),jorb=1,Norb)
           write(51,"(100(2F10.7,1x))")(Hk(iorb,jorb),jorb=1,Norb)
        enddo


        Hlocsite=Hlocsite + Hksite/dble(Nk)
        Hloc=Hloc+Hk/dble(Nk)
        do i=1,L
           w = dcmplx(wr(i),eps)+xmu
           fgr(i,:,:)=fgr(i,:,:) + inverse_g0k_afm(w,Hk)/dble(Nk)
           w = xi*wm(i)+xmu
           fg(i,:,:) =fg(i,:,:)  + inverse_g0k_afm(w,Hk)/dble(Nk)
        enddo
        count=count+1
        !if(mod(count,250))print*,count,"/",Nkx*Nkx
        call progress(count,Nk)
     enddo
  enddo
  call stop_progress

  open(10,file="bhz_afm_hloc.dat")
  do iorb=1,Norb
     write(10,"(90F21.12)")(dreal(Hloc(iorb,jorb)),jorb=1,Norb)
  enddo
  write(10,*)""
  do iorb=1,Norb
     write(10,"(90F21.12)")(dimag(Hloc(iorb,jorb)),jorb=1,Norb)
  enddo
  write(10,*)""
  do iorb=1,Norb
     write(10,"(90F21.12)")(dreal(Hlocsite(iorb,jorb)),jorb=1,Norb)
  enddo
  write(10,*)""
  do iorb=1,Norb
     write(10,"(90F21.12)")(dimag(Hlocsite(iorb,jorb)),jorb=1,Norb)
  enddo
  close(10)

  ik = 0
  open(10,file="Eigenbands.dat")
  !From \Gamma=(0,0) to X=(pi,0): 100 steps
  do ix=1,100
     ik=ik+1
     kx = 0.d0 + pi*real(ix-1,8)/100.d0
     ky = 0.d0
     Hk(:,:) = hk_bhz_afm(mh,lambda,kx,ky)
     eig = Eigk_afm(Hk)
     write(10,"(I,16F25.12)")ik,(eig(i),i=1,Norb)
  enddo
  !From X=(pi,0) to M=(pi,pi): 100 steps
  do iy=1,100
     ik=ik+1
     kx = pi
     ky = 0.d0 + pi*real(iy-1,8)/100.d0
     Hk(:,:) = hk_bhz_afm(mh,lambda,kx,ky)
     eig = Eigk_afm(Hk)
     write(10,"(I,16F25.12)")ik,(eig(i),i=1,Norb)
  enddo
  !From M=(pi,pi) to \Gamma=(0,0): 100 steps
  do ix=1,100
     ik=ik+1
     iy=ix
     kx = pi - pi*real(ix-1,8)/100.d0
     ky = pi - pi*real(iy-1,8)/100.d0
     Hk(:,:) = hk_bhz_afm(mh,lambda,kx,ky)
     eig = Eigk_afm(Hk)
     write(10,"(I,16F25.12)")ik,(eig(i),i=1,Norb)
  enddo
  close(10)

  open(10,file="DOS.dat")
  open(11,file="G_iw.dat")
  do i=1,L
     write(10,"(20F20.12)") wr(i),(-dimag(fgr(i,j,j))/pi,j=1,size(fgr,3))
     write(11,"(100F20.12)") wm(i),(dimag(fg(i,j,j)),j=1,size(fgr,3)),(dreal(fg(i,j,j)),j=1,size(fgr,3))
  enddo
  close(10)
  close(11)

  do iorb=1,Norb
     n(iorb) = -2.d0*sum(dimag(fgr(:,iorb,iorb))*fermi(wr(:),beta))*fmesh/pi
  enddo
  open(10,file="observables.dat")
  write(10,"(24F20.12)")mh,lambda,rh,delta,xmu,(n(iorb),iorb=1,Norb),sum(n)
  close(10)
  write(*,"(A,16F14.9)")"Occupations                   =",n


  write(*,*)"Remember to open the file:"//trim(file)
  write(*,*)"Parameters for "//trim(file)//" are in +paramaters4_"//trim(file)
  ! open(10,file="parameters_bhz_"//trim(file))
  ! write(10,nml=hkvars)
  ! close(10)
  close(50)




contains




  function hk_bhz_afm(mh,lambda,kx,ky) result(hk)
    real(8)                     :: mh
    real(8)                     :: lambda
    real(8)                     :: kx,ky
    complex(8),dimension(16,16) :: hk
    hk            = zero
    hk(1:8,1:8)   = hk_bhz_afm8x8(mh,lambda,kx,ky)
    hk(9:16,9:16) = conjg(hk_bhz_afm8x8(mh,lambda,-kx,-ky))
  end function hk_bhz_afm

  function hk_bhz_afm8x8(mh,lambda,kx,ky) result(hk)
    real(8)                     :: mh
    real(8)                     :: lambda
    real(8)                     :: kx,ky
    complex(8),dimension(8,8)   :: hk
    hk(1,1:8)=[one*mh,zero,-one/2.d0,xi*lambda/2.d0,-one/2.d0,-one*lambda/2.d0,zero,zero]
    hk(2,1:8)=[zero,-one*mh,xi*lambda/2.d0,one/2.d0,one*lambda/2.d0,one/2.d0,zero,zero]
    hk(3,1:8)=[-one/2.d0,-xi*lambda/2.d0,one*mh,zero,zero,zero,-one/2.d0,-one*lambda/2.d0]
    hk(4,1:8)=[-xi*lambda/2.d0,one/2.d0,zero,-one*mh,zero,zero,one*lambda/2.d0,one/2.d0]
    hk(5,1:8)=[-one/2.d0,one*lambda/2.d0,zero,zero,one*mh,zero,-one/2.d0,xi*lambda/2.d0]
    hk(6,1:8)=[-one*lambda/2.d0,one/2.d0,zero,zero,zero,-one*mh,xi*lambda/2.d0,one/2.d0]
    hk(7,1:8)=[zero,zero,-one/2.d0,one*lambda/2.d0,-one/2.d0,-xi*lambda/2.d0,one*mh,zero]
    hk(8,1:8)=[zero,zero,-one*lambda/2.d0,one/2.d0,-xi*lambda/2.d0,one/2.d0,zero,-one*mh]
    !
    hk(1,1:8)=hk(1,1:8)+[zero,zero,exp(-xi*2.d0*kx)*(-0.5d0),exp(-xi*2.d0*kx)*(-xi*lambda/2.d0),exp(xi*2.d0*ky)*(-0.5d0),exp(xi*2.d0*ky)*(lambda/2.d0),zero,zero]
    hk(2,1:8)=hk(2,1:8)+[zero,zero,exp(-xi*2.d0*kx)*(-xi*lambda/2.d0),exp(-xi*2.d0*kx)*(0.5d0),exp(xi*2.d0*ky)*(-lambda/2.d0),exp(xi*2.d0*ky)*(0.5d0),zero,zero]
    hk(3,1:8)=hk(3,1:8)+[exp(xi*2.d0*kx)*(-0.5d0),exp(xi*2.d0*kx)*(xi*lambda/2.d0),zero,zero,zero,zero,exp(xi*2.d0*ky)*(-0.5d0),exp(xi*2.d0*ky)*(lambda/2.d0)]
    hk(4,1:8)=hk(4,1:8)+[exp(xi*2.d0*kx)*(xi*lambda/2.d0),exp(xi*2.d0*kx)*(0.5d0),zero,zero,zero,zero,exp(xi*2.d0*ky)*(-lambda/2.d0),exp(xi*2.d0*ky)*(0.5d0)]
    hk(5,1:8)=hk(5,1:8)+[exp(-xi*2.d0*ky)*(-0.5d0),exp(-xi*2.d0*ky)*(-lambda/2.d0),zero,zero,zero,zero,exp(-xi*2.d0*kx)*(-0.5d0),exp(-xi*2.d0*kx)*(-xi*lambda/2.d0)]
    hk(6,1:8)=hk(6,1:8)+[exp(-xi*2.d0*ky)*(lambda/2.d0),exp(-xi*2.d0*ky)*(0.5d0),zero,zero,zero,zero,exp(-xi*2.d0*kx)*(-xi*lambda/2.d0),exp(-xi*2.d0*kx)*(0.5d0)]
    hk(7,1:8)=hk(7,1:8)+[zero,zero,exp(-xi*2.d0*ky)*(-0.5d0),exp(-xi*2.d0*ky)*(-lambda/2.d0),exp(xi*2.d0*kx)*(-0.5d0),exp(xi*2.d0*kx)*(xi*lambda/2.d0),zero,zero]
    hk(8,1:8)=hk(8,1:8)+[zero,zero,exp(-xi*2.d0*ky)*(lambda/2.d0),exp(-xi*2.d0*ky)*(0.5d0),exp(xi*2.d0*kx)*(xi*lambda/2.d0),exp(xi*2.d0*kx)*(0.5d0),zero,zero]
  end function hk_bhz_afm8x8







  function inverse_g0k_afm(zeta,hk,type) result(g0k)
    complex(8)                  :: zeta
    complex(8),dimension(16,16) :: hk
    complex(8),dimension(16,16) :: g0k,g0ktmp
    integer,optional            :: type
    integer                     :: type_
    type_=0 ; if(present(type))type_=type
    if(type_==0)then
       g0k = -hk
       forall(i=1:16)g0k(i,i)=zeta - hk(i,i)
       call matrix_inverse(g0k)
    else
       g0k=zero
       g0k(1:8,1:8)   = inverse_g0k_afm8x8(zeta,hk(1:8,1:8))
       g0k(9:16,9:16) = inverse_g0k_afm8x8(zeta,hk(9:16,9:16))
    endif
  end function inverse_g0k_afm

  function inverse_g0k_afm8x8(zeta,hk) result(g0k)
    complex(8)                :: zeta
    complex(8),dimension(8,8) :: hk
    complex(8),dimension(8,8) :: g0k
    integer :: i
    g0k = -hk
    forall(i=1:8)g0k(i,i)=zeta - hk(i,i)
    call matrix_inverse(g0k)
  end function inverse_g0k_afm8x8









  function Eigk_afm(hk) result(eig)
    complex(8),dimension(16,16) :: hk
    real(8),dimension(16)       :: eig
    call matrix_diagonalize(hk,eig)
  end function Eigk_afm





  function iso2site(hktmp) result(hk)
    complex(8),dimension(16,16) :: hktmp,hk
    integer                     :: icount,isite,ispin,iband,iorb
    integer                     :: jcount,jsite,jspin,jband,jorb
    integer                     :: row,col
    icount=0
    do isite=1,4
       do ispin=1,2
          do iorb=1,2
             icount=icount+1
             row = (ispin-1)*2*4 + (isite-1)*2 + iorb
             jcount=0
             do jsite=1,4
                do jspin=1,2
                   do jorb=1,2
                      jcount=jcount+1
                      col = (jspin-1)*2*4 + (jsite-1)*2 + jorb
                      hk(icount,jcount) = hktmp(row,col)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function iso2site

  function site2iso(hk) result(hktmp)
    complex(8),dimension(16,16) :: hktmp,hk
    integer                     :: icount,isite,ispin,iband,iorb
    integer                     :: jcount,jsite,jspin,jband,jorb
    integer                     :: row,col
    icount=0
    do isite=1,4
       do ispin=1,2
          do iorb=1,2
             icount=icount+1
             row = (ispin-1)*2*4 + (isite-1)*2 + iorb
             jcount=0
             do jsite=1,4
                do jspin=1,2
                   do jorb=1,2
                      jcount=jcount+1
                      col = (jspin-1)*2*4 + (jsite-1)*2 + jorb
                      hktmp(row,col) = hk(icount,jcount)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function site2iso





  elemental function fermi(x,beta)
    real(8),intent(in) :: x, beta 
    real(8)            :: fermi
    if(x*beta > 100.d0)then
       fermi=0.d0
       return
    endif
    fermi = 1.d0/(1.d0+exp(beta*x))
  end function fermi


end program bhz_afm







