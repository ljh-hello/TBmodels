! include 'dmft_get_gloc_dos.f90'
program hm_2b_flat
  USE SCIFOR
  USE DMFT_TOOLS
  ! USE DMFT_GET_GLOC_DOS
  implicit none

  integer,parameter                       :: Norb=2,Nspin=1,Nso=Nspin*Norb
  integer                                 :: Le,L
  integer                                 :: i,j,ie,iorb,jorb
  real(8),dimension(:,:),allocatable      :: Dbands
  real(8),dimension(:,:),allocatable      :: Ebands
  real(8)                                 :: Hloc(Nso)
  real(8)                                 :: Wband(Nso)
  real(8)                                 :: delta
  real(8)                                 :: xmu,beta,eps,Eout(2)
  real(8)                                 :: dens(Nso),de(Nso)


  complex(8),dimension(:,:,:,:,:),allocatable :: Gmats,Greal,Sfoo !(Nspin,Nspin,Norb,Norb,L)
  complex(8),dimension(:,:,:),allocatable :: Spoo !(Nso,Nso,L)
  character(len=20)                       :: finput

  call parse_cmd_variable(finput,"FINPUT",default="inputHM2b.conf")
  call parse_input_variable(Le,"LE",finput,default=500)
  call parse_input_variable(L,"L",finput,default=2048)
  call parse_input_variable(Wband,"WBAND",finput,default=[1d0,0.5d0])
  call parse_input_variable(delta,"DELTA",finput,default=0d0)
  call parse_input_variable(xmu,"XMU",finput,default=0.d0)
  call parse_input_variable(eps,"EPS",finput,default=1.d-2)
  call parse_input_variable(beta,"BETA",finput,default=1000.d0)
  call save_input_file(finput)
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-10d0,"wini")
  call add_ctrl_var(10d0,"wfin")
  call add_ctrl_var(eps,"eps")


  allocate(Ebands(Nso,Le))
  allocate(Dbands(Nso,Le))
  Ebands(1,:) = linspace(-Wband(1),Wband(1),Le,mesh=de(1))
  Ebands(2,:) = linspace(-Wband(2),Wband(2),Le,mesh=de(2))
  !
  Dbands(1,:) = dens_flat(Ebands(1,:),Wband(1))*de(1)
  Dbands(2,:) = dens_flat(Ebands(2,:),Wband(2))*de(2)


  Hloc=[-Delta/2,Delta/2]
  call TB_write_Hloc(one*diag(Hloc))


  !Build the local GF:
  allocate(Gmats(Nspin,Nspin,Norb,Norb,L))
  allocate(Greal(Nspin,Nspin,Norb,Norb,L))
  allocate(Sfoo(Nspin,Nspin,Norb,Norb,L))
  Gmats=zero
  Greal=zero
  Sfoo =zero
  ! call dmft_get_gloc_matsubara_normal_dos(Ebands,Dbands,Hloc,Gmats,Sfoo,iprint=1)
  call dmft_gloc_matsubara(Ebands,Dbands,Hloc,Gmats,Sfoo,iprint=1)
  Sfoo =zero
  ! call dmft_get_gloc_realaxis_normal_dos(Ebands,Dbands,Hloc,Greal,Sfoo,iprint=1)
  call dmft_gloc_realaxis(Ebands,Dbands,Hloc,Greal,Sfoo,iprint=1)
  do iorb=1,Nso
     dens(iorb) = fft_get_density(Gmats(1,1,iorb,iorb,:),beta)
  enddo
  !plot observables
  open(10,file="observables.nint")
  write(10,"(20F20.12)")(dens(iorb),iorb=1,Nso),sum(dens)
  close(10)
  write(*,"(A,20F14.9)")"Occupations =",(dens(iorb),iorb=1,Nso),sum(dens)

  allocate(Spoo(Nso,Nso,L));Spoo=zero
  Eout = dmft_kinetic_energy(Ebands,Dbands,Hloc,Spoo)
  print*,Eout




contains




  function dens_flat(ebands,wband) result(dens)
    real(8),dimension(:)            :: ebands
    real(8)                         :: wband
    real(8),dimension(size(ebands)) :: dens
    integer                         :: i
    real(8)                         :: e
    do i=1,size(ebands)
       e=ebands(i)
       dens(i)= step(wband-abs(e))/(2*wband)
    enddo
  end function dens_flat




  ! function dmft_kinetic_energy_normal_dos(Ebands,Dbands,Hloc,Sigma) result(Eout)
  !   real(8),dimension(:,:),intent(in)                           :: Ebands  ![Nspin*Norb][Lk]
  !   real(8),dimension(size(Ebands,1),size(Ebands,2)),intent(in) :: Dbands  ![Nspin*Norb][Lk]
  !   real(8),dimension(size(Ebands,1)),intent(in)                :: Hloc    ![Nspin*Norb]
  !   complex(8),dimension(:,:,:)                                 :: Sigma   ![Nspin*Norb][Nspin*Norb][L]
  !   !
  !   integer                                                     :: Lk,Nso,Liw
  !   integer                                                     :: i,ik,iso
  !   !
  !   integer                                                     :: Norb,Nporb
  !   integer                                                     :: Nspin  
  !   real(8)                                                     :: beta
  !   real(8)                                                     :: xmu
  !   !
  !   real(8),dimension(size(Ebands,1),size(Ebands,1))            :: Sigma_HF
  !   !
  !   complex(8)                                                  :: Ak,Bk,Ck,Dk
  !   complex(8)                                                  :: Gk,Tk
  !   real(8)                                                     :: Tail0,Tail1,Lail0,Lail1,spin_degeneracy
  !   !
  !   real(8)                                                     :: H0,Hl
  !   real(8)                                                     :: Ekin,Eloc
  !   real(8)                                                     :: Eout(2)
  !   !
  !   real(8),dimension(:),allocatable                            :: wm
  !   !Retrieve parameters:
  !   call get_ctrl_var(Norb,"NORB")
  !   call get_ctrl_var(Nspin,"NSPIN")
  !   call get_ctrl_var(beta,"BETA")
  !   call get_ctrl_var(xmu,"XMU")
  !   !
  !   Nso = size(Ebands,1)
  !   Lk  = size(Ebands,2)
  !   Liw = size(Sigma,3)
  !   !Testing:
  !   if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_normal_dos: Nso != Norb*Nspin [from Hk]"
  !   call assert_shape(Sigma,[Nso,Nso,Liw],"dmft_kinetic_energy_normal_main","Sigma")
  !   !
  !   !Allocate and setup the Matsubara freq.
  !   if(allocated(wm))deallocate(wm);allocate(wm(Liw))
  !   wm = pi/beta*dble(2*arange(1,Liw)-1)
  !   !
  !   !Get HF part of the self-energy
  !   Sigma_HF = dreal(Sigma(:,:,Liw))
  !   !
  !   !
  !   write(*,"(A)") "Kinetic energy computation"
  !   call start_timer()
  !   H0=0d0
  !   Hl=0d0
  !   !Get principal part: Tr[ Hk.(Gk-Tk) ]
  !   do ik=1,Lk
  !      do iso=1,Nso
  !         Ak = Ebands(iso,ik)
  !         Bk =-Ebands(iso,ik) - Sigma_HF(iso,iso) 
  !         do i=1,Liw
  !            Gk = (xi*wm(i)+xmu) - Sigma(iso,iso,i) - Ebands(iso,ik)
  !            Gk = 1d0/Gk
  !            Tk = 1d0/(xi*wm(i)) - Bk/(xi*wm(i))**2
  !            Ck = Ak*(Gk - Tk)
  !            Dk = Hloc(iso)*(Gk - Tk)
  !            H0 = H0 + Dbands(iso,ik)*Ck
  !            Hl = Hl + Dbands(iso,ik)*Dk
  !         enddo
  !      enddo
  !      call eta(ik,Lk)
  !   enddo
  !   call stop_timer()
  !   spin_degeneracy=3d0-Nspin     !2 if Nspin=1, 1 if Nspin=2
  !   H0=H0/beta*2*spin_degeneracy
  !   Hl=Hl/beta*2*spin_degeneracy
  !   !
  !   !get tail subtracted contribution: Tr[ Hk.Tk ]
  !   Tail0=0d0
  !   Tail1=0d0
  !   Lail0=0d0
  !   Lail1=0d0
  !   do ik=1,Lk
  !      do iso=1,Nso
  !         Ak = Ebands(iso,ik)
  !         Bk =-Ebands(iso,ik) - Sigma_HF(iso,iso)
  !         Ck= Ak*Bk
  !         Dk= Hloc(iso)*Bk
  !         Tail0 = Tail0 + 0.5d0*Dbands(iso,ik)*Ak
  !         Tail1 = Tail1 + 0.25d0*Dbands(iso,ik)*Ck
  !         Lail0 = Lail0 + 0.5d0*Dbands(iso,ik)*Hloc(iso)
  !         Lail1 = Lail1 + 0.25d0*Dbands(iso,ik)*Dk
  !      enddo
  !   enddo
  !   Tail0=Tail0*spin_degeneracy
  !   Tail1=Tail1*beta*spin_degeneracy
  !   Lail0=Lail0*spin_degeneracy
  !   Lail1=Lail1*beta*spin_degeneracy
  !   !
  !   Ekin=H0+Tail0+Tail1
  !   Eloc=Hl+Lail0+Lail1
  !   Eout = [Ekin,Eloc]
  !   !
  !   call write_kinetic_info()
  !   call write_kinetic_value(Eout)
  !   !
  !   deallocate(wm)
  ! end function dmft_kinetic_energy_normal_dos


  ! !+-------------------------------------------------------------------+
  ! !PURPOSE  : write legend, i.e. info about columns 
  ! !+-------------------------------------------------------------------+
  ! subroutine write_kinetic_info()
  !   integer :: unit
  !   unit = free_unit()
  !   open(unit,file="dmft_kinetic_energy.info")
  !   write(unit,"(A1,90(A14,1X))")"#",reg(txtfy(1))//"<K>",reg(txtfy(2))//"<Eloc>"
  !   close(unit)
  ! end subroutine write_kinetic_info


  ! !+-------------------------------------------------------------------+
  ! !PURPOSE  : Write energies to file
  ! !+-------------------------------------------------------------------+
  ! subroutine write_kinetic_value(Eout)
  !   real(8) :: Eout(2),Ekin,Eloc
  !   integer :: unit
  !   unit = free_unit()
  !   Ekin=Eout(1)
  !   Eloc=Eout(2)
  !   open(unit,file="dmft_kinetic_energy.dat")
  !   write(unit,"(90F15.9)")Ekin,Eloc
  !   close(unit)
  ! end subroutine write_kinetic_value


  ! function get_kinetic_energy(Hk,Liw) result(ed_Ekin)
  !   integer                                  :: Lk,No,Liw
  !   integer                                  :: i,ik,iorb
  !   complex(8),dimension(:,:,:)              :: Hk ![Nso][Nso][Nk]
  !   real(8),dimension(size(Hk,3))            :: Wtk
  !   !
  !   real(8),dimension(:),allocatable         :: wm
  !   complex(8),dimension(:,:),allocatable    :: Ak,Bk
  !   complex(8),dimension(:,:),allocatable    :: Ck,Zk
  !   complex(8),dimension(:,:),allocatable    :: Zeta,Gk,Tk
  !   real(8)                                  :: Tail0,Tail1,spin_degeneracy
  !   !
  !   real(8)                                  :: H0,ed_Ekin
  !   !
  !   No = size(Hk,1)
  !   Lk = size(Hk,3)
  !   if(No/=size(Hk,2))stop "get_kinetic_energy: size(Hk,1)!=size(Hk,2) [Norb_total]"
  !   if(Lk/=size(Wtk))stop "get_kinetic_energy: size(Wtk)!=size(Hk,3) [L_k]"
  !   !
  !   allocate(wm(Liw))
  !   allocate(Ak(No,No),Bk(No,No),Ck(No,No),Zk(No,No),Zeta(No,No),Gk(No,No),Tk(No,No))
  !   !
  !   wm = pi/beta*dble(2*arange(1,Liw)-1)
  !   Wtk=1d0/Lk
  !   !
  !   H0=0d0
  !   Zk=0d0 ; forall(i=1:No)Zk(i,i)=1d0
  !   do ik=1,Lk
  !      Bk=-Hk(:,:,ik)
  !      do i=1,Liw
  !         Gk = (xi*wm(i)+xmu)*Zk(:,:) - Hk(:,:,ik)
  !         select case(No)
  !         case default
  !            call matrix_inverse(Gk)
  !         case(1)
  !            Gk = 1d0/Gk
  !         end select
  !         Tk = Zk(:,:)/(xi*wm(i)) - Bk(:,:)/(xi*wm(i))**2
  !         Ck = matmul(Hk(:,:,ik),Gk - Tk)
  !         H0 = H0 + Wtk(ik)*trace_matrix(Ck,No)
  !      enddo
  !   enddo
  !   spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
  !   H0=H0/beta*2.d0*spin_degeneracy
  !   !
  !   Tail0=0d0
  !   Tail1=0d0
  !   do ik=1,Lk
  !      Ak= Hk(:,:,ik)
  !      Bk=-Hk(:,:,ik)
  !      Ck= matmul(Ak,Bk)
  !      Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak,No)
  !      Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck,No)
  !   enddo
  !   Tail0=spin_degeneracy*Tail0
  !   Tail1=spin_degeneracy*Tail1*beta
  !   ed_Ekin=H0+Tail0+Tail1
  !   deallocate(wm,Ak,Bk,Ck,Zk,Zeta,Gk,Tk)
  ! end function get_kinetic_energy


  ! function get_local_energy(Hk,Liw) result(ed_Eloc)
  !   integer                                  :: Lk,No,Liw
  !   integer                                  :: i,ik,iorb
  !   complex(8),dimension(:,:,:)              :: Hk ![Nso][Nso][Nk]
  !   real(8),dimension(size(Hk,3))            :: Wtk
  !   !
  !   real(8),dimension(:),allocatable         :: wm
  !   complex(8),dimension(:,:),allocatable    :: Ak,Bk,Hloc
  !   complex(8),dimension(:,:),allocatable    :: Ck,Zk
  !   complex(8),dimension(:,:),allocatable    :: Zeta,Gk,Tk
  !   real(8)                                  :: Tail0,Tail1,spin_degeneracy
  !   !
  !   real(8)                                  :: H0,ed_Eloc
  !   !
  !   No = size(Hk,1)
  !   Lk = size(Hk,3)
  !   if(No/=size(Hk,2))stop "get_kinetic_energy: size(Hk,1)!=size(Hk,2) [Norb_total]"
  !   if(Lk/=size(Wtk))stop "get_kinetic_energy: size(Wtk)!=size(Hk,3) [L_k]"
  !   !
  !   allocate(wm(Liw))
  !   allocate(Ak(No,No),Bk(No,No),Ck(No,No),Zk(No,No),Zeta(No,No),Gk(No,No),Tk(No,No),Hloc(No,No))
  !   !
  !   wm = pi/beta*dble(2*arange(1,Liw)-1)
  !   Wtk=1d0/Lk
  !   !
  !   Hloc=sum(Hk,3)/Lk
  !   H0=0d0
  !   Zk=0d0 ; forall(i=1:No)Zk(i,i)=1d0
  !   do ik=1,Lk
  !      Bk=-Hk(:,:,ik)
  !      do i=1,Liw
  !         Gk = (xi*wm(i)+xmu)*Zk(:,:) - Hk(:,:,ik)
  !         select case(No)
  !         case default
  !            call matrix_inverse(Gk)
  !         case(1)
  !            Gk = 1d0/Gk
  !         end select
  !         Tk = Zk(:,:)/(xi*wm(i)) - Bk(:,:)/(xi*wm(i))**2
  !         Ck = matmul(Hloc,Gk - Tk)
  !         H0 = H0 + Wtk(ik)*trace_matrix(Ck,No)
  !      enddo
  !   enddo
  !   spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
  !   H0=H0/beta*2.d0*spin_degeneracy
  !   !
  !   Tail0=0d0
  !   Tail1=0d0
  !   do ik=1,Lk
  !      Ak= Hloc
  !      Bk=-Hk(:,:,ik)
  !      Ck= matmul(Ak,Bk)
  !      Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak,No)
  !      Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck,No)
  !   enddo
  !   Tail0=spin_degeneracy*Tail0
  !   Tail1=spin_degeneracy*Tail1*beta
  !   ed_Eloc=H0+Tail0+Tail1
  !   deallocate(wm,Ak,Bk,Ck,Zk,Zeta,Gk,Tk,Hloc)
  ! end function get_local_energy


  ! function trace_matrix(M,dim) result(tr)
  !   integer                       :: dim
  !   complex(8),dimension(dim,dim) :: M
  !   complex(8) :: tr
  !   integer                       :: i
  !   tr=dcmplx(0d0,0d0)
  !   do i=1,dim
  !      tr=tr+M(i,i)
  !   enddo
  ! end function trace_matrix

end program hm_2b_flat

















