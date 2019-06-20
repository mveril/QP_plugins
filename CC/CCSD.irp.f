subroutine CCSD

! CCSD module

  implicit none

! Input variables (provided by IRP)

  integer                       :: maxscf
  double precision              :: thresh

  logical                       :: doCCSDT
  integer                       :: nBas,nEl
  double precision              :: ERHF
  double precision,allocatable  :: eHF(:)
  double precision,allocatable  :: ERI(:,:,:,:)

! Local variables

  integer                       :: p,q,r,s,i,j,a,b
  double precision              :: start_CCSDT,end_CCSDT,t_CCSDT
  integer                       :: nBas2
  integer                       :: nO
  integer                       :: nV
  integer                       :: nSCF
  double precision              :: Conv
  double precision              :: EcMP2
  double precision              :: ECCSD,EcCCSD
  double precision              :: EcCCT
  double precision              :: get_two_e_integral,u_dot_v 
  double precision,allocatable  :: seHF(:)
  double precision,allocatable  :: sERI(:,:,:,:)
  double precision,allocatable  :: dbERI(:,:,:,:)
  double precision,allocatable  :: delta_OV(:,:)
  double precision,allocatable  :: delta_OOVV(:,:,:,:)

  double precision,allocatable  :: OOOO(:,:,:,:)
  double precision,allocatable  :: OOOV(:,:,:,:)
  double precision,allocatable  :: OVOO(:,:,:,:)
  double precision,allocatable  :: VOOO(:,:,:,:)
  double precision,allocatable  :: OOVV(:,:,:,:)
  double precision,allocatable  :: OVVO(:,:,:,:)
  double precision,allocatable  :: OVVV(:,:,:,:)
  double precision,allocatable  :: VOVV(:,:,:,:)
  double precision,allocatable  :: VVVO(:,:,:,:)
  double precision,allocatable  :: VVVV(:,:,:,:)

  logical,allocatable           :: socc(:)
  double precision,allocatable  :: eO(:)
  double precision,allocatable  :: eV(:)
  double precision,allocatable  :: hvv(:,:)
  double precision,allocatable  :: hoo(:,:)
  double precision,allocatable  :: hvo(:,:)
  double precision,allocatable  :: gvv(:,:)
  double precision,allocatable  :: goo(:,:)
  double precision,allocatable  :: aoooo(:,:,:,:)
  double precision,allocatable  :: bvvvv(:,:,:,:)
  double precision,allocatable  :: hovvo(:,:,:,:)

  double precision,allocatable  :: r1(:,:)
  double precision,allocatable  :: r2(:,:,:,:)

  double precision,allocatable  :: t1(:,:)
  double precision,allocatable  :: t2(:,:,:,:)
  double precision,allocatable  :: tau(:,:,:,:)

! Hello world

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|         CCSD calculation           |'
  write(*,*)'**************************************'
  write(*,*)

! IRP init  
  provide cc_mode
  maxSCF=cc_n_it_max
  thresh=cc_thresh
  if (trim(cc_mode)=='CCSD(T)') then
    doCCSDT=.true.
  else
    doCCSDT=.false.
  end if
  nBas=mo_num
  nEl=elec_num
  ERHF=hf_energy
  allocate(eHF(nBas))
  eHF(:)=fock_matrix_diag_mo(:)
  provide mo_two_e_integrals_in_map
  allocate(ERI(nBas,nBas,nBas,nBas))
  do s=1,nBas
    do r=1,nBas
      do q=1,nBas
        do p=1,nBas 
          ERI(p,q,r,s)=get_two_e_integral(p,q,r,s,mo_two_e_integrals_in_map)
        end do
      end do
    end do
  end do

! Spatial to spin orbitals

  nBas2 = 2*nBas

  allocate(seHF(nBas2),sERI(nBas2,nBas2,nBas2,nBas2))

  call spatial_to_spin_MO_energy(nBas,eHF,nBas2,seHF)
  deallocate(eHF)
  call spatial_to_spin_ERI(nBas,ERI,nBas2,sERI)
  deallocate(ERI)

! Antysymmetrize ERIs

  allocate(dbERI(nBas2,nBas2,nBas2,nBas2))

  call antisymmetrize_ERI(2,nBas2,sERI,dbERI)

  deallocate(sERI)

! Define occupied and virtual spaces

  nO = nEl
  nV = nBas2 - nO

! Form energy denominator  

  allocate(socc(nbas2))
  call spatial_to_spin_occ(mo_num,mo_occ,mo_num*2,socc)
  allocate(eO(nO),eV(nV))
  allocate(delta_OV(nO,nV),delta_OOVV(nO,nO,nV,nV))

  do p=1,nbas2
    if(socc(p)) then
      eO(count(socc(1:p))) = seHF(p)
    else
      eV(count(.not. socc(1:p))) = seHF(p)
    endif
  enddo

  call form_delta_OV(nO,nV,eO,eV,delta_OV)
  call form_delta_OOVV(nO,nV,eO,eV,delta_OOVV)

  deallocate(seHF)

! Create integral batches

  allocate(OOOO(nO,nO,nO,nO),                                     & 
           OOOV(nO,nO,nO,nV),OVOO(nO,nV,nO,nO),VOOO(nV,nO,nO,nO), &
           OOVV(nO,nO,nV,nV),OVVO(nO,nV,nV,nO),                   & 
           OVVV(nO,nV,nV,nV),VOVV(nV,nO,nV,nV),VVVO(nV,nV,nV,nO), & 
           VVVV(nV,nV,nV,nV))

  do s=1,nbas2
    do r=1,nbas2
      do q=1,nbas2
        do p=1,nbas2
          if(socc(p) .AND. socc(q) .AND. socc(r) .AND. socc(s)) then
            OOOO(count(socc(1:p)),count(socc(1:q)),count(socc(1:r)):,count(socc(1:s))) = dbERI(p,q,r,s)
          endif
          if(socc(p) .AND. socc(q) .AND. socc(r) .AND. .NOT. socc(s)) then
            OOOV(count(socc(1:p)),count(socc(1:q)),count(socc(1:r)):,count(.NOT. socc(1:s))) = dbERI(p,q,r,s)
          endif
          if(socc(p) .AND. .NOT. socc(q) .AND. socc(r) .AND. socc(s)) then
            OVOO(count(socc(1:p)),count(.NOT. socc(1:q)),count(socc(1:r)):,count(socc(1:s))) = dbERI(p,q,r,s)
          endif
          if(.NOT. socc(p) .AND. socc(q) .AND. socc(r) .AND. socc(s)) then
            VOOO(count(.NOT. socc(1:p)),count(socc(1:q)),count(socc(1:r)),count(socc(1:s))) = dbERI(p,q,r,s)
          endif
          if(socc(p) .AND. socc(q) .AND. .NOT. socc(r) .AND. .NOT. socc(s)) then
            OOVV(count(socc(1:p)),count(socc(1:q)),count(.NOT. socc(1:r)),count( .NOT. socc(1:s))) = dbERI(p,q,r,s)
          endif
          if(socc(p) .AND. .NOT. socc(q) .AND. .NOT. socc(r) .AND. socc(s)) then
            OVVO(count(socc(1:p)),count(.NOT. socc(1:q)),count(.NOT. socc(1:r)),count(socc(1:s))) = dbERI(p,q,r,s)
          endif
          if(socc(p) .AND. .NOT. socc(q) .AND. .NOT. socc(r) .AND. .NOT. socc(s)) then
            OVVV(count(socc(1:p)),count(.NOT. socc(1:q)),count(.NOT. socc(1:r)), count(.NOT. socc(1:s))) = dbERI(p,q,r,s)
          endif
          if(.NOT. socc(p) .AND. socc(q) .AND. .NOT. socc(r) .AND. .NOT. socc(s)) then
            VOVV(count(.NOT. socc(1:p)),count(socc(1:q)),count(.NOT. socc(1:r)),count(.NOT. socc(1:s))) = dbERI(p,q,r,s)
          endif
          if(.NOT. socc(p) .AND. .NOT. socc(q) .AND. .NOT. socc(r) .AND. socc(s)) then
            VVVO(count(.NOT. socc(1:p)),count(.NOT. socc(1:q)),count(.NOT. socc(1:r)), count(socc(1:s))) = dbERI(p,q,r,s)
          endif
          if(.NOT. socc(p) .AND. .NOT. socc(q) .AND. .NOT. socc(r) .AND. .NOT. socc(s)) then
            VVVV(count(.NOT. socc(1:p)),count(.NOT. socc(1:q)),count(.NOT. socc(1:r)), count(.NOT. socc(1:s))) = dbERI(p,q,r,s)
          endif
        enddo
      enddo
    enddo
  enddo
  deallocate(socc)

  deallocate(dbERI)
 
! MP2 guess amplitudes

  allocate(t1(nO,nV),t2(nO,nO,nV,nV),tau(nO,nO,nV,nV))

  t1(:,:)     = 0d0
  do b=1,nV
    do a=1,nV
      do j=1,nO
        do i=1,nO
          if(delta_OOVV(i,j,a,b)==0d0) then
            t2(i,j,a,b) = 0d0
          else
            t2(i,j,a,b) = -OOVV(i,j,a,b)/delta_OOVV(i,j,a,b)
          end if
       enddo
     enddo
   enddo
 enddo
  call form_tau(nO,nV,t1,t2,tau)

  EcMP2 = 0.5d0*u_dot_v(pack(OOVV,.true.),pack(tau,.true.),size(OOVV))
  write(*,'(1X,A10,1X,F10.6)') 'Ec(MP2) = ',EcMP2

! Initialization

  allocate(hvv(nV,nV),hoo(nO,nO),hvo(nV,nO), &
           gvv(nV,nV),goo(nO,nO), & 
           aoooo(nO,nO,nO,nO),bvvvv(nV,nV,nV,nV),hovvo(nO,nV,nV,nO), &
           r1(nO,nV),r2(nO,nO,nV,nV))

  Conv = 1d0
  nSCF = 0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------
  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| CCSD calculation                                  |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(CCSD)','|','Ec(CCSD)','|','Conv','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!   Scuseria Eqs. (5), (6) and (7)

    call form_h(nO,nV,eO,eV,OOVV,t1,tau,hvv,hoo,hvo)

!   Scuseria Eqs. (9), (10), (11), (12) and (13)

    call form_g(nO,nV,hvv,hoo,VOVV,OOOV,t1,gvv,goo)

    call form_abh(nO,nV,OOOO,OVOO,OOVV,VVVV,VOVV,OVVO,OVVV,t1,tau,aoooo,bvvvv,hovvo)

!   Compute residuals

    call form_r1(nO,nV,OVVO,OVVV,OOOV,hvv,hoo,hvo,t1,t2,tau,r1)

    call form_r2(nO,nV,OOVV,OVOO,OVVV,OVVO,gvv,goo,aoooo,bvvvv,hovvo,t1,t2,tau,r2)

!   Check convergence 

    Conv = max(maxval(abs(r1(:,:))),maxval(abs(r2(:,:,:,:))))

!   Update 
    do a=1,nV
      do i=1,nO
        if(delta_OV(i,a) /= 0d0) then
          t1(i,a)      = t1(i,a)     - r1(i,a)    /delta_OV  (i,a)
        endif
      enddo
    enddo
    do b=1,nV
      do a=1,nV
        do j=1,nO
          do i=1,nO
            if(delta_OOVV(i,j,a,b) /= 0d0) then
              t2(i,j,a,b) =  t2(i,j,a,b) - r2(i,j,a,b)/delta_OOVV(i,j,a,b)
            end if
          end do
        end do
      end do
    end do

    call form_tau(nO,nV,t1,t2,tau)
 
!   Compute correlation energy

    EcCCSD = 0.5d0*u_dot_v(pack(OOVV,.true.),pack(tau,.true.),size(OOVV))

!   Dump results

    ECCSD = ERHF + EcCCSD

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nSCF,'|',ECCSD,'|',EcCCSD,'|',Conv,'|'

  end do
  write(*,*)'----------------------------------------------------'
!------------------------------------------------------------------------
! End of SCF loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  end if

! Deallocate memory

  deallocate(hvv,hoo,hvo,         &
             delta_OV,delta_OOVV, &
             gvv,goo,             &
             aoooo,bvvvv,hovvo,   &
             tau,                 &
             r1,r2)

!------------------------------------------------------------------------
! (T) correction
!------------------------------------------------------------------------
  if(doCCSDT) then
    write(*,*) "Starting (T) calculation"
!   call cpu_time(start_CCSDT)V
    call CCSDT(nO,nV,eO,eV,OOVV,VVVO,VOOO,t1,t2,EcCCT)
!   call cpu_time(end_CCSDT)
    call write_time(6)

!    t_CCSDT = end_CCSDT - start_CCSDT
!    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for (T) = ',t_CCSDT,' seconds'
     write(*,*)

    write(*,*)
    write(*,*)'----------------------------------------------------'
    write(*,*)'                 CCSDT(T) energy                    '
    write(*,*)'----------------------------------------------------'
    write(*,'(1X,A20,1X,F15.10)')' E(CCSD(T))  = ',ECCSD  + EcCCT
    write(*,'(1X,A20,1X,F10.6)') ' Ec(CCSD(T)) = ',EcCCSD + EcCCT
    write(*,*)'----------------------------------------------------'
    write(*,*)

    call save_energy(ECCSD  + EcCCT)

  else

    call save_energy(ECCSD)

  end if

end subroutine CCSD
