subroutine CCD

! CCD module

  implicit none

! Input variables (provided by IRP)

  integer                       :: maxSCF
  double precision              :: thresh

  integer                       :: nBas,nEl
  double precision              :: ENuc,ETHF
  double precision,allocatable  :: eHF(:)
  double precision,allocatable  :: ERI(:,:,:,:)

! Local variables

  integer                       :: p,q,r,s,i,j,a,b
  integer                       :: nBas2
  integer                       :: nO
  integer                       :: nV
  integer                       :: nSCF
  double precision              :: Conv
  double precision              :: EcMP2
  double precision              :: get_two_e_integral,u_dot_v
  double precision              :: ECCD,EcCCD
  logical,allocatable           :: socc(:)
  double precision,allocatable  :: seHF(:)
  double precision,allocatable  :: sERI(:,:,:,:)
  double precision,allocatable  :: dbERI(:,:,:,:)

  double precision,allocatable  :: eO(:)
  double precision,allocatable  :: eV(:)
  double precision,allocatable  :: delta_OOVV(:,:,:,:)

  double precision,allocatable  :: OOOO(:,:,:,:)
  double precision,allocatable  :: OOVV(:,:,:,:)
  double precision,allocatable  :: OVOV(:,:,:,:)
  double precision,allocatable  :: VVVV(:,:,:,:)

  double precision,allocatable  :: X1(:,:,:,:)
  double precision,allocatable  :: X2(:,:)
  double precision,allocatable  :: X3(:,:)
  double precision,allocatable  :: X4(:,:,:,:)

  double precision,allocatable  :: u(:,:,:,:)
  double precision,allocatable  :: v(:,:,:,:)

  double precision,allocatable  :: r2(:,:,:,:)
  double precision,allocatable  :: t2(:,:,:,:)

! IRP init  
  maxSCF=cc_n_it_max
  thresh=cc_thresh

  nBas=mo_num
  nEl=elec_num
! ENuc
  ETHF=hf_energy
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

! Hello world

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|          CCD calculation           |'
  write(*,*)'**************************************'
  write(*,*)

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

  allocate(eO(nO),eV(nV))
  allocate(delta_OOVV(nO,nO,nV,nV))

  allocate(socc(nbas2))

  call spatial_to_spin_occ(mo_num,mo_occ,mo_num*2,socc)
  do p=1,nbas2
    if(socc(p)) then
      eO(count(socc(1:p))) = seHF(p)
    else
      eV(count(.not. socc(1:p))) = seHF(p)
    endif
  enddo

  call form_delta_OOVV(nO,nV,eO,eV,delta_OOVV)

  deallocate(seHF)

! Create integral batches
  
  allocate(OOOO(nO,nO,nO,nO),OOVV(nO,nO,nV,nV),OVOV(nO,nV,nO,nV),VVVV(nV,nV,nV,nV))
  do s=1,nbas2
    do r=1,nbas2
      do q=1,nbas2
        do p=1,nbas2
          if(socc(p) .AND. socc(q).AND. socc(r) .AND. socc(s)) then
            OOOO(count(socc(1:p)),count(socc(1:q)),count(socc(1:r)):,count(socc(1:s))) = dbERI(p,q,r,s)
          endif
          if(socc(p) .AND. socc(q) .AND. .NOT. socc(r) .AND. .NOT. socc(s)) then
            OOVV(count(socc(1:p)),count(socc(1:q)),count(.NOT. socc(1:r)):,count( .NOT. socc(1:s))) = dbERI(p,q,r,s)
          endif
          if(socc(p) .AND. .NOT. socc(q) .AND. socc(r) .AND. .NOT. socc(s)) then
            OVOV(count(socc(1:p)), count(.NOT. socc(1:q)),count(socc(1:r)):,count( .NOT. socc(1:s))) = dbERI(p,q,r,s)
          endif
          if( .NOT. socc(p) .AND. .NOT. socc(q) .AND. .NOT. socc(r) .AND. .NOT. socc(s)) then
            VVVV(count(.NOT. socc(1:p)), count(.NOT. socc(1:q)),count(.NOT. socc(1:r)),count( .NOT. socc(1:s))) = dbERI(p,q,r,s)
          endif
        enddo
      enddo
    enddo
  enddo

  deallocate(socc)

  deallocate(dbERI)
 
! MP2 guess amplitudes

  allocate(t2(nO,nO,nV,nV))
  do b=1,nV
    do a=1,nV
      do j=1,nO
        do i=1,nO
          if(delta_OOVV(i,j,a,b) == 0d0) then
            t2(i,j,a,b) = 0d0
          else
            t2(i,j,a,b) = -OOVV(i,j,a,b)/delta_OOVV(i,j,a,b)
          end if
        end do
      end do
    end do
  end do

  EcMP2 = 0.25d0*u_dot_v(pack(OOVV,.true.),pack(t2,.true.),size(OOVV))
  write(*,'(1X,A10,1X,F10.6)') 'Ec(MP2) = ',EcMP2

! Initialization

  allocate(r2(nO,nO,nV,nV),u(nO,nO,nV,nV),v(nO,nO,nV,nV))
  allocate(X1(nO,nO,nO,nO),X2(nV,nV),X3(nO,nO),X4(nO,nO,nV,nV))

  Conv = 1d0
  nSCF = 0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------
  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| CCD calculation                                  |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(CCD)','|','Ec(CCD)','|','Conv','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!   Form linear array

    call form_u(nO,nV,OOOO,VVVV,OVOV,t2,u)

!   Form interemediate arrays

    call form_X(nO,nV,OOVV,t2,X1,X2,X3,X4)

!   Form quadratic array

    call form_v(nO,nV,X1,X2,X3,X4,t2,v)

!   Compute residual
    r2(:,:,:,:) = OOVV(:,:,:,:) + delta_OOVV(:,:,:,:)*t2(:,:,:,:) + u(:,:,:,:) + v(:,:,:,:)

!   Check convergence 

    Conv = maxval(abs(r2(:,:,:,:)))
  
!   Update amplitudes

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

!   Compute correlation energy

    EcCCD = 0.25d0*u_dot_v(pack(OOVV,.true.),pack(t2,.true.),size(OOVV))

!   Dump resultS

    ECCD = ETHF + EcCCD

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nSCF,'|',ECCD,'|',EcCCD,'|',Conv,'|'

  enddo
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

  else
    call save_energy(ECCD)
  endif

end subroutine CCD
