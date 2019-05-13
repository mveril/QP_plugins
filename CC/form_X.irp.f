subroutine form_X(nO,nV,OOVV,t2,X1,X2,X3,X4)

! Form intermediate arrays X's in CCD

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)
  double precision              :: rt2(nO**2,nV**2),rOOVV(nO**2,nV**2),rOOVVOV(nO*nV,nO*nV),rt2OV(nO*nV,nO*nV)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: ij,ab,kl,cd,bd,jl,jb,kd,ia,lc,lb,ld
  integer                       :: contract
  double precision              :: gt2T(nO**2,nO**2)
  double precision              :: gtt2(nv**2,nv**2)
  double precision              :: gt2OV(nO*nV,nO*nV)


! Output variables

  double precision,intent(out)  :: X1(nO,nO,nO,nO) 
  double precision,intent(out)  :: X2(nV,nV) 
  double precision,intent(out)  :: X3(nO,nO) 
  double precision,intent(out)  :: X4(nO,nO,nV,nV)

! Initialization

  do b=1,nV
    do a=1,nV
      ab=contract(a,b,nV)
      do j=1,nO
          jb=contract(j,b,nO)
        do i=1,nO
          ij=contract(i,j,nO)
          ia=contract(i,a,nO)
          rt2(ij,ab)=t2(i,j,a,b)
          rt2OV(ia,jb)=t2(i,j,a,b)
          rOOVV(ij,ab)=OOVV(i,j,a,b)
          rOOVVOV(ia,jb)=OOVV(i,j,a,b)
        end do
      end do
    end do
  end do

  X1(:,:,:,:) = 0d0
  X2(:,:)     = 0d0
  X3(:,:)     = 0d0
  X4(:,:,:,:) = 0d0

! Build X1
! pure intrinsic fortran equivalent
! gt2T = matmul(rOOVV,transpose(rt2))
! dgemm
  call dgemm('N','T',nO**2,nO**2,nV**2,1.d0,rOOVV,nO**2,rt2,nO**2,0.d0,gt2T,nO**2)
  do j=1,nO
    do i=1,nO
     ij=contract(i,j,nO)
      do l=1,nO
        do k=1,nO
          kl=contract(k,l,nO)
          X1(k,l,i,j) = X1(k,l,i,j) + gt2T(kl,ij)
        enddo
      enddo
    enddo
  enddo

! Build X2
! pure intrinsic fortran equivalent
! gt2OV = matmul(-rOOVVOV,-rt2OV))
! dgemm
  call dgemm('N','N',nO*nV,nO*nV,nO*nV,1.d0,-rOOVVOV,nO*nV,-rt2OV,nO*nV,0.d0,gt2OV,nO*nV)
  do l=1,nO
    do c=1,nV
      lc=contract(l,c,nO)
      do b=1,nV
        lb=contract(l,b,nO)
        X2(b,c) = X2(b,c) + gt2OV(lc,lb)
      enddo
    enddo
  enddo

! Build X3

  do l=1,nO
    do j=1,nO
      jl=contract(j,l,nO)
      do k=1,nO
        kl=contract(k,l,nO)
        X3(k,j) = X3(k,j) + gt2T(kl,jl)
      enddo
    enddo
  enddo

! Build X4

  do d=1,nV
    do l=1,nO
      ld=contract(l,d,nO)
      do i=1,nO
        do a=1,nO
          ia=contract(i,a,nO)
          do k=1,nO
            X4(i,l,a,d) = X4(i,l,a,d) + gt2OV(ld,ia) 
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine form_X
