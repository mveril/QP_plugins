subroutine form_v(nO,nV,X1,X2,X3,X4,t2,v)

! Form quadratic array in CCD

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: X1(nO,nO,nO,nO)
  double precision,intent(in)   :: X2(nV,nV)
  double precision,intent(in)   :: X3(nO,nO)
  double precision,intent(in)   :: X4(nO,nO,nV,nV)

! Local variables

  integer                       :: i,j,ij,k,l,kl
  integer                       :: a,b,ab,c,d
  integer                       :: contract 
  double precision              :: rt2(nO**2,nV**2),rX1(nO**2,nO**2)
  double precision              :: iqxX1Tt2(nO**2,nV**2)

! Output variables

  double precision,intent(out)  :: v(nO,nO,nV,nV)
  v(:,:,:,:) = 0d0

!TODO multithread
! multithread independent task1
  do b=1,nV
    do a=1,nV
      ab=contract(a,b,nV)
      do j=1,nO
        do i=1,nO
          ij=contract(i,j,nO)
          rt2(ij,ab)=t2(i,j,a,b)
        end do
      end do
    end do
  end do
! multithread independent task2
  
  do l=1,nO
    do k=1,nO
      kl=contract(k,l,nO)
      do j=1,nO
        do i=1,nO
          ij=contract(i,j,nO)
          rX1(ij,kl)=X1(i,j,k,l)
        end do
      end do
    end do
  end do
!End multithread

! pure intrinsic fortran equivalent
! X1Tt2 = 25d-2*matmul(transpose(rX1),rt2)
! dgemm
  call dgemm('T','N',nO**2,nV**2,nO**2,0.25d0,rX1,nO**2,rt2,nO**2,0.d0,iqxX1Tt2,nO**2)
  do b=1,nV
    do a=1,nV
      ab=contract(a,b,nV)
      do j=1,nO
        do i=1,nO
          ij=contract(i,j,nO)
          v(i,j,a,b) = v(i,j,a,b) + iqxX1Tt2(ij,ab)
        enddo
      enddo
    enddo
  enddo

  do b=1,nV
    do a=1,nV
      do j=1,nO
        do i=1,nO
          do c=1,nV
            v(i,j,a,b) = v(i,j,a,b) - 0.5d0*(X2(b,c)*t2(i,j,a,c) + X2(a,c)*t2(i,j,c,b))
          enddo
        enddo
      enddo
    enddo
  enddo

  do b=1,nV
    do a=1,nV
      do j=1,nO
        do i=1,nO
          do k=1,nO
            v(i,j,a,b) = v(i,j,a,b) - 0.5d0*(X3(k,j)*t2(i,k,a,b) + X3(k,i)*t2(k,j,a,b))
          enddo
        enddo
      enddo
    enddo
  enddo

  do b=1,nV
    do a=1,nV
      do j=1,nO
        do i=1,nO
          do c=1,nV
            do K=1,nO
              v(i,j,a,b) = v(i,j,a,b) + (X4(i,k,a,c)*t2(j,k,b,c) + X4(i,k,b,c)*t2(k,j,a,c))
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine form_v
