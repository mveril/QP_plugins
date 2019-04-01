subroutine form_X(nO,nV,OOVV,t2,X1,X2,X3,X4)

! Form intermediate arrays X's in CCD

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)
  double precision  :: rt2(nO**2,nV**2),rOOVV(nO**2,nV**2)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: ij,ab,kl,cd
  integer                       :: contract
  double precision  :: gt2TX1(nO**2,nO**2) 
  double precision  :: rX1(nO**2,nO**2) 

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
        do i=1,nO
          ij=contract(i,j,nO)
          rt2(ij,ab)=t2(i,j,a,b)
          rOOVV(ij,ab)=OOVV(i,j,a,b)
        end do
      end do
    end do
  end do

  X1(:,:,:,:) = 0d0
  X2(:,:)     = 0d0
  X3(:,:)     = 0d0
  X4(:,:,:,:) = 0d0

! Build X1
  gt2TX1 = matmul(rOOVV,transpose(rt2))
  do j=1,nO
    do i=1,nO
     ij=contract(i,j,nO)
      do l=1,nO
        do k=1,nO
          kl=contract(k,l,nO)
          X1(k,l,i,j) = X1(k,l,i,j) + gt2TX1(kl,ij)
        enddo
      enddo
    enddo
  enddo

! Build X2

  do c=1,nV
    do b=1,nV
      do l=1,nO
        do k=1,nO
          do d=1,nV
            X2(b,c) = X2(b,c) + OOVV(k,l,c,d)*t2(k,l,b,d)
          enddo
        enddo
      enddo
    enddo
  enddo

! Build X3

  do j=1,nO
    do k=1,nO
      do l=1,nO
        do d=1,nV
          do c=1,nV
            X3(k,j) = X3(k,j) + OOVV(k,l,c,d)*t2(j,l,c,d)
          enddo
        enddo
      enddo
    enddo
  enddo

! Build X4

  do d=1,nV
    do a=1,nV
      do l=1,nO
        do i=1,nO
          do k=1,nO
            do c=1,nV
              X4(i,l,a,d) = X4(i,l,a,d) + OOVV(k,l,c,d)*t2(i,k,a,c) 
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine form_X
