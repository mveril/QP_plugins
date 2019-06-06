subroutine form_u(nO,nV,OOOO,VVVV,OVOV,t2,u)

! Form linear array in CCD

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: OOOO(nO,nO,nO,nO)
  double precision,intent(in)   :: VVVV(nV,nV,nV,nV)
  double precision,intent(in)   :: OVOV(nO,nV,nO,nV)

! Local variables
  integer                       :: contract
  integer                       :: i,j,ij,k,l,kl
  integer                       :: a,b,ab,c,d,cd
  double precision              :: rt2(nO**2,nV**2),rVVVV(nV**2,nV**2),rOOOO(nO**2,nO**2)

  double precision              :: idgt2T(nV**2,nO**2)
  double precision              :: idgTt2(nO**2,nV**2)

! Output variables

  double precision,intent(out)  :: u(nO,nO,nV,nV)

  u(:,:,:,:) = 0d0
  do b=1,nV
    do a=1,nV
      ab=contract(a,b,nV)
      do j=1,nO
        do i=1,nO
          ij=contract(i,j,nO)
          rt2(ij,ab)=t2(i,j,a,b)
        enddo
      enddo
      do c=1,nV
        do d=1,nV
          cd=contract(c,d,nV)
          rVVVV(ab,cd)=VVVV(a,b,c,d)
        enddo
      enddo
    enddo
  enddo
  do l=1,nV
    do k=1,nV
      kl=contract(k,l,nO)
      do j=1,nO
        do i=1,nO
          ij=contract(i,j,nO)
          rOOOO(ij,kl)=OOOO(i,j,k,l)
        enddo
      enddo
    enddo
  enddo
 
! pure intrinsic fortran equivalent
! idgt2T = 0.5d0*matmul(rVVVV,transpose(rt2))
! dgemm
  call dgemm('N','T',nV**2,nO**2,nV**2,0.5d0,rVVVV,nV**2,rt2,nO**2,0.d0,idgt2T,nV**2)
  do b=1,nV
    do a=1,nV
      ab=contract(a,b,nV)
      do j=1,nO
        do i=1,nO
          ij=contract(i,j,nO)
          u(i,j,a,b) = u(i,j,a,b) +idgt2T(ab,ij)
        enddo
      enddo
    enddo
  enddo

! pure intrinsic fortran equivalent
! idgTt2 = 0.5d0*matmul(transpose(rOOOO),rt2)
! dgemm
  call dgemm('T','N',nO**2,nV**2,nO**2,0.5d0,rOOOO,nO**2,rt2,nO**2,0.d0,idgTt2,nO**2)
  do b=1,nV
    do a=1,nV
      ab=contract(a,b,nV)
      do j=1,nO
        do i=1,nO
          ij=contract(i,j,nO)
          u(i,j,a,b) = u(i,j,a,b) + idgTt2(ij,ab)
        enddo
      enddo
    enddo
  enddo

  do b=1,nV
    do a=1,nV
      do j=1,nO
        do i=1,nO
          do c=1,nV
            do k=1,nO
              u(i,j,a,b) = u(i,j,a,b) - OVOV(k,b,j,c)*t2(i,k,a,c) &
                                      + OVOV(k,a,j,c)*t2(i,k,b,c) &
                                      - OVOV(k,a,i,c)*t2(j,k,b,c) &
                                      + OVOV(k,b,i,c)*t2(j,k,a,c)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine form_u
