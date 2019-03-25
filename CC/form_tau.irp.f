subroutine form_tau(nO,nV,t1,t2,tau)

! Form tau in CCSD

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

! Output variables

  double precision,intent(out)  :: tau(nO,nO,nV,nV)

  do b=1,nV
    do a=1,nV
      do j=1,nO
        do i=1,nO

          tau(i,j,a,b) = 0.5d0*t2(i,j,a,b) + t1(i,a)*t1(j,b)

        enddo
      enddo
    enddo
  enddo

end subroutine form_tau
