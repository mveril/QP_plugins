subroutine form_T(nO,nV,delta_OOOVVV,ub,ubb,EcCCT)

! Compute (T) correction

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV

  double precision,intent(in)   :: delta_OOOVVV(nO,nO,nO,nV,nV,nV)
  double precision,intent(in)   :: ub(nO,nO,nO,nV,nV,nV)
  double precision,intent(in)   :: ubb(nO,nO,nO,nV,nV,nV)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

! Output variables

  double precision,intent(out)  :: EcCCT

  EcCCT = 0d0

  do c=1,nV
    do b=1,nV
      do a=1,nV
        do k=1,nO
          do j=1,nO
            do i=1,nO
              if(delta_OOOVVV(i,j,k,a,b,c) /= 0d0) then
                EcCCT = EcCCT                                &
                    + (ub(i,j,k,a,b,c) + ubb(i,j,k,a,b,c)) &
                    * ubb(i,j,k,a,b,c)/delta_OOOVVV(i,j,k,a,b,c) 
              endif
            end do
          end do
        end do
      end do
    end do
  end do

  EcCCT = - EcCCT/36d0

end subroutine form_T
