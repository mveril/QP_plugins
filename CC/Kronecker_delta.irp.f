function kronecker_delta(i,j)
  implicit none
  BEGIN_DOC
    ! Apply Kronecker_delta operator
  END_DOC
  integer,intent(in)  :: i,j
  integer :: kronecker_delta
  if (i==j) then
    Kronecker_delta=1
  else
    Kronecker_delta=0
  end if
end function kronecker_delta
