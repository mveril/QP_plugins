function contract(i,j,N) result(ij)
  integer,intent(in) :: i,j,N
  integer :: ij
      ij=i+N*(j-1)
end function
