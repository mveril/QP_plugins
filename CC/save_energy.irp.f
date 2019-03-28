subroutine save_energy(E)
  implicit none
  BEGIN_DOC
! Saves the energy in |EZFIO|.
  END_DOC
  double precision, intent(in) :: E
  call ezfio_set_cc_energy(E)
end
