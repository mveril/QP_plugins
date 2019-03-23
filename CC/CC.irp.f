program CC
  implicit none
  BEGIN_DOC
  ! CC
  ! Coupled cluster
  ! 
  ! This program perform coupled cluster calculation
  ! CCD
  ! CCSD
  ! or CCSD(T)
  END_DOC
  provide cc_mode
  select case ( trim(cc_mode) )
    case ( 'CCD' )
    call CCD
    call write_time(6)
    case ( 'CCSD', 'CCSDT' )
    call CCSD
    call write_time(6)
    case default
    print *, "Error method ", trim(cc_mode),"not available"
    print *, "Only CCD, CCSD and CCSDT are available"
  end select
end
