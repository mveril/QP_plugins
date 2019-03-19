program CC
  implicit none
  BEGIN_DOC
  ! CC
  ! Coupled cluster
  ! Coupled cluster
  ! This program perform
  END_DOC
  provide cc_mode
  provide mo_two_e_integrals_in_map
  select case ( trim(cc_mode) )
    case ( 'CCD' )
    call CCD
    case ( 'CCSD' )
    call CCSD
    case ( 'CCSDT' )
    call CCSDT
    case default
    print *, "Error noa", trim(cc_mode),"available"
  end select
end
