subroutine spatial_to_spin_occ(nBas,occ,nBas2,socc)

! Convert ERIs from spatial to spin orbitals

  implicit none

! Input variables

  integer,intent(in)            :: nBas,nBas2
  double precision,intent(in)   :: occ(nBas)

! Local variables

  integer                       :: p,sp,occval

! Output variables

  logical,intent(out)  :: socc(nBas2)
          sp=0
          do p=1,nBas
            sp+=1
            occval=nint(occ(p))
            socc(sp)=occval>0
            sp+=1
            socc(sp)=occval==2
          enddo

end subroutine spatial_to_spin_occ
