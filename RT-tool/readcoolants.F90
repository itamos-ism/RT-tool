subroutine readcoolants
use RTtool_module

allocate(coolant(1:coo))
do i = 1,coo
  open(unit=44,file=coolfile(i),status='old')
  read(44,'()') 
  read(44,*) coolant(i)%cname
  read(44,'()')
  read(44,*) coolant(i)%molweight
  read(44,'()')
  read(44,*) cur_nlev, cur_ntemp
  coolant(i)%cnlev = cur_nlev
  coolant(i)%cntemp = cur_ntemp
  close(44)
  allocate(coolant(i)%energies(1:cur_nlev))
  allocate(coolant(i)%weights(1:cur_nlev))
  allocate(coolant(i)%A_COEFFS(1:cur_nlev,1:cur_nlev))
  allocate(coolant(i)%B_COEFFS(1:cur_nlev,1:cur_nlev))
  allocate(coolant(i)%frequencies(1:cur_nlev,1:cur_nlev))
  allocate(coolant(i)%temperatures(1:7,1:cur_ntemp))
  allocate(coolant(i)%H_COL(1:cur_NLEV,1:cur_NLEV,1:cur_NTEMP))
  allocate(coolant(i)%HP_COL(1:cur_nlev,1:cur_nlev,1:cur_ntemp))
  allocate(coolant(i)%EL_COL(1:cur_NLEV,1:cur_NLEV,1:cur_NTEMP))
  allocate(coolant(i)%HE_COL(1:cur_NLEV,1:cur_NLEV,1:cur_NTEMP))
  allocate(coolant(i)%H2_COL(1:cur_NLEV,1:cur_NLEV,1:cur_NTEMP))
  allocate(coolant(i)%PH2_COL(1:cur_NLEV,1:cur_NLEV,1:cur_NTEMP))
  allocate(coolant(i)%OH2_COL(1:cur_NLEV,1:cur_NLEV,1:cur_NTEMP))
  coolant(i)%energies=0.0D0
  coolant(i)%weights=0.0D0
  coolant(i)%A_COEFFS=0.0D0
  coolant(i)%B_COEFFS=0.0D0
  coolant(i)%frequencies=0.0D0
  coolant(i)%temperatures=0.0D0
  coolant(i)%H_COL=0.0D0
  coolant(i)%HP_COL=0.0D0
  coolant(i)%EL_COL=0.0D0
  coolant(i)%HE_COL=0.0D0
  coolant(i)%H2_COL=0.0D0
  coolant(i)%PH2_COL=0.0D0
  coolant(i)%OH2_COL=0.0D0
  call readinput(coolfile(i),coolant(i)%cnlev,coolant(i)%cntemp,&
          coolant(i)%energies,coolant(i)%weights,&
          coolant(i)%A_COEFFS,coolant(i)%B_COEFFS,coolant(i)%frequencies,coolant(i)%temperatures,&
          coolant(i)%H_COL,coolant(i)%HP_COL,coolant(i)%EL_COL,coolant(i)%HE_COL,coolant(i)%H2_COL,&
          coolant(i)%PH2_COL,coolant(i)%OH2_COL)
enddo
return
end subroutine
