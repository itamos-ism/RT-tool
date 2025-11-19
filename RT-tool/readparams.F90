subroutine readparams
use RTtool_module
character(len=20)::cfile

close(1);open(unit=1,file='paramsRT.dat',status='old')

read(1,*) directory
read(1,*) prefix
filein = adjustl(trim(directory))//'/'//adjustl(trim(prefix))
open(unit=2,file=adjustl(trim(filein))//trim(adjustl('.RTspop.fin')),status='old')
nheader = 0
read(2,*) network
nheader = nheader + 1
if (network.eq.'REDUCED') then 
  nspec=33
  chemsuf = 'reduced.d'
  iH2=31
endif
if (network.eq.'MEDIUM') then
  nspec=77
  chemsuf = 'medium.d'
  iH2=7
endif
if (network.eq.'FULL') then 
  nspec=215
  chemsuf = 'full.d'
  iH2=213
endif
allocate(species(1:nspec))
write(6,*) 'Chemical network = ',network
read(1,*) vturb
vturb = vturb*1d5 !km/s to cm/s
read(1,*) z
Tcmb = 2.725*(1+z)
coo=0
do 
  read(2,'(A)') cfile
  cfile=trim(adjustl(cfile))
  if (cfile=='ENDCOOLFILES') exit
  coo=coo+1
  coolfile(coo)=cfile
  write(6,*) coolfile(coo)
enddo
100 continue
write(6,*) 'Coolants found = ',coo
close(2)
nheader = nheader + coo + 1 !1 for ENDCOOLFILES

return
end subroutine
