subroutine excitation
use RTtool_module

!Calculate excitation temperature & black body emission
do p=1,itot
   do j=1,coo
     do il=1,coolant(j)%cnlev
       do jl=1,coolant(j)%cnlev
         if (jl.ge.il) exit
         if (coolant(j)%A_COEFFS(il,jl).eq.0.0D0) cycle
         if (pdr(p)%coolant(j)%pop(jl).eq.0) then
            pdr(p)%coolant(j)%Tex(il,jl) = 0.
            pdr(p)%coolant(j)%Bnu(il,jl) = 0.
         else
            pdr(p)%coolant(j)%Tex(il,jl) = (HP*coolant(j)%FREQUENCIES(il,jl)/KB) / log(coolant(j)%WEIGHTS(il)*&
                    pdr(p)%coolant(j)%pop(jl)/pdr(p)%coolant(j)%pop(il)/coolant(j)%WEIGHTS(jl))
            pdr(p)%coolant(j)%Bnu(il,jl) = (2.*HP*coolant(j)%FREQUENCIES(il,jl)**3/C**2) / (dexp(HP*&
                    coolant(j)%FREQUENCIES(il,jl)/KB/pdr(p)%coolant(j)%Tex(il,jl))-1.0)
            Bckgr = (2.*HP*coolant(j)%FREQUENCIES(il,jl)**3/C**2) / (dexp(HP*coolant(j)%FREQUENCIES(il,jl)/KB/Tcmb)-1.0)

            pdr(p)%coolant(j)%Bnu(il,jl) = pdr(p)%coolant(j)%Bnu(il,jl) - Bckgr
          endif
        enddo
      enddo
    enddo
enddo

return
end subroutine
