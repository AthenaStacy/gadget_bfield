#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////             R A T E _ E Q                 \\\\\\\\\\
c
c=======================================================================
c
      subroutine rate_eq(neq, time, y, ydot, rpar, ipar)
#ifdef CHEMCOOL
      implicit NONE
#include "cool.h"
#include "non_eq.h"
      integer neq
c
      REAL time, y(nspec), ydot(nspec), rpar(nrpar)
      integer ipar(nipar)
c
      integer itemp, i
      REAL abundances(nabn), rates(nrates), 
     $     rates_chem(nrates_chem)

      REAL phrates(npr)
c
      REAL abe     , abhp    , abHI    , abh2    ,
     $     abcp    , abco
c
      REAL ylam    , temp    ,
     $     yn      , dl      , divv    , yntot   , gamma , 
     $     energy  , dtemp   , ncrinv  , h2var0  , h2var1, 
     $     h2var2  , ykdh    , ykdh2   , ykdis   ,
     $     ylam_chem
c
      REAL G_dust

#if CHEMISTRYNETWORK == 5
      REAL k0, k1, gamma_chx, gamma_co, beta
#endif
c
      REAL phi     , h_gr    , hgrvar1 , hgrvar2
c
      REAL tdust2  , temp2   , fa
c
      REAL cr1
c
      REAL ch1     , ch2     , ch3     , ch4     , ch5  ,
     $     ch6     , ch7     , ch11    , ch13    , ch25 ,
     $     ch34    , ch35
c
      REAL dtch1   , dtch2   , dtch3   , dtch4   , dtch5  ,
     $     dtch6   , dtch7   , dtch11  , dtch13  , dtch25 ,
     $     dtch34  , dtch35
c
      REAL mh
      parameter (mh  = PROTONMASS)
c
c Unpack parameters
c
      yn   = rpar(1)
      dl   = rpar(2)
      divv = rpar(3)
c
c Chemical species:
c
      call validate_input(y, ipar)
c
c Tracked, non-equilibrium species
c
      abhp = max(0d0, y(ihp))
      abh2 = max(0d0, y(ih2))
#if CHEMISTRYNETWORK == 5
      abco = max(0d0, y(ico))
#else
      abco = 0d0
#endif
c
c Non-equilibrium species computed using conservation laws
c
      abHI = max(1d0 - 2d0 * abh2 - abhp, 0d0)
#if CHEMISTRYNETWORK == 5
      abcp = max(0d0, abundc - abco)
      abe  = abhp + abcp + abundsi
#else
      UNUSED_PARAM(abcp)
      abe  = abhp + abundc + abundsi
#endif
c
c Internal energy (erg cm^-3)
c
      energy = y(itmp)
c
c Compute temperature
c
c [NB We ignore corrections to yntot & gamma arising from low abundance 
c     molecular species (e.g. H-, H2+) and metals]
c
      yntot = (1d0 + abhe - abh2 + abe) * yn
      gamma = (5d0 + 5d0 * abhe - 3d0 * abh2 + 5d0 * abe) /
     $        (3d0 + 3d0 * abhe - abh2 + 3d0 * abe)
      temp = (gamma - 1d0) * energy / (yntot * kboltz)
c
c Bounds checking
c
c Temp:
c
      if (temp .le. 0d0) then
        itemp   = 1
        dtemp   = 0d0
        temp    = tmin
        ipar(1) = 1
      elseif (temp .le. tmin) then
        itemp = 1
        dtemp = 0d0
      elseif (temp .ge. tmax) then
        itemp = nmd
        dtemp = temp - temptab(itemp)
      else
        itemp = int(dlog10(temp) / dtlog) + 1
        if (itemp .le. 0 .or. itemp .gt. nmd) then
          print*, 'Fatal error in rate_eq.F', itemp, temp
          ABORT(1)
        endif
        dtemp = temp - temptab(itemp)
      endif
c
c Read in coefficients
c
      dtch1  = dtchtab(1, itemp)
      dtch2  = dtchtab(2, itemp)
      dtch3  = dtchtab(3, itemp)
      dtch4  = dtchtab(4, itemp)
      dtch5  = dtchtab(5, itemp)
      dtch6  = dtchtab(6, itemp)
      dtch7  = dtchtab(7, itemp)
      dtch11 = dtchtab(11, itemp)
      dtch13 = dtchtab(13, itemp)
      dtch25 = dtchtab(25, itemp)
      dtch34 = dtchtab(34, itemp)
      dtch35 = dtchtab(35, itemp)
c      
      ch1  = chtab(1, itemp) + dtemp * dtch1
      ch2  = chtab(2, itemp) + dtemp * dtch2
      ch3  = chtab(3, itemp) + dtemp * dtch3
      ch4  = chtab(4, itemp) + dtemp * dtch4
      ch5  = chtab(5, itemp) + dtemp * dtch5
      ch6  = chtab(6, itemp) + dtemp * dtch6
      ch7  = chtab(7, itemp) + dtemp * dtch7
      ch11 = chtab(11, itemp) + dtemp * dtch11
      ch13 = chtab(13, itemp) + dtemp * dtch13
      ch25 = chtab(25, itemp) + dtemp * dtch25
      ch34 = chtab(34, itemp) + dtemp * dtch34
      ch35 = chtab(35, itemp) + dtemp * dtch35
c
c H2 collisional dissociation rates -- need special handling because of
c density dependence
c
      ncrinv   = (2d0 * abh2 * (ch6 - ch5) + ch5)
      h2var0   = 1d0 / ( 1d0 + yn * ncrinv)
      h2var1   = ch3**h2var0
      h2var2   = ch4**h2var0
      ykdh     = ch1 * h2var1
      ykdh2    = ch2 * h2var2
c
c Photodissociation:
c
      call calc_photo(temp, rpar, abh2, 0d0, abco, phrates, G_dust)
      ykdis = phrates(1)
#if CHEMISTRYNETWORK == 5
      gamma_co = phrates(10)
      gamma_chx = 5d0 * phrates(10)
#endif
c
c Rates for recombination on grain surfaces. These rates must be 
c computed here because of their dependence on the electron number
c density. 
c
      if (abe .eq. 0d0) then
c If the fractional ionization is zero, then there won't be any recombination,
c so the value we use for phi doesn't matter too much -- 1d20 is simply an 
c arbitrary large number
c 
        phi = 1d20
      else
        phi = G_dust * sqrt(temp) / (yn * abe)
      endif
c
c HII:
c
      if (phi .eq. 0d0) then
        h_gr = 1.225d-13 * dust_to_gas_ratio
      else
        hgrvar1  = 8.074d-6 * phi**1.378d0
        hgrvar2  = (1d0 + ch34 * phi**ch35)
        h_gr     = 1.225d-13 * dust_to_gas_ratio /
     $             (1d0 + hgrvar1 * hgrvar2)
      endif
c
c H2 formation rate on dust -- we calculate this here, as this allows us to
c vary tdust without needing to call cheminmo again
c                  
      temp2  = temp  * 1d-2
      tdust2 = tdust * 1d-2
      fa     = 1d0 / (1d0 + 1d4 * dexp(-6d2 / tdust))
      ch7    = ch7 * fa / (1d0 + 0.4d0 * (temp2 + tdust2)**0.5d0 +
     $                     0.2d0 * temp2 + 0.08d0 * temp2**2d0)
c
c Cosmic ray ionization rates: precomputed in cheminmo.
c
c HI:
      cr1 = crtab(1)
c
c Protons
c 
      ydot(ihp) = (ch11 * abe   * abHI          ! Collisional ion: H  + e-
     $          +  ch13 * abhp  * abe           ! Gas-phase H recomb.(case B)
     $          -  h_gr * abhp  * 1d0) * yn     ! Grain H recomb.
     $          +  cr1  * abHI                  ! Cosmic ray ion.:  H + CR
c
c Molecular hydrogen
c
      ydot(ih2) = (-ch25  * abe  * abh2  -       ! Collisional dissoc.: H2 + e- 
     $              ykdh  * abHI * abh2  -       ! Collisional dissoc.: H2 + H 
     $              ykdh2 * abh2 * abh2  +       ! Collisional dissoc.: H2 + H2
     $              ch7   * abHI * 1d0 ) * yn -  ! H2 formation on dust
     $              ykdis * abh2                 ! Photodissoc:         H2 + PH
c
#if CHEMISTRYNETWORK == 5
c
c CO -- from NL97
c
      k0 = 5d-16
      k1 = 5d-10
      beta = k1 * (abundo - abco) / (k1 * (abundo - abco) + 
     $       gamma_chx / (abh2 * yn))
c
      ydot(ico) = k0 * abcp * beta * yn**2 - gamma_co * abco * yn
#endif
c
c H2 collisional dissociation cooling:
c
      rates_chem(1) = 4.48d0 * eV * (ykdh  * abHI * abh2 + 
     $                               ykdh2 * abh2 * abh2 +
     $                               ch25  * abe  * abh2 ) * yn**2
c
c H2 photodissociation heating -- 0.4eV per photodissociation:
c
      rates_chem(2) = - 4d-1 * eV * abH2 * yn * ykdis
c
c Heating by UV pumping of H2 -- effective for n > n_crit
c
c We assume 2 eV per UV photon (Burton, Hollenbach & Tielens 1990),
c and adopt a pumping rate 8.5 times larger than the dissociation rate
c (Draine & Bertoldi 1996, table 2). 
c
      rates_chem(3) = - 2d0 * eV * abH2 * yn * ykdis * 8.5d0 * 
     $                  (1d0 - h2var0)
c
c H2 formation heating -- we assume a similar depedence on n / n_crit to 
c the collisional dissociation rates, although ideally one would like to
c do a more careful, level-dependent treatment.
c
c Note that energy input depends on formation mechanism:
c
c Dust -- 4.48eV 
c
      rates_chem(4) = - 4.48d0 * eV * ch7 * abHI * yn**2 
     $                * (1d0 - h2var0)
c
c HI collisional ionization cooling -- 13.6eV per ionization:
c
      rates_chem(5) = 13.6d0 * eV * ch11 * abe * abHI * yn**2
c
c HII recombination cooling -- we assume that kT ergs is radiated 
c per recombination and that this is independent of whether the
c recombination occurs in the gas phase or on a grain surface
c
      rates_chem(6) = kboltz * temp * abhp * (ch13 * abe + h_gr) * yn**2
c
      ylam_chem = rates_chem(1)  + rates_chem(2)  + rates_chem(3)  +
     $            rates_chem(4)  + rates_chem(5)  + rates_chem(6)
c
c Initialize all abundances to zero, then just set the non-zero ones
c
      do i = 1, nabn
        abundances(i) = 0d0
      enddo
c
      abundances(1)  = abh2
#if CHEMISTRYNETWORK == 5
      abundances(2)  = max(0d0, abundo - abco)
      abundances(5)  = abco
      abundances(7)  = abcp
#else
      abundances(2)  = abundo
      abundances(7)  = abundc
#endif
      abundances(9)  = abundsi
      abundances(10) = abe
      abundances(11) = abhp
      abundances(12) = abHI
c
#ifdef ADIABATIC_DENSITY_THRESHOLD
      if (yn .ge. yn_adiabatic) then
        ydot(itmp) = 0d0
        return
      endif
#endif
      call cool_func(temp, yn, dl, divv, G_dust, abundances, ylam,
     $               rates, 0)
c
      ydot(itmp) = -(ylam + ylam_chem)
c
#endif /* CHEMCOOL */
      return      
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////              R A T E _ E Q                \\\\\\\\\\
c
c=======================================================================
#endif /* CHEMISTRYNETWORK == 4 || 5 */
