#if CHEMISTRYNETWORK == 8
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////             R A T E _ E Q                 \\\\\\\\\\
c
c=======================================================================
c
c    Written by S. Glover (2007)
c
c  REFERENCES
c
c   BPA06  -- Banerjee, Pudritz & Anderson, 2006, astro-ph/0609428
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
      REAL    dtemp
      integer itemp, I
c
      REAL abundances(nabn), rates(nrates),
     $     rates_chem(nrates_chem), ylam, ylam_chem,
     $     CMB_temp, ylam_cmb, rates_cmb(nrates)
      REAL temp, yn, dl, divv, energy, gamma, yntot, ynsq, ekn
c
      REAL abh2, abhd, ab3M

      REAL gamma_cr, gamma_gd, lambda_dust, kappa, kappa_eff,
     $     sigma_dc, f_dc, L_jeans, rho, cs2

#ifdef FF_COLLAPSE_TEST
      REAL n0, tffc, ndot, compress, tff
#endif
      REAL abHI, abDI, abheI

      REAL ykdh, ykdh2, ykdhe,
     $     ncrinv, h2var0, h2var1, h2var2,
     $     fa, temp2, tdust2

      REAL  ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch102, ch103, 
     $     ch108, ch109, ch110, ch137, ch138, ch139

      REAL  dtch1, dtch2, dtch3, dtch4, dtch5,  dtch6,
     $      dtch7,  dtch102, dtch103, dtch108,
     $      dtch109, dtch110, dtch137, dtch138, dtch139
c
      REAL mh
      parameter (mh  = PROTONMASS)
      REAL G_n
      parameter (G_n = GRAVITY)
      REAL pi
      parameter (pi = PI)
c
c Unpack parameters
c
      yn    = rpar(1)
      dl    = rpar(2)
      divv  = rpar(3)
c
      ynsq = yn**2
c
      call validate_input(y, ipar)
      if (ipar(1) .eq. 1) then
        do I = 1, nspec
          ydot(I) = 0d0
        enddo
c We simply return here because we'll be redoing the step anyway
c (making any further work pointless) and because we may be doing 
c this on account of having negative chemical abundances, which 
c would otherwise cause serious problems later on.
        return
      endif
c
c From this point on, we know that our tracked abundances cannot be negative
c
      abh2   = y(ih2)
      abhd   = y(ihd)
c      tdust  = y(itd)
c
      abHI  = max(0d0, 1d0 - 2d0 * abH2 - abhd)
      abDI  = max(0d0, abundD - abhd)
      abheI = abhe
c
c For some three-body reactions, we only have a total rate, not individual
c rates for collisions with H, H2, He etc. To simplify our handling of 
c these reactions, we define an abundance ab3M that is simply the sum of
c the abundances of the three main colliders:
c
      ab3M = abHI + abh2 + abheI
      UNUSED_PARAM(ab3M)
c
c Internal energy (erg cm^-3)
c
      energy = y(itmp)
c
c Compute current adiabatic index and total number density
c
      yntot = (1d0 + abhe - abh2) * yn
      ekn = energy / (kboltz * yntot)
      call compute_gamma(abh2, ekn, gamma)
c
c Compute local Jeans length (for use in cool_func)
c
      rho = (1d0 + 4d0 * abhe) * mh * yn 
      cs2 = (gamma * (gamma - 1d0) * energy / rho)
      L_jeans  = dsqrt(pi * cs2 / (G_n * rho))
c
c Compute temperature
c
c [NB We ignore corrections to yntot & gamma arising from low abundance 
c     molecular species (e.g. H-, H2+), lithium, deuterium and metals]
c
      temp = (gamma - 1d0) * energy / (yntot * kboltz)
c
c Bounds checking
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
c Rates:
c
       dtch1  = dtchtab(1, itemp)
       dtch2  = dtchtab(2, itemp)
       dtch3  = dtchtab(3, itemp)
       dtch4  = dtchtab(4, itemp)
       dtch5  = dtchtab(5, itemp)
       dtch6  = dtchtab(6, itemp)
       dtch7  = dtchtab(7, itemp)
       dtch102 = dtchtab(102,itemp)
       dtch103 = dtchtab(103,itemp)
       dtch108 = dtchtab(108,itemp)
       dtch109 = dtchtab(109,itemp)
       dtch110 = dtchtab(110,itemp)
       dtch137 = dtchtab(137,itemp)
       dtch138 = dtchtab(138,itemp)
       dtch139 = dtchtab(139,itemp)
c
       ch1  = chtab(1, itemp) + dtemp * dtch1
       ch2  = chtab(2, itemp) + dtemp * dtch2
       ch3  = chtab(3, itemp) + dtemp * dtch3
       ch4  = chtab(4, itemp) + dtemp * dtch4
       ch5  = chtab(5, itemp) + dtemp * dtch5
       ch6  = chtab(6, itemp) + dtemp * dtch6
       ch7  = chtab(7, itemp) + dtemp * dtch7
       ch102 = chtab(102,itemp) + dtemp * dtch102
       ch103 = chtab(103,itemp) + dtemp * dtch103
       ch108 = chtab(108,itemp) + dtemp * dtch108
       ch109 = chtab(109,itemp) + dtemp * dtch109
       ch110 = chtab(110,itemp) + dtemp * dtch110
       ch137 = chtab(137,itemp) + dtemp * dtch137
       ch138 = chtab(138,itemp) + dtemp * dtch138
       ch139 = chtab(139,itemp) + dtemp * dtch139
c
c H2 collisional dissociation rates -- need special handling because of
c density dependence
c
      ncrinv   = (2d0 * abh2 * (ch6 - ch5) + ch5 + abheI * ch139)
      h2var0   = 1d0 / ( 1d0 + yn * ncrinv)
      h2var1   = ch3**h2var0
      h2var2   = ch4**h2var0
      ykdh     = ch1 * h2var1
      ykdh2    = ch2 * h2var2
      if (ch138 .lt. 1d-100) then
        ykdhe = 1d-100
      else
        ykdhe = ch138 * (ch137 / ch138)**h2var0
      endif
c
c H2 formation rate on dust -- we calculate this here, as this allows us to
c vary tdust without needing to call cheminmo again
c                  
      temp2  = temp  * 1d-2
      tdust2 = tdust * 1d-2
      fa     = 1d0 / (1d0 + 1d4 * dexp(-6d2 / tdust))
      ch7    = ch7 * fa / (1d0 + 0.4d0 * dsqrt(temp2 + tdust2) +
     $                     0.2d0 * temp2 + 0.08d0 * temp2 * temp2)
c
c Molecular hydrogen
c
      ydot(ih2) = (ch108  * abHI**3           ! 3-body H2 form: 3H -> H2 + H 
     $           + ch109  * abHI**2 * abh2    ! "  : 2H + H2 -> H2 + H2
     $           + ch110  * abHI**2 * abheI) * ynsq ! "  : 2H + He -> H2 + He 
     $           + (ch7   * abHI    * 1d0     ! H2 formation on dust
     $           -  ykdh  * abHI    * abh2    ! Collisional dissoc.: H2 + H 
     $           -  ykdh2 * abh2    * abh2    ! Collisional dissoc.: H2 + H2
     $           -  ykdhe * abh2    * abheI   ! Collisional dissoc.: H2 + He
     $           - ch102  * abh2    * abDI    ! HD formation:         D + H2
     $           + ch103  * abhd    * abHI) * yn  ! HD destruction:      HD + H
c
c  Hydrogen deuteride (HD)
c
c  NB We follow Flower & Harris (2007) and assume that the three-body rates for
c  HD formation are the same as those for H2 formation
c
      if (abundD .eq. 0d0) then
        ydot(ihd) = 0d0
      else
        ydot(ihd) = (ch108 * abHI
     $            +  ch109 * abh2
     $            +  ch110 * abhe) * abHI * abDI * ynsq 
     $            + (ch102 * abh2 * abDI            ! HD form: D + H2
     $            -  ch103 * abhd * abHI) * yn      ! HD dest: HD + H
      endif
c
c Chemical heating & cooling
c
c
c H2 collisional dissociation cooling:
c
      rates_chem(1) = 4.48d0 * eV * (ykdh  * abHI  * abh2 + 
     $                               ykdh2 * abh2  * abh2 +
     $                               ykdhe * abheI * abh2) * ynsq
c
c H2 formation heating -- we assume a similar dependence on n / n_crit to 
c the collisional dissociation rates, although ideally one would like to
c do a more careful, level-dependent treatment.
c
c Note that energy input depends on formation mechanism:
c
c Dust   -- 4.48eV 
c 3-body -- 4.48eV
c
      rates_chem(2) = - 4.48d0 * eV * (ch7 * abHI * ynsq)
     $                * (1d0 - h2var0)
c
      rates_chem(3) = - 4.48d0 * eV * (ch108 * abHI + 
     $                  ch109 * abh2 + 
     $                  ch110 * abheI) * abHI * abHI * ynsq * yn
      ylam_chem = 0d0
      do I = 1, nrates_chem
        ylam_chem = ylam_chem + rates_chem(I)
      enddo
c
c Initialize all abundances to zero, then just set the non-zero ones
c
      do i = 1, nabn
        abundances(i) = 0d0
      enddo
c
      abundances(1)  = abh2
      abundances(12) = abHI
      abundances(13) = abhd
      abundances(14) = abheI
c
#ifdef ADIABATIC_DENSITY_THRESHOLD
      if (yn .ge. yn_adiabatic) then
        ydot(itmp) = 0d0
#ifdef DEBUG_RATE_EQ
        call validate_output(ydot)
#endif
        return
      endif
#endif
      call cool_func(temp, yn, L_jeans, divv, 0d0, abundances, ylam, 
     $               rates, 0)
c
c Save gas->dust energy transfer rate for use below
c
      gamma_gd = rates(1)
c
c At high z, correct for effects of CMB. (In the local ISM, T_cmb
c is too low for this to be important: other sources of heating
c will dominate at v. low T)
c
      if (redshift .gt. 0d0) then
        CMB_temp  = 2.726d0 * (1d0 + redshift)
        call cool_func(CMB_temp, yn, L_jeans, divv, 0d0, abundances, 
     $                 ylam_cmb, rates_cmb, 1)
      else
        ylam_cmb = 0d0
      endif
c
      ydot(itmp) = -(ylam - ylam_cmb + ylam_chem)
c
#ifdef FF_COLLAPSE_TEST
      ndot = yn**1.5d0 * tffc
c
      compress = ((1.1d0 - abh2) * BOLTZMANN * temp) * ndot
      ydot(itmp) = ydot(itmp) + compress + 
     $             (ndot / yn) * y(itmp)
#endif

c
c Dust temperature
c
      if (dust_to_gas_ratio .gt. 0d0) then
c
c Dust heating by cosmic rays: from BPA06, based on GOL01
c
        gamma_cr = 1.95d-28 * yn
c
c Dust cooling: based on treatment in BPA06 
c
        f_dc     = 1.67d0
        sigma_dc = 6.85d-5

        if (tdust .lt. 2d2) then
          kappa = 1.65d-26 * (tdust / 18.24d0)**2 * yn
        elseif (tdust .lt. 1.5d3) then
          kappa = 1.984d-24 * yn
        else
          kappa = 1.984d-24 * (tdust / 1.5d3)**(-12) * yn
        endif
c
        kappa_eff = min(kappa, f_dc / L_jeans)

        lambda_dust = kappa_eff * sigma_dc * tdust**4

c        ydot(itd) = gamma_cr + gamma_gd - lambda_dust
      else
c        ydot(itd) = 0d0
      endif
c
#ifdef DEBUG_RATE_EQ
      call validate_output(ydot)
#endif

#ifdef FF_COLLAPSE_TEST
      if (itest .gt. 0) then
        tff = 1d0 / (tffc * dsqrt(yn)) 
        do I = 1, num_eqb_species
          print*, I, time_to_equilibrium(I) / tff
        enddo
        itest = 0
      endif
#endif
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
#endif /* CHEMISTRYNETWORK == 8 */
