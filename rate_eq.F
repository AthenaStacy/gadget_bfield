c Low metallicity chemistry network as described in Glover & Jappsen (2007)
#if CHEMISTRYNETWORK == 2
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////             R A T E _ E Q                 \\\\\\\\\\
c
c=======================================================================
c
c  Written by S. Glover (AMNH, 2004-2005; AIP, 2006-2008; ITA 2008)
c
c  REFERENCES
c
c    BHRB04 -- Bergin et al, 2004, ApJ, 612, 921
c
      subroutine rate_eq(neq, time, y, ydot, rpar, ipar)
#ifdef CHEMCOOL
      implicit NONE
#include "cool.h"
#include "non_eq.h"
#include "shield_data.h"
      integer neq
c
      REAL time, y(nspec), ydot(nspec), rpar(nrpar)
      integer ipar(nipar)
      REAL phrates(npr)
c
      integer niter
      parameter (niter = 10)
      REAL eqtol
      parameter (eqtol = 1d-3)
c
      REAL    dtemp
      integer itemp, I
c
      REAL abundances(nabn), rates(nrates), rates_cmb(nrates),
     $     rates_chem(nrates_chem)
c
      REAL ylam, ylam_chem, ylam_cmb, CMB_temp
c
#ifdef METALS_TG
      REAL metallicity, metal_cool
#endif
c
      REAL temp, yn, dl, divv, energy, gamma, yntot,
     $     ynsq, NH
c
      REAL abe     , abhp    , abHI    , abh2    , abhm  , 
     $     abh2p   , abcI    , abcII   , absiI   , absiII,
     $     aboI    , aboII   , abhd    , abdp    , abDI  ,
     $     abheI   , abhep   , abhepp  , absipp  , abhm_old,
     $     abh2p_old
c
      REAL num, denom
      REAL comp_equilib
      external comp_equilib

      REAL time_to_equilibrium(num_eqb_species)
c
      REAL ykdh, ykdh2, h_gr, ykdis, 
     $     ncrinv, h2var0, h2var1, h2var2, G_dust, phi,
     $     hgrvar1, hgrvar2, cii_gr, c2grvar1, c2grvar2, si2_gr,
     $     si2grvar1, si2grvar2, fa, temp2, tdust2,
     $     oii_gr, dii_gr, heii_gr, hegrvar1, hegrvar2
c
      REAL ph2     , ph3     , ph4     , ph5     , ph6  ,
     $     ph7     , ph8     , ph9     , ph35    , ph38 , 
     $     ph39    , ph40    , ph41    , ph42    , ph43 , 
     $     ph44    , ph45    , ph46    , ph47
c
      REAL cr1     , cr2     , cr3     , cr4     , cr5  ,
     $     cr6     , cr7     , cr1a
c
      REAL ch1     , ch2     , ch3     , ch4     , ch5  ,
     $     ch6     , ch7     , ch8     , ch9     , ch10 ,
     $     ch11    , ch13    , ch14    , ch15    , ch17 , 
     $     ch18    , ch20    , ch21    , ch22    , ch23 ,
     $     ch24    , ch25    , ch26    , ch27    , ch28 ,
     $     ch29    , ch30    , ch31    , ch32    , ch33 ,
     $     ch34    ,
     $     ch35    , ch36    , ch37    , ch38    , ch39 ,
     $     ch40    , ch41    , ch42    , ch43    , ch44 ,
     $     ch45    , ch46    , ch47    , ch48    , ch49 ,
     $     ch58    , ch91    , ch94    , ch97    ,
     $     ch95    , ch96    , ch98    , ch99    , ch100, 
     $     ch101   , ch102   , ch103   , ch108   , ch109,
     $     ch110   , ch121   , ch126   , ch127   , ch135, 
     $     ch136
c

      REAL dtch1   , dtch2   , dtch3   , dtch4   , dtch5  ,
     $     dtch6   , dtch7   , dtch8   , dtch9   , dtch10 ,
     $     dtch11  , dtch13  ,
     $     dtch14  , dtch15  , dtch17  , dtch18  , dtch20 ,
     $     dtch21  , dtch22  , dtch23  , dtch24  , dtch25 ,
     $     dtch26  , dtch27  , dtch28  , dtch29  , dtch30 ,
     $     dtch31  , dtch32  , dtch33  , dtch34  ,
     $     dtch35  , dtch36  , dtch37  , dtch38  , dtch39 ,
     $     dtch40  , dtch41  , dtch42  , dtch43  , dtch44 ,
     $     dtch45  , dtch46  , dtch47  , dtch48  , dtch49 ,
     $     dtch58  , dtch91  , dtch94  , dtch97  ,
     $     dtch95  , dtch96  , dtch98  , dtch99  , dtch100,
     $     dtch101 , dtch102 , dtch103 , dtch108 , dtch109,
     $     dtch110 , dtch121 , dtch126 , dtch127 , dtch135,
     $     dtch136
c
      REAL mh
      parameter (mh = PROTONMASS)
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
c Tracked, non-equilibrium species
c
      abhp   = max(0d0, y(ihp))
      abh2   = max(0d0, y(ih2))
      abhd   = max(0d0, y(ihd))
      abdp   = max(0d0, y(idp))
      abhep  = max(0d0, y(ihep))
      abhepp = max(0d0, y(ihepp))
      if (iflag_highT .eq. 1) then
        abcII  = 0d0
        absiII = 0d0
        absipp = 0d0
        aboII  = 0d0
      elseif (iflag_fixed_ion .eq. 1) then
        abcII  = 0d0
        absiII = 0d0
        absipp = 0d0
        aboII  = max(0d0, y(io))
      else
        abcII  = max(0d0, y(ic))
        absiII = max(0d0, y(isi))
        absipp = max(0d0, y(isipp))
        aboII  = max(0d0, y(io))
      endif
c
c Non-equilibrium species computed using conservation laws
c
      abHI  = max(1d0 - 2d0 * abh2 - abhp - abhd, 0d0)
      abHeI = max(abhe - abhep - abhepp, 0d0)
      abDI  = max(abundD - abhd - abdp, 0d0) 
      if (iflag_highT .eq. 1) then
        abcI  = 0d0
        absiI = 0d0
        aboI  = 0d0
        abe   = abhp + abdp + abhep
     $        + 2d0 * abhepp
      elseif (iflag_fixed_ion .eq. 1) then
        abcI  = 0d0
        absiI = 0d0
        aboI  = max(abundo  - aboII,  0d0)
        abe   = abhp + abundc + abundsi + aboII + abdp + abhep
     $        + 2d0 * abhepp
      else
        abcI  = max(abundc  - abcII,  0d0)
        absiI = max(abundsi - absiII - absipp, 0d0)
        aboI  = max(abundo  - aboII,  0d0)
        abe   = abhp + abcII + absiII + aboII + abdp + abhep
     $        + 2d0 * abhepp + 2d0 * absipp
      endif
c
c Internal energy (erg cm^-3)
c
      energy = y(itmp)
c
c Compute temperature
c
c [NB We ignore corrections to yntot & gamma arising from low abundance 
c     molecular species (e.g. H-, H2+), deuterium and metals]
c
      yntot = (1d0 + abhe - abh2 + abe) * yn
      gamma = (5d0 + 5d0 * abhe - 3d0 * abh2 + 5d0 * abe) /
     $        (3d0 + 3d0 * abhe - abh2 + 3d0 * abe)
      temp = (gamma - 1d0) * energy / (yntot * kboltz)
c
c Useful auxiliary functions
c
      ynsq = yn**2
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
      elseif (temp .gt. tmin .and. temp .lt. tmax) then
        itemp = int(dlog10(temp) / dtlog) + 1
        if (itemp .ne. itemp .or. itemp .le. 0
     $     .or. itemp .gt. nmd) then
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
      dtch8  = dtchtab(8, itemp)
      dtch9  = dtchtab(9, itemp)
      dtch10 = dtchtab(10, itemp)
      dtch11 = dtchtab(11, itemp)
c
      dtch13 = dtchtab(13, itemp)
      dtch14 = dtchtab(14, itemp)
      dtch15 = dtchtab(15, itemp)
c
      dtch17 = dtchtab(17, itemp)
      dtch18 = dtchtab(18, itemp)
c
      dtch20 = dtchtab(20, itemp)
      dtch21 = dtchtab(21, itemp)
      dtch22 = dtchtab(22, itemp)
      dtch23 = dtchtab(23, itemp)
      dtch24 = dtchtab(24, itemp)
      dtch25 = dtchtab(25, itemp)
      dtch26 = dtchtab(26, itemp)
      dtch27 = dtchtab(27, itemp)
      dtch28 = dtchtab(28, itemp)
      dtch29 = dtchtab(29, itemp)
      dtch30 = dtchtab(30, itemp)
      dtch31 = dtchtab(31, itemp)
      dtch32 = dtchtab(32, itemp)
      dtch33 = dtchtab(33, itemp)
      dtch34 = dtchtab(34, itemp)
      dtch35 = dtchtab(35, itemp)
      dtch36 = dtchtab(36, itemp)
      dtch37 = dtchtab(37, itemp)
      dtch38 = dtchtab(38, itemp)
      dtch39 = dtchtab(39, itemp)
      dtch40 = dtchtab(40, itemp)
      dtch41 = dtchtab(41, itemp)
      dtch42 = dtchtab(42, itemp)
      dtch43 = dtchtab(43, itemp)
      dtch44 = dtchtab(44, itemp)
      dtch45 = dtchtab(45, itemp)
      dtch46 = dtchtab(46, itemp)
      dtch47 = dtchtab(47, itemp)
      dtch48 = dtchtab(48, itemp)
      dtch49 = dtchtab(49, itemp)
c
      dtch58 = dtchtab(58, itemp)
c
      dtch91 = dtchtab(91, itemp)
c
      dtch94  = dtchtab(94, itemp)
      dtch95  = dtchtab(95, itemp)
      dtch96  = dtchtab(96, itemp)
      dtch97  = dtchtab(97, itemp)
      dtch98  = dtchtab(98, itemp)
      dtch99  = dtchtab(99, itemp)
      dtch100 = dtchtab(100, itemp)
      dtch101 = dtchtab(101, itemp)
      dtch102 = dtchtab(102, itemp)
      dtch103 = dtchtab(103, itemp)
c
      dtch108 = dtchtab(108, itemp)
      dtch109 = dtchtab(109, itemp)
      dtch110 = dtchtab(110, itemp)
c
      dtch121 = dtchtab(121, itemp)
c
      dtch126 = dtchtab(126, itemp)
      dtch127 = dtchtab(127, itemp)
c
      dtch135 = dtchtab(135, itemp)
      dtch136 = dtchtab(136, itemp)
c
      ch1  = chtab(1, itemp) + dtemp * dtch1
      ch2  = chtab(2, itemp) + dtemp * dtch2
      ch3  = chtab(3, itemp) + dtemp * dtch3
      ch4  = chtab(4, itemp) + dtemp * dtch4
      ch5  = chtab(5, itemp) + dtemp * dtch5
      ch6  = chtab(6, itemp) + dtemp * dtch6
      ch7  = chtab(7, itemp) + dtemp * dtch7
      ch8  = chtab(8, itemp) + dtemp * dtch8
      ch9  = chtab(9, itemp) + dtemp * dtch9
      ch10 = chtab(10, itemp) + dtemp * dtch10
      ch11 = chtab(11, itemp) + dtemp * dtch11
c
      ch13 = chtab(13, itemp) + dtemp * dtch13
      ch14 = chtab(14, itemp) + dtemp * dtch14
      ch15 = chtab(15, itemp) + dtemp * dtch15
c
      ch17 = chtab(17, itemp) + dtemp * dtch17
      ch18 = chtab(18, itemp) + dtemp * dtch18
c
      ch20 = chtab(20, itemp) + dtemp * dtch20
      ch21 = chtab(21, itemp) + dtemp * dtch21
      ch22 = chtab(22, itemp) + dtemp * dtch22
      ch23 = chtab(23, itemp) + dtemp * dtch23
      ch24 = chtab(24, itemp) + dtemp * dtch24
      ch25 = chtab(25, itemp) + dtemp * dtch25
      ch26 = chtab(26, itemp) + dtemp * dtch26
      ch27 = chtab(27, itemp) + dtemp * dtch27
      ch28 = chtab(28, itemp) + dtemp * dtch28
      ch29 = chtab(29, itemp) + dtemp * dtch29
      ch30 = chtab(30, itemp) + dtemp * dtch30
      ch31 = chtab(31, itemp) + dtemp * dtch31
      ch32 = chtab(32, itemp) + dtemp * dtch32
      ch33 = chtab(33, itemp) + dtemp * dtch33
      ch34 = chtab(34, itemp) + dtemp * dtch34
      ch35 = chtab(35, itemp) + dtemp * dtch35
      ch36 = chtab(36, itemp) + dtemp * dtch36
      ch37 = chtab(37, itemp) + dtemp * dtch37
      ch38 = chtab(38, itemp) + dtemp * dtch38
      ch39 = chtab(39, itemp) + dtemp * dtch39
      ch40 = chtab(40, itemp) + dtemp * dtch40
      ch41 = chtab(41, itemp) + dtemp * dtch41
      ch42 = chtab(42, itemp) + dtemp * dtch42
      ch43 = chtab(43, itemp) + dtemp * dtch43
      ch44 = chtab(44, itemp) + dtemp * dtch44
      ch45 = chtab(45, itemp) + dtemp * dtch45
      ch46 = chtab(46, itemp) + dtemp * dtch46
      ch47 = chtab(47, itemp) + dtemp * dtch47
      ch48 = chtab(48, itemp) + dtemp * dtch48
      ch49 = chtab(49, itemp) + dtemp * dtch49
c
      ch58 = chtab(58, itemp) + dtemp * dtch58
c
      ch91  = chtab(91, itemp)  + dtemp * dtch91
c
      ch94  = chtab(94, itemp)  + dtemp * dtch94
      ch95  = chtab(95, itemp)  + dtemp * dtch95
      ch96  = chtab(96, itemp)  + dtemp * dtch96
      ch97  = chtab(97, itemp)  + dtemp * dtch97
      ch98  = chtab(98, itemp)  + dtemp * dtch98
      ch99  = chtab(99, itemp)  + dtemp * dtch99
      ch100 = chtab(100, itemp) + dtemp * dtch100
      ch101 = chtab(101, itemp) + dtemp * dtch101
      ch102 = chtab(102, itemp) + dtemp * dtch102
      ch103 = chtab(103, itemp) + dtemp * dtch103
c
      ch108 = chtab(108, itemp) + dtemp * dtch108
      ch109 = chtab(109, itemp) + dtemp * dtch109
      ch110 = chtab(110, itemp) + dtemp * dtch110
c
      ch121 = chtab(121, itemp) + dtemp * dtch121
c
      ch126 = chtab(126, itemp) + dtemp * dtch126
      ch127 = chtab(127, itemp) + dtemp * dtch127
c
      ch135 = chtab(135, itemp) + dtemp * dtch135
      ch136 = chtab(136, itemp) + dtemp * dtch136
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
      cr1  = crtab(1)  ! HI
      cr1a = crtab(1)  ! DI, assumed same as HI
      cr2  = crtab(2)  ! H2
      cr3  = crtab(3)  ! CI
      cr4  = crtab(4)  ! OI
      cr5  = crtab(5)  ! SiI
      cr6  = crtab(6)  ! HeI
      cr7  = crtab(7)  ! Si+
c
      if (dust_to_gas_ratio .gt. 0d0 .and. iflag_highT .eq. 0) then
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
c H+:
c
        if (phi .eq. 0d0) then
          h_gr = 1.225d-13 * dust_to_gas_ratio
        else
          hgrvar1  = 8.074d-6 * phi**1.378d0
          hgrvar2  = (1d0 + ch34 * phi**ch35)
          h_gr     = 1.225d-13 * dust_to_gas_ratio /
     $               (1d0 + hgrvar1 * hgrvar2)
        endif
c
c He+
c
        if (phi .eq. 0d0) then
          heii_gr = 5.572d-14 * dust_to_gas_ratio
        else
          hegrvar1  = 3.185d-7 * phi**1.512d0
          hegrvar2  = (1d0 + ch135 * phi**ch136)
          heii_gr   = 5.572d-14 * dust_to_gas_ratio /
     $               (1d0 + hegrvar1 * hegrvar2)
        endif
c
c D+ and O+
c
        dii_gr = h_gr / dsqrt(2d0)
        oii_gr = h_gr / 4d0
c
c C+:
c
        if (phi .eq. 0d0) then
          cii_gr = 4.558d-13 * dust_to_gas_ratio
        else
          c2grvar1  = 6.089d-3 * phi**1.128d0
          c2grvar2  = (1d0 + ch36 * phi**ch37)
          cii_gr    = 4.558d-13 * dust_to_gas_ratio /
     $                (1d0 + c2grvar1 * c2grvar2)
        endif
c
c Si+:
c
        if (phi .eq. 0d0) then
          si2_gr = 2.166d-14 * dust_to_gas_ratio
        else
          si2grvar1 = 5.678d-8 * phi**1.874d0
          si2grvar2 = (1d0 + ch38 * phi**ch39)
          si2_gr    = 2.166d-14 * dust_to_gas_ratio /
     $                (1d0 + si2grvar1 * si2grvar2)
        endif
c
c H2 formation rate on dust -- we calculate this here, as this allows us to
c vary tdust without needing to call cheminmo again
c
        temp2  = temp  * 1d-2
        tdust2 = tdust * 1d-2
        fa     = 1d0 / (1d0 + 1d4 * dexp(-6d2 / tdust))
        ch7    = ch7 * fa / (1d0 + 0.4d0 * (temp2 + tdust2)**0.5d0 +
     $                       0.2d0 * temp2 + 0.08d0 * temp2**2d0)
c
      else
        h_gr    = 0d0
        heii_gr = 0d0
        dii_gr  = 0d0
        cii_gr  = 0d0
        oii_gr  = 0d0
        si2_gr  = 0d0
        ch7     = 0d0
      endif
c
      NH = 0.5d0 * dl * yn
c
      call calc_photo(temp, rpar, abh2, abhd, 0d0, phrates, G_dust)
c
      ykdis = phrates(1)
      ph2   = phrates(2)
      ph3   = phrates(3)
      ph4   = phrates(4)
      ph5   = phrates(5)
      ph6   = phrates(6)
      ph7   = phrates(7)
      ph8   = phrates(8)
      ph9   = phrates(9)
c
      ph35  = phrates(35)
c
      ph38  = phrates(38)
      ph39  = phrates(39)
      ph40  = phrates(40)
      ph41  = phrates(41)
      ph42  = phrates(42)
      ph43  = phrates(43)
      ph44  = phrates(44)
      ph45  = phrates(45)
      ph46  = phrates(46)
      ph47  = phrates(47)
c
#ifdef RAYTRACE_TG
      if(ray_flag_sun .eq. 0) then
        if(ray_H_coeff .gt. 0d0) then
          ph38  = 4.55d-7/ray_H_coeff/ray_H_coeff ! HI ion
          ph39  = 4.18d-7/ray_H_coeff/ray_H_coeff ! HeI ion
          ph42  = 4.04d-18/ray_H_coeff/ray_H_coeff ! HI heat
          ph43  = 4.14d-18/ray_H_coeff/ray_H_coeff ! HeI heat
        endif

        if(ray_He_coeff .gt. 0d0) then
          ph40  = 6.73d-9/ray_He_coeff/ray_He_coeff ! HeII ion
          ph44  = 7.21d-20/ray_He_coeff/ray_He_coeff ! HeII heat
        endif

        if(ray_LW_coeff .gt. 0d0) then
          ykdis = 1.28d-7/ray_LW_coeff/ray_LW_coeff ! H2  diss
          ph35 = 1.28d-7/3.0/ray_LW_coeff/ray_LW_coeff ! HD diss
        endif
      else if(ray_flag_sun .eq. 1) then
        if(ray_H_coeff .gt. 0d0) then
          ph38  = 1.32d-6/ray_H_coeff/ray_H_coeff ! HI ion
          ph39  = 1.43d-6/ray_H_coeff/ray_H_coeff ! HeI ion
          ph42  = 1.28d-17/ray_H_coeff/ray_H_coeff ! HI heat
          ph43  = 1.57d-17/ray_H_coeff/ray_H_coeff ! HeI heat
        endif

        if(ray_He_coeff .gt. 0d0) then
          ph40  = 3.72d-8/ray_He_coeff/ray_He_coeff ! HeII ion
          ph44  = 4.46d-19/ray_He_coeff/ray_He_coeff ! HeII heat
        endif

        if(ray_LW_coeff .gt. 0d0) then
          ykdis = 3.38d-7/ray_LW_coeff/ray_LW_coeff ! H2  diss
          ph35 = 3.38d-7/3.0/ray_LW_coeff/ray_LW_coeff ! HD diss
        endif
      else if(ray_flag_sun .eq. 2) then
        if(ray_H_coeff .gt. 0d0) then
          ph38  = 3.69d-6/ray_H_coeff/ray_H_coeff ! HI ion
          ph39  = 4.29d-6/ray_H_coeff/ray_H_coeff ! HeI ion
          ph42  = 3.74d-17/ray_H_coeff/ray_H_coeff ! HI heat
          ph43  = 4.94d-17/ray_H_coeff/ray_H_coeff ! HeI heat
        endif

        if(ray_He_coeff .gt. 0d0) then
          ph40  = 1.36d-7/ray_He_coeff/ray_He_coeff ! HeII ion
          ph44  = 1.71d-18/ray_He_coeff/ray_He_coeff ! HeII heat
        endif

        if(ray_LW_coeff .gt. 0d0) then
          ykdis = 9.07d-7/ray_LW_coeff/ray_LW_coeff ! H2  diss
          ph35 = 9.07d-7/3.0/ray_LW_coeff/ray_LW_coeff ! HD diss
        endif
      else
        print*, "Flag_sun not set! Aborting..."
        ABORT(1)
      endif
#endif
c
c Estimate whether we're optically thick to ionizing radiation;
c if we're not, then we need to set ch18 to zero.

      if (2d0 * NH * abhI .lt. 1.6d17) then
        ch18 = 0d0
      endif
c
c Equilibrium species
c
c Since the equilibrium abundances are coupled, we have to iterate to find 
c accpetable solutions. We iterate niter times (or until all abundances have
c converged).
c
      abhm_old  = 0d0
      abh2p_old = 0d0
c
      abhm  = 0d0
      abh2p = 0d0
#ifdef COLLECT_EQB_DATA
      do I = 1, num_eqb_species
        time_to_equilibrium(I) = 0d0
      enddo
#endif
c
      do I = 1, niter
c
c H-
c
        num   = ch21 * abHI * abe + ch126 * abh2 * abe
        denom = ch22 * abHI + ch26 * abe + ch27 * abHI + ch28 * abhp
     $        + ch29 * abhp + (ph4 / yn) + ch8 * abhep
     $        + ch31 * abheI + 2d0 * ch127 * abh2p
c
        abhm = comp_equilib(num, denom, 1d0, 1, yn, time_to_equilibrium)
        HM = abhm
c
        abe  = max(0d0, abe  - (abhm - abhm_old))
        abHI = max(0d0, abHI - (abhm - abhm_old))
c
c H2+ :
c
        num   = ch23 * abHI * abhp + ch24 * abh2 * abhp 
     $        + ch29 * abhp * abhm + (cr2 * abh2 / yn)
     $        + cst(75) * abhep * abh2 + (ph41 / yn) * abh2
        denom = cst(1) * abHI + ch30 * abe  + (ph5 / yn)
     $        + 2d0 * ch127 * abhm
c
        abh2p = comp_equilib(num, denom, 0.5d0, 2, yn, 
     $                       time_to_equilibrium)
        H2II = abh2p
c
        abe  = max(0d0, abe  + (abh2p - abh2p_old))
        abHI = max(0d0, abHI - 2d0*(abh2p - abh2p_old))
c
c Check for convergence:
c
        if (I .ge. 2) then
          if ((abs(abhm_old - abhm)   .lt. eqtol * abhm)  .and.
     $        (abs(abh2p_old - abh2p) .lt. eqtol * abh2p)) then
            goto 100
          endif
        endif
        abhm_old  = abhm
        abh2p_old = abh2p
      enddo
c
100   continue
c
c Protons
c 
      ydot(ihp) = (ch11 * abe   * abHI          ! Collisional ion: H  + e-
     $          + cst(1)* abh2p * abHI          ! Charge transfer: H  + H2+
     $          -  ch13 * abhp  * abe           ! Gas-phase H recomb.(case B)
     $          -  ch23 * abHI  * abhp          ! H2+ formation:   H  + H+
     $          -  ch24 * abh2  * abhp          ! Charge transfer: H2 + H+
     $          -  ch28 * abhp  * abhm          ! Neutralization:  H- + H+
     $          -  ch29 * abhp  * abhm          ! H2+ formation:   H- + H+
     $          -  h_gr * abhp  * 1d0           ! Grain H recomb.
     $          +  ch32 * abhep * abHI          ! CT: He+ + H
     $          -  ch33 * abheI * abhp          ! CT: He  + H+
     $          +  ch9  * abhep * abh2          ! Dissoc. CT: He+ + H2
     $          + ch18  * abhep * abe           ! Ion. from He+ rec: He+ + e
     $          +  ch46 * aboII * abHI          ! Charge transfer: H  + O+
     $          -  ch47 * aboI  * abhp          ! Charge transfer: H+ + O
     $          -  ch48 * abcI  * abhp          ! Charge transfer: H+ + C
     $          +  ch58 * abcII * abHI          ! Charge transfer: C+ + H
     $          -  ch49 * absiI * abhp          ! Charge transfer: H+ + Si
     $          -  ch96 * absiII * abhp         ! Charge transfer: Si+ + H+
     $          +  ch95 * absipp * abHI         ! Charge transfer: Si++ + H
     $          -  ch98 * abhp * abDI           ! Charge transfer: D + H+
     $          -  ch91 * abhp * abDI           ! HD+ formation:   H+ + D 
     $          +  ch99 * abdp * abHI           ! Charge transfer: H + D+
     $          +  ch100* abdp * abh2           ! HD form:         D+ + H2
     $          -  ch101* abhd * abhp) * yn     ! HD dest:         HD + H+
     $          +  ph5  * abh2p                 ! Photodissoc:    H2+ + PH
     $          +  cr1  * abHI                  ! Cosmic ray ion.:  H + CR
     $          +  ph38 * abHI                  ! Photoionization: H + PH
c
c Molecular hydrogen
c
      ydot(ih2) =  (ch108 * abHI**3           ! 3-body H2 form: 3H -> H2 + H 
     $           + ch109  * abHI**2 * abh2    ! "  : 2H + H2 -> H2 + H2
     $           + ch110  * abHI**2 * abheI) * ynsq ! "  : 2H + He -> H2 + He 
     $           + (ch22 * abHI  * abhm          ! H2 formation:        H- + H
     $          + cst(1) * abHI  * abh2p         ! Charge transfer:     H  + H2+
     $          -  ch24  * abhp  * abh2          ! Charge transfer:     H2 + H+
     $          -  ch25  * abe   * abh2          ! Collisional dissoc.: H2 + e- 
     $          -  ykdh  * abHI  * abh2          ! Collisional dissoc.: H2 + H 
     $          -  ykdh2 * abh2  * abh2          ! Collisional dissoc.: H2 + H2
     $          +  ch7   * abHI  * 1d0           ! H2 formation on dust
     $          -  ch100 * abh2  * abdp          ! HD form:             D+ + H2
     $          +  ch101 * abhd  * abhp          ! HD dest:             HD + H+
     $          -  ch102 * abh2  * abDI          ! HD form:             D + H2
     $          +  ch103 * abhd  * abHI          ! HD dest:             HD + H
     $          -  ch9   * abhep * abh2          ! Dissoc. CT: He+ + H2
     $          - cst(75) * abhep * abh2          ! CT: He+ + H2
     $          - ch126  * abh2  * abe           ! H2 + e -> H- + H
     $          + ch127  * abhm  * abh2p) * yn   ! H2+ + H- -> H2 + H
     $          -  ykdis * abh2                  ! Photodissoc:         H2 + PH
     $          -  cr2   * abh2                  ! Cosmic ray ion.:     H2 + CR
     $          -  ph41  * abh2                  ! Photoionization: H2 + PH -> H2+ + e
c
c Singly-ionized carbon (C+)
c
      if (abundc .eq. 0d0 .or. iflag_highT .eq. 1 .or.
     $    iflag_fixed_ion .eq. 1) then
        ydot(ic) = 0d0
      else
        ydot(ic)   = (-ch40 * abcII * abe          ! Gas-phase C  recomb.
     $             - cii_gr * abcII * 1d0          ! Grain C recomb.
     $             + ch43   * abcI  * abe          ! Collisional ion:  C  + e-
     $             + ch48   * abcI  * abhp         ! Charge transfer:  C  + H+
     $             - ch58   * abcII * abHI         ! Charge transfer:  C+ + H
     $             + ch94   * abcI  * abhep        ! Charge transfer:  C + He+
     $             - cst(2) * abcII * absiI) * yn  ! Charge transfer: C+ + Si 
     $             + cr3    * abcI                 ! Cosmic ray ion.:  C + CR
     $             + ph2    * abcI                 ! Photoion:         C + PH
      endif
c
c Singly-ionized silicon (Si+)
c
      if (abundsi .eq. 0d0 .or. iflag_highT .eq. 1 .or.
     $    iflag_fixed_ion .eq. 1) then
        ydot(isi) = 0d0
      else
        ydot(isi) = (-ch41 * absiII * abe          ! Gas-phase Si recomb.
     $            - si2_gr * absiII * 1d0          ! Grain Si recomb.
     $            + ch44   * absiI  * abe          ! Collisional ion:  Si + e-
     $            + cst(2) * abcII  * absiI        ! Charge transfer:  Si + C+
     $            + ch49   * absiI  * abhp         ! Charge transfer:  Si + H+ 
     $            + ch95   * absipp * abHI         ! Charge transfer: Si++ + H
     $            - ch96   * absiII * abhp         ! Charge transfer: Si+ + H+
     $            + cst(52)* absiI  * abhep        ! CT: Si + He+ 
     $            + ch121  * absipp * abe) * yn    ! Rad. recomb: Si++ + e-
     $            + cr5    * absiI                 ! Cosmic ray ion.:  Si + CR
     $            + ph3    * absiI                 ! Photoion:         Si + PH
     $            - cr7    * absiII                ! CR ion: Si+ + CR
     $            - ph47   * absiII                ! Photoion: Si+ + PH -> Si++
      endif
c
c Doubly-ionized silicon (Si++)
c
      if (abundsi .eq. 0d0 .or. iflag_highT .eq. 1 .or.
     $    iflag_fixed_ion .eq. 1) then
        ydot(isipp) = 0d0
      else
        ydot(isipp) = (ch96  * absiII * abhp
     $              -  ch95  * absipp * abHI
     $              -  ch121 * absipp * abe) * yn
     $              +  cr7   * absiII
     $              + ph47   * absiII
      endif
c
c Singly-ionized oxygen (O+) 
c
      if (abundo .eq. 0d0 .or. iflag_highT .eq. 1) then
        ydot(io) = 0d0
      else
        ydot(io) = (-ch42 * aboII * abe           ! Gas-phase O recomb.
     $           - oii_gr * aboII * 1d0           ! Grain surf. OII recomb
     $           -   ch46 * aboII * abHI          ! Charge transfer: H  + O+
     $           +   ch47 * aboI  * abhp          ! Charge transfer: H+ + O
     $           +   ch97 * aboI  * abhep         ! CT: He+ + O
     $           +   ch45 * aboI  * abe) * yn     ! Collisional ion: O + e-
     $           +   cr4  * aboI                  ! Cosmic ray ion:  O + CR
     $           +   ph46 * aboI                  ! Photoionization
      endif
c
c Singly-ionized helium (He+)
c
      if (abhe .eq. 0d0) then
        ydot(ihep) = 0d0
      else
        ydot(ihep) = (ch14    * abheI * abe
     $              - ch15    * abhep * abe
     $              - heii_gr * abhep * 1d0
     $              - ch17    * abhep * abe
     $              + ch20    * abhepp * abe
     $              - ch32    * abhep * abHI
     $              - ch32    * abhep * abDI
     $              + ch33    * abheI * abhp
     $              + ch33    * abheI * abdp
     $              -  ch8    * abhep * abhm
     $              -  ch9    * abhep * abh2
     $              - cst(75) * abhep * abh2
     $              - cst(52) * absiI * abhep
     $              - ch94    * abcI  * abhep
     $              - ch97    * aboI  * abhep) * yn
     $              + cr6     * abheI
     $              + ph39    * abheI
     $              - ph40    * abhep
      endif
c
c Doubly-ionized helium (He++)
c
      if (abhe .eq. 0d0) then
        ydot(ihepp)  = 0d0
      else
        ydot(ihepp) = (ch15 * abhep  * abe
     $              -  ch20 * abhepp * abe) * yn
     $              +  ph40 * abhep
      endif
c
c  Ionized deuterium (D+)
c
      if (abundD .eq. 0d0) then
        ydot(idp) = 0d0
      else
        ydot(idp) = (ch98 * abhp  * abDI            ! Charge transfer: D + H+
     $            -  ch99 * abdp  * abHI            ! Charge transfer: H + D+
     $            -  ch100* abdp  * abh2            ! HD form: D+ + H2
     $            +  ch101* abhd  * abhp            ! HD dest: HD + H+
     $            +  ch32 * abhep * abDI            ! CT: He+ + D
     $            -  ch33 * abheI * abdp            ! CT: He  + D+
     $            -  ch13 * abdp  * abe             ! Case B recomb.
     $            - dii_gr* abdp  * 1d0             ! Grain surf. recomb. 
     $            - ch91  * abHI  * abdp            ! HD+ formation: H + D+ 
     $            +  ch11 * abDI  * abe) * yn       ! Collisional ionization
     $            +  cr1a * abDI                    ! Cosmic ray ion: D + CR
     $            +  ph38 * abDI                    ! Photoionization: H + PH
      endif
c
c  Hydrogen deuteride (HD)
c
c  NB For simplicity, we assume that each photoionization leads to the destruction 
c  of the molecule; this is incorrect when the fractional ionization is small, but
c  we do not expect HD photoionization to be important in those conditions
c
      if (abundD .eq. 0d0) then
        ydot(ihd) = 0d0
      else
        ydot(ihd) = (ch100 * abh2 * abdp            ! HD form: D+ + H2
     $            -  ch101 * abhd * abhp            ! HD dest: HD + H+
     $            +  ch102 * abh2 * abDI            ! HD form: D + H2
     $            -  ch103 * abhd * abHI) * yn      ! HD dest: HD + H
     $            -  ph35  * abhd                   ! Photodissociation
     $            -  ph41  * abhd                   ! Photoionization
      endif
c
c Chemical heating & cooling
c
c
c H2 collisional dissociation cooling:
c
      rates_chem(1) = 4.48d0 * eV * (ykdh  * abHI * abh2 + 
     $                               ykdh2 * abh2 * abh2 +
     $                               ch25  * abe  * abh2 ) * ynsq
c
c H2 destruction by charge transfer with HII -- endothermic by 1.83eV 
c
      rates_chem(2) = 1.83d0 * eV * ch24 * abh2 * abhp * ynsq
c
c H2 photodissociation heating -- 0.4eV per photodissociation:
c
      rates_chem(3) = - 4d-1 * eV * abH2 * yn * ykdis
c
c Heating by UV pumping of H2 -- effective for n > n_crit
c
c We assume 2 eV per UV photon (Burton, Hollenbach & Tielens 1990),
c and adopt a pumping rate 8.5 times larger than the dissociation rate
c (Draine & Bertoldi 1996, table 2). 
c
      rates_chem(4) = - 2d0 * eV * abH2 * yn * ykdis * 8.5d0 * 
     $                  (1d0 - h2var0)
c
c H2 formation heating -- we assume a similar depedence on n / n_crit to 
c the collisional dissociation rates, although ideally one would like to
c do a more careful, level-dependent treatment.
c
c Note that energy input depends on formation mechanism:
c
c Dust -- 4.48eV 
c H-   -- 3.53eV
c H2+  -- 1.83eV
c 3-body -- 4.48eV
c
      rates_chem(5) = - 4.48d0 * eV * (ch7 * abHI * ynsq)
     $                * (1d0 - h2var0)
c
      rates_chem(6) = - 3.53d0 * eV * ch22 * abHI * abhm * ynsq
     $                * (1d0 - h2var0)
c
      rates_chem(7) = - 1.83d0 * eV * cst(1) * abHI * abh2p * 
     $                ynsq * (1d0 - h2var0)
c
      rates_chem(8) = - 4.48d0 * eV * (ch108 * abHI**3 + 
     $                  ch109 * abHI**2 * abh2 + 
     $                  ch110 * abHI**2 * abheI) * yn**3
c
c HI collisional ionization cooling -- 13.6eV per ionization:
c
      rates_chem(9) = 13.6d0 * eV * ch11 * abe * abHI * ynsq
c
c HeI collisional ionization cooling -- 24.6 eV per ionization
c
      rates_chem(10) = 24.6d0 * eV * ch14 * abe * abheI * ynsq
c
c HeII collisional ionization cooling -- 54.4 eV per ionization
c
      rates_chem(11) = 54.4d0 * eV * ch15 * abe * abhep * ynsq
c
c HII recombination cooling -- we assume that kT ergs is radiated 
c per recombination and that this is independent of whether the
c recombination occurs in the gas phase or on a grain surface
c
      rates_chem(12) = kboltz * temp * abhp * (ch13 * abe + h_gr) * 
     $                 ynsq
c
c HeII recombination cooling. Dielectronic recombination results
c in the loss of ~ 40.8eV per recombination, while for standard 
c radiative recombination, we make the same assumption as for HII
c
c Dielectronic:
      rates_chem(13) = 40.8d0 * eV * ch10 * abhep * abe * ynsq
c
c Radiative:
      rates_chem(14) = kboltz * temp * ((ch17 - ch10) * abe + heii_gr) *
     $                 abhep * ynsq
c
c HeIII recombination cooling -- see note for HII
c
      rates_chem(15) = kboltz * temp * abhepp * ch20 * abe * ynsq
c
      if (iflag_highT .eq. 1 .or. iflag_fixed_ion .eq. 1) then
        rates_chem(16) = 0d0
        rates_chem(17) = 0d0
      else
c
c CI photoionization heating
c
        rates_chem(16) = - ph6 * abcI * yn
c
c SiI photoionization heating
c
        rates_chem(17) = - ph7 * absiI * yn
      endif
c
c H- photodissociation heating
c
      rates_chem(18) = - ph8 * abhm * yn
c
c H2+ photodissociation heating
c
      rates_chem(19) = - ph9 * abh2p * yn
c
c HI photoionization heating
c 
      rates_chem(20) = - ph42 * (abHI + abDI) * yn
c
c HeI photoionization heating
c 
      rates_chem(21) = - ph43 * abheI * yn
c
c HeII photoionization heating
c 
      rates_chem(22) = - ph44 * abhepp * yn
c
c H2 photoionization heating
c
      rates_chem(23) = - ph45 * (abh2 + abhd) * yn
c
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
      abundances(2)  = aboI
      if (iflag_fixed_ion .eq. 1) then
        abundances(6)  = 0d0
        abundances(7)  = abundc
        abundances(8)  = 0d0
        abundances(9)  = abundsi
      else
        abundances(6)  = abcI
        abundances(7)  = abcII
        abundances(8)  = absiI
        abundances(9)  = absiII
      endif
      abundances(10) = abe
      abundances(11) = abhp
      abundances(12) = abHI
      abundances(13) = abhd
      abundances(14) = abheI
      abundances(15) = abhep
      abundances(16) = abhepp
c
      call cool_func(temp, yn, dl, divv, G_dust, abundances, ylam,
     $               rates, 0)
c
#ifdef METALS_TG
      if (abratio_o .gt. 0d0) then
        metallicity = abundo / abratio_o
      elseif (abratio_c .gt. 0d0) then
        metallicity = abundc / abratio_c
      elseif (abratio_si .gt. 0d0) then
        metallicity = abundsi / abratio_si
      else
        metallicity = 0d0
      endif
c
      if (iflag_highT .eq. 1) then
        call sdcool(metallicity, temp, metal_cool)

        if (metal_cool .ne. metal_cool .or. metal_cool .lt. 0d0) then
          print*, 'metal_cool = ', metal_cool
          ABORT(9753)
        endif

        metal_cool = metal_cool * yn**2

        ylam = ylam + metal_cool
      endif
#endif
c
c Now correct for effects of CMB:
c
      if (redshift .gt. 0d0) then
        CMB_temp  = 2.726d0 * (1d0 + redshift)
        call cool_func(CMB_temp, yn, dl, divv, G_dust, abundances, 
     $                 ylam_cmb, rates_cmb, 1)
      else
        ylam_cmb = 0d0
      endif
c
      ydot(itmp) = -(ylam - ylam_cmb + ylam_chem)
c
#ifdef DEBUG_RATE_EQ
      call validate_output(ydot)
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
c=======================================================================
c
c    \\\\\\\\\\        B E G I N   F U N C T I O N        //////////
c    //////////          C O M P _ E Q U I L I B          \\\\\\\\\\
c
c=======================================================================
c
#ifdef CHEMCOOL
      REAL function comp_equilib(num, denom, max_abundance, id, yn, 
     $                           time_to_equilibrium)
      implicit NONE
#include "cool.h"
#include "non_eq.h"
      REAL num, denom, max_abundance, yn
      REAL time_to_equilibrium(num_eqb_species)
      integer id
c
      if (num .eq. 0d0) then
        comp_equilib = 0d0
      elseif (denom .eq. 0d0) then
c If this happens, then our assumption that this species is in
c equilibrium has fatally broken down...
        print*, "Creation without destruction in comp_equilib"
        print*, "Species = ", id
        ABORT(9751)
      else
#ifdef COLLECT_EQB_DATA
        time_to_equilibrium(id) = 1d0 / (denom * yn)
#endif
        comp_equilib = num / denom
        if (comp_equilib .gt. max_abundance) then
          print*, "Error: impossible abundance!"
          print*, "Species ", id
          print*, num, denom, num/denom
          ABORT(9752)
        endif
      endif
c
      return
      end
#endif /* CHEMCOOL */
c=======================================================================
c
c    \\\\\\\\\\          E N D     F U N C T I O N        //////////
c    //////////          C O M P _ E Q U I L I B          \\\\\\\\\\
c
c=======================================================================
#endif /* CHEMISTRYNETWORK == 2 */
