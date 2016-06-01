c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                C O _ C O O L              \\\\\\\\\\
c
c=======================================================================
c
      subroutine co_cool(temp, N_co_eff, N_co_eff_vib, co_rot_L0, 
     $                   co_rot_lte, co_rot_alpha, co_rot_n05, 
     $                   co_vib_lte_out)
#ifdef CHEMCOOL
c
c     written by: Simon Glover, AMNH, 2004-2005, AIP, 2006
c
c  PURPOSE: Compute CO cooling rates based on tabulated Neufeld & 
c           Kaufman (1993) and Neufeld, Lepp & Melnick (1995) 
c           cooling functions [see coolinmo.F]
c
c  INPUT VARIABLES: temp, N_co_eff
c
c  OUTPUT VARIABLES: co_rot_L0, co_rot_lte, co_rot_alpha, co_rot_n05,
c                    co_vib_lte
c
c--------------------------------------------------------------------
c
      implicit none
#include "cool.h"
c
      REAL       temp, dtemp_co, dN_co, N_co_eff, N_co_eff_vib
      integer    itemp_co, iN_co
c
      REAL       co_rot_L0, co_rot_lte, co_rot_n05,
     $           co_rot_alpha, co_rot_lte_1, co_rot_lte_2,
     $           co_rot_alp_1, co_rot_alp_2, co_rot_n05_1,
     $           co_rot_n05_2, co_vib_lte_1, co_vib_lte_2,
     $           co_vib_lte_out
c
c CO rotational cooling
c
      if (temp .le. co_temptab(1)) then
        itemp_co = 1
        dtemp_co = 0d0
      elseif (temp .ge. co_temptab(nTco)) then
        itemp_co = nTco
        dtemp_co = 0d0
      else
        itemp_co = int(temp) - 9    ! XXX: Table currently starts at 10K
        dtemp_co = temp - int(temp)
      endif
c
c For column densities that do not lie within the region covered by the
c NK93 or NLM95 data, we use the smallest or largest of the tabulated 
c values, as appropriate.
c
      if (N_co_eff .le. co_colntab(1)) then
        iN_co = 1
        dN_co = 0d0
      elseif (N_co_eff .ge. co_colntab(ncdco)) then
        iN_co = ncdco
        dN_co = 0d0
      else
        iN_co = int((10 * N_co_eff) - 144)
        dN_co = (N_co_eff - co_colntab(iN_co)) / 0.1d0
      endif
c
      co_rot_L0 = co_L0(itemp_co) + dtemp_co * dTco_L0(itemp_co)
c
      co_rot_lte_1 =   co_lte(iN_co,itemp_co) + dtemp_co *
     $               dTco_lte(iN_co,itemp_co)
      co_rot_alp_1 =   co_alp(iN_co,itemp_co) + dtemp_co *
     $               dTco_alp(iN_co,itemp_co)
      co_rot_n05_1 =   co_n05(iN_co,itemp_co) + dtemp_co *
     $               dTco_n05(iN_co,itemp_co)
      if (iN_co .eq. ncdco) then
        co_rot_lte   = co_rot_lte_1
        co_rot_alpha = co_rot_alp_1
        co_rot_n05   = co_rot_n05_1
      else
        co_rot_lte_2 =   co_lte(iN_co+1,itemp_co) + dtemp_co *
     $                 dTco_lte(iN_co+1,itemp_co)
        co_rot_alp_2 =   co_alp(iN_co+1,itemp_co) + dtemp_co *
     $                 dTco_alp(iN_co+1,itemp_co)
        co_rot_n05_2 =   co_n05(iN_co+1,itemp_co) + dtemp_co *
     $                 dTco_n05(iN_co+1,itemp_co)
c
        co_rot_lte = co_rot_lte_1 + (co_rot_lte_2 - co_rot_lte_1) *
     $               dN_co
        co_rot_alpha = co_rot_alp_1 + (co_rot_alp_2 - co_rot_alp_1) *
     $                 dN_co
        co_rot_n05 = co_rot_n05_1 + (co_rot_n05_2 - co_rot_n05_1) *
     $               dN_co
      endif
c
c Do final conversion to correct units:
c
      co_rot_L0  = 10d0**(-co_rot_L0)
      co_rot_lte = 10d0**(-co_rot_lte)
      co_rot_n05 = 10d0**(co_rot_n05)
c
c CO vibrational cooling
c
      if (temp .le. co_vib_temptab(1)) then
        itemp_co = 1
        dtemp_co = 0d0
      elseif (temp .ge. co_vib_temptab(nTco_vib)) then
        itemp_co = nTco_vib
        dtemp_co = 0d0
      else
        itemp_co = int(temp) - 99    ! Table starts at 100K
        dtemp_co = temp - int(temp)
      endif
c
      if (N_co_eff_vib .le. co_vib_colntab(1)) then
        iN_co = 1
        dN_co = 0d0
      elseif (N_co_eff_vib .ge. co_vib_colntab(ncdco_vib)) then
        iN_co = ncdco_vib
        dN_co = 0d0
      else
        iN_co = int((10 * N_co_eff) - 129)
        dN_co = (N_co_eff_vib - co_vib_colntab(iN_co)) / 0.1d0
      endif

      co_vib_lte_1 =   co_vib_lte_final(iN_co,itemp_co) + dtemp_co *
     $               dTco_vib_lte(iN_co,itemp_co)

      if (iN_co .eq. ncdco_vib) then
        co_vib_lte_out = co_vib_lte_1
      else
        co_vib_lte_2 =   co_vib_lte_final(iN_co+1,itemp_co) + dtemp_co *
     $                 dTco_vib_lte(iN_co+1,itemp_co)
        co_vib_lte_out = co_vib_lte_1 + (co_vib_lte_2 - co_vib_lte_1) *
     $               dN_co
      endif
c
      co_vib_lte_out = 10d0**(-co_vib_lte_out) * exp(-3.08d3 / temp)
c
#endif /* CHEMCOOL */
c
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////               C O _ C O O L               \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////          H 2 O _ R O T _ C O O L          \\\\\\\\\\
c
c=======================================================================
c
      subroutine h2o_rot_cool(temp, N_h2o_eff, h2o_rot_L0, h2o_rot_lte, 
     $                        h2o_rot_alpha, h2o_rot_n05, iop_flag)
#ifdef CHEMCOOL
c
c     written by: Simon Glover, AMNH, 2004-2005, AIP 2006
c
c  PURPOSE: Compute H2O rotational cooling rates based on tabulated 
c           Neufeld & Kaufman (1993) and Neufeld, Lepp & Melnick (1995) 
c           cooling functions [see coolinmo.F]
c
c  INPUT VARIABLES: temp, N_h2o_eff, iop_flag (0 == para, 1 == ortho)
c
c  OUTPUT VARIABLES: h2o_rot_L0, h2o_rot_lte, h2o_rot_alpha, h2o_rot_n05
c
c--------------------------------------------------------------------
c
      implicit none
#include "cool.h"
c
      REAL     temp, dtemp_h2o, dN_h2o, N_h2o_eff
      integer  itemp_h2o, iN_h2o, iop_flag
c
      REAL     h2o_rot_L0, h2o_rot_lte, h2o_rot_n05,
     $         h2o_rot_alpha, h2o_rot_lte_1, h2o_rot_lte_2,
     $         h2o_rot_alp_1, h2o_rot_alp_2, h2o_rot_n05_1,
     $         h2o_rot_n05_2
c
      if (temp .le. h2o_temptab(1)) then
        itemp_h2o = 1
        dtemp_h2o = 0d0
      elseif (temp .ge. h2o_temptab(nTh2o)) then
        itemp_h2o = nTh2o
        dtemp_h2o = 0d0
      else
        itemp_h2o = int(temp) - 9    ! XXX: Table currently starts at 10K
        dtemp_h2o = temp - int(temp)
      endif
c
c For column densities that do not lie within the region covered by the
c NK93 or NLM95 data, we use the smallest or largest of the tabulated 
c values, as appropriate.
c
      if (N_h2o_eff .le. h2o_colntab(1)) then
        iN_h2o = 1
        dN_h2o = 0d0
      elseif (N_h2o_eff .ge. h2o_colntab(ncdh2o)) then
        iN_h2o = ncdh2o
        dN_h2o = 0d0
      else
        iN_h2o = int((10 * N_h2o_eff) - 144)
        dN_h2o = (N_h2o_eff - h2o_colntab(iN_h2o)) / 0.1d0
      endif
c
      if (iop_flag .eq. 0) then
        h2o_rot_L0 = h2o_L0_para(itemp_h2o) + dtemp_h2o * 
     $             dTh2o_L0_para(itemp_h2o)
        h2o_rot_lte_1 = h2o_lte_para(iN_h2o,itemp_h2o) + 
     $    dtemp_h2o * dTh2o_lte_para(iN_h2o,itemp_h2o)
        h2o_rot_alp_1 = h2o_alp_para(iN_h2o,itemp_h2o) + 
     $    dtemp_h2o * dTh2o_alp_para(iN_h2o,itemp_h2o)
        h2o_rot_n05_1 = h2o_n05_para(iN_h2o,itemp_h2o) + 
     $    dtemp_h2o * dTh2o_n05_para(iN_h2o,itemp_h2o)
        if (iN_h2o .eq. ncdh2o) then
          h2o_rot_lte   = h2o_rot_lte_1
          h2o_rot_alpha = h2o_rot_alp_1
          h2o_rot_n05   = h2o_rot_n05_1
        else
          h2o_rot_lte_2 = h2o_lte_para(iN_h2o+1,itemp_h2o) + 
     $      dtemp_h2o * dTh2o_lte_para(iN_h2o+1,itemp_h2o)
          h2o_rot_alp_2 = h2o_alp_para(iN_h2o+1,itemp_h2o) + 
     $      dtemp_h2o * dTh2o_alp_para(iN_h2o+1,itemp_h2o)
          h2o_rot_n05_2 = h2o_n05_para(iN_h2o+1,itemp_h2o) + 
     $      dtemp_h2o * dTh2o_n05_para(iN_h2o+1,itemp_h2o)
c
          h2o_rot_lte = h2o_rot_lte_1 + (h2o_rot_lte_2 - 
     $                  h2o_rot_lte_1) * dN_h2o
          h2o_rot_alpha = h2o_rot_alp_1 + (h2o_rot_alp_2 - 
     $                    h2o_rot_alp_1) * dN_h2o
          h2o_rot_n05 = h2o_rot_n05_1 + (h2o_rot_n05_2 - 
     $                  h2o_rot_n05_1) * dN_h2o
        endif
      else
        h2o_rot_L0 = h2o_L0_ortho(itemp_h2o) + dtemp_h2o * 
     $             dTh2o_L0_ortho(itemp_h2o)
        h2o_rot_lte_1 = h2o_lte_ortho(iN_h2o,itemp_h2o) + 
     $    dtemp_h2o * dTh2o_lte_ortho(iN_h2o,itemp_h2o)
        h2o_rot_alp_1 = h2o_alp_ortho(iN_h2o,itemp_h2o) + 
     $    dtemp_h2o * dTh2o_alp_ortho(iN_h2o,itemp_h2o)
        h2o_rot_n05_1 = h2o_n05_ortho(iN_h2o,itemp_h2o) + 
     $    dtemp_h2o * dTh2o_n05_ortho(iN_h2o,itemp_h2o)
        if (iN_h2o .eq. ncdh2o) then
          h2o_rot_lte   = h2o_rot_lte_1
          h2o_rot_alpha = h2o_rot_alp_1
          h2o_rot_n05   = h2o_rot_n05_1
        else 
          h2o_rot_lte_2 = h2o_lte_ortho(iN_h2o+1,itemp_h2o) + 
     $        dtemp_h2o * dTh2o_lte_ortho(iN_h2o+1,itemp_h2o)
          h2o_rot_alp_2 = h2o_alp_ortho(iN_h2o+1,itemp_h2o) + 
     $        dtemp_h2o * dTh2o_alp_ortho(iN_h2o+1,itemp_h2o)
          h2o_rot_n05_2 = h2o_n05_ortho(iN_h2o+1,itemp_h2o) + 
     $        dtemp_h2o * dTh2o_n05_ortho(iN_h2o+1,itemp_h2o)
c
          h2o_rot_lte = h2o_rot_lte_1 + (h2o_rot_lte_2 - 
     $                  h2o_rot_lte_1) * dN_h2o
          h2o_rot_alpha = h2o_rot_alp_1 + (h2o_rot_alp_2 - 
     $                    h2o_rot_alp_1) * dN_h2o
          h2o_rot_n05 = h2o_rot_n05_1 + (h2o_rot_n05_2 - 
     $                  h2o_rot_n05_1) * dN_h2o
        endif
      endif
c
c Do final conversion to correct units:
c
      h2o_rot_L0  = 10d0**(-h2o_rot_L0)
      h2o_rot_lte = 10d0**(-h2o_rot_lte)
      h2o_rot_n05 = 10d0**(h2o_rot_n05)
c
#endif /* CHEMCOOL */
c
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////          H 2 O _ R O T _ C O O L          \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////          H 2 O _ V I B _ C O O L          \\\\\\\\\\
c
c=======================================================================
c
      subroutine h2o_vib_cool(temp, N_h2o_eff, h2o_vib_lte_out)
#ifdef CHEMCOOL
c
c     written by: Simon Glover, AMNH, 2004-2005, AIP 2006
c
c  PURPOSE: Compute H2O vibrational cooling rates based on tabulated 
c           Neufeld & Kaufman (1993) and Neufeld, Lepp & Melnick (1995) 
c           cooling functions [see coolinmo.F]
c
c  INPUT VARIABLES: temp, N_h2o_eff
c
c  OUTPUT VARIABLES: h2o_vib_lte_out
c
c--------------------------------------------------------------------
c
      implicit none
#include "cool.h"
c
      REAL     temp, dtemp_h2o, dN_h2o, N_h2o_eff
      integer  itemp_h2o, iN_h2o
c
      REAL     h2o_vib_lte_1, h2o_vib_lte_2,
     $         h2o_vib_lte_out
c
c H2O vibrational cooling
c
      if (temp .le. h2o_vib_temptab(1)) then
        itemp_h2o = 1
        dtemp_h2o = 0d0
      elseif (temp .ge. h2o_vib_temptab(nTh2o_vib)) then
        itemp_h2o = nTh2o_vib
        dtemp_h2o = 0d0
      else
        itemp_h2o = int(temp) - 99    ! Table starts at 100K
        dtemp_h2o = temp - int(temp)
      endif
c
      if (N_h2o_eff .le. h2o_vib_colntab(1)) then
        iN_h2o = 1
        dN_h2o = 0d0
      elseif (N_h2o_eff .ge. h2o_vib_colntab(ncdh2o_vib)) then
        iN_h2o = ncdh2o_vib
        dN_h2o = 0d0
      else
        iN_h2o = int((10 * N_h2o_eff) - 129)
        dN_h2o = (N_h2o_eff - h2o_vib_colntab(iN_h2o)) / 0.1d0
      endif

      h2o_vib_lte_1 =   h2o_vib_lte_final(iN_h2o,itemp_h2o) + 
     $               dtemp_h2o * dTh2o_vib_lte(iN_h2o,itemp_h2o)

      if (iN_h2o .eq. ncdh2o_vib) then
        h2o_vib_lte_out = h2o_vib_lte_1
      else
        h2o_vib_lte_2 =   h2o_vib_lte_final(iN_h2o+1,itemp_h2o) + 
     $                 dtemp_h2o * dTh2o_vib_lte(iN_h2o+1,itemp_h2o)
        h2o_vib_lte_out = h2o_vib_lte_1 + (h2o_vib_lte_2 - 
     $                    h2o_vib_lte_1) * dN_h2o
      endif
c
      h2o_vib_lte_out = 10d0**(-h2o_vib_lte_out) * exp(-2.325d3 / temp)
c
#endif /* CHEMCOOL */
c
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////          H 2 O _ V I B _ C O O L          \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////      T H R E E _ L E V E L _ P O P S      \\\\\\\\\\
c
c=======================================================================
c
      subroutine three_level_pops(r01, r02, r12, r10, r20, r21, 
     $                            n0, n1, n2)
#ifdef CHEMCOOL
      implicit none
#include "cool.h"
      REAL r01, r02, r12, r10, r20, r21
      REAL n0 , n1 , n2
      REAL a1 , a2 , a3 , b1 , b2 , b3
c
c If excitation rates are negligibly small, then we assume that all
c of the atoms are in level 0:
c
      if (r01 .eq. 0d0 .and. r02 .eq. 0d0) then
        n0 = 1d0
        n1 = 0d0
        n2 = 0d0
        return
      endif

      a1 = r01 + r02
      a2 = -r10
      a3 = -r20
      b1 = r01
      b2 = -(r10 + r12)
      b3 = r21
c
      n2 = -a1 * (a1 * b2 - b1 * a2) / ((a1 - a2) * 
     $     (a1 * b3 - b1 * a3) - (a1 - a3) * 
     $     (a1 * b2 - b1 * a2))
c
      n1 = (a1 / (a1 - a2)) - ((a1 - a3) / (a1 - a2)) * n2
c
      n0 = 1d0 - n1 - n2
c
#endif /* CHEMCOOL */
c
      return
      end
c=======================================================================
c
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E      //////////
c    //////////      T H R E E _ L E V E L _ P O P S      \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////          C O M P U T E _ S T I M          \\\\\\\\\\
c
c=======================================================================
c
      subroutine compute_stim(a10, e10, rad_temp, b10)
#ifdef CHEMCOOL
      implicit none
#include "cool.h"
      REAL a10, e10, rad_temp, b10
      REAL x

      x = e10 / (kboltz * rad_temp)
      if (x .lt. 5d0) then
        b10 = a10 / (dexp(x) - 1d0)
      else
        b10 = 0d0
      endif

#endif /* CHEMCOOL */
c
      return
      end
c=======================================================================
c
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E      //////////
c    //////////          C O M P U T E _ S T I M          \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////        L O A D _ H 2 _ T A B L E          \\\\\\\\\\
c
c=======================================================================
c
      subroutine load_H2_table
#ifdef CHEMCOOL
      implicit none
#include "cool.h"
#include "mol_data.h"
      integer I, J

      open(12, file='H2-cooling-ratios.dat', status='old')
      do I = 1, nh2op
        do J = 1, nh2op
          read(12,*) h2_opac_temp(I), h2_opac_column(J), h2_opac(I,J)
        enddo
      enddo
      close (12, status='keep')

#endif /* CHEMCOOL */
c
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////         L O A D _ H 2 _ T A B L E         \\\\\\\\\\
c
c=======================================================================

c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////    C O M P U T E _ H 2 _ O P A C I T Y    \\\\\\\\\\
c
c=======================================================================
c
      subroutine compute_h2_opacity(temp, N_H2_eff, opac)
#ifdef CHEMCOOL
      implicit none
#include "cool.h"
#include "mol_data.h"
      REAL temp, N_H2_eff, opac
      REAL column_min, column_max, logN, logT, diff, dN, dT
      REAL opac_tmp(2)
      integer ii, jj

      column_min = 10**(h2_opac_column(1))
      column_max = 10**(h2_opac_column(nh2op))

      if (N_H2_eff .le. column_min) then
        opac = 1d0
        return
      elseif (N_H2_eff .ge. column_max) then
        opac = 0d0
        return
      else
        logN = log10(N_H2_eff)
        jj   = 1 + int(10 * (logN - 17.0))
        diff = h2_opac_column(jj+1) - h2_opac_column(jj)
        dN   = (logN - h2_opac_column(jj)) / diff
c
        logT = log10(temp)
        if (logT .le. h2_opac_temp(1)) then
          ii = 1
          dT = 0d0
        elseif (logT .ge. h2_opac_temp(nh2op)) then  
          ii = nh2op
          dT = 0d0
        else
          ii   = 1 + int((logT - 1.5) / 0.03)
          diff = h2_opac_temp(ii+1) - h2_opac_temp(ii)
          dT   = (logT - h2_opac_temp(ii)) / diff
        endif
c
        if (dT .gt. 0d0) then
          opac_tmp(1) = h2_opac(ii,jj) + dT * (h2_opac(ii+1,jj) 
     $                - h2_opac(ii,jj))
          opac_tmp(2) = h2_opac(ii,jj+1) + dT * (h2_opac(ii+1,jj+1) 
     $                - h2_opac(ii,jj+1))
        else
          opac_tmp(1) = h2_opac(ii,jj)
          opac_tmp(2) = h2_opac(ii,jj+1)
        endif
c
        opac = opac_tmp(1) + dN * (opac_tmp(2) - opac_tmp(1))
        opac = 1d1**opac
      endif
#endif /* CHEMCOOL */
c
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////    C O M P U T E _ H 2 _ O P A C I T Y    \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                S D C O O L                \\\\\\\\\\
c
c=======================================================================
c
      subroutine sdcool(Z, T, cool)
#ifdef CHEMCOOL
      implicit none
#include "cool.h"
#include "sd_metal.h"
      REAL Z, T, lZ, lT, cool
      REAL cool1, cool2, dZ, dR1, dR2, dtemp
      integer idx_temp, offset1, offset2, offset3
c
c Compute temperature index
c
      lT = log10(T)
      if (lT .lt. 4d0) then
c This shouldn't happen, but we want to be robust in case it does...
        cool = 0d0
        return
      elseif (lT .ge. 8.5d0) then
c Above 10**8.5 K, we assume a constant cooling rate.
        idx_temp = 91
      else
        idx_temp = int(1 + 20 * (lT - 4d0))
        dtemp    = (lT - 4d0) - ((idx_temp - 1) * 5d-2)
      endif 
c
c Compute metal index
c
      if (Z .eq. 0d0) then
c No metal cooling if no metals
        cool = 0d0
        return
      endif
c
      lZ = log10(Z)
c
      if (lZ .gt. 0.5) then
        if (idx_temp .eq. 91) then
          cool = sdrates(91)
        else
          cool1 = sdrates(idx_temp)
          cool2 = sdrates(idx_temp + 1)
          cool = cool1 + (cool2 - cool1) * (dtemp / 5d-2)
        endif
c Scale from log Z = +0.5 entry to true metallicity -- assume linear 
c dependence on Z
        cool = 1d1**(cool) * (Z / dsqrt(1d1))
      elseif (lZ .lt. -2d0) then
        if (lZ .lt. -3d0) then
          if (idx_temp .eq. 91) then
            cool = sdrates(637)
          else
            cool1 = sdrates(idx_temp + 546)
            cool2 = sdrates(idx_temp + 547)
            cool = cool1 + (cool2 - cool1) * (dtemp / 5d-2)
          endif
c Scale from log Z = -3.0 entry to true metallicity -- assume linear 
c dependence on Z
          cool = 1d1**(cool) * (Z / 1d-3)
        else
          dZ = lZ + 2d0
          if (idx_temp .eq. 91) then
            cool = sdrates(546) - dZ * (sdrates(637) - sdrates(546))
          else
            dR1   = sdrates(idx_temp + 546) - sdrates(idx_temp + 455)
            cool1 = sdrates(idx_temp + 455) - dZ * dR1
            dR2   = sdrates(idx_temp + 547) - sdrates(idx_temp + 456)
            cool2 = sdrates(idx_temp + 456) - dZ * dR2
            cool  = cool1 + (cool2 - cool1) * (dtemp / 5d-2)
          endif
          cool = 1d1**cool
        endif
      else
        offset1 = 91 * int(1 - lZ / 5d-1)
        offset2 = offset1 + 91
        offset3 = offset1 + 182
        dZ = (lZ - 5d-1 * (1d0 - (offset1 / 91)))
        if (idx_temp .eq. 91) then
          cool = sdrates(offset2)
     $         - (dZ / 5d-1) * (sdrates(offset3) - sdrates(offset2))
        else
          dR1   = sdrates(idx_temp + offset2) 
     $          - sdrates(idx_temp + offset1)
          cool1 = sdrates(idx_temp + offset1) - (dZ / 5d-1) * dR1
          dR2   = sdrates(idx_temp + offset2 + 1) 
     $          - sdrates(idx_temp + offset1 + 1)
          cool2 = sdrates(idx_temp + offset1 + 1) - (dZ / 5d-1) * dR2
          cool  = cool1 + (cool2 - cool1) * (dtemp / 5d-2)
        endif
        cool = 1d1**cool
      endif
#endif /* CHEMCOOL */
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////               S D C O O L                 \\\\\\\\\\
c
c=======================================================================