#ifdef CHEMCOOL
c
c Gamma as a function of e/kn, for (e/kn) = 1.5 dex -> 4.5 dex in 0.125 dex increments, for
c fh2 = 0.1, 0.2, 0.3, 0.4, and 0.5
c
      integer nhng
      parameter (nhng = 125)
      REAL hng(nhng)

      DATA hng  /1.62412, 1.62412, 1.62412, 1.62412, 1.62412, 1.62412, 
     $           1.62412, 1.62410, 1.62389, 1.62268, 1.61924, 1.61337, 
     $           1.60644, 1.60007, 1.59517, 1.59182, 1.58970, 1.58842, 
     $           1.58767, 1.58725, 1.58700, 1.58686, 1.58678, 1.58674, 
     $           1.58671, 1.57877, 1.57877, 1.57877, 1.57877, 1.57877, 
     $           1.57877, 1.57877, 1.57876, 1.57852, 1.57697, 1.57193,
     $           1.56243, 1.55031, 1.53864, 1.52916, 1.52237, 1.51793, 
     $           1.51518, 1.51355, 1.51261, 1.51206, 1.51176, 1.51158, 
     $           1.51148, 1.51143, 1.53033, 1.53033, 1.53033, 1.53033, 
     $           1.53033, 1.53033, 1.53033, 1.53033, 1.53016, 1.52879, 
     $           1.52352, 1.51249, 1.49734, 1.48157, 1.46813, 1.45799, 
     $           1.45106, 1.44664, 1.44397, 1.44239, 1.44149, 1.44097, 
     $           1.44067, 1.44051, 1.44041, 1.47848, 1.47848, 1.47848, 
     $           1.47848, 1.47848, 1.47848, 1.47848, 1.47847, 1.47841, 
     $           1.47758, 1.47332, 1.46257, 1.44631, 1.42811, 1.41149, 
     $           1.39833, 1.38885, 1.38257, 1.37865, 1.37631, 1.37495, 
     $           1.37417, 1.37372, 1.37347, 1.37333, 1.42282, 1.42282, 
     $           1.42282, 1.42282, 1.42282, 1.42282, 1.42282, 1.42282, 
     $           1.42281, 1.42246, 1.41990, 1.41128, 1.39590, 1.37711, 
     $           1.35854, 1.34296, 1.33105, 1.32276, 1.31741, 1.31414,
     $           1.31221, 1.31108, 1.31043, 1.31006, 1.30986/

#endif /* CHEMCOOL */
