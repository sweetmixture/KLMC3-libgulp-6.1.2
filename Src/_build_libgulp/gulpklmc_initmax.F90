subroutine gulpklmc_initmax

  use klmc
  use reallocate
  implicit none

  lklmcfreshrun = .true.

  lklmc_maxat            = .true.
  lklmc_maxatloc         = .true.
  lklmc_maxatot          = .true.
  lklmc_maxbond          = .true.
  lklmc_maxbondq         = .true.
  lklmc_maxccspec        = .true.
  lklmc_maxcfg           = .true.
  lklmc_maxconnect       = .true.
  lklmc_maxdef           = .true.
  lklmc_maxeamden        = .true.
  lklmc_maxeamfnspec     = .true.
  lklmc_maxeamspec       = .true.
  lklmc_maxedipspec      = .true.
  lklmc_maxfgrad         = .true.
  lklmc_maxfit           = .true.
  lklmc_maxfor           = .true.
  lklmc_maxfstrain       = .true.
  lklmc_maxgcmcmol       = .true.
  lklmc_maxlambda        = .true.
  lklmc_maxlib           = .true.
  lklmc_maxmcswaps       = .true.
  lklmc_maxmcswapspec    = .true.
  lklmc_maxmctrans       = .true.
  lklmc_maxmol           = .true.
  lklmc_maxnboa          = .true.
  lklmc_maxnboo          = .true.
  lklmc_maxnbopot        = .true.
  lklmc_maxnboq0         = .true.
  lklmc_maxnboq          = .true.
  lklmc_maxnbor          = .true.
  lklmc_maxnboz          = .true.
  lklmc_maxnebreplicatot = .true.
  lklmc_maxnppa          = .true.
  lklmc_maxnpts          = .true.
  lklmc_maxobs           = .true.
  lklmc_maxone           = .true.
  lklmc_maxplanepot      = .true.
  lklmc_maxpot           = .true.
  lklmc_maxqrange        = .true.
  lklmc_maxr1at          = .true.
  lklmc_maxreaxffspec    = .true.
  lklmc_maxreaxffval3    = .true.
  lklmc_maxregion        = .true.
  lklmc_maxsix           = .true.
  lklmc_maxspcellbo      = .true.
  lklmc_maxspcell        = .true.
  lklmc_maxspec          = .true.
  lklmc_maxtdfield       = .true.
  lklmc_maxtempramp      = .true.
  lklmc_maxthb           = .true.
  lklmc_maxtitle         = .true.
  lklmc_maxneighk        = .true.
  lklmc_maxpdfcfg        = .true.
!
! 08/23 wkjee: memory from module 'reallocate'
!
  Wordslo  = 0
  Wordsi2  = 0
  Wordsi4  = 0
  Wordsr4  = 0
  Wordsr8  = 0
  Wordsch  = 0
  Wordsc8  = 0
  Wordsc16 = 0
  PeakMemory = 0

  return
end subroutine
