  subroutine outfile
!
!  Output required restart and graphics files
!
!   8/95 Output for arcfiles added for noenergy calcs
!   8/98 FDF file format added
!   3/01 CIF file format added
!   2/02 Output for Marvin2 added
!   8/02 Output for DLV added
!  11/07 Output for Biograf file added
!  10/08 Length of filenames that can be output extended 
!   2/09 Modified to accommodate new version of FoX and gulp_cml
!   3/11 Output of lammps potentials added
!   3/12 Output of DCD file added
!  10/12 Output of Lammps files added
!   8/13 Ensure that this routine hasn't been called in parallel from
!        a non-io processor
!   4/18 Modified so that file names are not truncated
!   4/18 XSF file added
!   2/22 JDFTx output added
!
!  Conditions of use:
!
!  GULP is available free of charge to academic institutions
!  and non-commerical establishments only. Copies should be
!  obtained from the author only and should not be distributed
!  in any form by the user to a third party without the express
!  permission of the author. This notice applies to all parts
!  of the program, except any library routines which are
!  distributed with the code for completeness. All rights for
!  such routines remain with the original distributor.
!
!  No claim is made that this program is free from errors and
!  no liability will be accepted for any loss or damage that
!  may result. The user is responsible for checking the validity
!  of their results.
!
!  Copyright Curtin University 2022
!
!  Julian Gale, CIC, Curtin University, February 2022
!
  use control
  use dump
  use gulp_files
  use iochannels
  use gulp_cml,      only : lcml, lvcml, cmlfilename
  use parallel,      only : ioproc
  implicit none
!
!  Local variables
!
  logical             :: lout
!
!  If this is not the IO processor then return
!
  if (.not.ioproc) return
!
  lout = .false.
  write(ioout,'(/)')
!
!  Write out restart file
!
  if (idump.gt.0) then
    lout = .true.
    call dumpdur(idump,-1_i4)
    if (dfile(1:1).ne.' ') then
      if (ldumpnooverwrite) then
        write(ioout,'(''  Dump files written with root name '',a)') trim(dfile)
      else
        write(ioout,'(''  Dump file written as '',a)') trim(dfile)
      endif
    else
      write(ioout,'(''  Dump file written on channel '',i3)') idump
    endif
  endif
!
!  Write out additional output files
!
  if (lmarv) then
    lout = .true.
    if (lmarv2) then
      call outmarv2(13_i4)
      if (marvfile(1:1).ne.' ') then
        write(ioout,'(''  Marvin2 file written as '',a)') trim(marvfile)
      else
        write(ioout,'(''  Marvin2 file written on channel 13'')')
      endif
    else
      call outmarv(13_i4)
      if (marvfile(1:1).ne.' ') then
        write(ioout,'(''  Marvin file written as '',a)') trim(marvfile)
      else
        write(ioout,'(''  Marvin file written on channel 13'')')
      endif
    endif
  endif
  if (lthb) then
    lout = .true.
    call outthb(14_i4)
    if (thbfile(1:1).ne.' ') then
      write(ioout,'(''  THBREL file written as '',a)') trim(thbfile)
    else
      write(ioout,'(''  THBREL file written on channel 14'')')
    endif
  endif
  if (lxtl) call outxtl(15_i4)
  if (lxr) call outxr(17_i4)
  if (lcssr) call outcssr(19_i4)
  if (lbio) call outbiograf(23_i4)
  if (lxsf) call outxsf(38_i4)
  if (lfdf) then
    lout = .true.
    call outfdf(11_i4)
    if (fdffile(1:1).ne.' ') then
      write(ioout,'(''  FDF file written as '',a)') trim(fdffile)
    else
      write(ioout,'(''  FDF file written on channel 11'')')
    endif
  endif
  if (lcif) then
    lout = .true.
    call outcif(20_i4)
    if (ciffile(1:1).ne.' ') then
      write(ioout,'(''  CIF file written as '',a)') trim(ciffile)
    else
      write(ioout,'(''  CIF file written on channel 20'')')
    endif
  endif
  if (ldlv) call outdlv(21_i4)
  if (lcml) then
    lout = .true.
    if (cmlfilename(1:1).ne.' ') then
      write(ioout,'(''  XML file (CML format) written as '',a)') trim(cmlfilename)
      if (lvcml) write(ioout,'(''   (XML file was verbose) '')')
    else
      write(ioout,'(''  XML file (CML format) written '')')
      if (lvcml) write(ioout,'(''  (XML file was verbose) '')')
    endif
  endif
  if (lqbo) then
    lout = .true.
    if (qbofile(1:1).ne.' ') then
      write(ioout,'(''  QBO file written as '',a)') trim(qbofile)
    else
      write(ioout,'(''  QBO file written on channel 24'')')
    endif
  endif
  if (ldcd) then
    lout = .true.
    if (dcdfile(1:1).ne.' ') then
      write(ioout,'(''  DCD file written as '',a)') trim(dcdfile)
    else
      write(ioout,'(''  DCD file written on channel 36'')')
    endif
  endif
  if (ljdftx) then
    lout = .true.
    call outjdftx(27_i4)
    if (jdftxfile(1:1).ne.' ') then
      write(ioout,'(''  JDFTx file written as '',a)') trim(jdftxfile)
    else
      write(ioout,'(''  JDFTx file written on channel 27'')')
    endif
  endif
  if (llammpspots) then
    lout = .true.
    call outlammpspots(25_i4)
    if (lammpspotsfile(1:1).ne.' ') then
      write(ioout,'(''  LAMMPS potentials file written as '',a)') trim(lammpspotsfile)
    else
      write(ioout,'(''  LAMMPS potentials file written on channel 25'')')
    endif
  endif
  if (llammps) then
    lout = .true.
    call outlammps(26_i4,27_i4)
    if (lammpsfile(1:1).ne.' ') then
      write(ioout,'(''  LAMMPS file written as '',a)') trim(lammpsfile)
    else
      write(ioout,'(''  LAMMPS file written on channel 26'')')
    endif
  endif
  if (ltrj) then
    lout = .true.
    if (trjfile(1:1).ne.' ') then
      if (ltrjascii) then
        write(ioout,'(''  GULP ASCII trajectory file written as '',a)') trim(trjfile)
      else
        write(ioout,'(''  GULP trajectory file written as '',a)') trim(trjfile)
      endif
    else
      if (ltrjascii) then
        write(ioout,'(''  GULP ASCII trajectory file written on channel 31'')')
      else
        write(ioout,'(''  GULP trajectory file written on channel 31'')')
      endif
    endif
  endif
!
!  Output arc file only if noenergy calc is to be done as outarc is normally called elsewhere
!
  if (larc.and.lnoenergy) then
    call outarc(16_i4,.false.,.false.)
  endif
  if (lxyz.and.lnoenergy) then
    call outxyz(18_i4,.false.,.false.)
  endif
  if (lout) write(ioout,'(/)')
!
  return
  end
