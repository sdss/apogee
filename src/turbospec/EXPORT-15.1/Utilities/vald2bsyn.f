* 
      program vald2linelist
*
* Convert VALD "show all" (long format) output to bsyn line data lines
*
* Date: Thu, 24 Sep 1998 14:31:01 +0200 (MET DST)
* From: Bengt Edvardsson <be@astro.uu.se>
* To: plez@Ferrum.fysik.lu.se
* Subject: Artikel, ny jonabs, spectrum
* 
* Linjelistan faar nu 2 kolumner (foer ekvivalentbredd och fel) 
* efter "levels" och en full linjeID paa varje rad. 
*
      implicit none
      integer i,j,iel,ion,n,nu(92,2),natom,llower,lupper,nutotal,Nlmax
      integer ii
      parameter (Nlmax=100000)
      real chil,gflog,fdamp,gamrad,eqw,eqwerr,jupper
      real jlow,chiu,landelow,landemean,landeup
      doubleprecision w
      real sunabund(92)
      character*110 array(Nlmax,92,2)
      character*1 lelec(2),blip
      character*2 el,lele(92),cion(2)
      character filename*100,string*100,levels*30,species*5,string2*100
      character tenchar*10,lower*1,upper*1,desl*11,desu*11
*     logical spectrum
*
      data lele /
     &           'H ','HE','LI','BE','B ','C ','N ','O ','F ','NE',
     &           'NA','MG','AL','SI','P ','S ','CL','AR','K ','CA',
     &           'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',
     &           'GA','GE','AS','SE','BR','KR','RB','SR','Y ','ZR',
     &           'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN',
     &           'SB','TE','I ','XE','CS','BA','LA','CE','PR','ND',
     &           'PM','SM','EU','GD','TB','DY','HO','ER','TM','YB',
     &           'LU','HF','TA','W ','RE','OS','IR','PT','AU','HG',
     &           'TL','PB','BI','PO','AT','RN','FR','RA','AC','TH',
     &           'PA','U '/
*
* Solar abundances ref: Grevesse Sauval 1998, Space Science Rev, in press
      data sunabund /
     & 12.00, 10.93,  1.10,  1.40,  2.55,  8.52,  7.92,  8.83,  4.56,   |  1 -  9
     &  8.08,  6.33,  7.58,  6.47,  7.55,  5.45,  7.33,  5.50,  6.40,   | 10 - 18
     &  5.12,  6.36,  3.17,  5.02,  4.00,  5.67,  5.39,  7.50,  4.92,   | 19 - 27
     &  6.25,  4.21,  4.60,  2.88,  3.41,  2.37,  3.41,  2.63,  3.31,   | 28 - 36
     &  2.60,  2.97,  2.24,  2.60,  1.42,  1.92, -99.0,  1.84,  1.12,   | 37 - 45
     &  1.69,  0.94,  1.77,  1.66,  2.00,  1.00,  2.24,  1.51,  2.17,   | 46 - 54
     &  1.13,  2.13,  1.17,  1.58,  0.71,  1.50, -99.0,  1.01,  0.51,   | 55 - 63
     &  1.12, -0.10,  1.14,  0.26,  0.93,  0.00,  1.08,  0.06,  0.88,   | 64 - 72
     & -0.13,  1.11,  0.28,  1.45,  1.35,  1.80,  1.01,  1.13,  0.90,   | 73 - 81
     &  1.95,  0.71, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0,  0.09,   | 82 - 90
     & -99.0, -0.47 /                                                   | 91 - 92
*
*
* Solar abundances ref: Grevesse Noels Sauval 1996 ASP Conf 99, 117
*     data sunabund /
*    & 12.00, 10.99,  1.16,  1.15,  2.60,  8.55,  7.97,  8.87,  4.56,
*    &  8.08,  6.33,  7.58,  6.47,  7.55,  5.45,  7.33,  5.50,  6.52,
*    &  5.12,  6.36,  3.17,  5.02,  4.00,  5.67,  5.39,  7.50,  4.92,
*    &  6.25,  4.21,  4.60,  2.88,  3.41,  2.37,  3.38,  2.63,  3.23,
*    &  2.60,  2.97,  2.24,  2.60,  1.42,  1.92, -9.00,  1.84,  1.12,
*    &  1.69,  0.94,  1.77,  1.66,  2.00,  1.00,  2.24,  1.51,  2.23,
*    &  1.13,  2.13,  1.17,  1.58,  0.71,  1.50, -9.00,  1.01,  0.51,
*    &  1.12, -0.10,  1.14,  0.26,  0.93,  0.00,  1.08,  0.76,  0.88,
*    & -0.13,  1.11,  0.28,  1.45,  1.35,  1.80,  1.01,  1.17,  0.90,
*    &  1.95,  0.71, -9.00, -9.00, -9.00, -9.00, -9.00, -9.00,  0.09,
*    & -9.00, -0.47 /
*
      data eqw/0.0/, eqwerr/1.0/
      data cion/'I ','II'/
*
      equivalence (lelec,el)
      do i=1,2
        do j=1,92
          nu(j,i)=0
        enddo
      enddo
      natom=0.0
*
      print *, 'VALD ''show line'' (long format) file name?'
      read(*,'(a)') filename
      open(10,file=filename,form='formatted',status='old')
      print *, 'output metallic line file name?'
      read(*,'(a)') filename
      open(11,file=filename,form='formatted',status='unknown')
*     print *,'Metallicity [Fe/H]?'
*     read(*,*) feh
      nutotal=0
      do i=1,500000
        read(10,'(a)',end=99) string
        read(string,1000) blip
cccc 1000   format(a1,a2,i2,2x,f10.4,1x,f7.3,1x,f8.4,16x,f5.1,22x,f6.3)
cc change of format 13/07-2007  BPz
 1000   format(a1,a2,i2,2x,f12.4,1x,f7.3,1x,f8.4,16x,f5.1,22x,f6.3)
***************************************************************************************************
*                                                            Lande factors     Damping parameters
*Elm Ion  WL(A)     log(gf) Exc. lo   J lo Exc. up   J up  lower upper  mean   Rad.   Stark  Waals
*'O  1', 6363.7760,-10.303,  0.0200,  1.0,  1.9670,  2.0,99.000,99.000,99.000,-2.170, 0.000, 0.000,
*         '  2p4 3P    2p4 1D    NBS        2   2   2   2   2   2   2   2   2'
***************************************************************************************************
        if (blip.ne.'#') then
          if(blip.eq.'''') then
             
cc             read(string,1000) blip,el,ion,w,gflog,chil,jupper,gamrad
* change to free format, as wavelength may have different length
* BPz 02/04-2008
            read(string,1000) blip,el,ion
            read(string(8:len(string)),*) w,gflog,chil,jlow,chiu,
     &                  jupper,landelow,landemean,landeup,gamrad
*****            print1000, blip,el,ion,w,gflog,chil,jupper,gamrad
* try free format here also, as there are or 10 leading blanks depending on extract "all" or "stellar"
* BPz  09/12-2009
cc            read(10,'(10x,a30)') levels
            read(10,'(a)') string
            do ii=1,11
              if (string(ii:ii).eq.'''') then
                string2=string(ii+1:100)
                goto 999
              endif
            enddo
            stop 'blip not found in string !!!!'
 999        continue
            read(string2,'(a30)') levels
ccc            print*,levels
* Keep only neutral or singly ionized lines with CHIlow<15 eV:
            if(ion.le.2.and.chil.lt.15.) then 
              if(gamrad.gt.3.0) then
                gamrad=10.0**gamrad
              else
                gamrad=0.0
              endif
* every second line contains the energy level info for the first one
              read(levels,'(a10)') tenchar
              write(desl,'(1x,a10)') tenchar
              read(levels,'(10x,a10)') tenchar
              write(desu,'(1x,a10)') tenchar
              call define_transition(desl,desu,lower,upper,llower,
     &                               lupper)
* Convert element designation to upper case if necessary:
              call chcase(lelec(1),1)
              call chcase(lelec(2),0)
              write(species,'(a2,x,a2)') el,cion(ion)
              call chcase(lelec(2),1)
* Find element number:
              do j=1,92
                if(el.eq.lele(j)) iel=j
              enddo
* Apply FDAMP, for references See BDP (A&A 275,101) or notes below
              if(ion.eq.1) then
                if(iel.eq.11) then
                  fdamp=2.0
                else if(iel.eq.14) then
                  fdamp=1.3
                else if(iel.eq.20) then
                  fdamp=1.8
                else if(iel.eq.26) then
                  fdamp=1.4
                else
                  fdamp=2.5
                endif
              else if (ion.eq.2) then
                if(iel.eq.20) then
                  fdamp=1.4
* from fit of H&K in the HM model to the fts intensity spectrum
                else if(iel.eq.38) then
                  fdamp=1.8
* from fit of Sr II 4077.724 in the HM model to the fts intensity spectrum
                else if(iel.eq.56) then
                  fdamp=3.0
                else
                  fdamp=2.5
                endif
              endif
* OK, write this line to the relevant array in the speqw format:
              nu(iel,ion)=nu(iel,ion)+1
              n=nu(iel,ion)
              if(n.gt.Nlmax-1) then
                print*, 'More than ',Nlmax-1,' lines of a species'
                stop
              endif
              write(array(n,iel,ion),1020) w,chil,gflog,fdamp,
     &              (2.*jupper+1.),gamrad,lower,upper,eqw,eqwerr,
     &               species,levels
 1020         format(f10.3,x,f6.3,x,f7.3,x,f5.2,x,f6.1,x,1p,e9.2,0p,
     &               x,'''',a1,'''',x,'''',a1,'''',x,f5.1,x,f6.1,x,
     &                '''',a5,x,a30,'''')
            endif
          endif
        endif
      enddo
      stop 'More than 500000 lines in the file!'
   99 continue
      close(10)
* write line data file header lines:
*       write(11,3010)
*3010   format('*NMY=4 IDAMP=2')
*     write(11,4000) feh, sunabund(1), sunabund(2), sunabund(6)+feh,
*    &            sunabund(7)+feh, sunabund(8)+feh,sunabund(10)+feh,
*    &           sunabund(11)+feh,sunabund(12)+feh,sunabund(13)+feh,
*    &           sunabund(14)+feh,sunabund(16)+feh,sunabund(19)+feh,
*    &           sunabund(20)+feh,sunabund(24)+feh,sunabund(26)+feh,
*    &           sunabund(28)+feh
*4000  format('* H     He    C     N     O     Ne    Na    Mg '/,
*    & '* Al    Si    S     K     Ca    Cr    Fe    Ni * [M/H]=',
*    & f5.2,/,8f6.2,/,8f6.2)
*     write(11,2000)
*2000 format('IP=0 EPS=0.0010',/,'XITE5=0.00',/,'XLBOFF=0.0000')
      do iel=1,92
        do ion=1,2
          if(nu(iel,ion).gt.0) natom=natom+1
        enddo
      enddo
      print 3000, natom
 3000 format('number of metal line species=',i3)
      do iel=2,92
        do ion=1,2
          if(nu(iel,ion).gt.0) then
            nutotal=nutotal+nu(iel,ion)
            el=lele(iel)
            call chcase(lelec(2),0)
            write(11,4050) iel,ion,nu(iel,ion),el,cion(ion)
 4050       format('''',i2,'.000',14x,'''',i5,i10,/,'''',a2,x,a2,'''')
            do n=1,nu(iel,ion)
              write(11,'(a)') array(n,iel,ion)
            enddo
          endif
        enddo
      enddo
      close(11)
      print *, nutotal,' lines written to file ',filename
      stop 'normal end'
      end
c
c
c
      subroutine chcase(c,idir)
* changes upper to lower case if idir=0
*      or lower to upper case if idir<>0
* for one character c
      integer idir
      character*1 c
*
      if(idir.eq.0) then
        if(c.eq.'A') then
              c='a'
        else if(c.eq.'B') then
                   c='b'
        else if(c.eq.'C') then
                   c='c'
        else if(c.eq.'D') then
                   c='d'
        else if(c.eq.'E') then
                   c='e'
        else if(c.eq.'F') then
                   c='f'
        else if(c.eq.'G') then
                   c='g'
        else if(c.eq.'H') then
                   c='h'
        else if(c.eq.'I') then
                   c='i'
        else if(c.eq.'J') then
                   c='j'
        else if(c.eq.'K') then
                   c='k'
        else if(c.eq.'L') then
                   c='l'
        else if(c.eq.'M') then
                   c='m'
        else if(c.eq.'N') then
                   c='n'
        else if(c.eq.'O') then
                   c='o'
        else if(c.eq.'P') then
                   c='p'
        else if(c.eq.'Q') then
                   c='q'
        else if(c.eq.'R') then
                   c='r'
        else if(c.eq.'S') then
                   c='s'
        else if(c.eq.'T') then
                   c='t'
        else if(c.eq.'U') then
                   c='u'
        else if(c.eq.'V') then
                   c='v'
        else if(c.eq.'X') then
                   c='x'
        else if(c.eq.'Y') then
                   c='y'
        else if(c.eq.'Z') then
                   c='z'
        endif
      else
        if(c.eq.'a') then
              c='A'
        else if(c.eq.'b') then
                   c='B'
        else if(c.eq.'c') then
                   c='C'
        else if(c.eq.'d') then
                   c='D'
        else if(c.eq.'e') then
                   c='E'
        else if(c.eq.'f') then
                   c='F'
        else if(c.eq.'g') then
                   c='G'
        else if(c.eq.'h') then
                   c='H'
        else if(c.eq.'i') then
                   c='I'
        else if(c.eq.'j') then
                   c='J'
        else if(c.eq.'k') then
                   c='K'
        else if(c.eq.'l') then
                   c='L'
        else if(c.eq.'m') then
                   c='M'
        else if(c.eq.'n') then
                   c='N'
        else if(c.eq.'o') then
                   c='O'
        else if(c.eq.'p') then
                   c='P'
        else if(c.eq.'q') then
                   c='Q'
        else if(c.eq.'r') then
                   c='R'
        else if(c.eq.'s') then
                   c='S'
        else if(c.eq.'t') then
                   c='T'
        else if(c.eq.'u') then
                   c='U'
        else if(c.eq.'v') then
                   c='V'
        else if(c.eq.'x') then
                   c='X'
        else if(c.eq.'y') then
                   c='Y'
        else if(c.eq.'z') then
                   c='Z'
        endif
      endif
      return
      end
c
c
c
c      program test
c      character*11 desl,desu,summary
c      character*1 lower,upper
c      integer llower,lupper
c      common /dummary/ summary
c      open(27,file='gf22.dat',status='old')
c      n=0
c   13 continue
c        if(n.gt.1) then
c          do i=1,n
c            read(27,*,end=17)
c            nread=nread+1
c          enddo
c        else if(n.lt.0) then
c          stop
c        endif
c        read(27,1000,end=17) code,desl,desu
c        nread=nread+1
c 1000   format(51x,f9.2,a11,a11)
c        call define_transition(desl,desu,lower,upper,llower,lupper)
c      print *,
c     &code,desl,' ',desu,' ',summary,' =>"',lower,'-',upper,'"',nread
c*       read(*,*) n
c        n=0
c      goto 13
c   17 continue
c      end
c
c
c
      subroutine define_transition(desl,desu,lower,upper,llower,lupper)
* try to tell what kind of transition this is
* desl, desu are 11 character level designations e.g. " (1G)sp u3F"
* lower, upper designate "l" quantum design. eg s p d ..., "X" = unknown
* llower, lupper designate "l" quantum numb. eg 0 1 2 ...
* GUESS (the lowest possible) if only one level is undefined
*
      implicit none
      integer llorbit(2),luorbit(2),llower,lupper
      character*11 desl,desu
      character*1 lower,upper
*     character*11 summary
*     common /dummary/ summary
*
      call identify_level(desl,llorbit)
      call identify_level(desu,luorbit)
*     write(summary,'(4i2)') llorbit,luorbit
* define transition type if possibly possible
      llower=llorbit(1)
      lupper=luorbit(1)
* these are the defaults, the rest are for exceptions
      if(llorbit(1).ge.0 .and. luorbit(1).ge.0 .and.
     &   (iabs(llorbit(1)-luorbit(1)).ne.1)) then
* if Delta(l) .ne. +-1 then try with the second alternatives
        if(llorbit(2).lt.0) then
* which of the 2 luorbits goes best with llorbit(1)
          if(luorbit(2).ge.0) then
            if((luorbit(1).eq.llorbit(1)) .or.
     &         (iabs(llorbit(1)-luorbit(2)).eq.1)) then
              lupper=luorbit(2)
            endif
          endif
        else
          if(luorbit(2).lt.0) then
* which of the 2 llorbits goes best with luorbit(1)
            if((llorbit(1).eq.luorbit(1)) .or.
     &         (iabs(luorbit(1)-llorbit(2)).eq.1)) then
              llower=llorbit(2)
            endif
          else
* llorbit 1 and 2 and luorbit 1 and 2 are all >=0
* find out which of the 3 remaining combinations make Delta(l)=+-1
* the ordering of these 3 alternatives is arbitrary (I've found no rule)
            if(iabs(llorbit(1)-luorbit(2)).eq.1) then
              lupper=luorbit(2)
            else
              if(iabs(llorbit(2)-luorbit(1)).eq.1) then
                llower=llorbit(2)
              else
                if(iabs(llorbit(2)-luorbit(2)).eq.1) then
* last chance: take both alternatives
                  llower=llorbit(2)
                  lupper=luorbit(2)
                endif
              endif
            endif
          endif
        endif
      endif
* GUESS (the lowest possible l) if one level is undefined:
      if(llower.lt.0 .and. lupper.ge.0) then
        if(lupper.eq.0) then
          llower=1
        else if(lupper.eq.1) then
          llower=0
        else if(lupper.eq.2) then
          llower=1
        else if(lupper.eq.3) then
          llower=2
        else if(lupper.eq.4) then
          llower=3
        else if(lupper.eq.5) then
          llower=4
        else if(lupper.eq.6) then
          llower=5
        else if(lupper.eq.7) then
          llower=6
        endif
      endif
      if(llower.ge.0 .and. lupper.lt.0) then
        if(llower.eq.0) then
          lupper=1
        else if(llower.eq.1) then
          lupper=0
        else if(llower.eq.2) then
          lupper=1
        else if(llower.eq.3) then
          lupper=2
        else if(llower.eq.4) then
          lupper=3
        else if(llower.eq.5) then
          lupper=4
        else if(llower.eq.6) then
          lupper=5
        else if(llower.eq.7) then
          lupper=6
        endif
      endif
* convert to letter labels
      if(llower.lt.0) then
        lower='X'
      else if(llower.eq.0) then
        lower='s'
      else if(llower.eq.1) then
        lower='p'
      else if(llower.eq.2) then
        lower='d'
      else if(llower.eq.3) then
        lower='f'
      else if(llower.eq.4) then
        lower='g'
      else if(llower.eq.5) then
        lower='h'
      else if(llower.eq.6) then
        lower='i'
      else if(llower.eq.7) then
        lower='k'
      else
        print *,'why is llower=',llower,' ???'
        stop 'subr. define_transition'
      endif
      if(lupper.lt.0) then
        upper='X'
      else if(lupper.eq.0) then
        upper='s'
      else if(lupper.eq.1) then
        upper='p'
      else if(lupper.eq.2) then
        upper='d'
      else if(lupper.eq.3) then
        upper='f'
      else if(lupper.eq.4) then
        upper='g'
      else if(lupper.eq.5) then
        upper='h'
      else if(lupper.eq.6) then
        upper='i'
      else if(lupper.eq.7) then
        upper='k'
      else
        print *,'why is lupper=',lupper,' ???'
        stop 'subr. define_transition'
      endif
      return
      end

      subroutine identify_level(des,lorb)
* try to tell if this is an s, p, d, f, g, h or i orbital
* method:
* 1) separate into configuration and term (look for blank), if designation
*    does not consist of 2 or more parts separated by 1 blank: NO RESULT
* 2) throw away term, concentrate on configuration
* 3) find the last lower-case letter, that gives the orbital lorb(1)
*    if any character before is also a lower-case letter, that gives the lorb(2)
*    lorb < 0 means no result
*
      implicit none
      character*1 des(11)
      integer i,l,n,lorb(2)
*
      lorb(1)=-1
      lorb(2)=-1
      i=12
   11 continue
        i=i-1
      if(des(i).eq.' '.and.i.gt.1) goto 11
* OK: we've found the last character in the term
   12 continue
        i=i-1
      if(i.le.0) return
      if(des(i).ne.' '.and.i.gt.1) goto 12
* OK: we've found the (first) blank before the term
* now look for lower case letters:
      n=1
   13 continue
        i=i-1
        if(i.le.0) return
        if(des(i).eq.'s') then
          l=0
        else if(des(i).eq.'p') then
          l=1
        else if(des(i).eq.'d') then
          l=2
        else if(des(i).eq.'f') then
          l=3
        else if(des(i).eq.'g') then
          l=4
        else if(des(i).eq.'h') then
          l=5
        else if(des(i).eq.'i') then
          l=6
        else if(des(i).eq.'k') then
          l=7
        else
          l=-1
          if(i.gt.1) goto 13
        endif
        lorb(n)=l
      if(n.eq.2) then
        return
      else
        n=n+1
        if(i.gt.1) goto 13
      endif
      return
* lorb(1) and lorb(2) contain the two last l quantum numbers
      end
*
