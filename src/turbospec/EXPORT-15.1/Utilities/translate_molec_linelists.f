*To: plez@Ferrum.fysik.lu.se (Bertrand Plez)
*Subject: ny translatelist
*Date: Fri, 17 Jul 1998 16:32:48 +0200
*From: Nils Ryde <ryde@astro.uu.se>
* 
* Corrected error for C2. BPz 22/09-2003
* Added identification of lines.  BPz 24/03-2006
****   STILL TO BE DONE for some species 
*
* this version with selection of strongest lines only. BPz 21/10-2010
* this version with new VALD TiO, and replacement of laboratory wavenumber,
* and corrected delta wavenumbers. Bpz 9/12-2011
*

      program translatelinelists
*
* this program translates linelists from the TiO/linelist format 
* to the new bsyn ascii format including Fdamp.
*   BPz 10/07-1997
*
      implicit none
      character filein*100,fileout*100,ethop*22,isotopes*10,
     &          elname*7,transition*6,comment*78,formatin*1,
     &          elstring*21,system*1,whichvo*1,q1*9,bran*4,
     &          whichcah*1,whichoh*1,sys*3,molname*3,blabla*85,
     &          deltaN*1,deltaJ*1,e*1,mol*3,trans*2,branch2*7,
     &          identification*20,lab*3,levellow*2,levelup*2,
     &          translow*1,transup*1,sym*1,br*2,idwave*3,red*1,
     &          idvald*14,levup*1,levlow*1
      doubleprecision lambda,x,x2,refrind,wave,isodw,edmf,A,
     &                S,xx,isofact(3),sysfact(3),wavekotlar,
     &                massfact(3),dw(30),dwmax,dwmin,waveram,
     &                nwaveplez(4),nwavekot(4),nwaveram(4)
      real            chie,fdamp,gf,gflog,gup,raddamp,el,
     &                lambdamin,lambdamax,limite,chiup,
     &                planck,lumiere,electron,kayser2eV,
     &                strength,agam,sgam,n,d,xjup,xjlow,
     &                strongest,tempselect,weakratio,ngf(4),
     &                nexc(4),Eup,nup,nlow,eqw,eqwerr
      integer i,nline,iel,ion,id(6),kkk,jup,imol,ivlow,ivup,
     &        jlow,ibranch,isot,outfil,itrans,isowant,isys,
     &        mo,iso,ierf,iers,ierh,ireff,irefs,irefh,
     &        v1,v2,lev1,lev2,icode,niso,symup,symlow

      logical notend,readVO,readH2O,readCNugj,readCaH,readTiO,
     &        readCH12,readCH13,increasewave,readC2,readSiO,
     &        readFeH, readCNBPz, readCNBpz_old, readCO, readOH,
     &        readCH4, readHCN, readC2H2,readZrO,readNH,readSiH,
     &        readMgH,readAlH,readC2querci,readC2Kurucz,first,
     &        readSiS,readLaO,readMgHKurucz,cnred,readTiOVALD,
     &        readCNsinglefile
      real tioisotopf(5),fact
      data tioisotopf/0.08, 0.073, 0.738, 0.055, 0.054/

      isotopes='          '
      transition='      '
      iel=0
      nline=0
      strongest=1.e-30
      identification=' '
      planck=6.62620e-34
      lumiere=2.997925e10
      electron=1.602192e-19
      kayser2eV = planck * lumiere / electron
      xjup=0.
      outfil=20
      elname='      '
      first=.true.
      levup='X'
      levlow='X'
      eqw=0.0
      eqwerr=1.0
*
      print*,'input linelist to translate?'
*
cc      print*,'For H2O it is : /ste2/uffegj/scanh2o_vo7.dat'
      read(5,10) filein
      open(10,file=filein,status='old',pad='yes')
      print*,'Species? (e.g. TiO, or VO)'
      read(5,10) elname
*
      print*,'output linelist ?'
*
      read(5,10) fileout
      open(30,status='scratch')
*
      print*,'lower and upper limit for lambda to select'
*
      read(5,*) lambdamin,lambdamax
      if (lambdamin.gt.lambdamax) then
        chie=lambdamax
        lambdamax=lambdamin
        lambdamin=chie
      endif
*
      print*,'temperature for selection of lines '
      read(5,*) tempselect
      print*,'strength of weakest line to be selected',
     &       'relative to strongest'
      read(5,*) weakratio
      weakratio=log10(weakratio)
*
      readFeH=.false.
      readC2=.false.
      readC2Kurucz=.false.
      readC2querci=.false.
      readVO=.false.
      readTiO=.false.
      readTiOVALD=.false.
      readLaO=.false.
      readSiO=.false.
      readSiS=.false.
      readH2O=.false.
      readCNugj=.false.
      readCNBPz=.false.
      readCNBpz_old=.false.
      readCaH=.false.
      readCH12=.false.
      readCH13=.false.
      readCO=.false.  
      readOH=.false.
      readNH=.false.
      readMgH=.false.
      readMgHKurucz=.false.
      readAlH=.false.
      readSiH=.false.
      readCH4=.false.
      readHCN=.false.
      readC2H2=.false.
      readZrO=.false.
      readCNsinglefile=.false.!!
       
      increasewave=.false.
      limite=-30.
      raddamp=0.0
      ion=1
      fdamp=0.
      fact=1.0
      if (elname.eq.'TiO   ') then
        elstring='0822.000000          '
        readTiO=.true.
        fact=1.0
      else if (elname(3:6).eq.'TiOV') then
        readTiOVALD=.true.
        if (elname(1:2).eq.'46') then
          elstring='0822.000046          '
        else if (elname(1:2).eq.'47') then
          elstring='0822.000047          '
        else if (elname(1:2).eq.'48') then
          elstring='0822.000048          '
        else if (elname(1:2).eq.'49') then
          elstring='0822.000049          '
        else if (elname(1:2).eq.'50') then
          elstring='0822.000050          '
        endif
      else if (elname.eq.'46TiO ') then
        elstring='0822.000046          '
        readTiO=.true.
        fact=1.0/tioisotopf(1)
      else if (elname.eq.'47TiO ') then
        elstring='0822.000047          '
        readTiO=.true.
        fact=1.0/tioisotopf(2)
      else if (elname.eq.'48TiO ') then
        elstring='0822.000048          '
        readTiO=.true.
        fact=1.0/tioisotopf(3)
      else if (elname.eq.'49TiO ') then
        elstring='0822.000049          '
        readTiO=.true.
        fact=1.0/tioisotopf(4)
      else if (elname.eq.'50TiO ') then
        elstring='0822.000050          '
        readTiO=.true.
        fact=1.0/tioisotopf(5)
      else if (elname.eq.'VO    ') then
cc        print*,'VO from Brett (B) or VO from Plez (P)?'
cc        read(5,'(A)') whichvo
        whichvo='P'
        elstring='0823.000000          '
        readVO=.true.
        limite=-99.
      else if (elname.eq.'LaO   ') then
        elstring='0857.000000          '
        readLaO=.true.
        limite=-99.
      else if (elname.eq.'FeH   ') then
        elstring='0126.000000          '
        readFeH=.true.
        limite=-99.
      else if (elname.eq.'CaH   ') then
cc        print*,'CaH from Brett (B) or CaH from Plez (P)?'
cc        read(5,'(A)') whichcah
        whichcah='P'
        elstring='0120.000000          '
        readCaH=.true.
        limite=-99.
      else if (elname.eq.'H2O   ') then
        elstring='010108.000000000       '
        readH2O=.true.
        increasewave=.true.
        jup=1
        xjup=float(jup)
* we set jup=1 because we don't know what it is.
      else if (elname.eq.'28SiO ') then
        elstring='0814.000028          '
        readSiO=.true.
        isowant=28
        increasewave=.true.
      else if (elname.eq.'29SiO ') then
        elstring='0814.000029          '
        readSiO=.true.
        isowant=29
        increasewave=.true.
      else if (elname.eq.'30SiO ') then
        elstring='0814.000030          '
        readSiO=.true.
        isowant=30
        increasewave=.true.
      else if (elname.eq.'C1212K') then
* reads Kurucz linelist.
        elstring='0606.012012          '
        isowant=1212
        readC2Kurucz=.true.
      else if (elname.eq.'C1213K') then
* reads Kurucz linelist.
        isowant=1213
        elstring='0606.012013          '
        readC2Kurucz=.true.
      else if (elname.eq.'C1313K') then
* reads Kurucz linelist.
        isowant=1313
        elstring='0606.013013          '
        readC2Kurucz=.true.
      else if (elname.eq.'C1212Q') then
* reads Querci&Plez new linelist.
        elstring='0606.012012          '
        isowant=1
        readC2querci=.true.
      else if (elname.eq.'C1213Q') then
* reads Querci&Plez new linelist.
        isowant=2
        elstring='0606.012013          '
        readC2querci=.true.
      else if (elname.eq.'C1313Q') then
* reads Querci&Plez new linelist.
        isowant=3
        elstring='0606.013013          '
        readC2querci=.true.
      else if (elname.eq.'C12C12') then
* reads Rurik&Plez new linelist.
        elstring='0606.012012          '
        isowant=1
        readC2=.true.
      else if (elname.eq.'C12C13') then
* reads Rurik&Plez new linelist.
        isowant=2
        elstring='0606.012013          '
        readC2=.true.
      else if (elname.eq.'C13C13') then
* reads Rurik&Plez new linelist.
        isowant=3
        elstring='0606.013013          '
        readC2=.true.
      else if (elname.eq.'C12N14') then
* these first 4 are BPz violet CN or red (nils)
        elstring='0607.012014          '
        readCNBPz=.true.
      else if (elname.eq.'C12N15') then
        elstring='0607.012015          '
        readCNBPz=.true.
      else if (elname.eq.'C13N14') then
        elstring='0607.013014          '
        readCNBPz=.true.
      else if (elname.eq.'C13N15') then
        elstring='0607.013015          '
        readCNBPz=.true.
      else if (elname.eq.'C2N4sg') then
* The next four cases are for the single large CN file with 4 isotopes !!
* created by BPz and Ruben Pedro Hedrosa, on 01/06-2011 !!
        elstring='0607.012014          '
        isowant=1
        readCNsinglefile=.true.
      else if (elname.eq.'C2N5sg') then
        elstring='0607.012015          '
        isowant=2
        readCNsinglefile=.true.
      else if (elname.eq.'C3N4sg') then
        elstring='0607.013014          '
        isowant=3
        readCNsinglefile=.true.
      else if (elname.eq.'C3N5sg') then
        elstring='0607.013015          '
        isowant=4
        readCNsinglefile=.true.
      else if (elname.eq.'C12N  ') then
* these are for joergensen and Larsson list red CN.
        elstring='0607.012000          '
        readCNugj=.true.
* increasewave tells us that the line list is in increasing wavenumber order.
        increasewave=.true.
      else if (elname.eq.'C13N  ') then
        elstring='0607.013000          '
        readCNugj=.true.
        increasewave=.true.
      else if (elname.eq.'C12H  ') then
        elstring='0106.000012          '
        readCH12=.true.
      else if (elname.eq.'C13H  ') then
        elstring='0106.000013          '
        readCH13=.true.
        dwmin=1.d13
        dwmax=0.d0
        do i=1,30
          read(10,*) itrans,ivup,ivlow,jlow,ibranch,gf,wave,chie,dw(i)
          dwmax=max(dwmax,dw(i))
          dwmin=min(dwmin,dw(i))
        enddo
        if ((int(dwmax+0.1d0).eq.3).and.(int(dwmin+0.1).eq.3)) then
          print*,'WARNING! reading 13CH from wrong file!!'
          print*,'WARNING! Should be scan_CH.dat with mixed isotopes!!'
          stop
        else
          rewind(10)
        endif
      else if (elname.eq.'26CO  ') then
        elstring='0608.012016          '
        readCO=.true.
        isowant=26
      else if (elname.eq.'27CO  ') then
        elstring='0608.012017          '
        readCO=.true.
        isowant=27
      else if (elname.eq.'28CO  ') then
        elstring='0608.012018          '
        readCO=.true.
        isowant=28
      else if (elname.eq.'36CO  ') then
        elstring='0608.013016          '
        readCO=.true.
        isowant=36
      else if (elname.eq.'37CO  ') then
        elstring='0608.013017          '     
        readCO=.true.
        isowant=37
      else if (elname.eq.'38CO  ') then
        elstring='0608.013018          '
        readCO=.true.
        isowant=38
      else if (elname.eq.'46CO  ') then
        elstring='0608.014016          '
        readCO=.true.
        isowant=46
      else if (elname.eq.'14NH  ') then
        elstring='0107.000014         '
        readNH=.true.
        isowant=14
      else if (elname.eq.'15NH  ') then
        elstring='0107.000015         '
        readNH=.true.
        isowant=15
      else if (elname.eq.'OH-IR ') then
c IR OH from HItran.
        elstring='0108.000016         '
        readOH=.true.
        whichoh='G'
        print*,'Hitran vib-rot OH 16OH'
      else if (elname.eq.'OH-AX ') then
        elstring='0108.000000         '
        readOH=.true.
        whichoh='A'
        print*,'GOLDMAN 2001 OH A-X  or Hitran IR'
      else if (elname.eq.'16OH  ') then
        elstring='0108.000016         '
        readOH=.true.
        isowant=116
        whichoh='K'
        print*,'Kurucz A-X OH 16OH'
      else if (elname.eq.'18OH  ') then
        elstring='0108.000018         '
        readOH=.true.
        isowant=118
        whichoh='K'
        print*,'Kurucz A-X OH 18OH'
      else if (elname.eq.'AlH     ') then
        elstring='0113.000000          '
        readAlH=.true.
      else if (elname.eq.'MgH     ') then
        elstring='0112.000000          '
        readMgH=.true.
      else if (elname.eq.'24MgH   ') then
        elstring='0112.000024          '
        readMgHKurucz=.true.
        isowant=24
      else if (elname.eq.'25MgH   ') then
        elstring='0112.000025          '
        readMgHKurucz=.true.
        isowant=25
      else if (elname.eq.'26MgH   ') then
        elstring='0112.000026          '
        readMgHKurucz=.true.
        isowant=26
      else if (elname.eq.'28SiH   ') then
        elstring='0114.000028          '
        readSiH=.true.
        isowant=28
      else if (elname.eq.'29SiH   ') then
        elstring='0114.000029          '
        readSiH=.true.
        isowant=29
      else if (elname.eq.'30SiH   ') then
        elstring='0114.000030          '
        readSiH=.true.
        isowant=30
      else if (elname.eq.'90ZrO   ') then
        elstring='0840.000090          '
        readZrO=.true.
      else if (elname.eq.'91ZrO   ') then
        elstring='0840.000091          '
        readZrO=.true.
      else if (elname.eq.'92ZrO   ') then
        elstring='0840.000092          '
        readZrO=.true.
      else if (elname.eq.'94ZrO   ') then
        elstring='0840.000094          '
        readZrO=.true.
      else if (elname.eq.'96ZrO   ') then
        elstring='0840.000096          '
        readZrO=.true.
      else if (elname.eq.'HF      ') then
        elstring='0109.000000          '
      else if (elname.eq.'SiS     ') then
        elstring='1416.000000          '
        readSiS=.true.
        increasewave=.true.
      else
* put if (...eq.other molecules...) etc here
        do 1966 kkk=1,6
          id(kkk)=ichar(elname(kkk:kkk))-48
1966    continue
        iel=id(6)+id(5)*10+id(4)*100+id(3)*1000+id(2)*10000
        ion=iel-int(iel/100.)*100
        iel=(iel-ion)/100
        ion=ion+1
        fdamp =2.0
        el=iel
        write(elstring,1000) iel
1000    format(i2,'.000               ')
      endif
**
* check if CNred, with extra wavenumbers
      if (readCNBPz) then
        read(10,1081) mol,ivup,ivlow,jlow,branch2,gf,
     &       wave,chie,trans
        if (trans.eq.'R') then
          cnred=.true.
        else
          cnred=.false.
        endif
        rewind(10)
      endif
**
      print 1234,elname,isotopes,transition,iel,ion,nline
      print*,elstring
      write(comment,10) filein(1:78)
cc      write(30,19) elstring,ion,nline
cc      write(30,18) comment

      if (readC2querci) then
*
* isofact contains a correction factor to bring all isotopes up to the
* original gf value. Querci have included factors = 1, 1/2, and 1/16 
* into 12C2, 12C13C, and 13C2 respectively. 
* I remove the factor 2phi that QQK introduced for 12C2 and 13C2.
* This is the additional 0.5 factor in isofact.
* Following Goorvitch  1990, ApJS 74, 769, the partition function for
* 12C2 is half the full partition function computed summing all
* levels of all parities. The full partition function is
* twice that of Irwin. I checked it by summing (2J+1)exp(-E/kT)
* for the X1sigma and a3pi levels of QQK's linelist (for 12C13C where
* all levels appear), and I get about twice Irwin's partition function. 
* therefore, for 12C2, the partition function to be used is Irwin's.
* For 12C13C the partition function should be 4 times larger, and
* the statistical weights of the levels be 2 or 2, depending on J and
* their parity. As Lambda doubling is not resolved, we can assume
* that all levels have a weight of 4 in the Swan bands (we should see
* an intensity alternation in the other bands?). The factor 4
* cancels out in the population computation.
* For 13C2 the situation is the same, with weights of 3 and 1.
*
* See also below (computation of gf), for an extra factor 2 due to 
* lambda-doubling in Swan 12C13C transitions
*
*  BPz 02/10-2003

        isofact(1)=1. *0.5
        isofact(2)=2. *1.0 
        isofact(3)=16.*0.5 
*****
*****  ADDITIONAL FACTOR 2 IN ISOFACT(2) AND FACTOR 4 IN ISOFACT(3)
*****  ADDED ON MAY 16 2005, AS I CANNOT FIT 13C2 BANDS IN SPECTRA
*****       REMOVED AGAIN RIGHT AWAY. LOOKS AS WRONG IN SPECTRA!
*****
********        isofact(2)=4.      ! typo error! changed 22/09-2003
* sysfact contains the numerical factor to scale the Querci gf values
* to modern lifetime measurements. The gf values in the line list are from
* Querci, Querci and Kunde, 1971, AA 15, 256.
*
        sysfact(1)=1.68
        sysfact(2)=0.243
        sysfact(3)=1.27
*
* massfact contains the factors M(XY)/M(C2_solar_isotopic_mix) that
* allows to put all the gf values on the same scale. In this way,
* the opacity will be per g of C2 (solar_iso_mix) for all isotopes.
* this is needed in OS files, as we later multiply by the mass of C12 (solar)
* in MARCS.
*
        massfact(1)=24.*1.66056d-24
        massfact(2)=25.*1.66056d-24
        massfact(3)=26.*1.66056d-24
      endif

      limite=10.**limite
      notend=.true.
      do while (notend)
          if (readVO) then
            if (whichvo.eq.'B') then
              read(10,1002,end=99) jup,gf,wave,chie
              xjup=float(jup)
* we shift the zero of energies from Brett definition (bottom of 
* potential) to the definition of Sauval and Tatum for partition 
* functions (lowest energy level).
* BEWARE!!! VO linelists contain WAVELENGTH and not wavenumber!!!!
              chie=chie-0.063
            else if (whichvo.eq.'P') then
              read(10,1005,end=99) jup,gf,wave,chie,system
              xjup=float(jup)
            else
              stop 'unknown VO format'
            endif
          else if (readSiO) then
123         read(10,*,end=99) wave,chie,gf,ivup,ivlow,jup,jlow,isot
            xjup=float(jup)
            if (isot.ne.isowant) goto 123
            chie=chie*kayser2eV
          else if (readCaH) then
            if (whichcah.eq.'B') then
              read(10,1002,end=99) jup,gf,wave,chie
              xjup=float(jup)
* Same comment as for VO.
              chie=chie-0.080
            else if (whichcah.eq.'P') then
              read(10,1005,end=99) jup,gf,wave,chie,system
              xjup=float(jup)
            else
              stop 'unknown CaH format'
            endif
          else if (readH2O) then
234         read(10,*,end=99) wave,gf,chie
* BPz 14/08-1997: inclusion d'un cut a 10-4 de la raie la plus forte 
*                 a 4000K. Le facteur numerique est pour cm-1 -> eV.
*                 la raie la plus forte a 6.096, la plus faible -1.15
            strength=log10(gf)-chie*5040./4000.*1.24e-4
            if (strength.lt.(6.096-4.)) goto 234
            chie=chie*kayser2eV
* ?????????????????????????????????????????????????????????
* Transformation S (km/mol) --> gf / ( 1.02704/(B^2*A)^0.5 )
*  car Q = Qvib * Qrot = Qvib * T^(3/2)*1.02704/(B^2*A)^0.5
            gf=1.87572e-12*gf*1.e5
* ?????????????????????????????????????????????????????????
          else if (readCNugj) then
* CN Larsson + Joergensen
11          read(10,1003,end=99) imol,isot,ivup,ivlow,jlow,ibranch,gf,
     &               wave,chie
1003        format(i3,i2,i3,i3,i5,i3,e13.0,f11.0,f10.0)
            chie = chie*kayser2eV
            gf=gf*0.734
            if (jlow.lt.0) goto 11
            xjup=float(jlow)-0.5
* should be jlow-0.5 + f(pqr), i.e. should depend of the branch. Here we
* assume jup = jlow.
            
* CN Bertrand Plez           
* edited on 28/04-2011 to account for wavenumbers from Kotlar, and Ram et al., 
* now added at end of data line.
*
          else if (readCNBPz) then
12          continue
            wavekotlar=-1.0
            waveram=-1.0
            if (cnred) then
              read(10,2000,end=99) mol,ivup,ivlow,jlow,branch2,gf,
     &                         wave,chie,trans,wavekotlar,waveram
            else
              read(10,1081,end=99) mol,ivup,ivlow,jlow,branch2,gf,
     &                         wave,chie,trans
            endif
* By default, use Plez wavenumber
            idwave='   '
            if (mol(3:3).eq.'3') then
* this is for 13C14N, where we assume that Ram et al. is more accurate than Kotlar,
* because Ram proviede specific constants for 13C14N, whereas our Kotlar is a scaled
* 12C14N
              if (waveram.gt.-1.0) then
                wave=waveram
                idwave='Ram'
              else if (wavekotlar.gt.-1.0) then
                wave=wavekotlar
                idwave='Kot'
              endif
            else
* other isotopes (we have Ram data only for 12C14N)
              if (wavekotlar.gt.-1.0) then
* if Kotlar data is available, we use it, ... unless
                wave=wavekotlar
                idwave='Kot'
              endif
* ... unless Ram et al. data exists
              if (waveram.gt.-1.0) then
* But don't do it if they seem not to be reliable (12C14N)
                if (ivup.eq.8.and.jlow.gt.25) then
                else if (ivup.eq.7.and.jlow.gt.70) then
                else if (ivlow.eq.11.and.jlow.gt.30) then
                else if (ivlow.eq.12.and.jlow.gt.30) then
                else
                  wave=waveram
                  idwave='Ram'
                endif
              endif
            endif

1081        format(a3,i3,i3,i4,a7,e12.0,f11.0,f7.0,1x,a1)
2000        format(a3,i3,i3,i4,a7,e12.0,f11.0,f7.0,1x,a1,1x,
     &             f10.0,1x,f10.0)
            write(identification,1082) ivup,ivlow,branch2(2:2),
     &                                 branch2(4:7),jlow,idwave(1:3)
1082        format('''',i2,i3,1x,a1,a4,i4,1x,a2,'''')
ccc            if (jlow.lt.0) goto 12
            xjlow=float(jlow)+0.5
            trans=branch2(5:5)
            if (trans .eq. 'P') then
              xjup=xjlow-1
            else if (trans .eq. 'Q') then
              xjup=xjlow
            else if (trans .eq. 'R') then
              xjup=xjlow+1
            endif
* Added on 01/06-2011 !!
          else if (readCNsinglefile) then
            READ(10,2020,end=99) mol,ivup,ivlow,jlow,branch2,
     &                  (ngf(niso),niso=1,4),
     &                  (nwaveplez(niso),niso=1,4),
     &                  (nwavekot(niso),niso=1,4),
     &                  (nwaveram(niso),niso=1,4),
     &                  (nexc(niso),niso=1,4),trans
2020    format(a3,i3,i3,i4,a7,4(1x,e11.4),3(1x,4(1x,f10.3)),
     &       1x,4(1x,f6.3),1x,a1)
!!
            if (isowant.eq.3) then
! This is 13C14N. We have :
!  - wavenumbers computed by Plez 
!  - Kotlar's 12C14N wavenumbers+isotopicshift from Plez
!  - Ram et al's wavenumbers
! we chose 1) Ram 2) Kotlar 3) Plez
              wave=nwaveplez(isowant)
              chie=nexc(isowant)
              gf=ngf(isowant)
              idwave='   '
              if (nwaveram(isowant).gt.-1.0) then
                wave=nwaveram(isowant)
                idwave='Ram'
              else if (nwavekot(isowant).gt.-1.0) then
                wave=nwavekot(isowant)
                idwave='Kot'
              endif
            else
              wave=nwaveplez(isowant)
              chie=nexc(isowant)
              gf=ngf(isowant)
              idwave='   '
* other isotopes (we have Ram data only for 12C14N)
              if (nwavekot(isowant).gt.-1.0) then
* if Kotlar data is available, we use it, ... unless
                wave=nwavekot(isowant)
                idwave='Kot'
              endif
* ... unless Ram et al. data exists
              if (nwaveram(isowant).gt.-1.0) then
* But don't do it if they seem not to be reliable (12C14N)
                if (ivup.eq.8.and.jlow.gt.25) then
                else if (ivup.eq.7.and.jlow.gt.70) then
                else if (ivlow.eq.11.and.jlow.gt.30) then
                else if (ivlow.eq.12.and.jlow.gt.30) then
                else
                  wave=nwaveram(isowant)
                  idwave='Ram'
                endif
              endif
            endif
            if (isowant.eq.2.or.isowant.eq.4) then
! This is 12C15N or 13C15N. We have :
!  - wavenumbers computed by Plez 
!  - Kotlar's 12C14N wavenumbers+isotopicshift from Plez
! We don't have Ram's data. In order to have better line positions than 
! those calculated by Plez, we try to use Ram's data for 12C14N. For that,
! we compute the correction ram-plez for isotope 1314. We multiply this
! by the ratio of isotopic shifts 1215/1314 (from Plez). This number is 
! close to 0.74. We then add this correction to Plez wavenumber.
! This was verified to give good results on lines of 12C15N for which
! we had observed line positions from Reginald Colin. 
! we do the same for 13C15N.
!
              if (nwaveram(3).gt.-1.0) then
cc                if (ivup.EQ.2.AND.ivlow.EQ.0.and.jlow.LT.20) then
                  wave=nwaveplez(isowant)+(nwaveram(3)-nwaveplez(3))*
     &                 (nwaveplez(isowant)-nwaveplez(1))/
     &                 (nwaveplez(3)-nwaveplez(1))
                  idwave='PzR'
cc                endif
              endif
            endif
!!
            write(identification,2082) ivup,ivlow,
     &                                 branch2(4:7),jlow,idwave(1:3)
!!
2082        format('''',i2,i3,1x,a4,i4,1x,a3,'''')!!
            xjlow=float(jlow)+0.5
            trans=branch2(5:5) ! maybe better to use directly branch(5:5)
            if (trans .eq. 'P') then
              xjup=xjlow-1
            else if (trans .eq. 'Q') then
              xjup=xjlow
            else if (trans .eq. 'R') then
              xjup=xjlow+1
            endif
!!
          else if (readCH12) then
* CH Joergensen et al 1997 or LIFBASE data (12/09-2000, BPz + AJ)
            read(10,*,end=99) itrans,ivup,ivlow,jlow,ibranch,gf,wave,
     &                        chie,isodw
            chie = chie*kayser2eV
            xjup=float(jlow)-0.5
* should be jlow-0.5 + f(pqr), i.e. should depend of the branch. Here we
* assume jup = jlow.
          else if (readCH13) then
* CH Joergensen et al 1997 or LIFBASE data (12/09-2000, BPz + AJ)
            read(10,*,end=99) itrans,ivup,ivlow,jlow,ibranch,gf,wave,
     &                        chie,isodw
            wave=wave+isodw
            chie = chie*kayser2eV
            xjup=float(jlow)-0.5
* should be jlow-0.5 + f(pqr), i.e. should depend of the branch. Here we
* assume jup = jlow.
          else if (readC2querci) then
C* Old C2
C            read(10,*,end=99) isot, gf,wave,chie
C            if (isot.eq.1) then
C              gf=gf*0.9890*0.9890
C            else if (isot.eq.2) then
C              gf=gf*0.9890*0.0110
C            else if (isot.eq.3) then
C              gf=gf*0.0110*0.0110
C            else
C              print*, 'wrong isotope for C2: ',isot
C              stop
C            endif       
* New Querci C2. added 10/03-1999 BPz
* Note that in fact "jup" is K" in Querci's file
*
126         read(10,125,end=99) imol,isys,iso,ivup,ivlow,bran,jup,
     &                   chie,wave,gf
            xjup=float(jup)
 125    format(3i2,2i3,a4,i4,e15.0,f12.0,e15.0)
* Format in file: (after reformatting by reformat_C2_Querci.f and sorting
*
* C2  iso
*   sys  v' v" bran K"    exc (cm-1)    lambda (A)     gf (cgs)
*  2 2 1  0  0  R    0 0.92402356E+03  12086.3190 0.42069549E+09
*
*  sys =1: Ballik-Ramsay        iso =1 12C12C
*      =2: Phillips                 =2 12C13C
*      =3: Swan                     =3 13C13C

            if (iso.ne.isowant) goto 126
            if(isys.lt.1.or.isys.gt.3) stop 'bad system!'
            chie=(chie-900.)*kayser2eV
* the lowest levels are at 888cm-1 for 13C2, 906 for 12C13C, and 921 for 12C2.
            gf=gf*sysfact(isys)*isofact(iso)/8.8528d-13*
     &         massfact(iso)
* Now, for the Swan bands, lambda-doubling appears, but is not included
* by QQK. However, because of intensity alternation, it does not
* show up in 12C2 and 13C2. It does affect 12C13C, as all lines should be 
* present. An extra factor of 2 in the gf-values is necessary!
            if (isys.eq.3.and.iso.eq.2) then
              gf=gf*2.
            endif
* 8.8528d-13 is = pie2/mec2  ; * 0.5 is to bring the gf values in accordance
* with Irwin's partition function which is twice as small as Quercis (and
* corresponds to 12C12C, or 13C13C, whereas Quercis is for 12C13C).

          else if (readCO) then
124         read(10,*,end=99) wave,edmf,A,chie,gf,S,ivup,ivlow,
     &      transition,jlow,isot
            if (isot.ne.isowant) goto 124
            chie=chie*kayser2eV
            xjup=float(jlow)
*  we use here jlow in stead of jup, doesn't make a lot of difference 
          else if (readNH) then
* NH Kurucz  14 and 15NH
322           read(10,*,end=99) molname,iso,wave,gf,chie,sys
              if (iso.ne.isowant) goto 322
              gf=10.**(gf-0.807)
* gf correction estimated from Meyer and Roth 1991, ApJ376, L49 (2 lines).
              chie=chie*kayser2eV
          else if (readOH) then
c         voor de lijnlijsten van Jacques Sauval
c            read(10,25,end=99) wave,codeO,codeH,codekur1,codekur2,chie,
c     &      flog,factor,s_j,code,num0,molline2,branch,spec,jlow
c            chie=chie*kayser2eV
c            f=10**(flog)
c            gf=s_j*f 
c            if (branch .le. 'P') then
c              xjup=float(jlow)-1.
c              goto 500
c             else if (branch .le. 'Q') then
c              xjup=float(jlow)
c              goto 500
c             else if (branch .le. 'R') then
c              xjup=float(jlow)+1.
c              goto 500
c             endif 
c         voor de OH lijnlijsten van D. Schwenke
c            read(10,*,end=99) code,wave,gf,edmf,chie,jup,jlow,vnodup,
c     &      xup,fineup,vnodlow,xlow,finelow,rootup,rootlow,parity
c            chie=chie*kayser2eV
c            xjup=float(jup)
            if (whichoh.eq.'G'.or. whichoh.eq.'A') then
              if (whichoh.eq.'G') then
c         voor de OH lijnlijst van Goldman (Hitran IR)
                read(10,28,end=99) mo,iso,wave,S,A,agam,sgam,chie,n,d,
     &          v1,v2,q1,lev1,lev2,deltaN,deltaJ,xjlow,e,ierf,iers,ierh,
     &           ireff,irefs,irefh
              else if (whichoh.eq.'A') then
c Goldman OH A-X 2001 or Hitran IR
                read(10,29,end=99) mo,iso,wave,S,A,agam,sgam,chie,n,d,
     &           v1,v2,deltaN,deltaJ,br,xjlow,sym
29              format(i2,i1,f12.6,e10.3,e10.3,f5.4,f5.4,f10.4,f4.2,
     &                 f8.6,i3,i3,8x,a1,a1,a2,f5.1,a1)
C13131333.712295 9.478E-21 5.592E+03.0086.0000 9601.92040.660.000000  2  2        RQ21 12.5E442 1 1 1 1
              endif
              chie=chie*kayser2eV
              if (deltaJ .eq. 'P') then
                xjup=xjlow-1
                gf=1.499*(2.*xjup+1.)*A/wave/wave
              else if (deltaJ .eq. 'Q') then
                xjup=xjlow
                gf=1.499*(2.*xjup+1.)*A/wave/wave
              else if (deltaJ .eq. 'R') then
                xjup=xjlow+1
                gf=1.499*(2.*xjup+1.)*A/wave/wave
              endif
***********************
              write(identification,1085) v1,v2,deltaN,deltaJ,br,xjlow,
     &          sym
1085       format('''',i2,i2,1x,a1,a1,a2,f5.1,a1,'''')
***********************
            else if (whichoh.eq.'K') then
* OH A-X Kurucz  16 and 18OH
321           read(10,*,end=99) molname,iso,wave,gf,chie,sys
              if (iso.ne.isowant) goto 321
              gf=10.**(gf)
              chie=chie*kayser2eV
            else
              stop 'Bad choice for OH !'
            endif
          else if (readC2Kurucz) then
* C2 Kurucz  1212, 1213, 1313
329           read(10,*,end=99) molname,iso,wave,gf,chie,sys
              system=sys(1:1)
              if (iso.ne.isowant) goto 329
              gf=10.**(gf)
              chie=chie*kayser2eV
          else if (readMgH) then
* MgH from Skory, Stancil, Weck et al.
*'MgH',vup,vlow,jup,jlow,gf,wave,E(jlow,vlow,1),'AX'
*'MgH' 13  0   5   4  0.44591E-12   32039.738    112.872 'AX'
* implemented 30/09-2009 BPz
            read(10,*,end=99) molname,ivup,ivlow,jup,jlow,gf,wave,chie,
     &                        trans
            if (molname.ne.'MgH') stop 'ERROR: linelist not for MgH!'
            chie=chie*kayser2ev
            write(identification,1090) ivup,ivlow,jup,jlow,trans
1090        format('''',i2,1x,i2,1x,i3,1x,i3,1x,a2,'''')
          else if (readMgHKurucz) then
* MgH Kurucz  24, 25 and 26MgH
323         read(10,*,end=99) molname,iso,wave,gf,chie,sys
            if (iso.ne.isowant) goto 323
            gf=10.**(gf)
            chie=chie*kayser2eV
          else if (readSiH) then
* SiH Kurucz  28, 29 and 30SiH
cc OLD Kurucz line list removed on 25/08-2006 BPz
cc324           read(10,*,end=99) molname,iso,wave,gf,chie,sys
cc              if (iso.ne.isowant) goto 324
cc              gf=10.**(gf)
cc              chie=chie*kayser2eV
cc              wave=1.d8/wave
            if (first) then
              read(10,'(a)') blabla
              do while (blabla(1:1).eq.'*')
                read(10,'(a)') blabla
              enddo
              backspace(10)
              first=.false.
            endif
324         read(10,3244,end=99) wave,gf,xjlow,chie,xjup,chiup,
     &                      icode,translow,
     &                      ivlow,levellow,transup,ivup,levelup,iso
3244        format(f9.4,f7.3,f5.1,f10.3,f5.1,f11.3,i4,a1,i2,a2,3x,
     &             a1,i2,a2,3x,i2)
            if (iso.ne.isowant) goto 324
            gf=10.**(gf)
            chie=abs(chie)*kayser2eV
            wave=wave*10.d0
            write(identification,1084) transup,translow,ivup,ivlow,
     &                                 levelup,levellow,xjlow
1084        format('''',a1,a1,1x,i2,1x,i2,1x,a2,1x,a2,f4.1,'''')
          else if (readC2) then
            read(10,1006,end=99) ivup,ivlow,jup,gf,wave,chie,system
            xjup=float(jup)
          else if (readTiO.or.readZrO.or.readLaO) then
c ZrO, LaO have no lab wavenumbers, but the identification is similar to TiO
C LaO 12 15 111 + OP12  0.8330E-02   9629.819  1.960 A
            lab='   '
            read(10,1008,end=99,advance='yes') ivup,ivlow,
     &               jlow,branch2,gf,wave,chie,system,lab
            if (branch2(2:2).eq.'P'.or.branch2(3:3).eq.'P') then
              jup=jlow-1
            else if (branch2(2:2).eq.'Q'.or.branch2(3:3).eq.'Q') then
              jup=jlow
            else if (branch2(2:2).eq.'R'.or.branch2(3:3).eq.'R') then
              jup=jlow+1
            else
              print*,branch2,'problem!'
              stop
            endif
            xjup=float(jup)
            gf=gf*fact
           write(identification,1083) ivup,ivlow,branch2(1:4),jlow,
     &                                system,lab(1:1)
1083       format('''',i2,i3,1x,a4,i4,1x,a1,1x,a1,'''')
1006  format(3x,i3,i3,i4,8x,e11.0,x,f10.0,x,f6.0,x,A)
1007  format(3x,i3,i3,i4,3x,a4,1x,e11.0,x,f10.0,x,f6.0,x,a)
1008  format(3x,i3,i3,i4,3x,a4,1x,e11.0,x,f10.0,x,f6.0,x,a,1x,a3)
          else if (readSiS) then
c SiS from Cami et al. (2009) ApJ 690, L122
            read(10,1009,end=99) wave,chie,A,ivlow,ivup,jlow,jup
1009        format(f12.6,f14.6,f15.8,i4,i4,i4,i4)
            xjup=float(jup)
c A is the Einstein transition probability is s-1
            gf=1.499*(2.*xjup+1.)*A/wave/wave
            chie=chie*kayser2eV
            write(identification,1086) ivup,ivlow,jup,jlow
1086       format('''',i3,i3,i4,i4,'''')
          else if (readTiOVALD) then
* newer TiO for VALD, with also files with laboratory wavelengths
* the format is very different from the original TiO format, and includes radiative damping
*     3164.8974 5.488800E-12   3473.5603   0   7.0   7.0 -1  35061.0197  15   8.0   8.0  1 1.816113E+07 'TiO a f  R    '
            lab='   '
            read(10,*,end=99) wave,gf,chie,ivlow,xjlow,nlow,symlow,Eup,
     &                        ivup,xjup,nup,symup,raddamp,idvald
            jlow=int(xjlow)
            chie=chie*kayser2eV
            branch2(1:4)=idvald(9:12)
            lab(1:1)=idvald(14:14)
            if (idvald(5:7).eq.'a b') then
              system='D'
* delta
            else if (idvald(5:7).eq.'a c') then
              system='B'
* beta
            else if (idvald(5:7).eq.'a f') then
              system='F'
* a-f
            else if (idvald(5:7).eq.'d b') then
              system='P'
* phi
            else if (idvald(5:7).eq.'X A') then
              system='G'
* gamma
            else if (idvald(5:7).eq.'X B') then
              system='H'
* gamma prime
            else if (idvald(5:7).eq.'X C') then
              system='A'
* alpha
            else if (idvald(5:7).eq.'X E') then
              system='E'
* epsilon
            else if (idvald(5:7).eq.'E B') then
              system='C'
* E-B
            endif
            write(identification,1083) ivup,ivlow,branch2(1:4),jlow,
     &                                 system,lab(1:1)
          else
            read(10,1001,end=99) jup,gf,wave,chie
            xjup=float(jup)
          endif
          if (readC2querci.or.readSiH.or.readTiOVALD) then
* already air_lambda
            lambda=wave
          else 
* must convert to air_lambda
            if ((readVO.and.whichvo.eq.'B').or.
     &          (readCaH.and.whichcah.eq.'B')) then
              lambda=wave
            else
              lambda=1.d8/wave
            endif
            x=lambda*1.e-4
            x2=x*x
            refrind=1.+1.e-6* ( 64.328 + 29498.1/(146.-(1./x2)) + 
     &                          255.4/(41.-(1./x2)) )
            lambda=lambda/refrind
          endif
          if(gf.gt.limite) then

          gflog=log10(gf)
          gup=2.*xjup+1.
* maybe not jup, but jlow. The difference is in any case small.
          if (lambda.ge.lambdamin) then
            if (lambda.le.lambdamax) then
             if (gflog.gt.-9.999) then
               write(30,20) lambda,chie,gflog,fdamp,gup,raddamp,
     &                      levup,levlow,eqw,eqwerr,
     &                      identification
             else
               write(30,21) lambda,chie,gflog,fdamp,gup,raddamp,
     &                      levup,levlow,eqw,eqwerr,
     &                      identification
             endif
cc            else if (.not.increasewave) then
cc              goto 99
            endif
          endif
*
* find out strongest line
          strongest=max(strongest,gflog+log10(lambda)-
     &              chie*5040./tempselect)
*

          endif
        goto 100
99      notend=.false.
100     continue
      enddo
*
* count lines in output and replace '123' by correct number
*    BPz 4/6-99
      nline=0
      rewind(30)
      do while (.true.)
        read(30,*,end=998) lambda,chie,gflog
        if ((gflog+log10(lambda)-chie*5040./tempselect).ge.
     &        (weakratio+strongest)) then
           nline=nline+1
        endif
      enddo
998   continue

      if (nline.ne.0) then
      rewind(30)
      open(outfil,file=fileout,status='unknown')
      write(outfil,19) elstring,ion,nline
      write(outfil,18) comment
      do while (.true.)
        read(30,10,end=999) blabla
        read(blabla,*) lambda,chie,gflog
        if ((gflog+log10(lambda)-chie*5040./tempselect).ge.
     &        (weakratio+strongest)) then
          write(outfil,10) blabla
        endif
      enddo
999   continue
      endif
      print*,comment
      print*,nline,' lines written'
      close(30)
      stop

10    format(a)
18    format('''',a,'''')
19    format('''',a21,'''',x,i2,x,i7)
20    format(f10.3,x,f6.3,x,f6.3,x,f5.2,x,f6.1,x,e8.2,1x,
     & '''',a,'''',1x,'''',a,'''',1x,f4.1,1x,f4.1,1x,a20)
21    format(f10.3,x,f6.3,x,f6.2,x,f5.2,x,f6.1,x,e8.2,1x,
     & '''',a,'''',1x,'''',a,'''',1x,f4.1,1x,f4.1,1x,a20)
28    format(i2,i1,f12.6,e10.3,e10.3,f5.4,f5.4,f10.4,f4.2,f8.6,
     +i3,i3,a9,
     +i1,i1,a1,a1,f4.1,a1,i1,i1,i1,i2,i2,i2)
1001  format(10x,i3,8x,e11.0,x,f10.0,x,f6.0)
1002  format(10x,i3,8x,e15.0,x,f10.0,x,f6.0)
1005  format(10x,i3,8x,e11.0,x,f10.0,x,f6.0,x,A)
1234  format(a6,1x,a10,1x,a6,' iel ',i4,' ion ',i2,1x,i9,' lines')

      end
