*
      subroutine CIAh2h(omega,t,propac)
*
* BPz 6/10/95. interpolates H2-H CIAs from A. Borysow.
* iw is the index of the wavenumber just larger than omega
* it is the index of the temperature just smaller than T
* Mod for new data 26 Feb. 2003 /BE
*
      implicit none

      integer j,iw,it,maxwave,maxtemp,i,nwave,ntemp,ii
      parameter (maxwave=4500,maxtemp=7)
      real omega,t,propac,omegalast,qt,qw,h2hop(maxtemp),x,exp10
      real wave(maxwave),temp(maxtemp),h2h(maxtemp, maxwave)
      character ident*10
      logical first

      data omegalast /0.00/
      data iw        /1/
      data first     /.true./

*      save omegalast,iw,h2hop,first,wave,temp,h2h
*      save ident,nwave,ntemp

* ATTENTION EXTRAPOLATION AUX OMEGAS EN DEHORS DE LA TABLE
* DEHORS=outside of
* DE BORYSOW!!  ET extrapolation a T < 1000K et >7000K!!
*
      if (first) then
        open(33,file='DATA/CIA.H2-H.dat', status='old')
        read(33,*) ident
        read(33,*) nwave
        read(33,*) ntemp
        if (ntemp.gt.maxtemp) stop 'CIAh2h: maxtemp too small'
        if (nwave.gt.maxwave) stop 'CIAh2h: maxwave too small'
        read(33,*) (temp(i),i=1,ntemp)
        read(33,*)
        read(33,*)
        do i=nwave,1,-1
          read(33,*) wave(i),(h2h(j,i),j=1,ntemp)
        enddo
        first=.false.
        close(33)
      endif

* no extrapolation to higher wavenumbers!!
      if (omega.gt.wave(1)) then
        propac=0.
        return
      endif

      if (omega.gt.omegalast) then
        iw=1
      endif
      if (omega.ne.omegalast) then
        do while ((omega.le.wave(iw+1)).and.(iw.lt.nwave))
          iw=iw+1
        enddo
        iw=min(iw,nwave-1)
        qw=(omega-wave(iw))/(wave(iw+1)-wave(iw))
        do j=1,ntemp
          h2hop(j)=qw*h2h(j,iw+1)+(1.-qw)*h2h(j,iw)
        enddo
      endif
      it=1
      do ii=1,ntemp-1
        if (t.ge.temp(ii)) it=ii
      enddo
      it=min(it,ntemp-1)
      qt=(t-temp(it))/(temp(it+1)-temp(it))
      propac=qt*log(h2hop(it+1))+(1.-qt)*log(h2hop(it))
      if(propac.gt.-70) then
        propac=exp(propac)
      else
        propac=1.e-30
      endif
      omegalast=omega

      return
      end
