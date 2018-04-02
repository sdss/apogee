program faltbo

  implicit none

  integer                                          :: iprf,ii,jj,ios,isize,isize2
  integer                                          :: resample,itype,ipad,jbeg,jend
  real(kind=kind(1.d0))                            :: FWHM,FWHM2,somme,dlam,dlam2
  real(kind=kind(1.d0)), dimension(600000)         :: lambda,Fnorm
  real(kind=kind(1.d0)), allocatable, dimension(:) :: Fconv1,lambda2,F2,dl,prof,Fconv2
  real(kind=kind(1.d0)), parameter                 :: pi = 4.*atan(1.)
  character(len=800)                                :: inspec,outfil,command
  logical                                          :: velocity

  resample=1
  ! read inputs and open files 
  print*,'input file ?'
  read(5,' (A) ') inspec
  print*,'output file ?'
  read(5,' (A) ') outfil 

  open(unit=10,file=inspec,status='old',iostat=ios,action='read')
  if (ios==0) print*,'input file opened'
  if (ios/=0) stop "problem in file opening ..."
  open(unit=20,file=outfil,status='unknown')

  print*,'FWHM (in km/s or mA, <0 if velocity) ?'
  read(5,*) FWHM
  print*,'profile type ? (1 = exp - 2 = gauss - 3 = rad.-tan. - 4 = rot.)'
  read(5,*) iprf         ! profile type
!
! we don't give the choice anymore.
! we suppose that spectra are always 2 column: 1st is lambda, the next 1 is Fnorm or Fabs
! BPz 30/03-2012
!

  velocity=.false.
  if (FWHM.LT.0.) velocity=.true.

! assumes no headers in inspec ... 

  ! read input spectrum 
  ii=1
  do while (ios == 0)
    read(unit=10,iostat=ios,fmt=*) lambda(ii),Fnorm(ii)
    ii=ii+1
    if (ii.gt.600000) stop 'Increase dimension'
  end do

  isize=ii-2
  print*,isize,' lines read in input file'
  allocate(Fconv1(isize))
  allocate(Fconv2(isize))

  ! define how many 0s necessary before and after and extend input spectrum 
  dlam=lambda(2)-lambda(1)
  dlam2=lambda(isize)-lambda(isize-1)
  if (velocity) then
     ipad=int(-50.*FWHM*lambda(1)/3.e5/dlam)+1
  else
     ipad=int(50.*FWHM/1.e3/dlam)+1
  endif
     print*,'ipad, isize ',ipad,isize
  isize2=isize+2*ipad
  print*,'ipad, isize2 = ',ipad,isize2 
  allocate(lambda2(isize2),F2(isize2))
  do ii=1,ipad
     lambda2(ii)=lambda(1)-dlam*(ipad-ii+1)
     F2(ii)=0.
     lambda2(isize+ipad+ii)=lambda(isize)+dlam2*ii
     F2(isize+ipad+ii)=0.
  end do

! loop to convolve 2 spectra (if needed)
    do ii=1,isize
     lambda2(ipad+ii)=lambda(ii)
     F2(ipad+ii)=Fnorm(ii)
    end do

    FWHM2=FWHM/1.e3
    do ii=1,isize,resample
     ! at current wavelength: check step, and allocate conv. profile accordingly
     if (velocity) FWHM2=-FWHM*lambda(ii)/3.e5 ! km/s -> A 
     jbeg=1
     jend=2*ipad+1
     allocate(dl(jend))
     allocate(prof(jend))
     do jj=jbeg,jend
        dl(jj)=lambda2(ii+jj-1)-lambda2(ii+ipad)
     end do

     ! call relevant procedure for convolution profile
     if (iprf==1 .or. iprf==2) then
       call gauss(ipad,dl,iprf,FWHM2,prof,jbeg,jend)
     else if (iprf==3) then
       call radtan(ipad,dl,FWHM2,prof,jbeg,jend)
     else if (iprf==4) then
       call rota(ipad,dl,FWHM2,lambda(ii),prof,jbeg,jend)
     else 
       stop 'undefined profile'
     endif
     ! if (iprf==5) call FTS(FWHM)

     ! normalization of convolution profile 
     somme=0.
!     print*,'jbeg jend',jbeg, jend
     do jj=jbeg,jend
        ! somme=somme+(prof(jj)+prof(jj+1))*abs(dl(jj))/2.
        somme=somme+prof(jj)
     end do
     if (somme.le.0) stop 'problem with profile normalization ...'
     somme=1./somme
     do jj=jbeg,jend
        prof(jj)=prof(jj)*somme
     end do

     ! convolution : f'(ii) = sum (C(j)*f(j))    j = ii-ipad -> ii+ipad 
     somme=0.
     do jj=jbeg,jend
        somme=somme+prof(jj)*F2(ii+jj-1) !
     end do
     Fconv2(ii)=somme
     deallocate(dl)
     deallocate(prof)
    end do

!  print*,'convolution done ! Now writing output ...'

  ! output, with resampling 
  do ii=1,isize,resample
    write(20,200) lambda(ii),Fconv2(ii)
  end do

200  format(f13.4,2(1x,1pe12.5))


end program faltbo

! ***************************************************************************

subroutine gauss(ipad,dl,iprf,FWHM,prof,jbeg,jend)

  ! exponential and gaussian profiles

  implicit none

  integer                                   :: ipad,iprf,i,jbeg,jend
  real(kind=kind(1.d0)),dimension(2*ipad+1) :: prof,dl
  real(kind=kind(1.d0))                     :: FWHM,const,y

  if (iprf==1) const=1.38629/FWHM
  if (iprf==2) const=1.66511/FWHM

! we require at least the 3 central points of the profile.
  jbeg=ipad
  jend=ipad+2

  do i=1,2*ipad+1
     prof(i)=0.
     y=dl(i)*const
     if (iprf==1) then
       if (exp(-abs(y)).gt.1.e-6) then
         prof(i)=exp(-abs(y))
         jbeg=min(jbeg,i)
         jend=max(jend,i)
       endif
     else if (iprf==2) then
       if (exp(-y*y).gt.1.e-6) then
         prof(i)=exp(-y*y)
         jbeg=min(jbeg,i)
         jend=max(jend,i)
       endif
     endif
  end do

  return

end subroutine gauss

! ***************************************************************************

subroutine radtan(ipad,dl,FWHM,prof,jbeg,jend)

  ! radial-tangential profile; macroturbulent velocity = 1.433*FWHM
  ! cf. Gray 1978, Solar Phys. 59, 193

  implicit none 

  integer                                       :: i,j,ipad,jbeg,jend
  real(kind=kind(1.d0)),dimension(2*ipad+1)     :: prof,dl
  real(kind=kind(1.d0))                         :: FWHM,width,delta,di,pp
  real(kind=kind(1.d0)),dimension(40),parameter ::  rtf = & 
       (/1.128,.939,.773,.628,.504,.399,.312,.240,.182,.133,   &
       .101,.070,.052,.037,.024,.017,.012,.010,.009,.007,.006, &
       .005,.004,.004,.003,.003,.002,.002,.002,.002,.001,.001, & 
       .001,.001,.001,.001,.000,.000,.000,.000 /)

! delta is the wavelength distance between given RTF points. 
  width=FWHM*1.433
  delta=width/10.

  do j=1,2*ipad+1
     prof(j)=0.
  end do

  prof(ipad+1)=rtf(1)

  j=ipad+2
  di=delta
  do i=2,35
     di=delta*float(i-1)
     do while (dl(j).le.di)
        pp=log(rtf(i))+(log(rtf(i-1))-log(rtf(i)))*(di-dl(j))/delta
        prof(j)=exp(pp)
        jend=j
        j=j+1
     end do
  end do

  j=ipad
  do i=2,36
     di=delta*float(i-1)
     do while (abs(dl(j)).le.di)
        pp=log(rtf(i))+(log(rtf(i-1))-log(rtf(i)))*(di-abs(dl(j)))/delta
        prof(j)=exp(pp)
        jbeg=j
        j=j-1
     end do
  end do

  return

end subroutine radtan

! ***************************************************************************

subroutine rota(ipad,dl,FWHM,lambda,prof,jbeg,jend)

  ! rotational profile; eps is wavelength dependant linear (in mu)
  ! limb-darkening coefficient found for F V -K IV stars (improve !!!)
  ! FWHM is v*sin i in wavelength units 

  implicit none 

  integer                                   :: i,ipad,jprof,nmax,jbeg,jend
  real(kind=kind(1.d0))                     :: FWHM,dlam,dlambda,lambda,eps
  real(kind=kind(1.d0)),dimension(2*ipad+1) :: prof,x
  real(kind=kind(1.d0)),dimension(2*ipad)   :: dl
  real(kind=kind(1.d0)),parameter           :: pi = 4.*atan(1.)

  eps=1.-0.3*lambda/5000.
  dlambda=FWHM
  jbeg=2*ipad+1
  jend=1
  do i=1,2*ipad+1
     if (abs(dl(i)).le.dlambda) then
        jbeg=min(jbeg,i)
        jend=max(jend,i)
        prof(i)=(2.*(1.-eps)*sqrt(1.-(dl(i)/dlambda)**2)+pi*0.5*eps* &
          (1.-(dl(i)/dlambda)**2))/(pi*dlambda*(1.-eps/3.))
     endif
  end do

  return

end subroutine rota

! ***************************************************************************

subroutine FTS(FWHM)

  ! FTS profile (sin(x)/x) - from Kurucz- program broaden.f
! obsviously not ready for use.!!

  implicit none 

  real(kind=kind(1.d0)) :: FWHM ! assumed in km/s ...

  ! x=(i-1.)*vstep/FWHM*2.*1.8954942
  ! red(i)=sin(x)/x*exp(-0.06*x**2)

  return

end subroutine FTS
