import numpy as np
import matplotlib.pyplot as plt

import pdb
import astroplan
from astroplan import Observer
from astroplan import FixedTarget
from astroplan.plots import plot_airmass
from astroplan.plots import plot_parallactic
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.coordinates import get_moon, get_sun
from astropy import units

def airmass(ra=0., dec=0., obs='apo', date='2019-10-01',name='object',plot=False,tz='US/Mountain') :

    """  Get airmass table for specified object position, observatory, date
    """

    # set the site
    site=Observer.at_site(obs,timezone=tz)

    # set the objects
    if type(ra) is float :
        obj=FixedTarget(name=name,coord=SkyCoord(str(ra)+'d',str(dec)+'d'))
    else :
        obj=FixedTarget(name=name,coord=SkyCoord(ra+'h',dec+'d'))

    time = Time('{:s} 00:00:00'.format(date),scale='utc',
           location=(site.location.lon,site.location.lat),precision=0)
    sunset =site.sun_set_time(time)
    civil = site.twilight_evening_civil(time)
    nautical = site.twilight_evening_nautical(time)
    astronomical = site.twilight_evening_astronomical(time)
    for t in [sunset,civil,nautical,astronomical] :
        t.format='isot'
        t.precision=0
    print('Observatory: ',obs)
    print('Sunset: ',sunset)
    print('Civil twilight: ',civil)
    print('Nautical twilight: ',nautical)
    print('Astronomical twilight: ',astronomical)


    # loop over all UTC hours for this date (would prefer local!)
    print('{:8s}{:8s}{:8s}{:8s}{:8s}{:9s}{:8s}{:8s} {:16s}{:20s}'.format(
          'Local','UT','LST','HA','Airmass','ParAng','Phase','Moon Alt','Moon RA', 'Moon DEC'))
    for hr in np.arange(24) :
        time = Time('{:s} {:d}:00:00'.format(date,hr),scale='utc',
               location=(site.location.lon,site.location.lat),precision=0)
        sun=get_sun(time)
        if site.sun_altaz(time).alt.value > 10 : continue

        moon=get_moon(time)

        val=site.altaz(time,obj).secz
        if val < 0 :
            airmass = '     ...'
        else :
            airmass = '{:8.2f}'.format(val)
        val=site.moon_altaz(time).alt.value 
        if val < 0 :
            moonalt = '     ...'
        else :
            moonalt = '{:8.2f}'.format(val)
        lst=time.sidereal_time('mean').hms
        ha=site.target_hour_angle(time,obj)
        ha.wrap_angle=180 *units.deg
        local=site.astropy_time_to_datetime(time)
        print('{:02d}:{:02d}  {:02d}:{:02d}  {:02d}:{:02d} {:3d}:{:02d} {:8s} {:8.2f} {:8.2f} {:8s} {:s} {:s}'.format(
               local.hour,local.minute,
               time.datetime.hour,time.datetime.minute,
               int(round(lst[0])),int(round(lst[1])),
               int(round(ha.hms[0])),int(abs(round(ha.hms[1]))),
               airmass,
               #site.altaz(time,obj).secz,
               site.parallactic_angle(time,obj).deg,
               astroplan.moon_illumination(time),moonalt,
               str(moon.ra.to_string(units.hour)), str(moon.dec)))

    if plot: 
        fig,ax = plt.subplots(2,1)
        fig.subplots_adjust(hspace=0.001)

        plot_airmass(obj,site,time,ax=ax[0])
        plot_parallactic(obj,site,time,ax=ax[1])
