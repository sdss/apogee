#!/usr/bin/env python
#
# COORDS.PY - coordinate utility functions.
#

from __future__ import print_function

__authors__ = 'David Nidever <dnidever@noao.edu>'
__version__ = '20190723'  # yyyymmdd

import numpy as np
from scipy.spatial import cKDTree
from . import utils

def rotsph(lon,lat,clon,clat,anode=None,reverse=False,original=False):
    '''
    This rotates a spherical coordinate system to a new pole

    I got the equations for this from the paper
    Calabretta et al. 2002, A&A, 395, 1077
    Equation 5.

    Also, see rotate_lb.pro that uses a matrix method
    and for which you can specify the equator you'd like.
    rotsph.pro is faster than rotate_lb.pro
    By default, the origin is the point where the two equators
    cross (unless =ANODE is set).
    This should give you the same result (within ~1E-10")
    rotate_lb,lon,lat,[clon,clat],[clon+90,0.],nlon,nlat

    Parameters
    ----------
    lon       Array of longitudes to be rotated
    lat       Array of latitudes to be rotated
    clon      Longitude of the new NORTH POLE in the old coordinate system
    clat      Latitude of the new NORTH POLE in the old coordinate system
    =anode    The "Ascending Node" which is the longitude (in the new
             system) of the first point where the old equator cross
             the new one.  This sets the zero-point of the new
             longitude.  By default the zero-point of the new
             coordinate system is the point where the two equators
             cross.
    /original Set the new longitude zero-point to be clon (if clat>0)
             and clon+180 (if clat<0).  This is the way it was
             originally done.  DON'T USE WITH "ANODE"
    /stp      Stop at the end of the program
    /reverse  The reverse operation.  In that case (nlon,nlat) should be input
           as (lon,lat). E.g.

           rotsph,ra,dec,cra,cdec,nlon,nlat
           rotsph,nlon,nlat,cra,cdec,nra,ndec,/reverse
           
           (ra,dec) and (nra,ndec) should be identical to 1E-10.

    Returns
    -------
    nlon  Array of rotated longitudes
    nlat  Array of rotated latitudes

    '''

    radeg = 180.0/np.pi

    alphap = np.array(clon/radeg)
    deltap = np.array(clat/radeg)
    phip = np.array(90.0/radeg)
    if original: phip = np.array(180.0/radeg)   # original way
    thetap = np.array(90.0/radeg)

    # By default the origin of the new coordinate system is the point
    # where the two equators cross for the first time
    #  Unless /original is set.  Then the zero-point is at clon
    #   (if clat>0) and clon+180 (if clat<0)

    # NORMAL
    if reverse is False:
        alpha = np.array(lon/radeg)
        delta = np.array(lat/radeg)

        # arg(x,y) but atan(y,x)
        phi = phip + np.arctan2( -np.cos(delta)*np.sin(alpha-alphap), np.sin(delta)*np.cos(deltap)- \
                                 np.cos(delta)*np.sin(deltap)*np.cos(alpha-alphap) )

        theta = np.arcsin( utils.limit((np.sin(delta)*np.sin(deltap)+np.cos(delta)*np.cos(deltap)*np.cos(alpha-alphap)),-1,1) )

        # Preparing the output
        nlon = phi*radeg
        nlat = theta*radeg

        # Ascending Node
        #  By default the origin of nlon is the point where the two equators
        #  cross the first time
        if anode is not None: nlon += anode

    # REVERSE
    else:
        phi = np.array(lon/radeg)
        theta = np.array(lat/radeg)

        # Ascending Node
        if anode is not None: phi = (lon-anode)/radeg

        # arg(x,y) but atan(y,x)
        alpha = alphap + np.arctan2( -np.cos(theta)*np.sin(phi-phip), np.sin(theta)*np.cos(deltap) - \
                                     np.cos(theta)*np.sin(deltap)*np.cos(phi-phip))
        delta = np.arcsin( np.sin(theta)*np.sin(deltap) + np.cos(theta)*np.cos(deltap)*np.cos(phi-phip) )

        # Preparing the output
        nlon = alpha*radeg
        nlat = delta*radeg

    # Want everything less than 360.0
    nlon = nlon % 360.0

    # Make negative points positive
    bd, = np.where(nlon < 0.0)
    if len(bd)>0:
        nlon[bd] = nlon[bd]+360.0

    return nlon, nlat


def rotsphcen(lon,lat,clon,clat,polar=False,gnomic=False,reverse=False):
    '''
    This is very similar to rotsph.pro except that the coordinates
    input are not for the north pole but for the new equator.
    Everything is in DEGREES.
    
    Parameters
    ----------
    lon       Array of longitudes to be rotated
    lat       Array of latitudes to be rotated
    clon      Longitude of the new EQUATOR in the old coordinate system
    clat      Latitude of the new EQUATOR in the old coordinate system
    /polar    Return polar coordinates (rad,phi) instead of LON/LAT.
    phi starts at North.
    /gnomic   Also do a gnomic (tangent plane) projection.
    /reverse  The reverse operation.  In that case (nlon,nlat) should be input
    as (lon,lat). E.g.

           rotsphcen,ra,dec,cra,cdec,nlon,nlat
           rotsphcen,nlon,nlat,cra,cdec,nra,ndec,/reverse
           
           (ra,dec) and (nra,ndec) should be identical to 1E-10.

    Returns
    -------
    nlon  Array of rotated longitudes.  If /polar then this is PHI
       the polar angle (measured from N toward E).
       
    nlat  Array of rotated latitudes.  If /polar then this is RAD
       the polar radial distance.
    '''

    radeg = 180.0/np.pi

    # NOT polar coordinates
    if (polar is False) and (gnomic is False):

        # Get coordinates for the north pole
        np_lon = np.array(clon)
        np_lat = np.array(clat+90.0)
        if (np_lat > 90.0):
            np_lon = np.array(clon+180.0)
            np_lon = np_lon % 360.0
            np_lat = 90.0-clat

        # Run rotsph.pro
        # NORMAL
        if reverse is False:
            nlon, nlat = rotsph(lon,lat,np_lon,np_lat,original=True)
        # REVERSE
        else:
            nlon, nlat = rotsph(lon,lat,np_lon,np_lat,reverse=True,original=True)

            # need to flip them around by 180 deg b/c the zero-point
            #  is set by the NP lon
            if ((clat+90.0) > 90.0):
                nlon = (nlon+180.0) % 360.0
                nlat = -nlat

        # Make the longitudes continuous
        nlon = (nlon+180.0) % 360.0
        nlon = nlon-180.0


    # POLAR or GNOMIC
    else:

        # Making polar coordinates

        #-----------------
        # NORMAL
        #------------------
        if reverse is False:
            # Run rotsph.pro and specify the center of the field (the origin) as the
            #  the Npole
            phi, theta = rotsph(lon,lat,clon,clat,original=True)
            # phi is now going clockwise and from South
            orig_phi = phi
            phi = -phi+180.0      # counterclockwise
            phi = phi % 360.0
            rad = 90.0-theta

            # Making gnomic projection
            if gnomic:
                # Scale the radius
                rad = radeg * np.cos(theta/radeg)/np.sin(theta/radeg)
                # Now convert from gnomic polar to X/Y
                # phi is from N toward E
                # x = R*sin(phi)
                # y = R*cos(phi)
                nlon = rad*np.sin(phi/radeg)
                nlat = rad*np.cos(phi/radeg)

            # Output polar coordinates
            if polar:
                nlon = phi
                nlat = rad

        #-----------------
        # REVERSE
        #-----------------
        else:

            # Polar
            if polar:
                phi = lon
                rad = lat
                theta = 90.0-rad

            # Removing gnomic projection
            if gnomic:
                # Now convert from X/Y to gnomic polar
                # phi is from N toward E
                # x = R*sin(phi)
                # y = R*cos(phi)
                #nlon = rad*sin(phi/radeg)
                #nlat = rad*cos(phi/radeg)
                rad = np.sqrt(lon**2.0+lat**2.0)
                phi = radeg*np.arctan2(lon,lat)      # in degrees
                # Scale the radius
                #rad = radeg * cos(theta/radeg)/sin(theta/radeg)
                theta = radeg*np.arctan(radeg/rad)   # in degrees

            # phi is now going clockwise and from South
            phi = -phi+180.0       # reverse phi

            #Run rotsph.pro and specify the center of the field (the origin) as the
            #  the Npole
            nlon, nlat = rotsph(phi,theta,clon,clat,reverse=True,original=True)

    return nlon, nlat


def doPolygonsOverlap(xPolygon1, yPolygon1, xPolygon2, yPolygon2):
    """Returns True if two polygons are overlapping."""

    # How to determine if two polygons overlap.
    # If a vertex of one of the polygons is inside the other polygon
    # then they overlap.
    
    n1 = len(xPolygon1)
    n2 = len(xPolygon2)
    isin = False

    # Loop through all vertices of second polygon
    for i in range(n2):
        # perform iterative boolean OR
        # if any point is inside the polygon then they overlap   
        isin = isin or isPointInPolygon(xPolygon1, yPolygon1, xPolygon2[i], yPolygon2[i])

    # Need to do the reverse as well, not the same
    for i in range(n1):
        isin = isin or isPointInPolygon(xPolygon2, yPolygon2, xPolygon1[i], yPolygon1[i])

    return isin

def isPointInPolygon(xPolygon, yPolygon, xPt, yPt):
    """Returns boolean if a point is inside a polygon of vertices."""
    
    # How to tell if a point is inside a polygon:
    # Determine the change in angle made by the point and the vertices
    # of the polygon.  Add up the delta(angle)'s from the first (include
    # the first point again at the end).  If the point is inside the
    # polygon, then the total angle will be +/-360 deg.  If the point is
    # outside, then the total angle will be 0 deg.  Points on the edge will
    # outside.
    # This is called the Winding Algorithm
    # http://geomalgorithms.com/a03-_inclusion.html

    n = len(xPolygon)
    # Array for the angles
    angle = np.zeros(n)

    # add first vertex to the end
    xPolygon1 = np.append( xPolygon, xPolygon[0] )
    yPolygon1 = np.append( yPolygon, yPolygon[0] )

    wn = 0   # winding number counter

    # Loop through the edges of the polygon
    for i in range(n):
        # if edge crosses upward (includes its starting endpoint, and excludes its final endpoint)
        if yPolygon1[i] <= yPt and yPolygon1[i+1] > yPt:
            # if (P is  strictly left of E[i])    // Rule #4
            if isLeft(xPolygon1[i], yPolygon1[i], xPolygon1[i+1], yPolygon1[i+1], xPt, yPt) > 0: 
                 wn += 1   # a valid up intersect right of P.x

        # if edge crosses downward (excludes its starting endpoint, and includes its final endpoint)
        if yPolygon1[i] > yPt and yPolygon1[i+1] <= yPt:
            # if (P is  strictly right of E[i])    // Rule #4
            if isLeft(xPolygon1[i], yPolygon1[i], xPolygon1[i+1], yPolygon1[i+1], xPt, yPt) < 0: 
                 wn -= 1   # a valid up intersect right of P.x

    # wn = 0 only when P is outside the polygon
    if wn == 0:
        return False
    else:
        return True

def isLeft(x1, y1, x2, y2, x3, y3):
    # isLeft(): test if a point is Left|On|Right of an infinite 2D line.
    #   From http://geomalgorithms.com/a01-_area.html
    # Input:  three points P1, P2, and P3
    # Return: >0 for P3 left of the line through P1 to P2
    # =0 for P3 on the line
    # <0 for P3 right of the line
    return ( (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1) )


# from astroML
def crossmatch(X1, X2, max_distance=np.inf,k=1):
    """Cross-match the values between X1 and X2

    By default, this uses a KD Tree for speed.

    Parameters
    ----------
    X1 : array_like
        first dataset, shape(N1, D)
    X2 : array_like
        second dataset, shape(N2, D)
    max_distance : float (optional)
        maximum radius of search.  If no point is within the given radius,
        then inf will be returned.

    Returns
    -------
    dist, ind: ndarrays
        The distance and index of the closest point in X2 to each point in X1
        Both arrays are length N1.
        Locations with no match are indicated by
        dist[i] = inf, ind[i] = N2
    """
    X1 = np.asarray(X1, dtype=float)
    X2 = np.asarray(X2, dtype=float)

    N1, D = X1.shape
    N2, D2 = X2.shape

    if D != D2:
        raise ValueError('Arrays must have the same second dimension')

    kdt = cKDTree(X2)

    dist, ind = kdt.query(X1, k=k, distance_upper_bound=max_distance)

    return dist, ind

# from astroML, modified by D. Nidever
def xmatch(ra1, dec1, ra2, dec2, dcr=2.0,unique=False):
    """Cross-match angular values between RA1/DEC1 and RA2/DEC2

    Find the closest match in the second list for each element
    in the first list and within the maximum distance.

    By default, this uses a KD Tree for speed.  Because the
    KD Tree only handles cartesian distances, the angles
    are projected onto a 3D sphere.

    This can return duplicate matches if there is an element
    in the second list that is the closest match to two elements
    of the first list.

    Parameters
    ----------
    ra1/dec1 : array_like
        first dataset, arrays of RA and DEC
        both measured in degrees
    ra2/dec2 : array_like
        second dataset, arrays of RA and DEC
        both measured in degrees
    dcr : float (optional)
        maximum radius of search, measured in arcsec.
        This can be an array of the same size as ra1/dec1.

    Returns
    -------
    ind1, ind2, dist: ndarrays
        The indices for RA1/DEC1 (ind1) and for RA2/DEC2 (ind2) of the
        matches, and the distances (in arcsec).
    """
    X1 = np.vstack((ra1,dec1)).T
    X2 = np.vstack((ra2,dec2)).T
    
    X1 = X1 * (np.pi / 180.)
    X2 = X2 * (np.pi / 180.)
    if utils.size(dcr)>1:
        max_distance = (np.max(dcr) / 3600) * (np.pi / 180.)
    else:
        max_distance = (dcr / 3600) * (np.pi / 180.)

    # Convert 2D RA/DEC to 3D cartesian coordinates
    Y1 = np.transpose(np.vstack([np.cos(X1[:, 0]) * np.cos(X1[:, 1]),
                                 np.sin(X1[:, 0]) * np.cos(X1[:, 1]),
                                 np.sin(X1[:, 1])]))
    Y2 = np.transpose(np.vstack([np.cos(X2[:, 0]) * np.cos(X2[:, 1]),
                                 np.sin(X2[:, 0]) * np.cos(X2[:, 1]),
                                 np.sin(X2[:, 1])]))

    # law of cosines to compute 3D distance
    max_y = np.sqrt(2 - 2 * np.cos(max_distance))
    k = 1 if unique is False else 10 
    dist, ind = crossmatch(Y1, Y2, max_y, k=k)
    
    # convert distances back to angles using the law of tangents
    not_inf = ~np.isinf(dist)
    x = 0.5 * dist[not_inf]
    dist[not_inf] = (180. / np.pi * 2 * np.arctan2(x,
                                  np.sqrt(np.maximum(0, 1 - x ** 2))))
    dist[not_inf] *= 3600.0      # in arcsec
    
    # Allow duplicates
    if unique is False:

        # no matches
        if np.sum(not_inf)==0:
            return [], [], [np.inf]
        
        # If DCR is an array then impose the max limits for each element
        if utils.size(dcr)>1:
            bd,nbd = utils.where(dist > dcr)
            if nbd>0:
                dist[bd] = np.inf
                not_inf = ~np.isinf(dist)
    
        # Change to the output that I want
        ind1 = np.arange(len(ra1))[not_inf]
        ind2 = ind[not_inf]
        mindist = dist[not_inf]

    # Return unique one-to-one matches
    else:

        # no matches
        if np.sum(~np.isinf(dist[:,0]))==0:
            return [], [], [np.inf]
        
        done = 0
        niter = 1
        # Loop until we converge
        while (done==0):

            # If DCR is an array then impose the max limits for each element
            if utils.size(dcr)>1:
                bd,nbd = utils.where(dist[:,0] > dcr)
                if nbd>0:
                    for i in range(nbd):
                        dist[bd[i],:] = np.inf

            # no matches
            if np.sum(~np.isinf(dist[:,0]))==0:
                return [], [], [np.inf]

            # closest matches
            not_inf1 = ~np.isinf(dist[:,0])
            ind1 = np.arange(len(ra1))[not_inf1]
            ind2 = ind[:,0][not_inf1]
            mindist = dist[:,0][not_inf1]
            if len(ind2)==0:
                return [], [], [np.inf]
                #print('no elements')
                #import pdb; pdb.set_trace()
            index = utils.create_index(ind2)
            # some duplicates to deal with
            bd,nbd = utils.where(index['num']>1)
            if nbd>0:
                torem = []            
                for i in range(nbd):
                    indx = index['index'][index['lo'][bd[i]]:index['hi'][bd[i]]+1]
                    # keep the one with the smallest minimum distance
                    si = np.argsort(mindist[indx])
                    if index['num'][bd[i]]>2:
                        torem += list(indx[si[1:]])    # add list
                    else:
                        torem.append(indx[si[1:]][0])  # add single element
                ntorem = utils.size(torem)
                # For each object that was "removed" and is now unmatched, check the next possible
                # match and move it up in the dist/ind list if it isn't INF
                for i in range(ntorem):
                    # There is a next possible match 
                    if ~np.isinf(dist[torem[i],niter-1]):
                        ind[torem[i],:] = np.hstack( (ind[torem[i],niter:].squeeze(), np.repeat(-1,niter)) )
                        dist[torem[i],:] = np.hstack( (dist[torem[i],niter:].squeeze(), np.repeat(np.inf,niter)) )
                    # All INFs
                    else:
                        ind[torem[i],:] = -1
                        dist[torem[i],:] = np.inf
            else:
                ntorem = 0

            niter += 1
            # Are we done, no duplicates or hit the maximum 10
            if (ntorem==0) | (niter>=10): done=1
                                
    return ind1, ind2, mindist

def dist(x1, y1, x2, y2):
    """ Calculate Euclidian distance between two sets of points."""

    if (utils.size(x1) != utils.size(y1)):
        raise ValueError('x1/y1 must have same number of elements')
    if (utils.size(x2) != utils.size(y2)):
        raise ValueError('x2/y2 must have same number of elements')    

    return np.sqrt( (x1-x2)**2 + (y1-y2)**2 )
    

def sphdist(lon1, lat1, lon2, lat2):
    """Calculate the angular distance between two sets of points.

    Parameters
    ----------
    lon1/lat1 : scalar or array_like
        first dataset, arrays of LON/LAT or RA/DEC
        both measured in degrees
    lon2/lat2 : scalar array_like
        second dataset, arrays of LON/LAT or RA/DEC
        both measured in degrees

    Returns
    -------
    dist: ndarrays
        The angular distance in degrees.
    """

    if (utils.size(lon1) != utils.size(lat1)):
        raise ValueError('lon1/lat1 must have same number of elements')
    if (utils.size(lon2) != utils.size(lat2)):
        raise ValueError('lon2/lat2 must have same number of elements')    
    
    # From this website;
    # http://www2.sjsu.edu/faculty/watkins/sphere.htm

    cosa = np.sin(np.deg2rad(lat1))*np.sin(np.deg2rad(lat2))+np.cos(np.deg2rad(lat1))*np.cos(np.deg2rad(lat2))*np.cos(np.deg2rad(lon1-lon2))
    dist = np.rad2deg(np.arccos(cosa))

    return dist
