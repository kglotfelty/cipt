#
#
#
#
#

from ciao_contrib.cipt.enhanced_region import *

"""
Convert ds9 shapes to CIAO stacks of regions


annulus(3938.5052,4158.5084,0,31.085122,62.170244,93.255366,124.34049,155.42561)
   : annulus( x, y, r0, r1, r2, r3, ... rN)

ellipse(4270.4871,4030.506,0,0,36,25.999999,72,51.999999,108,77.999998,144,104,42.150903)
   : ellipse( x, y, m0, n0, m1, n1, m2, n2, ... , mN, nN, angle)

box(4656.5148,4106.5044,0,0,62.664,74.664,125.328,149.328,187.992,223.992,44.993503)
   : box( x, y, m0, n0, m1, n1, m2, n2, ... , mN, nN, angle)

panda(3874.5018,3550.4989,3.4160865e-09,359.9995,4,73.913415,147.82683,2)
   : panda( x, y, angle_min, angle_max, #angle, rmin, rmax, #rad)

epanda(4222.4993,3518.4984,3.4160865e-09,359.9995,3,0,0,120,196,6,3.4160865e-09)
   : epanda( x, y, angle_min, angle_max, #angle, rmin_major, rmin_minr, rmax_maj, rmax_min, #rad, angle)

bpanda(4686.4929,3384.499,3.4160865e-09,359.9995,2,71.998,105.998,143.996,211.996,2,29.993503)
   : bpanda( x, y, angle_min, angle_max, #angle, rmin_major, rmin_minr, rmax_maj, rmax_min, #rad, angle)

"""


def ds9_annulus( *args ):
    """
    annulus(3938.5052,4158.5084,0,31.085122,62.170244,93.255366,124.34049,155.42561)
        : annulus( x, y, r0, r1, r2, r3, ... rN)
    """

    if (len(args)<4):
        raise RuntimeError( "not engough parameters" )

    xx = args[0]
    yy = args[1]
    retval = []
    for i in range( 3, len(args),1 ):
        retval.append( annulus( xx, yy, args[i-1], args[i] ) )

    return retval


def ds9_ellipse( *args ):
    """
    ellipse(4270.4871,4030.506,0,0,36,25.999999,72,51.999999,108,77.999998,144,104,42.150903)
      : ellipse( x, y, m0, n0, m1, n1, m2, n2, ... , mN, nN, angle)
    """

    if (len(args)<5):
        raise RuntimeError( "not engough parameters" )

    xx = args[0]
    yy = args[1]
    angle = args[-1] # Last value is angle

    retval = []
    for i in range( 4, len(args)-1, 2 ):
        inner = ellipse( xx, yy, args[i-2], args[i-1], angle )
        outer = ellipse( xx, yy, args[i], args[i+1], angle )
        retval.append( outer-inner )

    return retval


def ds9_box( *args ):
    """
    box(4656.5148,4106.5044,0,0,62.664,74.664,125.328,149.328,187.992,223.992,44.993503)
       : box( x, y, m0, n0, m1, n1, m2, n2, ... , mN, nN, angle)

    """
    if (len(args)<5):
        raise RuntimeError( "not engough parameters" )

    xx = args[0]
    yy = args[1]
    angle = args[-1] # Last value is angle

    retval = []
    for i in range( 4, len(args)-1, 2 ):
        inner = box( xx, yy, args[i-2], args[i-1], angle )
        outer = box( xx, yy, args[i], args[i+1], angle )
        retval.append( outer-inner )

    return retval


def ds9_panda( xx, yy, start_angle, stop_angle, num_angle, inner, outer, num_rad):
    """
    panda(3874.5018,3550.4989,3.4160865e-09,359.9995,4,73.913415,147.82683,2)
       : panda( x, y, angle_min, angle_max, #angle, rmin, rmax, #rad)
    """

    if start_angle > stop_angle:
        stop_angle = stop_angle + 360

    da = (stop_angle-start_angle)/float(num_angle)

    dr = (outer-inner)/float(num_rad)
    
    retval = []
    for aa in range( num_angle ):        
        ss = sector( xx, yy, start_angle+(aa*da), start_angle+((aa+1)*da) )
        for rr in range( num_rad ):
            anl = annulus( xx, yy, inner+(rr*dr), inner+((rr+1)*dr) )
            retval.append( anl*ss )

    return retval


def ds9_epanda( xx, yy, start_angle, stop_angle, num_angle, mjr_inner, mnr_inner, mjr_outer, mnr_outer, num_rad, ang):
    """
    epanda(4222.4993,3518.4984,3.4160865e-09,359.9995,3,0,0,120,196,6,3.4160865e-09)
       : epanda( x, y, angle_min, angle_max, #angle, rmin_major, rmin_minr, rmax_maj, rmax_min, #rad, angle)
    """
    if start_angle > stop_angle:
        stop_angle = stop_angle + 360

    da = (stop_angle-start_angle)/float(num_angle)
    dm = (mjr_outer-mjr_inner)/float(num_rad)
    dj = (mnr_outer-mnr_inner)/float(num_rad)
    
    retval = []
    for aa in range( num_angle ):        
        ss = sector( xx, yy, start_angle+(aa*da), start_angle+((aa+1)*da) )
        for rr in range( num_rad ):
            inner = ellipse( xx, yy, mjr_inner+(rr*dm), mnr_inner+(rr*dj), ang )
            outer = ellipse( xx, yy, mjr_inner+((rr+1)*dm), mnr_inner+((rr+1)*dj), ang )
            retval.append( (outer-inner)*ss)
    
    return retval
    

def ds9_bpanda( xx, yy, start_angle, stop_angle, num_angle, mjr_inner, mnr_inner, mjr_outer, mnr_outer, num_rad, ang):
    """
    bpanda(4686.4929,3384.499,3.4160865e-09,359.9995,2,71.998,105.998,143.996,211.996,2,29.993503)
        : bpanda( x, y, angle_min, angle_max, #angle, rmin_major, rmin_minr, rmax_maj, rmax_min, #rad, angle)
    """
    if start_angle > stop_angle:
        stop_angle = stop_angle + 360

    da = (stop_angle-start_angle)/float(num_angle)
    dm = (mjr_outer-mjr_inner)/float(num_rad)
    dj = (mnr_outer-mnr_inner)/float(num_rad)
    
    retval = []
    for aa in range( num_angle ):        
        ss = sector( xx, yy, start_angle+(aa*da), start_angle+((aa+1)*da) )
        for rr in range( num_rad ):
            inner = box( xx, yy, mjr_inner+(rr*dm), mnr_inner+(rr*dj), ang )
            outer = box( xx, yy, mjr_inner+((rr+1)*dm), mnr_inner+((rr+1)*dj), ang )
            retval.append( (outer-inner)*ss)
    
    return retval
    

