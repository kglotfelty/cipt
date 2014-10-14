
__all__ = [ 'Box', 'Circle', 'Ellipse', 'Field', 'Pie', 'Point', 'Polygon', 'Region', 'Sector',
            'box', 'circle', 'ellipse', 'field', 'pie', 'point', 'polygon', 'region', 'sector' ]


import region as _region

class FilterMask(object):
    """
    The base class for the FilterMask objects
    
    The only piece of information stored here is the string 
    representation of the shape.  When shapes are combined
    that is all that remains. 
    
    The individual shape components have separate
    derived classes below to handle their construction.   
    
    """
    def __init__(self, try_string = None ):

        if None == try_string:
            try_string = self.__repr__()
        try:
            self.region_obj = _region.regParse( try_string )
            if self.region_obj is None:
                raise Exception("Region parse error")
            self.as_string = _region.regRegionString(self.region_obj).replace(")-",")*!").replace(")|", ")+").replace(")&",")*").lower() 
        except:
            raise RuntimeError("Illegal region syntax '{}'".format(try_string))
    

    def area(self):
        """
        
        """
        return _region.regArea( self.region_obj )
        

    def inside(self, xx, yy ):
        """
        """
        return _region.regInsideRegion( self.region_obj, xx, yy )
    
    def write(self, filename):        
        with open(filename, "w") as fp:
            fp.write( self.__repr__() )
            fp.write("\n")
    
    def __add__(self,other):        
        cpts = self.as_string
        cpto = other.as_string
        ret = "+".join([cpts,cpto])
        return FilterMask(ret)

    def __mul__(self,other):
        # (a+b)*c = a*c + b*c , not a + b*c
        cpts = self.as_string.split("+")
        cpto = other.as_string.split("+")
        ret = "+".join( [ "{}*{}".format(x, o) for x in cpts for o in cpto])
        return FilterMask(ret)

    def __sub__(self,other):
        # Inversion is hard.  We only do the simple case
        # where all the shapes are either anded or ored together.

        cpts = self.as_string.split("+")
        cpto = other.as_string

        # DeMorgan's Law:  !(A*B) = !A + !B
        # This is easy when just * or just + are present.  If both
        # then have to unroll expression to get down to atomic
        # shapes: !(A*B+C) = !(A*B)*!C = ((!A+!B)*!C) = !A*!C + !B*!C
        if '*' in cpto and '+' in cpto:
            raise NotImplementedError("Negation of complex regions is not implemented")

        if '+' in cpto:
            cpto = cpto.replace("+","*!")
            ret = "+".join( [ "{}*!{}".format(x, cpto) for x in cpts ])
        else:
            cpto = cpto.split("*")
            ret = "+".join( [ "{}*!{}".format(x, o) for x in cpts for o in cpto])
        ret = ret.replace("!!","") # two negatives make a postive

        return FilterMask(ret)

    def __repr__(self):
        return self.as_string

    def __str__(self):
        return self.__repr__()
        
    def plot(self):
        raise NotImplementedError("Plotting not supported at this level")
        

class FilterMaskDefbyRadiusAngle(FilterMask):
    """
    Most shapes are defined by some combination of center, lengths, angles

    This class deals with those an how to construct the region string
    from them.
    """
    def __init__(self, xlen, ylen=None, center=None, angle=None, angle2=None ):
        self.xlen = xlen
        self.ylen = ylen
        self.center = (0,0) if center == None else center
        if len(self.center) != 2:
            raise RuntimeError("Illegal center value, must have 2 elements")
        self.xcenter = self.center[0]
        self.ycenter = self.center[1]
        self.angle = angle
        self.angle2 = angle2

        FilterMask.__init__(self)        
    
    def __repr__( self ):
        vals = filter( lambda x: x is not None, [ self.xcenter, self.ycenter, self.xlen, self.ylen, self.angle, self.angle2 ] )
        inner = ",".join( [str(x) for x in vals ] )
        retval = "{}({})".format( self.isa, inner )
        return retval


class Box( FilterMaskDefbyRadiusAngle ):
    def __init__(self, xlen, ylen=None, center=None, angle=None ):        
        # todo (x,(cx,cy)) should also be workable
        ylen = xlen if ylen is None else ylen
        self.isa = "box"
        FilterMaskDefbyRadiusAngle.__init__( self, xlen, ylen, center, angle)
        
    def plot(self):
        _plot_box( self.xcenter, self.ycenter, self.xlen, self.ylen, self.angle )
        

class Ellipse( FilterMaskDefbyRadiusAngle ):
    def __init__(self, xlen, ylen, center=None, angle=0 ):        
        self.isa = "ellipse"
        FilterMaskDefbyRadiusAngle.__init__( self, xlen, ylen, center, angle)
    def plot( self ):
        _plot_ellipse( self.xcenter, self.ycenter, self.xlen, self.ylen, self.angle )

class Annulus( FilterMaskDefbyRadiusAngle ):
    def __init__(self, inner, outer, center=None ):        
        self.isa = "annulus"
        FilterMaskDefbyRadiusAngle.__init__( self, inner, outer, center)

class Circle( FilterMaskDefbyRadiusAngle ):
    def __init__(self, inner, center=None ):        
        self.isa = "circle"
        FilterMaskDefbyRadiusAngle.__init__( self, inner, center=center)
    def plot(self):
        _plot_ellipse( self.xcenter, self.ycenter, self.xlen, self.xlen, 0 )
            
class Point( FilterMaskDefbyRadiusAngle ):
    def __init__(self, center=None ):        
        self.isa = "point"
        FilterMaskDefbyRadiusAngle.__init__( self, None, center=center)
    def plot(self):
        pass # points not plotted

class Pie( FilterMaskDefbyRadiusAngle ):
    def __init__(self, xlen, ylen, start,stop, center=None ):        
        self.isa = "pie"
        FilterMaskDefbyRadiusAngle.__init__( self, xlen, ylen, center=center, angle=start, angle2=stop)

class Sector(FilterMaskDefbyRadiusAngle ):
    def __init__(self, start, stop, center=None ):        
        self.isa = "sector"
        FilterMaskDefbyRadiusAngle.__init__( self, None, center=center, angle=start, angle2=stop)

class Rectangle( FilterMaskDefbyRadiusAngle ):
    def __init__(self, xmin, ymin, xmax, ymax ):        
        self.isa = "rectangle"
        FilterMaskDefbyRadiusAngle.__init__( self, xmax, ylen=ymax, center=(xmin,ymin) )

    def plot(self):
        dx = self.xmax-self.xmin
        dy = self.ymax-self.ymin
        xc = (self.xmax+self.xmin)/2.0
        yc = (self.ymax+self.ymin)/2.0
        _plot_box( xc, yc, dx, dy, 0.0 )

class Field( FilterMask ):
    def __init__(self):
        self.isa = "field"
        FilterMask.__init__(self)

    def __repr__(self):
        return "field()"
        
    def plot(self):
        pass  # fields are not plotted
    
    
class Polygon(FilterMask):
    """
    Polygon( (x1,y1), (x2,y2), ..., (xn,yn) )
    
    """
    def __init__(self, *args):
        self.isa = "polygon"

        ll = len(args)
        if ll <3 :
            raise RuntimeError("Must supply at least 3 x,y verticies")
        if not all(map( lambda x: len(x) == 2, args )):
            raise RuntimeError("Each element must have an x,y coordinate")

        self.polypoints = [ x for xy in args for x in xy ]
        FilterMask.__init__(self)
    
    def __repr__(self):
        coords = ",".join(map(str, self.polypoints))        
        return "{}({})".format( self.isa, coords )
        
    def plot(self):
        _plot_polygon(self.polypoints)

    
class Region(FilterMask):
    """
    handle region(filename)
    """
    def __init__( self, filename ):
        #reg = _region.regParse( "region({})".format(filename))
        #ss = _region.regRegionString( reg )
        FilterMask.__init__( self, try_string="region({})".format(filename) )


def box( xc, yc, xlen, ylen ):
    """
    """
    return Box( xlen, ylen, center=(xc,yc))
    
def ellipse( xc, yc, xlen, ylen, angle=0):
    """
    """
    return Ellipse( xlen, ylen, center=(xc,yc), angle=angle)

def annulus( xc, yc, inner, outer ):
    """
    """
    return Annulus( inner, outer, center=(xc,yc))

def circle( xc, yc, rad ):
    """
    """
    return Circle( rad, center=(xc,yc))
    
def point( xc, yc ):
    """
    """
    return Point( center=(xc,yc))
    
def pie( xc, yc, inner, outer, start, stop ):
    """
    """
    return Pie( inner, outer, start, stop, center=(xc,yc))
    

def sector( xc, yc, start, stop):
    """
    """
    return Sector( start, stop, center=(xc,yc))

def rectangle( xmin, ymin, xmax, ymax ):
    """
    """
    return Rectangle(xmin,xmas,ymin,ymax)


def field():
    """
    """
    return Field()

def polygon( *args ):
    """
    """
    xs = [ args[i] for i in xrange( 0,len(args),2)]
    ys = [ args[i] for i in xrange( 1,len(args),2)]
    return Polygon( *zip( xs, ys ))

def region( filename ):
    """
    """
    return Region(filename)
    
# ---------------------------
# --------------------------------------
#

try:
    from pychips import add_region
except:
    def add_region( *arg, **kwargs):
        raise NotImplementedError("Chips is not installed.  Plotting is not supported")
        
import numpy as np


def _plot_box( x0, y0, xl, yl, rotang, style="" ):
    """
    """
    
    dx=xl/2.0
    dy=yl/2.0
    
    xx = np.array([ x0-dx, x0+dx, x0+dx, x0-dx])
    yy = np.array([ y0-dy, y0-dy, y0+dy, y0+dy])
    
    dx = xx-x0
    dy = yy-y0    
    cosa = np.cos( np.deg2rad( -rotang) )
    sina = np.sin( np.deg2rad( -rotang) )    
    rx = ( dx * cosa + dy * sina ) + x0
    ry = (-dx * sina + dy * cosa ) + y0

    # plot it
    add_region(rx,ry, style)


def _plot_polygon( points, style="" ):
    """
    """
    x = points[::2]
    y = points[1::2]
    add_region(x,y, style)


def _plot_ellipse( x0, y0, r1, r2, rotang, style="", nbins=100 ):
    """
    Okay, very long way of doing this, not terrible efficent but it works.

    """
    np.seterr(all='ignore')
    if nbins < 1:
        raise RuntimeError("Cannot plot this shape")
    
    r12 = r1*r1
    r22 = r2*r2
    mr = max([r1,r2])
    dx = 2.0*max([r1,r2])/nbins

    # Top of ellipse
    xp = np.arange( -r1-dx, r1+dx, dx )+x0
    yp = y0 + np.sqrt( r22 * ( 1 - ((xp-x0)**2)/r12))
    
    # Bottom of ellipse
    xm = xp[::-1]
    ym = y0 - np.sqrt( r22 * ( 1 - ((xm-x0)**2)/r12))
    
    # Combine them together, omit 1st point to avoid dup, degenerate polys
    xp = np.append(xp, xm[1:])
    yp = np.append(yp, ym[1:])
    
    # Remove where Y value is outside of ellipse
    xy = filter( lambda i: np.isfinite( i[1]), zip(xp,yp))
    x = np.array([i[0] for i in xy])
    y = np.array([i[1] for i in xy])

    # rotation matrix
    dx = x-x0
    dy = y-y0    
    cosa = np.cos( np.deg2rad( -rotang) )
    sina = np.sin( np.deg2rad( -rotang) )    
    rx = ( dx * cosa + dy * sina ) + x0
    ry = (-dx * sina + dy * cosa ) + y0

    # plot it

    try:
        add_region(rx,ry, style)
    except:
        # may be too many poits, try with fewer
        ibins = int(0.9*nbins)
        _plot_ellipse( x0, y0, r1, r2, rotang, style=style, nbins=ibins ):
