"""
Routines to plot various geometric shapes

Shapes are plotted with chips add_region() routine which allows them 
to be filled.

Note:  everything is assuming a cartesian coordinate system with 
angles measure CCW from the +X axis.  This is consistent with
CIAO regions stored in Physical coordinates.

"""

__all__ = [ "plot_ellipse", "plot_box", "plot_circle", "plot_rectangle", 
            "plot_polygon", 
            "plot_annulus", "plot_pie", "plot_region", 
            "set_plot_with_lines"] 


import numpy as np


class MatplotlibPlotEngine():
    _UseLines = False
    
    def __init__(self):
        import matplotlib.pylab as plt
        self.plt = plt
    
    def begin_plot(self):
        fig,ax = self.plt.subplots()
        ax.set_aspect(1)
        
    def plot(self, xx, yy, **style):
        if self._UseLines:
            # TODO?: append first point onto end to close shape
            self.plt.plot(xx,yy, **style)
        else:
            self.plt.fill(xx,yy, **style)

    def end_plot(self):
        self.plt.show()

_Plotter = MatplotlibPlotEngine()


def set_plot_with_lines( use_lines ):
    """
    Switch between using lines or regions to draw shapes.
    
    use_lines = True  -> will use lines
    use_lines = False -> will use regions
    """
    _Plotter._UseLines = True


def simplify_polygon( xx, yy, delta ):
    """
    Stolen from Doug's simplify_polygon routine
    """
   
    def filter_point(x, y, delta2, ix):
        """Filter the data around point ix"""

        # We don't want all points near d2, just the "connected" points.
        #
        d2 = (x[ix+1:] - x[ix])**2 + (y[ix+1:] - y[ix])**2
        idx = [ix]
        p = 0
        while p < d2.size and d2[p] <= delta2:
            idx.append(p + ix + 1)
            p += 1

        return (x[idx].mean(), y[idx].mean(), ix + len(idx))      
    

    def clean_points(x, y, delta):
        """Return a set of 'cleaned' points."""

        ox = []
        oy = []
        i = 0
        delta2 = delta * delta

        while i < x.size:
            (nx,ny,ni) = filter_point(x, y, delta2, i)
            ox.append(nx)
            oy.append(ny)
            i = ni

        return (np.asarray(ox), np.asarray(oy))


    xp, yp = clean_points( np.asarray(xx), np.asarray(yy), delta )
    
    return xp, yp


def add_region( xx, yy, delta=0, **style ):
    """
    Actually draw the region.  All regions are drawn as polygons.
    Uses the _Plotter object (wrapper around matplotlib).
    
    If the plot fails, then it simplifies the polygon
    until it either works or max delta is met.
    """
    
    if delta > 4:
        raise RuntimeError("Problem plotting shape")

    try:
        _Plotter.plot(xx,yy, **style)
    except Exception as e:
        delta = delta + 0.1
        xp, yp = simplify_polygon( xx, yy, delta=delta)
        add_region( xp, yp, delta=delta, **style )
    

def plot_ellipse( x0, y0, r1, r2, rotang=None, **style ):
    """
    plot a rotated ellipse
    
    >>> plot_ellipse( 10, 10, 20, 15, 45)
    >>> plot_ellipse( 10, 10, 20, 15, 45, {'fill.style' : 'none' })

    """

    # Okay, very long way of doing this, not terrible efficent but it works.
    # 
    # We generate a grid of X values and solve the ellipse eqn for y.
    # First for the top half (+sqrt()) and then in reverse order
    # for the bottom half (-sqrt()).
    # 
    # There will be x-values that are invalid, those will return np.nan and
    # will be filtered out.  (numpy error disable): TODO: re-enable at end
    # 
    # Once the ellipse is created, the x,y pairs are then rotated.
    # 
    # This requires a fine grid for the x-axis, which is normally fine;
    # but if the ellipse is small, it may generate a "degenerate polygon"
    # that seems to have overlapping/duplicate data points.

    if None == rotang:
        rotang = 0.0

    nbins = 200

    orig_err = np.seterr( all='ignore' )
    
    # Figure out the X grid based on the lenghts of the radii
    r12 = r1*r1  # r^2
    r22 = r2*r2

    #
    # Generate a grid of X values.  Compute a grid of angles
    # in the range from 0 to 180 degrees, and then take cos()
    # of the values to give us an X grid that is dense at the ends
    # to capture the curvature of the ellipse.
    #
    angls = np.arange( 0, nbins+1)
    angls = (180.0/nbins)*angls
    angls_rad = np.deg2rad( angls )
    unit_x = np.cos( angls_rad )  # values from -1 to 1

    # now scale the unit X values based on the radius
    xp = unit_x * r1


    # Make ellipse around 0,0 and then shift when done.
    # Top of ellipse
    yp = np.sqrt( r22 * ( 1 - ((xp)**2)/r12))
    
    # Bottom of ellipse
    xm = xp[::-1]
    ym = -np.sqrt( r22 * ( 1 - ((xm)**2)/r12))

    # Restore numpy error handling
    np.seterr(**orig_err)
    
    # Combine them together, omit 1st point to avoid dup, degenerate polys
    xp = np.append(xp, xm[1:])
    yp = np.append(yp, ym[1:])
    
    # Remove where Y value is outside of ellipse
    xy = [ ii for ii in zip(xp,yp) if np.isfinite(ii[1]) ]
    x = np.array([i[0] for i in xy])
    y = np.array([i[1] for i in xy])

    # rotation matrix
    cosa = np.cos( np.deg2rad( -rotang) )
    sina = np.sin( np.deg2rad( -rotang) )    
    rx = ( x * cosa + y * sina ) + x0
    ry = (-x * sina + y * cosa ) + y0

    # plot it
    add_region(rx,ry, **style)


def plot_circle( x0, y0, r, **style):
    """
    Plot a circular region
    
    >>> plot_circle( 0,0, 10 )
    >>> plot_circle( 0,0, 10, "edge.style=longdash")
    >>> plot_circle( 4096, 4096, 500, { 'fill.style' : 'userfill1' })

    """
    plot_ellipse( x0, y0, r, r, 0, **style)


def plot_box( x0, y0, xl, yl=None, rotang=None, **style ):
    """
    Plot a box region
    
    A box is described by the center (x0,y0), the total length of the
    X and Y axis, and a rotation angle

    >>> plot_box( 100, 100, 50, 20, 10 )
    >>> plot_box( 100, 100, 30, 30, 0, {'fill.style' : chips_wave} )  # a square

    """
    
    if None == rotang:
        rotang = 0
    
    if None == yl:
        yl = xl
    
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
    add_region(rx,ry, **style)


def plot_rectangle( xl, yl, xh, yh, **style):
    """
    Plot a rectangle
    
    A rectangle is described by a lower-left corner (xl,yl) 
    and an upper right corner (xh, yh).  No rotation.
    
    >>> plot_rectangle( 0,0, 10, 10 ) # a square
    >>> plot_rectangle( 0,0, 100, 10 , "edge.style=shortdash")
    
    """
    x0 = (xh+xl)/2.0
    y0 = (yh+yl)/2.0
    lx = (xh-xl)
    ly = (yh-yl)
    plot_box( x0, y0, lx, ly, 0, **style )


def plot_polygon( xvals, yvals, **style):
    """
    Plot a polygon
    
    A polygon is described by a series of X and Y values.

    >>> x = [ 1, 2, 3 ]
    >>> y = [ 0, 1, 0 ]
    >>> plot_polygon( x, y )  # a triangle
    
    >>> import numpy as np
    >>> xx = np.arange(yy)
    >>> yy = np.random.rand(100)
    >>> plus = filter( lambda x: x[1]>0.5, zip(xx,yy))
    >>> minus = filter( lambda x: x[1]<=0.5, zip(xx,yy))
    >>> plus.extend(minus[::-1])
    >>> x = [i[0] for i in plus] 
    >>> y = [i[1] for i in plus]    
    >>> plot_polygon(x,y, 'fill.color=firebrick')

    """    
    add_region(xvals, yvals, **style)

    
def plot_pie( x0, y0, rad_in, rad_out, ang_start, ang_stop , **style ):
    """
    Plot a pie shaped region
    
    A pie is defined by a center point, an inner and out radii, and 
    a start and stop angle.
    
    >>> plot_pie( 0,0, 1,4,-45, 45)
    """
    
    nbins = 100

    while (ang_start > ang_stop) : ang_start = ang_start - 360 

    da = (ang_stop - ang_start)/(1.0*nbins)

    angles = np.arange( ang_start, ang_stop+da, da )
    
    angles_rad = np.deg2rad( angles )

    if rad_in > 0:
        x_in = rad_in * np.cos(angles_rad)
        y_in = rad_in * np.sin(angles_rad)
    else:
        # special case of inner radius=0
        x_in = np.array( [0.0] )
        y_in = np.array( [0.0] )
    
    x_out = rad_out * np.cos( angles_rad[::-1])
    y_out = rad_out * np.sin( angles_rad[::-1])
    
    xx = np.append( x_in, x_out )+x0
    yy = np.append( y_in, y_out )+y0
    
    add_region( xx, yy, **style )



def plot_annulus( x0, y0, rad_in, rad_out, **style ):
    """
    Plot a pie shaped region
    
    A pie is defined by a center point, an inner and out radii


    >>> plot_annulus(10,10,10,20)

    Due to the way the annulus is drawn, there will be a line
    connecting the outer and inner radii along the +X axis, zooming
    into that location will reveal a small gap.  The line ca be
    removed by setting the region opacity to 1 or by setting
    the edge line style to none
    
    >>> set_region("edge.style=none")

    """
    
    plot_pie( x0, y0, rad_in, rad_out, 0.1, 359.9, **style )


def plot_point(x0, y0):
    raise NotImplementedError("Please use add_point()")    
    




def plot_region( regRegion, **style ):
    """
    Loop over each shape in a CXC region and plot it.
    
    Note: no distinction between include and exclude shapes.    
    """
    try:
        tst = len( regRegion )
        tst = regRegion[0].shapes
    except:
        raise TypeError("Input parameter does not appear to be a valid region object")

    _Plotter.begin_plot()
 
    for rr in regRegion:
        ss = rr.shapes[0]

        if hasattr(ss,"name"):
            # Translate to CXCRegion
            ss.shape = ss.name
            ss.rad = ss.radii
            ss.xx = ss.xpoints
            ss.yy = ss.ypoints
            ss.ang = ss.angles

        shape = ss.shape
        if 'annulus' == shape:
            plot_annulus( ss.xx[0], ss.yy[0], ss.rad[0], ss.rad[1], **style )
        elif 'box' == shape:
            plot_box( ss.xx[0], ss.yy[0], ss.rad[0], ss.rad[1], **style )
        elif 'circle' == shape:
            plot_circle( ss.xx[0], ss.yy[0], ss.rad[0], **style )
        elif 'ellipse' == shape:
            plot_ellipse( ss.xx[0], ss.yy[0], ss.rad[0], ss.rad[1], ss.ang[0], **style )
        elif 'field' == shape:
            pass # no fields please
        elif 'pie' == shape:
            plot_pie( ss.xx[0], ss.yy[0], ss.rad[0], ss.rad[1], ss.ang[0], ss.ang[1], **style )
        elif 'point' == shape:
            pass # no points please
        elif 'polygon' == shape:
            plot_polygon( ss.xx, ss.yy, **style )
        elif 'rectangle' == shape:
            plot_rectangle( ss.xx[0], ss.yy[0], ss.xx[1], ss.yy[1], **style )
        elif 'rotbox' == shape:
            plot_box( ss.xx[0], ss.yy[0], ss.rad[0], ss.rad[1], ss.ang[0], **style )
        elif 'sector' == shape:
            plot_pie( ss.xx[0], ss.yy[0], 0, 99999, ss.ang[0], ss.ang[1], **style )
        else:
            raise TypeError("Unknown shape type : {}".format(shape))

    _Plotter.end_plot()

