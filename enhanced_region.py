#
# Copyright (C) 2013-2015 Smithsonian Astrophysical Observatory
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

"""
Enhanced CIAO Region Package


The standard CIAO region package has routines to parse regions, and determine
whether a point is inside or outside that region, along with things like
area and extent calculations.

However, once a region is created it cannot be modified or used in combination
with other regions.

This Enhanced CIAO Region Package allows users to create regions,
use them in logical operations, and modify them

>>> from ciao_contrib.cipt.enhanced_region import *
>>> a = circle(1000,1000,50)
>>> b = box(1010,1010,50,100)
>>> c = a+b
>>> print c
circle(1000,1000,50)+box(1010,1010,50,100)
>>> print c.area()
8574.0
>>> d = a*b
>>> print d
circle(1000,1000,50)*box(1010,1010,50,100)
>>> e = a-b
>>> print e
circle(1000,1000,50)*!box(1010,1010,50,100)


"""



__all__ = [ 'annulus', 'box', 'circle', 'ellipse', 'field', 'pie', 'point', 'polygon', 'polygon_vec', 'region', 'rectangle', 'sector' ]

    
_AND_ = "*"
_OR_ = "+"
_NOT_ = "!"
_BLANK_ = ""


from ctypes import *
try:
    dll = cdll.LoadLibrary('libregion.so')
    dml = cdll.LoadLibrary('libascdm.so')
except:
    try: # OSX
        dll = cdll.LoadLibrary('libregion.dylib')
        dml = cdll.LoadLibrary('libascdm.dylib')
    except:
        raise ImportError("Cannot load region library")


# Set return types.  Probably don't need void's, int's, and long's, but
# do it anyway for completeness
dll.regArea.restype = c_double
dll.regCreateEmptyRegion.restype = c_void_p
dll.regInsideRegion.restype = c_int
dll.regUnionRegion.restype = c_void_p
dll.regIntersectRegion.restype = c_void_p
dll.regInvert.restype = c_void_p
dll.regGetNoShapes.restype = c_long
dll.regGetShapeNo.restype = c_void_p
dll.regShapeRadii.restype = c_long
dll.regShapeAngles.restype = c_long
dll.regShapeGetName.restype = c_int
dll.regShapeGetNoPoints.restype = c_long
dll.regShapeGetComponent.restype = c_long
dll.regShapeGetPoints.restype = c_long # unused
dll.regShapeGetAngles.restype = c_long # unused
dll.regShapeGetRadii.restype = c_long # unused
dml.dmRegParse.restype = c_void_p




# Adds plotting capability to each shape
try:
    import chips_contrib.plot_shapes as plt
    def with_plotting(plot_fun, engine=plt):
        """Decorator to check for plotting engine"""
        def plot_wrap(self):
            if not engine:
                return
            plot_fun(self, engine)
        return plot_wrap
except:
    import warnings as warn
    warn.warn("Cannot locate chips_contrib.plot_shape module.  Region plotting will not be available.")
    plt = None
    def with_plotting(plot_fun):
        """No plotting available"""
        def plot_wrap(self):
            pass
        return plot_wrap



def wrap_vals( vv ):
    """Utility to wrap arrays to be pased in byreference"""
    if vv:
        retval = byref(( c_double * len(vv) ) (*vv))
    else:
        retval = POINTER(c_int)()  # NULL pointer
    return retval


class EnhancedShape( object ):
    """
    All regions are composed of simple geometric shapes (even 1 region with only
    1 shape is still a EnhancedRegion() ). 
    
    This object isn't meant to be exposed to users    
    """    

    precise_string = "{:.12g}"  # Format string for double values when printed

    def __init__(self, shape_ptr ):
        """The input is a c_void_p regShape pointer.  Uses the region lib
        API to extract the shape parameters"""
        
        self._ptr = shape_ptr
        self.get_name()
        self.get_points()
        self.get_radii()
        self.get_angle()
        self.get_component()


    def get_name(self):
        """
        Get the name.  This also sets whether the region is inclusive
        or exclusive.
        """
        shape_name = create_string_buffer(100)
        iflag = dll.regShapeGetName(self._ptr, shape_name, 99)
        self.shape = shape_name.value.lower()
        self.include = _NOT_ if 0 == iflag else _BLANK_


    def get_points(self):
        """
        Gets the x,y values.
        
        The string version of the values are stored now too.  
        TODO: wrap w/ a setter that does the string conversion so that
        if class is exposed then string is kept in sync with double values
        """
        npts = dll.regShapeGetNoPoints( self._ptr )
        ox = (c_double * npts)()
        oy = (c_double * npts)()
        dll.regShapeGetPoints( self._ptr, byref(ox), byref(oy), c_long( npts ))
        self.xx = [x for x in ox]
        self.yy = [y for y in oy]


    def get_radii(self):
        """
        Gets the radius/radii 
        
        Also stores the string version (see TODO above)
        """
        nrad = dll.regShapeRadii( self._ptr )
        outr = (c_double * nrad)()
        dll.regShapeGetRadii( self._ptr, byref(outr) )
        self.rad = [r for r in outr]
    

    def get_angle(self):
        """
        Gets the angle/angles
        
        Also stores the string version (see TODO above)
        """
        nang = dll.regShapeAngles( self._ptr )
        oa = (c_double*nang)()
        dll.regShapeGetAngles( self._ptr, byref(oa))
        self.ang = [ a for a in oa ]


    def get_component(self):
        """
        Stores the component value.  The actual value is not used.  All 
        the shapes with the component value are grouped together and AND'ed.
        Diff components are OR'ed.
        """
        self.component = dll.regShapeGetComponent( self._ptr )


    def store_string(self):
        """
        Convert the values to strings using the desired precision/format
        """
        self.xx_str = [ self.precise_string.format(x) for x in self.xx ]
        self.yy_str = [ self.precise_string.format(y) for y in self.yy ]
        self.rad_str = [ self.precise_string.format(y) for y in self.rad ]
        self.ang_str = [ self.precise_string.format(y) for y in self.ang ]
        

    def __str__(self):
        raise NotImplementedError("Implement this in the derived classes")


    @with_plotting
    def plot(self, engine):
        raise NotImplementedError("Implement this in the derived classes")


class Annulus(EnhancedShape):
    """An annulus is defined by x_center, y_center, inner_radius, outer_radius"""
    def __str__( self ):
        self.store_string()
        return "{i}{n}({x0},{y0},{r0},{r1})".format(i=self.include,
            n=self.shape,x0=self.xx_str[0],y0=self.yy_str[0],r0=self.rad_str[0],r1=self.rad_str[1])


    @with_plotting
    def plot(self, engine):
        engine.plot_annulus( self.xx[0], self.yy[0], self.rad[0], self.rad[1] ) 


class Box(EnhancedShape):
    """A box is defined by x_center, y_center, x_length, y_length"""
    def __str__( self ):
        self.store_string()
        return "{i}{n}({x0},{y0},{r0},{r1})".format(i=self.include,
            n=self.shape,x0=self.xx_str[0],y0=self.yy_str[0],r0=self.rad_str[0],r1=self.rad_str[1])


    @with_plotting
    def plot(self,engine):
        engine.plot_box( self.xx[0], self.yy[0], self.rad[0], self.rad[1] ) 


class Circle(EnhancedShape):
    """A circle is defined by x_center, y_center, radius"""
    def __str__( self ):
        self.store_string()
        return "{i}{n}({x0},{y0},{r0})".format(i=self.include,
            n=self.shape,x0=self.xx_str[0],y0=self.yy_str[0],r0=self.rad_str[0])


    @with_plotting
    def plot(self,engine):
        engine.plot_circle( self.xx[0], self.yy[0], self.rad[0])


class Ellipse(EnhancedShape):
    """An ellipse is defined by x_center, y_center, major_axis, minor_axis, and angle"""
    def __str__( self ):
        self.store_string()
        return "{i}{n}({x0},{y0},{r0},{r1},{a0})".format(i=self.include,
            n=self.shape,x0=self.xx_str[0],y0=self.yy_str[0],r0=self.rad_str[0],r1=self.rad_str[1],a0=self.ang_str[0])


    @with_plotting
    def plot(self,engine):
        engine.plot_ellipse( self.xx[0], self.yy[0], self.rad[0], self.rad[1], self.ang[0]) 


class Field(EnhancedShape):
    """A field is defined to be the entire R^2 dataspace"""
    def __str__( self ):
        self.store_string()
        return "{i}{n}()".format( i=self.include, n=self.shape )


    @with_plotting
    def plot(self,engine):
        pass


class Pie(EnhancedShape):
    """A Pie is defined by x_center, y_center, inner_radius, outer_radius, start_angle, stop_angle"""
    def __str__( self ):
        self.store_string()
        return "{i}{n}({x0},{y0},{r0},{r1},{a0},{a1})".format(i=self.include,
            n=self.shape,x0=self.xx_str[0],y0=self.yy_str[0],r0=self.rad_str[0],r1=self.rad_str[1],a0=self.ang_str[0],a1=self.ang_str[1],)


    @with_plotting
    def plot(self,engine):
        engine.plot_pie( self.xx[0], self.yy[0], self.rad[0], self.rad[1], self.ang[0], self.ang[1]) 


class Point(EnhancedShape):
    """A point is defined by x_center, y_center"""
    def __str__( self ):
        self.store_string()
        return "{i}{n}({x0},{y0})".format(i=self.include,
            n=self.shape,x0=self.xx_str[0],y0=self.yy_str[0])


    @with_plotting
    def plot(self,engine):
        pass


class Polygon(EnhancedShape):
    """A Polygon is defined as x1,y1, x2,y2,x3,y3,..."""
    def __str__( self ):
        self.store_string()
        xyxy = ",".join([ "{},{}".format(x,y) for x,y in zip( self.xx_str, self.yy_str ) ])
        return "{i}{n}({xy})".format(i=self.include,n=self.shape,xy=xyxy)


    @with_plotting
    def plot(self,engine):
        engine.plot_polygon(self.xx, self.yy)


class Rectangle(EnhancedShape):
    """A Polygon is defined as lower_left_x,lower_left_y,upper_right_x,upper_right_y"""
    def __str__( self ):
        self.store_string()
        return "{i}{n}({x0},{y0},{x1},{y1})".format(i=self.include,
            n=self.shape,x0=self.xx_str[0],y0=self.yy_str[0],x1=self.xx_str[1],y1=self.yy_str[1])


    @with_plotting
    def plot(self,engine):
        engine.plot_rectangle( self.xx[0], self.yy[0], self.xx[1], self.yy[1] )


class Rotbox(EnhancedShape):
    """A Rotbox is defined as x_center, y_center, x_length, y_length, angle """
    def __str__(self):
        self.store_string()
        return "{i}{n}({x0},{y0},{r0},{r1},{a0})".format(i=self.include,
            n=self.shape,x0=self.xx_str[0],y0=self.yy_str[0],r0=self.rad_str[0],r1=self.rad_str[1],a0=self.ang_str[0])


    @with_plotting
    def plot(self,engine):
        engine.plot_box( self.xx[0], self.yy[0], self.rad[0], self.rad[1], self.ang[0]) 


class Sector(EnhancedShape):
    """A Sector is defined by x_center, y_center, start_angle, stop_angle"""
    def __str__( self ):
        self.store_string()
        return "{i}{n}({x0},{y0},{a0},{a1})".format(i=self.include,
            n=self.shape,x0=self.xx_str[0],y0=self.yy_str[0],a0=self.ang_str[0],a1=self.ang_str[1],)


    @with_plotting
    def plot(self,engine):
        engine.plot_pie( self.xx[0], self.yy[0], 0.0, 999999.9, self.ang[0], self.ang[1]) 



class EnhancedRegion( object ):
    """
    
    A region is made up of a collection of EnhancedShapes that are combined
    with AND and OR logical operators.

    Only the shape_classes shapes are implemented.
    """
    shape_classes = { 'annulus': Annulus,
               'box' : Box,
               'circle' : Circle,
               'ellipse' : Ellipse,
               'field' : Field,
               'pie' : Pie,
               'point' : Point,
               'polygon' : Polygon,
               'rectangle' : Rectangle,
               'rotbox' : Rotbox,
               'sector' : Sector }
               

    def __init__( self, ptr  ):
        """
        EnhancedRegion Object
        
        This is not meant to be exposed to user.  It takes in a c_void_p
        regREGION pointer         
        """
        self._ptr = ptr
        self.load_shapes()
        self.dll = dll # keep a ref to this so it can be used to free shapes during garbage collection


    def __str__( self ):
        """Construct the string value from the shapes & logic"""
        retval = _BLANK_.join( [logic+str(shape) for logic,shape in zip(self.logic,self.shapes) ] )
        return(retval)


    def __repr__(self):
        return self.__str__()

                
    def load_shapes(self):
        """
        Load the shapes into the region.
        """
        self.classify_shapes()
        self.determine_logic()


    @staticmethod
    def get_shape_name( shape_ptr):
        """Get the name of a shape"""
        shape_name = create_string_buffer(100)
        dll.regShapeGetName(shape_ptr, shape_name, 99)
        return shape_name.value.lower()


    def classify_shapes( self ):
        """Create the correct EnhancedShape derived class based on the shape's name."""
        nshapes = dll.regGetNoShapes( self._ptr )
        self.shapes = []
        for i in range(1,nshapes+1):
            shape_ptr = dll.regGetShapeNo( self._ptr, c_long(i) )
            shape_name = self.get_shape_name( shape_ptr)
            if shape_name in self.shape_classes:
                self.shapes.append( self.shape_classes[shape_name](shape_ptr) )
            else:
                raise ValueError("Unknown shape name {}".format(shape_name))


    def determine_logic(self):
        """Determine the logic used to combine shapes.
        
        The same component values are AND'ed together,  Different 
        are OR'd.  The standard convention is that the same component values
        are stored together, so we only need to check for when the value
        changes."""        

        cpt_vals = [s.component for s in self.shapes ]

        self.logic = [ _BLANK_ ] # First shape has no logic
        for i in range( 1, len(cpt_vals)):
            # If component values are equal then &, else |
            if cpt_vals[i]==cpt_vals[i-1]:
                self.logic.append(_AND_)
            else:
                self.logic.append(_OR_)


    def area( self ):
        """Determine the region area. 
        
        TODO: make can pass in field bound box and bin-size parameters.
        There is also a routine that forces pixelated area we could 
        optionally call
        """
        import sys as sys
        DBL_MAX = sys.float_info.max
        fld = wrap_vals( [-DBL_MAX, DBL_MAX] )
        return dll.regArea( self._ptr, fld, fld, c_double(1.0) )


    def write(self, filename):
        """
        Write the region to an ascii file.  Filename is always clobbered!
        
        TODO: optionally replace _OR_ with "\n"
        
        TODO: Write as FITS file.  Only hard thing is the NaN-pad the 
        arrays to max size.
        """
        with open(filename, "w") as fp:
            fp.write( self.__repr__() )
            fp.write("\n")


    def inside( self, x, y ):
        """
        Determine if x,y pair is inside the region
        """
        return (1 == dll.regInsideRegion(self._ptr, c_double(x), c_double(y)))


    def plot( self ):
        """
        Plot the region by plotting each shape
        """
        bad = 0
        for s in self.shapes:
            try:
                s.plot()
            except:
                bad = bad+1

        if bad != 0:
            print "Warning: {} shapes could not be plotted".format(bad)


    def __del__(self):
        """
        Use the self.dll which should exist during garbage collection.
        """
        self.dll.regFree( self._ptr )


    def __add__(self,other):
        """
        Logically OR (Union) two regions together
        """
        if not isinstance( other, EnhancedRegion):
            raise NotImplementedError("Cannot perform region logic with {} type".format(type(other)) )
            
        cpts = self._ptr
        cpto = other._ptr
        return EnhancedRegion( dll.regUnionRegion( cpts, cpto ) )


    def __mul__(self,other):
        """
        Logically AND (Intersect) two regions together
        """
        if not isinstance( other, EnhancedRegion):
            raise NotImplementedError("Cannot perform region logic with {} type".format(type(other)) )

        cpts = self._ptr
        cpto = other._ptr
        return EnhancedRegion( dll.regIntersectRegion( cpts, cpto ) )


    def __sub__(self,other):
        """
        Subtraction is actually a short-hand or AND NOT.
        Other is inverted and then AND'ed with self.
        """
        if not isinstance( other, EnhancedRegion):
            raise NotImplementedError("Cannot perform region logic with {} type".format(type(other)) )

        cpts = self._ptr
        cpto = other._ptr
        invrt = dll.regInvert( cpto )
        return EnhancedRegion( dll.regIntersectRegion( cpts, invrt ) )
        

    def __copy__(self):
        return EnhancedRegion( dll.regCopyRegion( self._ptr))
        

    def __eq__(self, other):
        if not isinstance( other, EnhancedRegion):
            raise NotImplementedError("Cannot perform region logic with {} type".format(type(other)) )

        return dll.regCompareRegion(self._ptr, other._ptr)


    def tweak( self, dx=0, dy=0, stretch=1, rotate=0 ):
        """
        Create a copy of the current region, shape-by-shape.
        
        Then apply tweaks to the x, y, r, angle parameters before
        creating and returning a new region.
        
        Each tweak is applied to each shape separately
        """
        copy_ptr = dll.regCreateEmptyRegion()
        for ii in range( len(self.shapes) ):
            shape = self.shapes[ii]

            reg_math = c_int(0) if _AND_ == self.logic[ii] else c_int(1)        
            reg_inc = c_int(0) if _NOT_ == shape.include else c_int(1)

            copy_xx = [ x+dx for x in shape.xx]
            copy_yy = [ y+dy for y in shape.yy]
            copy_rr = [ r*stretch for r in shape.rad]
            copy_aa = [ a+rotate for a in shape.ang ]

            dll.regAppendShape( copy_ptr,
                                shape.shape,
                                reg_inc, reg_math,
                                wrap_vals( copy_xx ),
                                wrap_vals( copy_yy ),
                                c_long( len(copy_xx) ),
                                wrap_vals( copy_rr ),
                                wrap_vals( copy_aa ),
                                c_int(0), c_int(0) )
        return EnhancedRegion(copy_ptr)
                

def SimpleGeometricShapes( name_as_string, xx, yy, radius, angle ):
    """
    For single shape regions crated with the routines below
    """
    new_ptr = dll.regCreateEmptyRegion()
    dll.regAppendShape( new_ptr,
                        name_as_string ,
                        c_int(1), c_int(1),
                        wrap_vals(xx),
                        wrap_vals(yy),
                        c_long( len(xx) ),
                        wrap_vals( radius ),
                        wrap_vals( angle ),
                        c_int(0), c_int(0) )
    if 0 == dll.regGetNoShapes( new_ptr ):
        raise RuntimeError("Bad region parameters")

    return EnhancedRegion(new_ptr)

# ------------------
# User Interface 

def circle( x, y, r ):
    """circle( x_center, y_center, radius )"""
    return SimpleGeometricShapes( "circle", [x], [y], [r], [] )


def ellipse( x, y, xlen, ylen, angle=0 ):
    """ellipse( x_center, y_center, major_axis, minor_axis[, rotation_angle])
    
    angle is measured in degrees from the +X axis
    """    
    return SimpleGeometricShapes( "ellipse", [x], [y], [xlen, ylen], [angle] )


def box( x,y, xlen, ylen, angle=0 ):
    """box( x_center, y_center, x_length, y_length[, rotation_angle])
    
    angle is measured in degrees from the +X axis
    """
    return SimpleGeometricShapes( "rotbox", [x], [y], [xlen, ylen], [angle] )


def annulus( x, y, inner, outer ):
    """annulus( x_center, y_center, inner_radius, outer_radius )"""
    return SimpleGeometricShapes( "annulus", [x], [y], [inner, outer], [] )


def field():
    """field()"""
    return SimpleGeometricShapes("field", [], [], [], [] )


def pie( x, y, inner, outer, start, stop ):
    """pie( x_center, y_center, inner_radius, outer_radius, start_angle, stop_angle )
    
    angles are measured in degrees from the +X axis
    """
    return SimpleGeometricShapes( "pie", [x], [y], [inner, outer], [start, stop] )


def point(x,y):
    """point( x_center, y_center)"""
    return SimpleGeometricShapes( "point", [x], [y], [], [] )


def rectangle(llx, lly, urx, ury ):
    """rectangle( lower_left_x, lower_left_y, upper_right_x, upper_right_y)"""
    return SimpleGeometricShapes( "rectangle", [llx,urx], [lly,ury], [], [] )


def sector( x, y, start, stop ):
    """sector( x_center, y_center, start_angle, stop_angle
    
    angles are measured in degrees from +X axis
    """
    return SimpleGeometricShapes( "sector", [x], [y], [], [start, stop] )


def polygon( *args ):
    """polygon(x1,y1,x2,y2,x3,y3,...)
    
    Must contain at least 3 points, and same number of x-y pairs.
    An extra point is added if last point is not equal to 1st point.
    """
    x = [args[i] for i in range(0,len(args),2) ]
    y = [args[i] for i in range(1,len(args),2) ]
    if (len(x) != len(y)):
        raise RuntimeError("must have equal number of x,y pairs")
    if (len(x)<3):
        raise RuntimeError("must have at least 3 points to make a polygon")
    return SimpleGeometricShapes( "polygon", x, y, [], [] )


def polygon_vec( x, y ):
    """polygon( x_values, y_values )
    
    Must contain at least 3 points, and same number of x-y pairs.
    An extra point is added if last point is not equal to 1st point.
    """    
    if (len(x) != len(y)):
        raise RuntimeError("must have equal number of x,y pairs")
    if (len(x)<3):
        raise RuntimeError("must have at least 3 points to make a polygon")
    return SimpleGeometricShapes( "polygon", x, y, [], [] )


def region( filename ):
    """region(filename)
    
    Load regions from file    
    """

    retval = None
    try:
        retval = dml.dmRegParse( filename )
        if not retval:
            raise IOError("try with a region()")
        return EnhancedRegion(retval)
    except:
        try:
            retval = dml.dmRegParse( "region({})".format(filename) )
            return EnhancedRegion(retval)
        except:
            raise IOError("Cannot parse region file")


def test():

    print circle(10,10,+123)+box(10,10,1,1)
    print circle(10,10,100)*box(10,10,1,1)
    print circle(10,10,100)-box(10,10,1,1)
    z = region("/lenin2.real/export/byobsid/repro/ds9.reg")
    print z
    #z.write("goo.reg")
    #z.plot()

    cc = region("/lenin2.real/Test/4.8b1/Regression/ciao4.8b1_linux64/dmcontour/04/dmcontour_4.fits")
    #print cc
    #cc.write("cntr.reg")
    
    a = pie(0,0, 1, 2, -45, 56)
    d = pie(0,0, 1, 2, -45, 56)
    b = a.__copy__()
    bb = box(5,5,10,20)
    print a
    print b
    print a == b
    print a == d
    print a == cc
    
    p = a+z+bb 
    q = p.tweak(dx=5).tweak(stretch=5).tweak(rotate=+45)
    print p
    print q

#__test__()


