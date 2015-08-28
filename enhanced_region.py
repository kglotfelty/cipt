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


__all__ = [ 'annulus', 'box', 'circle', 'ellipse', 'field', 'pie', 'point', 'polygon', 'region', 'rectangle', 'sector', 'dss' ]


_AND_ = "*"
_OR_ = "+"
_NOT_ = "!"
_BLANK_ = ""


from ctypes import *
try:
    region_lib = cdll.LoadLibrary('libregion.so')
    cxcdm_lib = cdll.LoadLibrary('libascdm.so')
except:
    try: # OSX
        region_lib = cdll.LoadLibrary('libregion.dylib')
        cxcdm_lib = cdll.LoadLibrary('libascdm.dylib')
    except:
        raise ImportError("Cannot load region library")


# Set return types.  Probably don't need void's, int's, and long's, but
# do it anyway for completeness
region_lib.regArea.restype = c_double
region_lib.regCreateEmptyRegion.restype = c_void_p
region_lib.regInsideRegion.restype = c_int
region_lib.regUnionRegion.restype = c_void_p
region_lib.regIntersectRegion.restype = c_void_p
region_lib.regInvert.restype = c_void_p
region_lib.regGetNoShapes.restype = c_long
region_lib.regGetShapeNo.restype = c_void_p
region_lib.regShapeRadii.restype = c_long
region_lib.regShapeAngles.restype = c_long
region_lib.regShapeGetName.restype = c_int
region_lib.regShapeGetNoPoints.restype = c_long
region_lib.regShapeGetComponent.restype = c_long
region_lib.regShapeGetPoints.restype = c_long # unused
region_lib.regShapeGetAngles.restype = c_long # unused
region_lib.regShapeGetRadii.restype = c_long # unused
cxcdm_lib.dmRegParse.restype = c_void_p




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

    _precise_string = "{:.12g}"  # Format string for double values when printed

    def __init__(self, shape_ptr ):
        """The input is a c_void_p regShape pointer.  Uses the region lib
        API to extract the shape parameters"""

        self._ptr = shape_ptr
        self._get_name()
        self._get_points()
        self._get_radii()
        self._get_angle()
        self._get_component()


    @property
    def xx(self):
        """Get the X values for the shape"""
        return self._xx
    
    @property
    def yy(self):
        """Get the Y values for the shape"""
        return self._yy
        
    @property
    def rad(self):
        """Get the radii of the shape"""
        return self._rad

    @property
    def ang(self):
        """Get the angle of the shape"""
        return self._ang

    @property
    def shape(self):
        """Get the shape name"""
        return self._shape

    @property
    def include(self):
        """Get the include flag"""
        return self._include


    def _get_name(self):
        """
        Get the name.  This also sets whether the region is inclusive
        or exclusive.
        """
        shape_name = create_string_buffer(100)
        iflag = region_lib.regShapeGetName(self._ptr, shape_name, 99)
        self._shape = shape_name.value.lower()
        self._include = _NOT_ if 0 == iflag else _BLANK_


    def _get_points(self):
        """
        Gets the x,y values.

        The string version of the values are stored now too.
        TODO: wrap w/ a setter that does the string conversion so that
        if class is exposed then string is kept in sync with double values
        """
        npts = region_lib.regShapeGetNoPoints( self._ptr )
        ox = (c_double * npts)()
        oy = (c_double * npts)()
        region_lib.regShapeGetPoints( self._ptr, byref(ox), byref(oy), c_long( npts ))
        self._xx = [x for x in ox]
        self._yy = [y for y in oy]


    def _get_radii(self):
        """
        Gets the radius/radii

        Also stores the string version (see TODO above)
        """
        nrad = region_lib.regShapeRadii( self._ptr )
        outr = (c_double * nrad)()
        region_lib.regShapeGetRadii( self._ptr, byref(outr) )
        self._rad = [r for r in outr]


    def _get_angle(self):
        """
        Gets the angle/angles

        Also stores the string version (see TODO above)
        """
        nang = region_lib.regShapeAngles( self._ptr )
        oa = (c_double*nang)()
        region_lib.regShapeGetAngles( self._ptr, byref(oa))
        self._ang = [ a for a in oa ]


    def _get_component(self):
        """
        Stores the component value.  The actual value is not used.  All
        the shapes with the component value are grouped together and AND'ed.
        Diff components are OR'ed.
        """
        self.component = region_lib.regShapeGetComponent( self._ptr )


    def _store_string(self, fmt_str):
        """
        Convert the values to strings using the desired precision/format
        """
        self.xx_str = [ self._precise_string.format(x) for x in self.xx ]
        self.yy_str = [ self._precise_string.format(y) for y in self.yy ]
        self.rad_str = [ self._precise_string.format(y) for y in self.rad ]
        self.ang_str = [ self._precise_string.format(y) for y in self.ang ]

        self.fmt_vals = {
            'x0' : self.xx_str[0] if len(self.xx_str) > 0 else None,
            'x1' : self.xx_str[1] if len(self.xx_str) > 1 else None,
            'y0' : self.yy_str[0] if len(self.yy_str) > 0 else None,
            'y1' : self.yy_str[1] if len(self.yy_str) > 1 else None,
            'r0' : self.rad_str[0] if len(self.rad_str) > 0 else None,
            'r1' : self.rad_str[1] if len(self.rad_str) > 1 else None,
            'a0' : self.ang_str[0] if len(self.ang_str) > 0 else None,
            'a1' : self.ang_str[1] if len(self.ang_str) > 1 else None,
            'i'  : self.include,
            'n'  : self.shape,
            'xy' : ",".join([ "{},{}".format(x,y) for x,y in zip( self.xx_str, self.yy_str ) ])
            }
            
        return fmt_str.format( **self.fmt_vals )


    def __str__(self):
        raise NotImplementedError("Implement this in the derived classes")

    def __repr__(self):
        return str(self)

    @with_plotting
    def plot(self, engine):
        raise NotImplementedError("Implement this in the derived classes")


class Annulus(EnhancedShape):
    """An annulus is defined by x_center, y_center, inner_radius, outer_radius"""
    def __str__( self ):
        return self._store_string("{i}{n}({x0},{y0},{r0},{r1})")


    @with_plotting
    def plot(self, engine):
        engine.plot_annulus( self.xx[0], self.yy[0], self.rad[0], self.rad[1] )


class Box(EnhancedShape):
    """A box is defined by x_center, y_center, x_length, y_length"""
    def __str__( self ):
        return self._store_string("{i}{n}({x0},{y0},{r0},{r1})")


    @with_plotting
    def plot(self,engine):
        engine.plot_box( self.xx[0], self.yy[0], self.rad[0], self.rad[1] )


class Circle(EnhancedShape):
    """A circle is defined by x_center, y_center, radius"""
    def __str__( self ):
        return self._store_string("{i}{n}({x0},{y0},{r0})")


    @with_plotting
    def plot(self,engine):
        engine.plot_circle( self.xx[0], self.yy[0], self.rad[0])


class Ellipse(EnhancedShape):
    """An ellipse is defined by x_center, y_center, major_axis, minor_axis, and angle"""
    def __str__( self ):
        return self._store_string("{i}{n}({x0},{y0},{r0},{r1},{a0})")
        
        
    @with_plotting
    def plot(self,engine):
        engine.plot_ellipse( self.xx[0], self.yy[0], self.rad[0], self.rad[1], self.ang[0])


class Field(EnhancedShape):
    """A field is defined to be the entire R^2 dataspace"""
    def __str__( self ):
        return self._store_string("{i}{n}()")


    @with_plotting
    def plot(self,engine):
        pass


class Pie(EnhancedShape):
    """A Pie is defined by x_center, y_center, inner_radius, outer_radius, start_angle, stop_angle"""
    def __str__( self ):
        return self._store_string("{i}{n}({x0},{y0},{r0},{r1},{a0},{a1})")
        

    @with_plotting
    def plot(self,engine):
        engine.plot_pie( self.xx[0], self.yy[0], self.rad[0], self.rad[1], self.ang[0], self.ang[1])


class Point(EnhancedShape):
    """A point is defined by x_center, y_center"""
    def __str__( self ):
        return self._store_string("{i}{n}({x0},{y0})")


    @with_plotting
    def plot(self,engine):
        pass


class Polygon(EnhancedShape):
    """A Polygon is defined as x1,y1, x2,y2,x3,y3,..."""
    def __str__( self ):
        return self._store_string("{i}{n}({xy})")


    @with_plotting
    def plot(self,engine):
        engine.plot_polygon(self.xx, self.yy)


class Rectangle(EnhancedShape):
    """A Polygon is defined as lower_left_x,lower_left_y,upper_right_x,upper_right_y"""
    def __str__( self ):
        return self._store_string("{i}{n}({x0},{y0},{x1},{y1})")


    @with_plotting
    def plot(self,engine):
        engine.plot_rectangle( self.xx[0], self.yy[0], self.xx[1], self.yy[1] )


class Rotbox(EnhancedShape):
    """A Rotbox is defined as x_center, y_center, x_length, y_length, angle """
    def __str__(self):
        return self._store_string("{i}{n}({x0},{y0},{r0},{r1},{a0})")


    @with_plotting
    def plot(self,engine):
        engine.plot_box( self.xx[0], self.yy[0], self.rad[0], self.rad[1], self.ang[0])


class Sector(EnhancedShape):
    """A Sector is defined by x_center, y_center, start_angle, stop_angle"""
    def __str__( self ):
        return self._store_string("{i}{n}({x0},{y0},{a0},{a1})")


    @with_plotting
    def plot(self,engine):
        engine.plot_pie( self.xx[0], self.yy[0], 0.0, 999999.9, self.ang[0], self.ang[1])



class EnhancedRegion( object ):
    """

    A region is made up of a collection of EnhancedShapes that are combined
    with AND and OR logical operators.

    Only the shape_classes shapes are implemented.
    """
    _shape_classes = { 'annulus': Annulus,
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
        self._load_shapes()
        self._region_lib = region_lib # keep a ref to this so it can be used to free shapes during garbage collection


    def __str__( self ):
        """Construct the string value from the shapes & logic"""
        retval = _BLANK_.join( [logic+str(shape) for logic,shape in zip(self.logic,self.shapes) ] )
        return(retval)


    def __repr__(self):
        return self.__str__()


    def _load_shapes(self):
        """
        Load the shapes into the region.
        """
        self._classify_shapes()
        self._determine_logic()


    @staticmethod
    def _get_shape_name( shape_ptr):
        """Get the name of a shape"""
        shape_name = create_string_buffer(100)
        region_lib.regShapeGetName(shape_ptr, shape_name, 99)
        return shape_name.value.lower()


    @property
    def shapes(self):
        """A list of the shapes in the region"""
        return self._shapes


    def _classify_shapes( self ):
        """Create the correct EnhancedShape derived class based on the shape's name."""
        nshapes = region_lib.regGetNoShapes( self._ptr )
        self._shapes = []
        for i in range(1,nshapes+1):
            shape_ptr = region_lib.regGetShapeNo( self._ptr, c_long(i) )
            shape_name = self._get_shape_name( shape_ptr)
            if shape_name in self._shape_classes:
                self._shapes.append( self._shape_classes[shape_name](shape_ptr) )
            else:
                raise ValueError("Unknown shape name {}".format(shape_name))


    @property
    def logic(self):
        """Return a list of logic operators, how are shapes combined"""
        return self._logic


    def _determine_logic(self):
        """Determine the logic used to combine shapes.

        The same component values are AND'ed together,  Different
        are OR'd.  The standard convention is that the same component values
        are stored together, so we only need to check for when the value
        changes."""

        cpt_vals = [s.component for s in self.shapes ]

        self._logic = [ _BLANK_ ] # First shape has no logic
        for i in range( 1, len(cpt_vals)):
            # If component values are equal then &, else |
            if cpt_vals[i]==cpt_vals[i-1]:
                self._logic.append(_AND_)
            else:
                self._logic.append(_OR_)


    def area( self ):
        """Determine the region area.

        TODO: make can pass in field bound box and bin-size parameters.
        There is also a routine that forces pixelated area we could
        optionally call
        """
        import sys as sys
        DBL_MAX = sys.float_info.max
        fld = wrap_vals( [-DBL_MAX, DBL_MAX] )
        return region_lib.regArea( self._ptr, fld, fld, c_double(1.0) )


    def write(self, filename, newline=False, fits=False):
        """
        Write the region to a file.  Filename is always clobbered!

        """
        import os as os
        if os.path.exists( filename ):
            os.unlink(filename)
        
        if fits:
            # Some error checking here  would be nice :)
            oDs = cxcdm_lib.dmDatasetCreate( filename )
            blk = cxcdm_lib.dmTableWriteRegion( oDs, "REGION", None, self._ptr)
            cxcdm_lib.dmDatasetClose(oDs)
            
        else:
            as_str = str(self)
            if newline:  # Replace OR (+) with new lines.  Have to keep *'s
                as_str = as_str.replace(_OR_, "\n")
                
            with open(filename, "w") as fp:
                fp.write( as_str )
                fp.write("\n")


    def inside( self, x, y ):
        """
        Determine if x,y pair is inside the region
        """
        return (1 == region_lib.regInsideRegion(self._ptr, c_double(x), c_double(y)))


    def plot( self ):
        """
        Plot the region by plotting each shape, if plotting is
        available.
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
        self._region_lib.regFree( self._ptr )


    def __add__(self,other):
        """
        Logically OR (Union) two regions together
        """
        if not isinstance( other, type(self)): 
            raise TypeError("Cannot perform region logic with {} type".format(type(other)) )

        cpts = self._ptr
        cpto = other._ptr
        return EnhancedRegion( region_lib.regUnionRegion( cpts, cpto ) )


    def __mul__(self,other):
        """
        Logically AND (Intersect) two regions together
        """
        if not isinstance( other, type(self)):
            raise TypeError("Cannot perform region logic with {} type".format(type(other)) )

        cpts = self._ptr
        cpto = other._ptr
        return EnhancedRegion( region_lib.regIntersectRegion( cpts, cpto ) )


    def __sub__(self,other):
        """
        Subtraction is actually a short-hand or AND NOT.
        Other is inverted and then AND'ed with self.
        """
        if not isinstance( other, type(self)):
            raise TypeError("Cannot perform region logic with {} type".format(type(other)) )

        cpts = self._ptr
        cpto = other._ptr
        invrt = region_lib.regInvert( cpto )
        return EnhancedRegion( region_lib.regIntersectRegion( cpts, invrt ) )


    def __neg__(self):
        """
        Invert region
        """
        return EnhancedRegion( region_lib.regInvert( self._ptr ))


    def __copy__(self):
        """Return a copy of the region"""
        return EnhancedRegion( region_lib.regCopyRegion( self._ptr))


    def __eq__(self, other):
        """Is this region equal to another (shape-by-shape comparison)"""
        if not isinstance( other, type(self)):
            raise TypeError("Cannot perform region logic with {} type".format(type(other)) )

        return region_lib.regCompareRegion(self._ptr, other._ptr)

    def __and__(self, other ):
        """shape & shape ==> shape * shape """
        return self*other

    def __or__(self, other):
        """ shape | shape => shape + shape """
        return self+other

    def __xor__(self,other):
        """ shape ^ shape ==> (a-b)+(b-a)"""
        return (self-other)+(other-self)
        
    def __invert__(self):
        """~shape == !shape"""
        return -self
    

    def __getitem__(self, idx ):
        """Returns the logic and shape as index"""
        try:
            shape = self.shapes[idx]
            logic = self.logic[idx]
        except IndexError, e:
            raise IndexError("Invalid index value for this region")
        except:
            raise

        copy_ptr = region_lib.regCreateEmptyRegion()

        reg_inc = c_int(0) if _NOT_ == shape.include else c_int(1)
        region_lib.regAppendShape( copy_ptr,
                            shape.shape,
                            reg_inc, c_int(0),
                            wrap_vals( shape.xx ),
                            wrap_vals( shape.yy ),
                            c_long( len(shape.xx) ),
                            wrap_vals( shape.rad ),
                            wrap_vals( shape.ang ),
                            c_int(0), c_int(0) )
        return ( logic, EnhancedRegion(copy_ptr) )
        

    def __len__(self):
        """Returns the number of shapes in the region"""
        return len(self.shapes)


    def __contains__(self, other):
        """Checks to see if a shape is contained within a region
        
        >>> circle(1,1,2) in box(10,1,3,3)+circle(1,1,2)*sector(1,1,34,44)
        True
        """

        if not isinstance( other, type(self)):
            raise TypeError("Cannot perform region logic with {} type".format(type(other)) )
    
        if len(other) != 1:
            raise IndexError("Other region must contain a single shape")
        
        oth_shape_ptr = other.shapes[0]._ptr

        for ss in self.shapes:
            slf_shape_ptr = ss._ptr
            
            if region_lib.regCompareShape( slf_shape_ptr, oth_shape_ptr ):
                return True
        
        return False
        

    def index(self, other):
        """Checks to see if a shape is contained within a region, returns
        a list of indexes
        
        >>> b = circle(1,1,2) 
        >>> a = box(10,1,3,3)+circle(1,1,2)*sector(1,1,34,44)
        >>> a.index(b)
        [1]
        """

        if not isinstance( other, type(self)):
            raise TypeError("Cannot perform region logic with {} type".format(type(other)) )
    
        if len(other) != 1:
            raise IndexError("Other region must contain a single shape")
        
        oth_shape_ptr = other.shapes[0]._ptr
        retval = []

        for ii in range(len(self)):
            slf_shape_ptr = self.shapes[ii]._ptr
            
            if region_lib.regCompareShape( slf_shape_ptr, oth_shape_ptr ):
                retval.append(ii)
        
        return retval
        

    def tweak( self, dx=0, dy=0, stretch=1, rotate=0 ):
        """
        Create a copy of the current region, shape-by-shape.

        Then apply tweaks to the x, y, r, angle parameters before
        creating and returning a new region.

        Each tweak is applied to each shape separately
        """
        copy_ptr = region_lib.regCreateEmptyRegion()
        for ii in range( len(self) ):
            logic,s = self[ii]
            shape = self.shapes[ii]

            reg_math = c_int(0) if _AND_ == logic else c_int(1)
            reg_inc = c_int(0) if _NOT_ == shape.include else c_int(1)

            copy_xx = [ x+dx for x in shape.xx]
            copy_yy = [ y+dy for y in shape.yy]
            copy_rr = [ r*stretch for r in shape.rad]
            copy_aa = [ a+rotate for a in shape.ang ]

            region_lib.regAppendShape( copy_ptr,
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
    new_ptr = region_lib.regCreateEmptyRegion()
    region_lib.regAppendShape( new_ptr,
                        name_as_string ,
                        c_int(1), c_int(1),
                        wrap_vals(xx),
                        wrap_vals(yy),
                        c_long( len(xx) ),
                        wrap_vals( radius ),
                        wrap_vals( angle ),
                        c_int(0), c_int(0) )
    if 0 == region_lib.regGetNoShapes( new_ptr ):
        raise RuntimeError("Bad region parameters")

    return EnhancedRegion(new_ptr)

# ------------------
# User Interface

def circle( x, y, r ):
    """circle( x_center, y_center, radius )"""
    return SimpleGeometricShapes( "circle", [x], [y], [r], [] )


def ellipse( x, y, major, minor, angle=0 ):
    """ellipse( x_center, y_center, semi_major_axis, semi_minor_axis[, rotation_angle])

    angle is measured in degrees from the +X axis
    """
    return SimpleGeometricShapes( "ellipse", [x], [y], [major, minor], [angle] )


def box( x,y, xlen, ylen, angle=0 ):
    """box( x_center, y_center, x_length, y_length[, rotation_angle])

    angle is measured in degrees from the +X axis
    """
    return SimpleGeometricShapes( "rotbox", [x], [y], [xlen, ylen], [angle] )


def rotbox( x,y, xlen, ylen, angle=0 ):
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
    or
    polygon( x_array, y_array ) 
    
    
    Must contain at least 3 points, and same number of x-y pairs.
    An extra point is added if last point is not equal to 1st point.

    """

    if len(args) == 2:
        x = args[0]
        y = args[1]
    else:
        x = [args[i] for i in range(0,len(args),2) ]
        y = [args[i] for i in range(1,len(args),2) ]

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
        retval = cxcdm_lib.dmRegParse( str(filename) )
        if not retval:
            raise IOError("try with a region()")
        return EnhancedRegion(retval)
    except Exception, E:
        try:
            retval = cxcdm_lib.dmRegParse( "region({})".format(str(filename)) )
            if not retval:
                raise
            return EnhancedRegion(retval)
        except:
            raise IOError("Cannot parse region file")



cxcdm_lib.dmBlockOpen.restype = c_void_p
cxcdm_lib.dmSubspaceColOpen.restype = c_void_p
cxcdm_lib.dmSubspaceColGetRegion.restype = c_void_p
cxcdm_lib.dmBlockGetDataset.restype = c_void_p

def dss( filename, colname="sky", component=1 ):
    """Okay, this is treading on thin ice ... but I know it gives
    me a regRegion *ptr """
    
    retval = None
    
    blk = cxcdm_lib.dmBlockOpen( None, filename )
    if 0 == blk:
        raise IOError("Unable to open file '{}'".format(filename))

    cxcdm_lib.dmBlockSetSubspaceCpt( blk, c_long( component ) )
    col = cxcdm_lib.dmSubspaceColOpen( c_void_p( blk), colname )
    ptr = cxcdm_lib.dmSubspaceColGetRegion( col )     

    retval = EnhancedRegion( ptr )
    ds = cxcdm_lib.dmBlockGetDataset( blk )
    cxcdm_lib.dmBlockClose( blk )
    cxcdm_lib.dmDatasetClose( ds )
    return retval
    

#
#
#
#
#

"""

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
    

    
def test():
    print annulus(1,2,3,4)
    print box(1,2,3,4)
    print box(1,2,3,4,5)
    print circle(1,2,3)
    print ellipse(1,2,3,4,5)
    print ellipse(1,2,3,4)
    print field()
    print pie(1,2,3,4,5,6)
    print point(1,2)
    print polygon( 1,2,3,4,5,6)
    print polygon( [1,2,3],[4,5,6])
    print rectangle(1,2,3,4)
    print rotbox(1,2,3,4,5)
    print sector(1,2,3,4)

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



