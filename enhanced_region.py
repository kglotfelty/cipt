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

>>> a in c
True
>>> b in c
True
>>> c.index(a)
[0]
>>> c.index(b)
[1]

>>> c[0]
circle(1000,1000,50)
>>> c[1]
+box(1010,1010,50,100)


Users can also get access to the individual shape properties
(radius, centers, etc)

>>> myreg = circle(10,-10,100)
>>> mycir = myreg.shapes[0]
>>> mycir.xx
10
>>> mycir.yy
-10
>>> mycir.rad
100
>>> mycir.include
ShapeInclusion(val=1, str='')


"""


__all__ = [ 'annulus', 'box', 'circle', 'ellipse', 'field', 'pie', 'point', 'polygon', 'region', 'rectangle', 'sector', 'dss' ]


_AND_ = "*"
_OR_ = "+"
_NOT_ = "!"
_BLANK_ = ""


import os
lib = os.environ["ASCDS_INSTALL"]+"/lib/"

from ctypes import *
try:
    region_lib = cdll.LoadLibrary(lib+'libregion.so')
    cxcdm_lib = cdll.LoadLibrary(lib+'libascdm.so')
except:
    try: # OSX
        region_lib = cdll.LoadLibrary(lib+'libregion.dylib')
        cxcdm_lib = cdll.LoadLibrary(lib+'libascdm.dylib')
    except:
        raise ImportError("Cannot load region library")


# Set return types.  OSX is picky about these

region_lib.regArea.restype = c_double
region_lib.regCreateEmptyRegion.restype = c_void_p
region_lib.regCompareRegion.restype = c_int
region_lib.regCopyRegion.restype = c_void_p
region_lib.regExtent.restype = c_int
region_lib.regGetNoShapes.restype = c_long
region_lib.regGetShapeNo.restype = c_void_p
region_lib.regInsideRegion.restype = c_int
region_lib.regIntersectRegion.restype = c_void_p
region_lib.regInvert.restype = c_void_p
#region_lib.regShapeAngles.restype = c_long
#region_lib.regShapeRadii.restype = c_long
region_lib.regShapeGetAngles.restype = c_long 
region_lib.regShapeGetComponent.restype = c_long
region_lib.regShapeGetName.restype = c_int
region_lib.regShapeGetNoPoints.restype = c_long
region_lib.regShapeGetPoints.restype = c_long 
region_lib.regShapeGetRadii.restype = c_long 
region_lib.regUnionRegion.restype = c_void_p

cxcdm_lib.dmDatasetCreate.restype = c_void_p
cxcdm_lib.dmBlockGetDataset.restype = c_void_p
cxcdm_lib.dmBlockOpen.restype = c_void_p
cxcdm_lib.dmRegParse.restype = c_void_p
cxcdm_lib.dmSubspaceColGetRegion.restype = c_void_p
cxcdm_lib.dmSubspaceColOpen.restype = c_void_p
cxcdm_lib.dmTableWriteRegion.restype = c_void_p


def wrap_vals( vv ):
    """Utility to wrap arrays to be pased in byreference"""
    if vv:
        retval = byref(( c_double * len(vv) ) (*vv))
    else:
        retval = POINTER(c_int)()  # NULL pointer
    return retval


from collections import namedtuple
# These provide the enumerated values in the region lib private header file.

# The 'val' is the value passed to the region_lib, the 'str' is the
# value that is printed
ShapeOperation = namedtuple("ShapeOperation", ['val', 'str'] )
opAND = ShapeOperation( 0, _AND_ )
opOR = ShapeOperation( 1, _OR_ )
opNOOP = ShapeOperation( 1, _BLANK_ )

ShapeInclusion = namedtuple("ShapeInclusion", ['val', 'str'] )
incINCLUDE = ShapeInclusion( 1, _BLANK_ )
incEXCLUDE = ShapeInclusion( 0, _NOT_ )



class EnhancedShape( object ):
    """
    All regions are composed of simple geometric shapes (even 1 region with only
    1 shape is still a EnhancedRegion() ).

    Users cannot manipulate an individual shape's properties.  They can 'tweak' 
    the properties and create a new region.

    """

    _precise_string = "{:.12g}"  # Format string for double values when printed

    def __init__(self, shape_ptr ):
        """The input is a c_void_p regShape pointer.  Uses the region lib
        API to extract the shape parameters"""

        self._ptr = c_void_p(shape_ptr)
        self._get_name()
        self._get_points()
        self._get_radii()
        self._get_angle()
        self._get_component()
        self._logic = None


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


    @property
    def logic(self):
        """Get the logic operation"""
        return self._logic


    def _get_name(self):
        """
        Get the name.  This also sets whether the region is inclusive
        or exclusive.
        """
        shape_name = create_string_buffer(100)
        iflag = region_lib.regShapeGetName(self._ptr, shape_name, 99)
        self._shape = shape_name.value.lower()
        self._include = incEXCLUDE if incEXCLUDE.val == iflag else incINCLUDE


    def _get_points(self):
        """
        Gets the x,y values from the regShape pointer

        """
        npts = region_lib.regShapeGetNoPoints( self._ptr )
        ox = (c_double * npts)()
        oy = (c_double * npts)()
        region_lib.regShapeGetPoints( self._ptr, byref(ox), byref(oy), c_long( npts ))
        self._xx = tuple([x for x in ox])  # tuples are immutable 
        self._yy = tuple([y for y in oy])


    def _get_radii(self):
        """
        Gets the radius/radii from the regShape pointer
        """
        #nrad = region_lib.regShapeRadii( self._ptr )
        outr = (c_double * self._nrad)()
        region_lib.regShapeGetRadii( self._ptr, byref(outr) )
        self._rad = tuple([r for r in outr])


    def _get_angle(self):
        """
        Gets the angle/angles from the regShape pointer
        """
        #nang = region_lib.regShapeAngles( self._ptr )
        oa = (c_double*self._nang)()
        region_lib.regShapeGetAngles( self._ptr, byref(oa))
        self._ang = tuple([ a for a in oa ])


    def _get_component(self):
        """
        Stores the component value.  The actual value is not used.  All
        the shapes with the component value are grouped together and AND'ed.
        Diff components are OR'ed.
        """
        self.component = region_lib.regShapeGetComponent( self._ptr )


    def _format_string(self, fmt_str):
        """
        Convert the values to strings using the desired precision/format
        """
        self.xx_str = [ self._precise_string.format(x) for x in self.xx ]
        self.yy_str = [ self._precise_string.format(y) for y in self.yy ]
        self.rad_str = [ self._precise_string.format(r) for r in self.rad ]
        self.ang_str = [ self._precise_string.format(a) for a in self.ang ]

        self.fmt_vals = {
            'x0' : self.xx_str[0] if len(self.xx_str) > 0 else None,
            'x1' : self.xx_str[1] if len(self.xx_str) > 1 else None,
            'y0' : self.yy_str[0] if len(self.yy_str) > 0 else None,
            'y1' : self.yy_str[1] if len(self.yy_str) > 1 else None,
            'r0' : self.rad_str[0] if len(self.rad_str) > 0 else None,
            'r1' : self.rad_str[1] if len(self.rad_str) > 1 else None,
            'a0' : self.ang_str[0] if len(self.ang_str) > 0 else None,
            'a1' : self.ang_str[1] if len(self.ang_str) > 1 else None,
            'i'  : self.include.str,
            'n'  : self.shape,
            'xy' : ",".join([ "{},{}".format(x,y) for x,y in zip( self.xx_str, self.yy_str ) ])
            }
            
        return fmt_str.format( **self.fmt_vals )


    def _set_logic(self, val ):
        self._logic = val


    def __str__(self):
        raise NotImplementedError("Implement this in the derived classes")


    def __repr__(self):
        return str(self)


class Annulus(EnhancedShape):
    """An annulus is defined by x_center, y_center, inner_radius, outer_radius"""

    _nang = 0
    _nrad = 2

    def __str__( self ):
        return self._format_string("{i}{n}({x0},{y0},{r0},{r1})")


class Box(EnhancedShape):
    """A box is defined by x_center, y_center, x_length, y_length"""

    _nang = 0
    _nrad = 2
    
    def __str__( self ):
        return self._format_string("{i}{n}({x0},{y0},{r0},{r1})")


class Circle(EnhancedShape):
    """A circle is defined by x_center, y_center, radius"""

    _nang = 0
    _nrad = 1

    def __str__( self ):
        return self._format_string("{i}{n}({x0},{y0},{r0})")


class Ellipse(EnhancedShape):
    """An ellipse is defined by x_center, y_center, major_axis, minor_axis, and angle"""

    _nang = 1
    _nrad = 2

    def __str__( self ):
        return self._format_string("{i}{n}({x0},{y0},{r0},{r1},{a0})")
        

class Field(EnhancedShape):
    """A field is defined to be the entire R^2 dataspace"""

    _nang = 0
    _nrad = 0

    def __str__( self ):
        return self._format_string("{i}{n}()")


class Pie(EnhancedShape):
    """A Pie is defined by x_center, y_center, inner_radius, outer_radius, start_angle, stop_angle"""

    _nang = 2
    _nrad = 2

    def __str__( self ):
        return self._format_string("{i}{n}({x0},{y0},{r0},{r1},{a0},{a1})")
        

class Point(EnhancedShape):
    """A point is defined by x_center, y_center"""

    _nang = 0
    _nrad = 0
    
    def __str__( self ):
        return self._format_string("{i}{n}({x0},{y0})")


class Polygon(EnhancedShape):
    """A Polygon is defined as x1,y1, x2,y2,x3,y3,..."""

    _nang = 0
    _nrad = 0

    def __str__( self ):
        return self._format_string("{i}{n}({xy})")


class Rectangle(EnhancedShape):
    """A Polygon is defined as lower_left_x,lower_left_y,upper_right_x,upper_right_y"""

    _nang = 0
    _nrad = 0

    def __str__( self ):
        return self._format_string("{i}{n}({x0},{y0},{x1},{y1})")


class Rotbox(EnhancedShape):
    """A Rotbox is defined as x_center, y_center, x_length, y_length, angle """

    _nang = 1
    _nrad = 2

    def __str__(self):
        return self._format_string("{i}{n}({x0},{y0},{r0},{r1},{a0})")


class Sector(EnhancedShape):
    """A Sector is defined by x_center, y_center, start_angle, stop_angle"""

    _nang = 2
    _nrad = 0

    def __str__( self ):
        return self._format_string("{i}{n}({x0},{y0},{a0},{a1})")

# --------------------------------------

class EnhancedRegion( object ):
    """

    A region is made up of a collection of EnhancedShapes that are combined
    with AND and OR logical operators.

    """

    # This provides a mapping between shape name and the class of the object to create
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
        self._ptr = c_void_p(ptr)
        self._load_shapes()
        self._region_lib = region_lib # keep a ref to this so it can be used to free shapes during garbage collection


    def __str__( self ):
        """Construct the string value from the shapes & logic"""
        retval = _BLANK_.join( [shape.logic.str+str(shape) for shape in self.shapes ] )
        return(retval)


    def __repr__(self):
        return self.__str__()


    def _load_shapes(self):
        """
        Load the shapes into the region.
        """
        self._classify_shapes()
        logic = self._determine_logic()
        if self._shapes:
            map( lambda s,l: s._set_logic(l), self._shapes, logic)


    @staticmethod
    def _get_shape_name( shape_ptr):
        """Get the name of a shape.
        
        This is kind of a look-ahead thing so that we can know the shape
        name so that we know which EhanceShape sub-class to create.
        """
        shape_name = create_string_buffer(100)
        region_lib.regShapeGetName(shape_ptr, shape_name, c_long(99))
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
            vptr = c_void_p(shape_ptr)
            shape_name = self._get_shape_name(vptr)
            if shape_name in self._shape_classes:
                # Create a new shape object of the correct type
                new_shape_class = self._shape_classes[shape_name]
                new_shape = new_shape_class(shape_ptr)
                self._shapes.append( new_shape )
            else:
                raise ValueError("Unknown shape name {}".format(shape_name))
        self._shapes = tuple(self._shapes)


    def _determine_logic(self):
        """Determine the logic used to combine shapes.

        The same component values are AND'ed together,  Different
        are OR'd.  The standard convention is that the same component values
        are stored together, so we only need to check for when the value
        changes."""

        cpt_vals = [s.component for s in self.shapes ]

        logic = [ opNOOP ] # First shape has no logic
        for i in range( 1, len(cpt_vals)):
            # If component values are equal then &, else |
            if cpt_vals[i]==cpt_vals[i-1]:
                logic.append( opAND )
            else:
                logic.append( opOR )
        return logic


    def area( self, res=1.0 ):
        """Determine the region area.

        TODO: make can pass in field bound box and bin-size parameters.
        There is also a routine that forces pixelated area we could
        optionally call
        """
        ###import sys as sys
        ###DBL_MAX = sys.float_info.max
        ###fld = wrap_vals( [-DBL_MAX, DBL_MAX] )

        xtnd = self.extent()
        fldx = wrap_vals( [ xtnd["xlo"], xtnd["xhi"] ] )
        fldy = wrap_vals( [ xtnd["ylo"], xtnd["yhi"] ] )

        return region_lib.regArea( self._ptr, fldx, fldy, c_double(res) )


    def extent( self ):
        """
        Determine the size of a bounding box around the region.
        
        Returns 2 list:  x limits and y limits
        """
        import sys as sys
        DBL_MAX = sys.float_info.max
        fld = wrap_vals( [-DBL_MAX, DBL_MAX] )
        x_range = (c_double * 2)()
        y_range = (c_double * 2)()
        region_lib.regExtent( self._ptr, fld, fld, byref(x_range), byref( y_range) )

        retval = { 'xlo' : x_range[0], 
                   'xhi' : x_range[1],
                   'ylo' : y_range[0],
                   'yhi' : y_range[1] }

        return retval


    def write(self, filename, newline=False, fits=False):
        """
        Write the region to a file.  Filename is always clobbered!

        """
        import os as os
        if os.path.exists( filename ):
            os.unlink(filename)
        
        if fits:
            # Some error checking here  would be nice :)
            oDs = cxcdm_lib.dmDatasetCreate( c_char_p(filename) )
            if not oDs:
                raise IOError("Cannot create filename {}".format(filename))
            blk = cxcdm_lib.dmTableWriteRegion( c_void_p(oDs), c_char_p("REGION"), None, self._ptr)
            cxcdm_lib.dmDatasetClose(c_void_p(oDs))
            
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


    def __del__(self):
        """
        Use the self.dll which should exist during garbage collection.
        """
        if hasattr( self, "_region_lib" ):
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
        return EnhancedRegion( region_lib.regIntersectRegion( cpts, c_void_p(invrt) ))


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

        retval = region_lib.regCompareRegion(self._ptr, other._ptr)

        return (retval != 0)

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
        """Returns a new region with the shape at index"""
        try:
            shape = self.shapes[idx]
        except IndexError, e:
            raise IndexError("Invalid index value for this region")
        except:
            raise

        copy_ptr = region_lib.regCreateEmptyRegion()
        vptr = c_void_p(copy_ptr)


        reg_inc = shape.include.val
        region_lib.regAppendShape( vptr,
                            c_char_p(shape.shape),
                            c_int(reg_inc), c_int(opNOOP.val),
                            wrap_vals( shape.xx ),
                            wrap_vals( shape.yy ),
                            c_long( len(shape.xx) ),
                            wrap_vals( shape.rad ),
                            wrap_vals( shape.ang ),
                            None, None )

        retval = EnhancedRegion(copy_ptr)
        retval.shapes[0]._set_logic( shape.logic )
        return retval
        

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
        

    def tweak( self, dx=0, dy=0, stretch=1, pad=0, rotate=0 ):
        """
        Create a copy of the current region, shape-by-shape and 
        apply tweaks to the x, y, r, angle parameters before
        creating and returning a new region.

        Each tweak is applied to each shape separately

        dx : each x coordinate is shifted by dx amount
        dy : each y coordinate is shifted by dy amount
        stretch: each radii is multiplied by the stretch amount
        pad: the pad amount is added to each radii
        rotate: the rotation angle 
        
        Note: rotate only applies to shapes with a rotation angle
        parameter.  For example, a rectangle, unlike a box, does 
        not have a rotation parameter.  
        
        Both strecth and pad are applied to the radii; first the
        radii are stretched, then the pad is applied. 

        When a 'box' is rotated, it will automatically be converted
        into a 'rotbox'

        Examples:
        
        >>> circle(0,0,10).tweak(dx=5,dy=-7,stretch=2,pad=3,rotate=45)
        circle(5,-7,23)
        
        >>> ellipse(0,0,5,10,0).tweak(rotate=90)
        ellipse(0,0,5,10,90)
        
        >>> box(0,0,5,5).tweak(rotate=45)
        rotbox(0,0,5,5,45)

        >>> polygon( 0,0, 0,10, 10,0 ).tweak(dx=5)
        polygon(5,0,5,10,15,0,5,0)
        """
        copy_ptr = region_lib.regCreateEmptyRegion()
        vptr = c_void_p(copy_ptr)
        for ii in range( len(self) ):
            shape = self.shapes[ii]

            reg_math =shape.logic.val
            reg_inc = shape.include.val

            copy_xx = [ x+dx for x in shape.xx]
            copy_yy = [ y+dy for y in shape.yy]
            copy_rr = [ r*stretch+pad for r in shape.rad]
            copy_aa = [ a+rotate for a in shape.ang ]

            if 0 != rotate and "box" == shape.shape.lower():
                use_shape = "rotbox"
                copy_aa = [ rotate ]
            else:
                use_shape = shape.shape

            region_lib.regAppendShape( vptr,
                                c_char_p(use_shape),
                                c_int(reg_inc), c_int(reg_math),
                                wrap_vals( copy_xx ),
                                wrap_vals( copy_yy ),
                                c_long( len(copy_xx) ),
                                wrap_vals( copy_rr ),
                                wrap_vals( copy_aa ),
                                None, None)
            if region_lib.regGetNoShapes( vptr ) != ii+1:
                raise ValueError("Problem creating region")

        return EnhancedRegion(copy_ptr)


    def xform( self, xfunc ):
        """
        Create a copy of the current region, shape-by-shape and 
        apply the transformation function to the parameters.
        
        The x,y are transformed by applying the xfunc
        directly to the values
        
        Radii and angles are transformed by computing the scaling function
        for a slightly shifted location and finding the shift and
        rotation in the transformed coordinates.  This isn't
        technically correct but without wanting to pass around 
        the scale separately this is pretty close.
        The error here is probably less than the assumption that a 
        circle in tangent plane is also a circle in ra/dec.
        
        """
        from math import atan2, hypot, degrees
        _dy = 1/3600.0

        copy_ptr = region_lib.regCreateEmptyRegion()
        vptr = c_void_p(copy_ptr)
        for ii in range( len(self) ):
            shape = self.shapes[ii]

            reg_math =shape.logic.val
            reg_inc = shape.include.val

            copy_xx = [ x for x in shape.xx]
            copy_yy = [ y for y in shape.yy]

            copy_xy = zip(copy_xx, copy_yy)
            
            # Apply the function to the x,y values directly
            xfrm_xy = xfunc(copy_xy)
            xfrm_xx, xfrm_yy = zip(*xfrm_xy)

            # Get the scale and the angle
            copy_yy = [ y+_dy for x in shape.yy]
            copy_xy = zip(copy_xx, copy_yy)
            
            xfrm_xy = xfunc(copy_xy)
            delt_xx, delt_yy = zip(*xfrm_xy)
            stretch = hypot((delt_xx[0]-xfrm_xx[0]),(delt_yy[0]-xfrm_yy[0]))/_dy
            xfrm_rr = [ r*stretch for r in shape.rad]

            rotate = degrees(atan2((delt_xx[0]-xfrm_xx[0]),(delt_yy[0]-xfrm_yy[0])))
            xfrm_aa = [ a+rotate for a in shape.ang ]

            if 0 != rotate and "box" == shape.shape.lower():
                use_shape = "rotbox"
                copy_aa = [ rotate ]
            else:
                use_shape = shape.shape

            region_lib.regAppendShape( vptr,
                                c_char_p(use_shape),
                                c_int(reg_inc), c_int(reg_math),
                                wrap_vals( xfrm_xx ),
                                wrap_vals( xfrm_yy ),
                                c_long( len(xfrm_xx) ),
                                wrap_vals( xfrm_rr ),
                                wrap_vals( xfrm_aa ),
                                None, None)
            if region_lib.regGetNoShapes( vptr ) != ii+1:
                raise ValueError("Problem creating region")

        return EnhancedRegion(copy_ptr)



def SimpleGeometricShapes( name_as_string, xx, yy, radius, angle ):
    """
    For single shape regions crated with the routines below
    """
    new_ptr = region_lib.regCreateEmptyRegion()
    vptr = c_void_p(new_ptr)
    region_lib.regAppendShape( vptr,
                        c_char_p(name_as_string) ,
                        c_int(incINCLUDE.val), 
                        c_int(opNOOP.val),
                        wrap_vals(xx),
                        wrap_vals(yy),
                        c_long( len(xx) ),
                        wrap_vals( radius ),
                        wrap_vals( angle ),
                        None, None)
    if 0 == region_lib.regGetNoShapes( vptr ):
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

    Load regions from file or parse a region string
    
    >>> reg = region("foo.fits")
    >>> reg = region("region(foo.fits)"
    >>> reg = region("bounds(region(foo.fits))")
    >>> reg = region("circle(1,2,4)")    
    """

    retval = None
    try:
        retval = cxcdm_lib.dmRegParse( c_char_p(str(filename) ))
        if not retval:
            raise IOError("try with a region()")
        return EnhancedRegion(retval)
    except IOError, E:
        try:
            retval = cxcdm_lib.dmRegParse( c_char_p("region({})".format(str(filename)) ))
            if not retval:
                raise
            return EnhancedRegion(retval)
        except:
            raise IOError("Cannot parse region file")


def dss( filename, colname="sky", component=1 ):
    """Extract the region from a DM files's data subspace.  
    The colname is the subspace column name (usually 'sky').
    The component value is the DM subspace component (1 to N); usually
    for different CCDs -- how component numbers map to 
    meaninful sets of data is beyond this scope.
    
    Examples:

    >>> reg = dss("foo.fits")
    >>> chip = dss("imap.fits", colname="chip")
    >>> ccd0 = dss("image.fits", component=5)
    """
    
    retval = None
    
    blk = cxcdm_lib.dmBlockOpen( None, c_char_p(filename) )
    if 0 == blk:
        raise IOError("Unable to open file '{}'".format(filename))

    cxcdm_lib.dmBlockSetSubspaceCpt( c_void_p(blk), c_long( component ) )
    col = cxcdm_lib.dmSubspaceColOpen( c_void_p( blk), c_char_p(colname) )
    ptr = cxcdm_lib.dmSubspaceColGetRegion( c_void_p(col) )     

    retval = EnhancedRegion( ptr )
    ds = cxcdm_lib.dmBlockGetDataset( c_void_p(blk) )
    cxcdm_lib.dmBlockClose( c_void_p(blk ))
    cxcdm_lib.dmDatasetClose(c_void_p(ds) )
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
    print -circle(50,50,100)
    z = region("/data/lenin2/export/byobsid/repro/ds9.reg")
    print z
    z.write("goo.reg")

    cc = region("/lenin2.real/Projects/ImproveRegression/Test/ciaox_20160125/dmcontour/04/dmcontour_4.fits")
    print len(cc)
    cc.write("cntr.reg")

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
    q = p.tweak(dx=5).tweak(stretch=5).tweak(rotate=+45).tweak(pad=3)
    print p
    print q

    q = dss('/data/lenin2/Projects/PIMMS/Doc/9768/repro/todetect/goo.fits')
    print len(q)
    print q[0]
    print q.index(q[0])
    print q[0] in q


    from pycrates import read_file
    c = circle(4274.5,3954.5,5)
    img = read_file("img.fits")
    wcs = img.get_transform("eqpos")
    print c
    print c.xform(wcs.apply)
    print c.xform(wcs.apply).xform(wcs.invert)
    

    


__all__.append("test")



