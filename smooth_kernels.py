


class CIAOImageKernel():
    """
    Base class for aconvolve convolution kernels 
    
    The aconvolve CIAO tool is used to convolve an image with a specified
    smoothing kernel.  The syntax it expects takes the form:
    
        lib:function(D,M,P,L1,L2)
        
    where 
        'lib:' is the token to signal that one of the built in
        convolution kernels should be used.  It also recognizes 
        'txt:' for a string serielization of the kernel values, and
        'file:' to specify an external file.  The later two are not
        implemented here.
    
        'function' is the shape of the convolution kernel.  Available
        functions include:

          Gaussian : gaus    - 2D non-rotated, Gaussian
          Tophat   : tophat  - a flat, elliptical shape 
          Boxcar   : box     - a flat box, or step, shape
          Mexhat   : mexhat  - a 2D mexican hat wavelet
          Exp      : exp     - a 2D exponential decay
          Power    : power   - a 2D powerlaw decay
          Beta      : beta    - a 2D lorentizan
          Sinc      : sinc    - sin(x)/x
          Cone      : cone    - a 2D cone shape
          Pyramid   : pyramid - a 2D pyramid
          Sphere    : sphere  - a hemisphere 1-(x^2) shape
          Hanning   : hanning - a Hann filter (from family of cosine filters)
          RaisedCosine : racos   - Generic raised cosine
    
        'D' is the number of dimensions, always 2 for this object.
    
        'M' is magnitude, or max value.  This value is ignored when
        the kernel normalization is not "none" (more below).  This 
        value is not used for the sphere function.
        
        'P' is an optional parameter.  For the quasi-infinite shapes
        (gaus,exp,power,beta,sinc), P scales the lenghts to determine
        the extent of the kernel.  For racos it is the parameter
        of the eqn.
        
        L1 and L2 are the X and Y lengths.  Some shapes only take a single
        length.  These are:  beta, sinc, cone, pyramid, hanning, sphere, racos
        
        The kernel can also accept a norm, normalization, parameter.  The
        available options are

          "area" : normalize the kernel to have sum(pixels) = 1.0.  This
          is required if flux is to be preserved; however, for kernels
          with negative pixel values (mexhat, sinc) this can be bad
          since the integral is ~0.
          
          "max" : normalize the kernel such that the pixel with
          the max is 1.0:  pixels / max(pixels).
          
          "none": allow the 'M' paramter or other as appropriate to 
          control the absolute normalization.
    
    Examples:
    
    >>> from smooth_kernels import *
    >>> g = Gaussian(3)
    >>> print g
    lib:gaus(2,5,1.0,3.0,3.0)
    >>> g = Gaussian(3,7)
    >>> print g
    lib:gaus(2,5,1.0,3.0,7.0)
    >>> t = Tophat(5,norm="area")
    >>> t.get_spec()
    'lib:tophat(2,1.0,5.0,5.0)'
    >>> t.get_norm()
    'area'
    
    """
    
    def __init__( self ):
        """
        Initialize sensible defaults
        """
        self.xlen = None
        self.ylen = None
        self.norm = "area"
        self.nsig = None
        self.ktype = "lib"
        self.kname = None
        self.maxval = 1.0  # maybe make optional arg
        self.dampen = None
        
            
    def _parse_kwargs( self, kwargs ):
        """
        Currently only arg is norm.
        
        TBD: could implement the edge treatment here too.
        They go with aconvolve method=slide so doesn't quite fit
        
        TBD: could also add the override to the kernel center here.
        """
        
        if len(kwargs) > 1:
            raise ValueError("Too many named parameters specified")
        if len(kwargs) == 0:
            return
        if 'norm' not in kwargs:
            raise ValueError("'norm' is the only recongize optional parameter")
        
        if kwargs['norm'] not in ["area", "none", "max"]:
            raise ValueError("")

        self.norm = kwargs['norm']


    
    def __repr__(self):
        """
        Create the representation of the kernel in aconvovle syntax
        """
        retval = "{}:{}(2".format( self.ktype, self.kname )

        for val in [ self.nsig, self.maxval, self.dampen, self.xlen, self.ylen ]:
            if val:
                retval += ",{}".format(val)
        retval += ")"

        return retval

    
    def get_spec(self):
        """
        Get the kernelspec representation of the kernel
        """
        return self.__repr__()

    
    def get_norm(self):
        """
        Get the kernel normalization
        """
        return self.norm



class CartesianKernel( CIAOImageKernel ):
    """fs
    A base class for a kernel that takes separate X and Y lengths

    The aconvolve kernels are either circularly symmetric or
    separable in their X and Y components.
    
    This class implements the later.
    
    When given a single paramter it is used for both the
    x and y values.  When given a pair, it sets them
    seperately.
    
    Example:

    >>> g = Gaussian(3)
    >>> print g
    lib:gaus(2,5,1.0,3.0,3.0)
    >>> g = Gaussian(3,7)
    >>> print g
    lib:gaus(2,5,1.0,3.0,7.0)
        
    """
    def __init__( self, *args, **kwargs ):
        """
        Add checks to set the lenght based on the number of
        parameters
        """
        CIAOImageKernel.__init__(self)
        
        if len(args) == 0:
            raise ValueError("Must have 1 value")
        elif len(args) == 1:
            self.xlen = float(args[0])
            self.ylen = float(args[0])
        elif len(args) == 2:
            self.xlen = float(args[0])
            self.ylen = float(args[1])
        else:
            raise ValueError("Too many sigma given")
        
        self._parse_kwargs( kwargs )

  
class PolarKernel( CIAOImageKernel ):
    """
    A base class for a kernel that takes separate X and Y lenghts

    The aconvolve kernels are either circularly symmetric or
    separable in their X and Y components.
    
    This class implements the former.
    
    It only accepts a single length parameter.
    
    Example:

    >>> c = Cone(10)
    >>> print c
    lib:cone(2,1.0,10.0)
    
    """
    def __init__( self, *args, **kwargs ):
        """
        Add checks to set the lenght based on the number of
        parameters
        """
        CIAOImageKernel.__init__(self)
        
        if len(args) == 0:
            raise ValueError("Must have 1 value")
        elif len(args) == 1:
            self.xlen = float(args[0])
            self.ylen = None
        else:
            raise ValueError("Too many sigma given")
        
        self._parse_kwargs( kwargs )


class Gaussian( CartesianKernel ):
    """
    A 2D Gaussian

    A 2D, non-rotated, seperable Gaussian 
    
        N*exp(-((x^2)/L1)+((y^2)/L2))
    
    where L1 and L2 are the 1-sigma sizes in the x and y
    directions.

    If only 1 parameter is supplied, it is used for both
    the X and Y axis.
    
    The kernel is clipped at +/-5*sigma in each direction.
    
    Example:
    
    >>> g = Gaussian(3)
    >>> g = Gaussian(3,5)
    >>> g = Gaussian(3, norm="max")
    >>> g = Gaussian(3,5, norm="none")

    The optional norm parameter must be one of "area" (default), "max", or "none"

    """
    def __init__(self, *args, **kwargs):
        CartesianKernel.__init__(self, *args, **kwargs )
        self.kname = "gaus"
        self.nsig = 5

class Boxcar( CartesianKernel ):    
    """
    A 2D Box

    A 2D, non-rotated, rectangular box with fixed height
    
        Box(x,y) = N where x between xc +/- L1 & y between yc +/- L2 
    
    If only 1 parameter is supplied, it is used for both
    the X and Y axis.  Thus making a square.
        
    Example:
    
    >>> g = Boxcar(3)
    >>> g = Boxcar(3,5)
    >>> g = Boxcar(3, norm="max")
    >>> g = Boxcar(3,5, norm="none")

    The optional norm parameter must be one of "area" (default), "max", or "none"

    """
    def __init__(self, *args, **kwargs):
        CartesianKernel.__init__(self, *args, **kwargs )
        self.kname = "box"


class Tophat( CartesianKernel ):    
    """
    A 2D Elliptical Tophat

    A 2D, non-rotated, elliptical tophat
    
        Tophat = N where (x-xc)^2 < L1^2 and (y-yc)^2 < L2^2
    
    where L1 and L2 are the radii in the x and y
    directions.

    If only 1 parameter is supplied, it is used for both
    the X and Y axis.
        
    Example:
    
    >>> g = Tophat(3)
    >>> g = Tophat(3,5)
    >>> g = Tophat(3, norm="max")
    >>> g = Tophat(3,5, norm="none")

    The optional norm parameter must be one of "area" (default), "max", or "none"

    """
    def __init__(self, *args, **kwargs):
        CartesianKernel.__init__(self, *args, **kwargs )
        self.kname = "tophat"

class Mexhat( CartesianKernel ):    
    """
    A 2D Mexican Hat Wavelet

    A 2D, non-rotated, seperable Mexican hat wavelet
    
    Mexhat =   Norm * exp( -0.5 * (x/L1)^2 + (y/L2)^2 )
    
    where L1 and L2 are the 1-sigma sizes in the x and y
    directions.

    If only 1 parameter is supplied, it is used for both
    the X and Y axis.
    
    The kernel is clipped at +/-5*sigma in each direction.
    
    Example:
    
    >>> g = Mexhat(3)
    >>> g = Mexhat(3,5)
    >>> g = Mexhat(3, norm="max")
    >>> g = Mexhat(3,5, norm="none")

    The optional norm parameter should be either "max" or "none".  
    Since the integral of the wavelet is 0, using "area" may result
    in numeric instability (divide by zero).

    """
    def __init__(self, *args, **kwargs):
        CartesianKernel.__init__(self, *args, **kwargs )
        self.kname = "mexhat"
        self.nsig = 5.0

class Exp( CartesianKernel ):  
    """
    A 2D Exponential decay

    A 2D, non-rotated, seperable exponential decay

     exp = Norm * exp(-abs(x/L1)) * exp(-abs(y/L2)) 
    
    where L1 and L2 are the 1-sigma sizes in the x and y
    directions.

    If only 1 parameter is supplied, it is used for both
    the X and Y axis.
    
    The kernel is clipped at +/-5*sigma in each direction.
    
    Example:
    
    >>> g = Exp(3)
    >>> g = Exp(3,5)
    >>> g = Exp(3, norm="max")
    >>> g = Exp(3,5, norm="none")

    The optional norm parameter must be one of "area" (default), "max", or "none"

    """
    def __init__(self, *args, **kwargs):
        CartesianKernel.__init__(self, *args, **kwargs )
        self.kname = "exp"
        self.nsig = 5.0


class Power( CartesianKernel ):
    """
    A 2D powerlaw

    A 2D, non-rotated, seperable powerlaw
    
     power = Norm * abs(x)^L1 * abs(y)^L2 

    where L1 and L2 are the powerlaw index.  Typically these
    will be negative values.

    If only 1 parameter is supplied, it is used for both
    the X and Y axis.
    
    The kernel is clipped at +/-5*sigma in each direction.
    
    Example:
    
    >>> g = Power(3)
    >>> g = Power(3,5)
    >>> g = Power(3, norm="max")
    >>> g = Power(3,5, norm="none")

    The optional norm parameter must be one of "area" (default), "max", or "none"

    """
    def __init__(self, *args, **kwargs):
        CartesianKernel.__init__(self, *args, **kwargs )
        self.kname = "power"
        self.nsig = 5.0

class Beta( PolarKernel ):    
    """
    A 2D beta (or lorentizan) function

    A 2D, radially symmetric beta function
    
        beta = Norm / ((x^2 +y^2)/ L1^2)

    where L1 is the beta function index
    
    The kernel is clipped at +/-5*index in each direction.
    
    Example:
    
    >>> g = Beta(3)
    >>> g = Beta(3, norm="max")

    The optional norm parameter must be one of "area" (default), "max", or "none"

    """
    def __init__(self, *args, **kwargs):
        PolarKernel.__init__(self, *args, **kwargs )
        self.kname = "beta"
        self.nsig = 5.0

class Sinc( PolarKernel ):    
    """
    A 2D Sinc function

    A 2D, radially symmetric sinc function 
        
        d = sqrt(x^2+y^2)
        sinc = Norm * sin(d/L1)/(d/L1)

    where L1 is related to the period of the sin function.
    
    The kernel is clipped at +/-5*index in each direction.
    
    Example:
    
    >>> g = Sinc(3)
    >>> g = Sinc(3, norm="max")

    The optional norm parameter must be one of "area" (default), "max", or "none"

    """
    def __init__(self, *args, **kwargs):
        PolarKernel.__init__(self, *args, **kwargs )
        self.kname = "sinc"
        self.nsig = 5.0

class Cone( PolarKernel ):    
    """
    A 2D Cone

    A 2D, radially symmetric cone
    
        cone =  Norm * 1 - sqrt(x^2+y^2+...)  / s1

    where L1 is radius of the base of the cone.
    
    Example:
    
    >>> g = Cone(3)
    >>> g = Cone(3, norm="max")

    The optional norm parameter must be one of "area" (default), "max", or "none"

    """
    def __init__(self, *args, **kwargs):
        PolarKernel.__init__(self, *args, **kwargs )
        self.kname = "cone"

class Pyramid( PolarKernel ):    
    """
    A 2D Pyramid

    A 2D, radially symmetric pyramid shape
    
    pyramid = Norm * 1 - MIN(x,y,...)  / L1

    where L1 is the radius of the base of the pyramid
    
    Example:
    
    >>> g = Pyramid(3)
    >>> g = Pyramid(3, norm="max")

    The optional norm parameter must be one of "area" (default), "max", or "none"

    """
    def __init__(self, *args, **kwargs):
        PolarKernel.__init__(self, *args, **kwargs )
        self.kname = "pyramid"


class Sphere( PolarKernel ):    
    """
    A 2D Hemisphere

    A 2D Hemisphere
    
    sphere = sqrt( (x - L1)^2 + ( y - L1 )^2)

    where L1 is radius of the base of the hemisphere.
    
    Example:
    
    >>> g = Sphere(3)
    >>> g = Sphere(3, norm="max")

    The optional norm parameter must be one of "area" (default), "max", or "none"

    """
    def __init__(self, *args, **kwargs):
        PolarKernel.__init__(self, *args, **kwargs )
        self.kname = "sphere"
        self.maxval = None

class Hanning( PolarKernel ):    
    """
    A 2D Hanning function

    A 2D, radially symmetric hanning function.  This is a specialized 
    instance of the RaisedCosine() class with the parameter value set to 0.5.

    hanning = Norm * ( 0.5 - 0.5 * cos( 2*pi*sqrt( x^2 + y^2 + ...) / (L1-1) ))

    where L1 is radius of the base of the function.
    
    Example:
    
    >>> g = Cone(3)
    >>> g = Cone(3, norm="max")

    The optional norm parameter must be one of "area" (default), "max", or "none"
    """
    def __init__(self, *args, **kwargs):
        PolarKernel.__init__(self, *args, **kwargs )
        self.kname = "hanning"

class RaisedCosine( PolarKernel ):    
    """
    A 2D Raised Cosine function

    A 2D, radially symmetric raised cosine function.  (Cosine()+offset)
    
    racos = Norm * ( alpha - ( 1 - alpha) * cos( 2*pi*sqrt( x^2 + y^2 ) / (L1-1) ))

    where L1 is radius of the base of the function.
    
    Example:
    
    >>> g = RaisedCosine(0.5, 3)
    >>> g = RaisedCosine(0.54, 3, norm="max")

    The optional norm parameter must be one of "area" (default), "max", or "none"

    """

    def __init__(self, dampen, *args, **kwargs):
        PolarKernel.__init__(self, *args, **kwargs )
        self.kname = "racos"
        self.dampen = dampen

#-----------------

class AdaptiveKernel( object ):
    
    def __init__(self, minrad=1, maxrad=100, numrad=100, radscale="linear"):
        self.minrad = minrad
        self.maxrad = maxrad
        self.numrad = numrad
        self.radscale = radscale
            
class AdaptGaussian( AdaptiveKernel):
    def __init__(self, *args, **kwargs):
        self.function="gaussian"
        AdaptiveKernel.__init__( self, *args, **kwargs )        

class AdaptTophat( AdaptiveKernel):
    def __init__(self, *args, **kwargs):
        self.function="tophat"
        AdaptiveKernel.__init__( self, *args, **kwargs )

class AdaptBoxcar( AdaptiveKernel):
    def __init__(self, *args, **kwargs):
        self.function="box"
        AdaptiveKernel.__init__( self, *args, **kwargs )

class AdaptCone( AdaptiveKernel):
    def __init__(self, *args, **kwargs):
        self.function="cone"
        AdaptiveKernel.__init__( self, *args, **kwargs )

class AdaptPyramid( AdaptiveKernel):
    def __init__(self, *args, **kwargs):
        self.function="pyramid"
        AdaptiveKernel.__init__( self, *args, **kwargs )

class AdaptSphere( AdaptiveKernel):
    def __init__(self, *args, **kwargs):
        self.function="hemisphere"
        AdaptiveKernel.__init__( self, *args, **kwargs )

class AdaptQuad( AdaptiveKernel):
    def __init__(self, *args, **kwargs):
        self.function="quad"
        AdaptiveKernel.__init__( self, *args, **kwargs )

class AdaptExp( AdaptiveKernel):
    def __init__(self, *args, **kwargs):
        self.function="exp"
        AdaptiveKernel.__init__( self, *args, **kwargs )







