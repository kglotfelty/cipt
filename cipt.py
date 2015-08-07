# CIAO Image Processing Toolkit


#
# TODO:  NamedTempFiles() use os.environ["ASCDS_WORK_PATH"]
#

__all__ = [ "CIAOImage" ]


from crateify import Crateify
from smooth_kernels import *
from filter_mask import *
from history_crate import *
from ciao_contrib.runtool import make_tool




class CIAOImage( HistoryIMAGECrate ):
    """
    CIAO Image 

    The CIAOImage class extends the standard CIAO IMAGECrate and adds
    methods to perform various image processing tasks.
    
    All the methods return a new instance of the CIAOImage object
    so they can be chained together.
    
    >>> img=CIAOImage("img.fits")
    >>> img.smooth(Gaussian(3)).median(Annulus(20,22)).write("out.fits")
    
    or if using Chips
    
    chips>  add_image( img.smooth(Cone(3)).filter(Region("ciao.reg")))

    All the currently implemented method return an image the
    same size as the original. New methods may change that.
    
    The CIAOImage extends the IMAGECrate to preserve the
    history records so that when the crate is written
    using the CIAOImage.write() method it will contain
    all the processing commands used.
    
    There is support for simple arithmetic operations between
    CIAOImage objects of the same size.
    
    >>> img1 = CIAOImage("img1.fits")
    >>> img2 = CIAOImage("img2.fits")
    >>> img1+img2
    >>> img1*img2
    >>> img1/img2
    >>> img1-img2
    >>> 3*img1
    >>> 1-img1-1
    >>> img/img.get_key_value("ONTIME")    
    >>> img**2

    Users can also do simple arithmetic with NumPy arrays that are the 
    same size as the CIAOImage.

    >>> import numpy as np
    >>> shape = img.get_image().values.shape
    >>> noise = np.random.poisson( 0.1, size=shape )
    >>> img2 = img + noise
            
    The following methods are available.  They are broken down
    into several categories.  They are built using various CIAO
    tools

    Smoothing:


        Perform linear smoothing with a fixed convolution kernel
        using the 'aconvolve' tool.
        
        The available smoothing kernels are:

        : gaus( x [,y [, norm="area"|"max"|"none" [,**kwargs]]] )
            A Gaussian profile
        : tophat( x [,y [, norm="area"|"max"|"none"[,**kwargs]]] )
            A circular (or elliptical) tophat
        : boxcar( x [,y [, norm="area"|"max"|"none"[,**kwargs]]] )
            A square (or rectangular) boxcar
        : mexhat( x [,y [, norm="max"|"none"[,**kwargs]]] )
            A Mexian Hat wavelet
        : expno( x [,y [, norm="area"|"max"|"none"[,**kwargs]]] )
            A exponential fallout
        : power( x [,y [, norm="area"|"max"|"none"[,**kwargs]]] )
            A powerlaw fallout
        : beta( x [, norm="area"|"max"|"none"[,**kwargs]]] )
            A beta (or Lorentizen) function
        : sinc( x [, norm="max"|"none"[,**kwargs]]] )
            A Sinc (sin(x)/x) function
        : cone( x [, norm="area"|"max"|"none"[,**kwargs]]] )
            A circularly symmetric cone
        : pyramid( x [, norm="area"|"max"|"none"[,**kwargs]]] )
            A square pyramid shape
        : sphere( x [, norm="area"|"max"|"none"[,**kwargs]]] )
            A hemisphere/quadratic shape
        : hanning( x [, norm="area"|"max"|"none"[,**kwargs]]] )
            A hanning filter
        : racos( coef, x [, norm="area"|"max"|"none"[,**kwargs]]] )
            A raised cosine filter with coefficient=coef.

        If the Y value is not specified it is taken to be same as
        the X value.
        
        The normalization parameter varies by the kernel type.
        
        The kwargs are optional parameter that are passed to the
        aconvolve tool.  The most useful will be the 'method'
        parameter to control whether to use FFT convolution
        (default) or to use sliding cell (slower, but 
        has control over the edges).

        All units are in logical/image pixel sizes

        Examples:
        
        >>> img.gaus(3)
        >>> img.gaus(3,6)
        >>> img.gaus(3,6, norm="max")
        >>> img.gaus(3, method="slide", edge="mirror")
        >>> img.tophat(5)
        >>> img.boxcar(3,3)
        >>> img.mexhat(20)-img.mexhat(2)
        >>> img.expno(-1)  # negative exponents are good thing!
        >>> img.power(-2) # 1/n^2
        >>> img.beta(-1) 
        >>> img.sinc(100)
        >>> img.cone(3)
        >>> img.pyramid(10)
        >>> img.sphere(5)
        >>> img.hanning(10)
        >>> img.racos(0.5, 10) - img.hanning(10)
        
        Users can also convolve two CIAOImage's or a CIAOImage with 
        a NumPY array.
        
        >>> img.convolve( img2 )
        >>> img.convolve( np.array( [[1, -1], [-1, 1]] ), norm="none")


    Adaptive Smoothing
            
        Perform an adaptive smooth.  Increase the kernel size
        until the min_counts value is reached.  This 
        uses the dmimgadapt tool.
        
        The available routines are for each smoothing kernel.
        
        : agaus( counts, **kwargs ) :
            Adaptive smooth with a Gaussian smoothing kernel

        : atophat( counts, **kwargs ) :
            Adaptive smooth with a circular tophat smoothing kernel

        : abox( counts, **kwargs ) :
            Adaptive smooth with a square box smoothing kernel

        : acone( counts, **kwargs ) :
            Adaptive smooth with a circular cone smoothing kernel

        : apyramid( counts, **kwargs ) :
            Adaptive smooth with a square pyramid smoothing kernel

        : asphere( counts, **kwargs ) :
            Adaptive smooth with a hemisphere smoothing kernel

        : aquad( counts, **kwargs ) :
            Adaptive smooth with a quadratic (1/n^2) smoothing kernel

        : aexp( counts, **kwargs ) :
            Adaptive smooth with an exponential/beta smoothing kernel
            

        The kwargs are the dmimgadapt parameter.  The following
        defaults are used:

            minrad=1
            maxrad=100
            numrad=100
            radscale="linear"
            
        Examples:
        
        >>> img.agaus(10)
        >>> img.agaus(10, maxrad=20)
        >>> img.agaus(100, minrad=1, maxrad=100, numrad=20)
        >>> img.atophat(400) 
        >>> img.abox(10)
        >>> img.acone(36)
        >>> img.apyramid(1.0e-3)
        >>> img.asphere(100)
        >>> img.aquad(25)
        >>> img.aexp(40)

        : csmooth( ... )
            
            An adaptive smoothing algorithm.

        >>> img.csmooth( )


        : adaptive_bin( snr )
            
            Though not often presented as such, binning can be considered
            a type of smoothing.  In the traditional binning schemes
            there is also a decimation (down-sampling) of the data to
            a new grid but it can also be done such that the
            original grid is preserved.
            
            Examples:
            
            >>> img.adaptive_bin( 5 )

    Non Linear Smoothing:

        Non-linear filtering or smoothing is the process of taking
        a set of pixel values, applying some function to them and
        creating a new image whose pixel value is the result of that
        function call.  The function need not do a linear 
        transformation of the value (such as a sum) but can do
        anything (like count the number of good pixels).  As such
        the units of the output image are not the same as the input
        and they definitely do no conserve flux.
        
        All the non-linear filters accept a FilterMask object
        that describes the pixels to include in the search.  The
        center of the mask is usually at (0,0) so that as it is moved
        from image pixel to pixel the surrounding pixels are included
        in the function.
        
        For example:  min(Box(3)) will create a 3x3 box centered at 0,0.
        That box will be moved across the image.  At each pixel the
        3x3 pixel around the current pixel location are inspected and
        the miniumum value of those 9 pixel values is returned.

        : min(FilterMask) : find minimum pixel value of those in masked region 

        : max(FilterMask) : find maximum pixel value of those in masked region 

        : mean(FilterMask) : find mean pixel value of those in masked region 

        : median(FilterMask) : find median(*) pixel value of those in 
                              masked region 

        : mode(FilterMask) : 3 * median() - 2*mean().  A measure of how 
                            skewed the pixel values are

        : nmode(FilterMask) : Normalized mode, ie mode()/mean().
        
        : mid(FilterMask) : The value mid way between the min and max

        : sigma(FilterMask) : The standard deviation of the pixel values in 
                            the mask

        : extreme(FilterMask) : If the mean is closer to min value, use the 
                            min; otherwise use the max.

        : locheq(FilterMask) : local histogram equalization.  The center 
                            pixel (0,0) is scaled based on the distribution 
                            of pixels in the masked region.

        : kuwahara(FilterMaskd) : An edge preserving smoothing.  The mask is 
                            divided into quadrants.  The mean value in
                            the quadrant with the smallest standard deviation
                            is used.

        : unsharp(FilterMask): center pixel minus the mean of the pixels in the mask.

        : range(FilterMask) : The range of pixel values (max-min).

        : variance(FilterMask) : The variance of the pixel values compared to the
                            center pixel

        : q25(FilterMask) : 25% of the pixels in the mask are below this value.

        : q33(FilterMask) : 33% of the pixels in the mask are below this value.

        : q67(FilterMask) : 67% of the pixels in the mask are below this value.

        : q75(FilterMask) : 75% of the pixels in the mask are below this value.

        : mcv(FilterMask) : most common value.  A coarse histogram of the pixel
                        values is created and the location of the peak is returned.

        : sum(FilterMask) : sum of the values in the mask

        : rclip(FilterMask) : If the center pixel value is less than the min 
                        pixel in the mask, replace with the min value. If it is 
                        greater than the max, replace with the max. Otherwise, 
                        use original value.  This should be used with FilterMasks
                        that exclude the center pixel (eg Annulus)

        : peak(FilterMask) : If the center pixel is the max value of the pixels
                        in the mask return it.  Otherwise return NaN.  Used
                        to identify local maximums.

        : valley(FilterMask) : If the center pixel is the min value of he pixels
                        in the mask return it.  Otherwise return NaN.  Used
                        to identify local minimum.

        : count(FilterMask) : Count the number of valid pixels in the mask.
                        Different at the edge of image and at/around
                        any NaN values.

        : olympic(FilterMask) : Same as mean() except the lowest and
                        largest values are excluded.

        : pmean(FilterMask) : Poisson mean computed as (#Pixels with 
                        value = 0 or 1 divided by #Pixels with value = 
                        0 ) minus 1. If #pixels_0 is 0 then use median.
                        Should only be used for integer images

        : mu3(FilterMask) : The 3rd moment of the pixel values

        : mu4(FilterMask) : The 4th moment of the pixel values 
 
        : jitter(FilterMask) : Randomly select one of the pixel values in the
                        mask

        : rms(FilterMask) : The root-mean-square

        : nslope(FilterMask) : The minimum difference between pixel values
                        in the mask.

        : sig3mean(FilterMask) : The mean, after 5 iterations where data 
                        outside +/-3 sigma are excluded

        : sig3median(FilterMask) : The median after 5 iterations where data 
                        outside +/-3 sigma are excluded

        Examples:
        
        >>> img.min(Box(3))
        >>> img.smooth(Gaussian(3)).rclip(Annulus(20,22))
        >>> img.median( Box(10,10)-Circle(3) )
        
        (*) The quantile based algorithms use a slightly modified definition.
        If the number of points is even, rather than compute the 
        average of the two value, the lower value is returned.

        The filters can be defined using the following objects
        
          : Box(xlen [,ylen[, center=(xc,yc) [,angle]]])
          : Ellipse(xlen [,ylen[, center=(xc,yc) [,angle]]])
          : Annulus(xlen ,ylen[, center=(xc,yc) ] )
          : Circle(xlen [, center=(xc,yc) ])
          : Point([center=(xc,yc)] )
          : Pie(inner, outer, start, stop[, center=(xc,yc)] )
          : Sector(start, stop [,center=(xc,yc)])
          : Rectangle(xmin, ymin, xmax, ymax)
          : Polygon( (x1,y1), (x2,y2), ..., (xn,yn))
          : Field()
          : Region(filename)

        Examples:
        
        >>> Box(3)
        >>> Box(3,4)
        >>> Ellipse(50,100,angle=45)
        >>> Circle(1,center=(50,50))
        >>> Point()
        >>> Rectangle( 10, 10, 100, 100)
        >>> Polygon( (1,1), (1,2), (2,2), (2,1) )
        >>> Region("ciao.reg")
        >>> Region("ciao.fits[#row=1:10]")

        The atomic filter shapes can be combined with the "+", "*",
        and "-" operators.  The "-" is limited in the complexity of the
        regions it can successful negate.

        Examples:

        >>> Box(10)-Box(4)  # a box annulus
        >>> Circle(10)+Circle(1,center=(-8,8))+Circle(1,center=(8,8))  # Micky Mouse

        >>> el1 = Ellipse( 50, 200, center=(300,300), angle=45)
        >>> el2 = Ellipse( 50, 200, center=(300,300), angle= 135)        
        >>> elxor = el1-el2 + el2-el1    # Exclusive OR, XOR
        
        >>> Field()-Region("ciao.reg")
    

    Scale
    
        The pixel values in the CIAOImage can be scaled with various
        filters.  These may be useful especially when displaying
        the images.  The most common scales are:  log(), sqrt(), power(),
        growing in popularity asinh().  The rest are provided more for
        completeness.

        : cos()
        : sin()
        : tan()
        : acos()
        : asin()
        : atan()
        : cosh()
        : sinh()
        : tanh()
        : acosh()
        : asinh()
        : atanh()
        : exp()
        : log()
        : ln()
        : sqrt()
        : fabs()
            absolute value      
    
        Examples:
        
        >>> img.ashin()
        >>> img.log()


    Cast

        The datatype of the image can also be modified by calling
        one of the datatype casting methods.  

        : byte()
        : short()
        : ushort()
        : int()
        : long()
        : ulong()
        : float()
        : double()
        
        Examples:
        
        >>> img.short()
        >>> img.smooth(Gaussian(3)).double()
        
        
    Filter

        : filter(region[,nullval=0[, coord="sky"]])            
            Applies a CIAO style region string to the image.  
            The NULL value (pixels outside the image) are set to
            the nullval value.  If nullval=None, the values are
            set to NaN.
            
            The coord parameter controls the coordinate names to
            filter on.  If not a Chandra image, the other common
            names are pos, (x,y), (#1,#2) or various grating
            (tg_r,tg_d) kind of values.

            The region string may also be constructed from
            the same objects used in the adaptive_filter
            section above.
            
            Examples:
            
            >>> img.filter("circle(4096,4096,10)")
            >>> img.filter(box(4096,4096,10,10)+ellipse(3290,4332,10,50,45))
            >>> img.filter("region(ds9.reg)")
            >>> img.filter(Circle(10,center=(4096,4096)))
            >>> img.filter(Field()-Region(ds9.reg))
            >>> img.filter(Box(100), coord="(#1,#2)")
            >>> img.double().filter("box(4000,4500,10,12)", nullval=None)
            
            The image size remains the same before and after 
            filtering.
            
        : thresh( [lo=0[,hi=None[, val=0]]] )
            Applies a lower and/or upper threshold to the pixel values.
            
            Examples:

            >>> img.thresh()   # values < 0 are set = 0
            >>> img.thresh(10)  # values < 10 are set 0
            >>> img.thresh(None, 10) # values > 10 are set to 0
            >>> img.thresh(None,10, 10) #  values > 10 are set to 10
            >>> img.thresh(1,val=1).thresh(hi=10,val=10)


         : fill( src[,bkg=None[,method="poisson"]] )
            Replaces the pixels in the source region with values
            derived by different methods in the background region.
            
            >>> img.fill( circle(4096,4096,10), annulus(4096,4096,20,30) )
            >>> img.fill( circle(3900,4200,20), method="global")
            
            The method parameter may take one of the following values:
    
            'poly' - a 2D polynomial is fit to the pixels in the background,
                    those fit parameters are used to fill in source region pixels.
                    Background must full enclose and include source.
            
            'dist' - a coarse histogram of the pixel values in the background
                    region is created. This is used a PDF for the src pixels values

            'global' - same as 'dist' but values are taken from the entire
                    image (except the source pixels).  No background region.

            'poisson' - Uses mean counts/pixel in background as mean
                    Poisson random variable for source region

            'bilint' - bilinear interpolation of pixels in the background
                    region over the source.  Background must fully
                    enclose, but not include, source.
            

    Translation

        : rotate(angle [,center=(xc,yc) [,**kwargs]])
            Rotate the image, counter clock wise from the +X axis,
            around the center (middle of the image if not specified).

        : scale(sx [,sy, [,center=(xc,yc) [,**kwargs]]])
            Scale the pixel size.  sx=2 means the pixels are 2*2=4 times
            as large so the image "shrinks".  A scale of 0.5 will
            double the pixel scale.
            
            If sy is not specified, it is set to the same as sx.

        : shift(dx, dy, [,**kwargs]])

            Apply a linear shift, in image pixel, to the data.  Values
            may be non-integer.
            
        The kwargs are passed onto the dmregrid2 routine.  The 
        most useful parameter is
        
            method="sum"|"average"
            
        If scaling a counts image you want to use the default "sum"
        method.  If you are scaling an exposure map then the "average"
        method is correct.

        Examples:
        
        >>> img.shift(10, 33.333)
        >>> img.rotate(45)
        >>> img.rotate(45, center=(0,0))
        
        These routines do not update the WCS of the image.  These
        routines also work in the image pixel size as the input image
        so if you
        
        >>> img.rotate(45).rotate(-45)
        >>> img.scale(0.5).scale(2)

        you may lose pixel value at the edge of the pixels that are
        not recoverable.


        : match(other [, **kwargs]) 

        The match method runs the reproject_image routine to match
        the coordinate systems for the current object with the
        WCS in the other object.  That is the current 
        CIAOImage is transformed via its WCS to match the other
        image.  
        
        >>> img1 = CIAOImage("chandra.fits")
        >>> img2 = CIAOImage("hst.fits")
        >>> overlay = img1.match(img2)
        
        The pixel values in the overlay CIAO image will now match
        the WCS in the img2 datase.

                    
    Transform

        The methods in this section are a collection of routines
        that do not otherwise fit into the above categories.
        They transform the input crate into something
        else entirely.

        : powerspectrum()
            Computes the powerspectrum of the image pixel values
            (magnitude of the FFT).
            
            >>> img.powerspectrum()

        : blob( [threshold=~0, srconly=True)
            The blob method identifies contiguous groups of pixels
            that are above threhold.  The returned
            object represents a map that identifies which pixels belong 
            together in which group.
            
            If no threshold is given, a value ~0 is used.
            
            Examples:
            
            >>> img.blob()
            >>> img.smooth(Gaussian(3)).blob(0.2)  # a crude source detect

        : distance( edgevalue )
            Using 'edgevalue' is the value below which is to be considered
            the edge of the image, how far away, ie number of pixels,
            is the current pixel away from that edge.  distance()
            use the city-block style distance computation such
            that ever pixel is an integer number of pixels away
            from the edge.
            
            Examples:

            >>> img.distance(0)
            >>> img.smooth(Gaussian(3)).distance(0.01).thresh(None, hi=1,val=None)
            
            The last example provides a way of generating a contour around
            the edgevalue, value. 

        : correlate( other=None )
            Correlation is the same as convolution with the convolution
            kernel flipped (so two are identical for a symmetrical kernel).
            The leaving the other parameter blank yields the autocorrelation.

            >>> img.correlate()
            >>> img.correlate(img2)
            >>> img.correlate( np.array( [[1,2,3],[1,2,3],[1,2,3]]) )

        : deconvolve( ... )
            Blind deconvolution is always a tricky subject, but 
            users can use the Richardson-Lucy deconvolution method
            available in CIAO's arestore         
        
            >>> img.deconvolve( img2)
            >>> img.deconvolve( np.array( [[1,1,1], [1,3,1], [1,1,1]]) )
                        

    The following method return a TABLECrate, not an Image
    
      : histogram( [grid="1"] )
        Returns a histogram of the image pixel values by running the
        dmimghist tool
        
        
            >>> tab = img.histogram()
            >>> tab = img.histogram("1:100:1")
            
      : project( "x"|"y")
      
        Returns projections of the image onto either the X or Y axis

            >>> tab = img.project("x")
            >>> tab = img.project("y")
            
    The following methods return a Region object

      : contour(levels)
      
         Create polygons that follow the specified contour levels.
         
             >>> reg = img.contour("1")
             >>> reg = img.contour("1,2,4,8,16")
    
      : ellipse( fraction, **kwargs )
      
          Returns an ellipse that encloses the specified fraction of
          the data in the image.  Additional keyword arguments are 
          passed to the dmellipse tool.
          
            >>> ell = img.ellipse( 0.5 )
            >>> bb = img.ellipse(0.5, shape="rotbox")
            >>> cir = img.ellipse(0.9, ellipticity=1, fix_ellipticity=True)

      : hull( [tolerance] )
        
          Returns the convex hull around the image pixel value above
          the specified tolerance value
          
            >>> hull = img.hull()
            >>> h1 = img.hull(1)

    The following methods return Python dictionaries
    
      Coordinates
      
      The coordinate routines below use the CIAO tool dmcoords
      to convert from the input coordinate systems using the
      CIAOImage's transform and Chandra's pixlib library.
      Non Chandra images can generally convert between
      physical, logical, and usually world coorindates.  The
      det, msc, and chip coordinates are exclusively for Chandra
      datasets.      
      
      Each of these methods returns a dictionary containing all the
      coordinate values.  Note:  Ra/Dec are always returned in decimal
      degrees.
      
      : world( ra, dec ) 
            
            World coordinates, ie Right Ascension  and Declination
            (for Chandra in the J2000, fk5,system).

            >>> c = img.world( 34.5569, -14.566 )
            >>> c = img.world( "23 45 12.232", "-43:32:19.6" ) # TODO 

      : physical( x, y )
            
            The physical image pixel system.  Generally this represents 
            the data within entire image size without and
            cropping.

            >>> c = img.physical( 4096, 4096 )

      : logical( lx, ly )

            The logical pixel system going from 1 to (M,N) in the
            Cartesian X and Y system.   [Note: Python using a 0
            based indexing system, so the values in the array
            are indexed from 0 to (M-1, N-1).]
            
            >>> c = img.logical(10,20)

      : det( detx, dety )

            The Chandra detector coordinate system, a 
            coordinate system defined in detector plane of the
            Chandra telescope.

            >>> c = img.det( 16385, 16385 )

      : msc( theta, phi)

            MSC: Mirror Spherical Coordinates, are the pair
            of angles needed to describe the location relative to
            the Chandra optical axis.  Theta is the 
            angular distance, measured in arcmin, away from the 
            optical axis (ie the 'off axis angle') and phi is
            the azimuthal angle in degrees around the optical axis.

            >>> c = img.msc( 5, 45 )

      : chip( id, chipx, chipy )
       
            Chandra chip coordinates for the detector specified 
            in the CIAImage's header/meta-data.  
            
            >>> c = c.chip(7, 512, 512 )
            >>> print img.chip(7,512,512)
            {'chip_id': 7.0, 'phi': 358.4371404326323, 'tdetx': 4429.0, 
            'tdety': 2214.0, 'logicaly': 183.1776032233852, 'logicalx': 
            182.7295184280113, 'dety': 4087.794395102537, 'detx': 
            4415.575824689214, 'ra': 198.8516702952409, 'chipy': 512.0, 
            'chipx': 512.0, 'y': 3823.210412893541, 'x': 
            4261.418073712045, 'dec': -16.40350088354444, 'theta': 
            2.61739491534039}
      
      Statistics
      
      : stats()
            
            Uses dmstat to compute various statistics about the 
            CIAOImage.  dmstat knows when pixels have been filtered out
            either via the image subspace information or due to the
            presense of Null or NaN pixel values.
      
            >>> stat = img.filter(circle(4096,4096,100)).stats()
            >>> print stat
            {'max_loc': [4170.5, 4042.5], 'numnull': 131256, 'max': 45.0, 
            'lcentroid': [144.55766115, 247.85339395], 'min': 0.0, 'sum': 
            22216.0, 'median': 10.0, 'pixelarea': 1968, 'centroid': 
            [4108.7306446, 4081.9135758], 'stddev': 5.7151233298, 
            'min_loc': [4022.5, 4162.5], 'mean': 11.288617886}

            The centroid, location of minimum pixel (min_loc),
            and location of maximum pixle (max_loc) are given in
            the physical coordinate system.  lcentroid is
            the centroid in the logical coordinate system.
      
      : moments()

            Uses imgmoment to compute the 2nd order moments matrix
            for the image pixel values.  The most useful application
            is then to use these to compute the centroid, as well
            as the rotation angle and the major and minor axes of 
            the moments of inertia -- commonly used as the elliptical
            region that encloses the pixel values.
      
            >>> mom = img.filter(reg).moments
            >>> print mom
            {'angle': 54.6766949999169, 'major': 102.4510379626248, 'centroid': [4281.081534650766, 3950.348949610924], 
            'matrix': array([[  1.28381000e+05,   3.45607987e-10,   3.20896463e+08],
                            [ -4.79485607e-09,  -1.13257789e+07,   1.17630918e+08],
                            [  3.28852636e+08,  -5.38440130e+08,   6.48365063e+11]]), 
            'minor': 98.73289942300798}


    Finally, if users have the pyds9 package installed, they can use the
    ds9 method to display the CIAOImage
    
      : ds9()

            If ds9 is running it will display the CIAOImage in 
            the current frame (overwriting anything that's already there); 
            ds9 will be launched if one is not running.
            
            It returns the pyds9 object which can be used to interact
            with the data being displayed

            >>> d = img.ds9()
            >>> d.set( "cmap bb")
            >>> d.set( "scale log")
      
            Basically the same XPA commands available on the 
            command line are available via the pyds9 module.


            
    """
    class _Callable():
        """
        Callable methods
        
        These are methods that implement a callable object
        """
        def __init__( self, crate, fun ):
            self.crate = crate
            self.fun = fun
        def __doc__(self):
            return "{0} : Apply {0} function to image pixel values".format(self.fun)

    class _NLFilter(_Callable):
        """
        Non Linear Filters
        
        Various forms of non-linear filters applied to the image
        """
        def __call__( self, mask ):
            dmimgfilt = make_tool("dmimgfilt")
            return _run_cmd( dmimgfilt, self.crate, function=self.fun, mask=mask.__str__())
        
    class _Scale(_Callable):
        """
        Mathematical scaling
        
        Include trigonometric, log, exponent, and absolute value
        """
        def __call__(self):
            dmimgcalc=make_tool("dmimgcalc")
            return _run_cmd( dmimgcalc, self.crate, op="imgout={}(img1*1.0)".format(self.fun), infile2="" )

    class _Cast(_Callable):
        """
        Cast to specific datatype
        """
        def __call__(self):
            dmimgcalc=make_tool("dmimgcalc")
            return _run_cmd( dmimgcalc, self.crate, op="imgout=({})(img1)".format(self.fun), infile2="" )


    def asmd( self, other, op, right ):
        """
        Add, Subtract, Multiply, Divide
        """
        import numbers as numbers
        dmimgcalc=make_tool("dmimgcalc")
        opmap = { '+' : 'add', '-' : 'sub', '*' : 'mul', '/' : 'div' }

        if isinstance(other, CIAOImage) or isinstance(other, np.ndarray):
            with serialize_temp_crate( other ) as infile2:
                return _run_cmd( dmimgcalc, self, infile2=infile2.name, op=opmap[op] )
        elif isinstance(other, numbers.Real):
            if right:            
                return _run_cmd( dmimgcalc, self, infile2=None, op="imgout=({1}{0}img1)".format(op,str(other)))
            else:
                return _run_cmd( dmimgcalc, self, infile2=None, op="imgout=(img1{0}{1})".format(op,str(other)))
        else:
            # Todo, add numpy array, must be same size as self
            raise NotImplementedError("Cannot combine class")

    def __add__(self,other):        
        return self.asmd( other, "+" ,False)    
    def __radd__(self, other):
        return self.asmd( other, "+", True)
    def __sub__(self,other):        
        return self.asmd( other, "-", False)    
    def __rsub__(self, other):
        return self.asmd( other, "-", True)
    def __mul__(self,other):        
        return self.asmd( other, "*", False)    
    def __rmul__(self, other):
        return self.asmd( other, "*", True)
    def __div__(self,other):        
        return self.asmd( other, "/", False)    
    def __rdiv__(self, other):
        return self.asmd( other, "/", True)

    def __pow__( self, other ):
        """        
        """
        import numbers as numbers
        dmimgcalc=make_tool("dmimgcalc")
        if isinstance(other, numbers.Real):
            return _run_cmd( dmimgcalc, self, op="imgout=(img1^{})".format(other), infile2="" )        
        raise NotImplementedError("Must use a real number for power")


    def __init__( self, filename ):        
        HistoryIMAGECrate.__init__(self,filename)
        self.get_history( filename)

        for nl in ["min","max","mean","median","mode","mid","sigma","extreme",
            "locheq","kuwahara","unsharp","range","variance","nmode",
            "q25","q33","q67","q75","mcv","sum","rclip","peak","valley",
            "count","olympic","pmean","mu3","mu4","jitter","rms","nslope"]:
                setattr(self, nl, self._NLFilter( self, nl ) )
                setattr(getattr(self,nl),"__doc__", "Apply {} funtion to all pixels".format(nl))

        # Cannot have methods that start with numbers, so handle this separately 
        setattr( self, "sig3mean", self._NLFilter( self, "3sigmean") )
        setattr( self, "sig3median", self._NLFilter(self, "3sigmedian"))

        for scl in  ["cos", "sin", "tan", "acos", "asin", "atan", "cosh", 
            "sinh", "tanh", "exp", "log", "ln", "sqrt", "fabs",
            "asinh", "acosh", "atanh"]:
                setattr(self, scl, self._Scale( self, scl ))

        for cst in ["byte", "short", "long", "ushort", "ulong", "float", "double", "int" ]:
                setattr(self, cst, self._Cast( self,cst ))

    
    def _smooth( self, kernel, method="fft", **kwargs ):
        """
        Don't call this one, call the ones below
        
        : _smooth( CIAOImageKernel(), **kwargs )
            Gaussian( x [,y [, norm="area"|"max"|"none"]] )
            Tophat( x [,y [, norm="area"|"max"|"none"]] )
            Boxcar( x [,y [, norm="area"|"max"|"none"]] )
            Mexhat( x [,y [, norm="max"|"none"]] )
            Exp( x [,y [, norm="area"|"max"|"none"]] )
            Power( x [,y [, norm="area"|"max"|"none"]] )
            Beta( x [, norm="area"|"max"|"none"]] )
            Sinc( x [, norm="max"|"none"]] )
            Cone( x [, norm="area"|"max"|"none"]] )
            Pyramid( x [, norm="area"|"max"|"none"]] )
            Sphere( x [, norm="area"|"max"|"none"]] )
            Hanning( x [, norm="area"|"max"|"none"]] )
            RaisedCosine( coef, x [, norm="area"|"max"|"none"]] )

            Any additional keyword arguments are passed to
            aconvolve.  These may include the method parameter
            to select FFT vs. sliding cell convolution and the
            associated parameters.
            
            Examples:
            
            >>> img.smooth(Gaussian(3))
            >>> img.smooth(Tophat(3,5))
            >>> img.smooth(Cone(7, norm="max"))
            >>> img.smooth(Boxcar(4,4), method="slide", edge="constant")
        """        
        aconvolve = make_tool("aconvolve")
        try:
            spec = kernel.get_spec()
            norm = kernel.get_norm()
        except Exception, e:
            print "Need to supply a Kernel object"
            raise e        
        out = _run_cmd( aconvolve, self, kernelsp=spec, norm=norm, method=method, **kwargs)
        return out

    def gaus( self, x, y=None, norm="area", **kwargs):
        """        
        """
        y = y if y else x
        k = Gaussian( x, y, norm=norm)
        return self._smooth( k, **kwargs )

    def tophat( self, x, y=None, norm="area", **kwargs):
        """        
        """
        y = y if y else x
        k = Tophat( x, y, norm=norm)
        return self._smooth( k, **kwargs )

    def box( self, x, y=None, norm="area", **kwargs):
        """        
        """
        y = y if y else x
        k = Boxcar( x, y, norm=norm)
        return self._smooth( k, **kwargs )

    def mexhat( self, x, y=None, norm="none", **kwargs):
        """        
        """
        y = y if y else x
        k = Mexhat( x, y, norm=norm)
        return self._smooth( k, **kwargs )

    def expo( self, x, y=None, norm="area", **kwargs):
        """        
        """
        y = y if y else x
        k = Exp( x, y, norm=norm)
        return self._smooth( k, **kwargs )

    def power( self, x, y=None, norm="area", **kwargs):
        """        
        """
        y = y if y else x
        k = Power( x, y, norm=norm)
        return self._smooth( k, **kwargs )

    def beta( self, x,  norm="area", **kwargs):
        """        
        """
        k = Beta( x,norm=norm)
        return self._smooth( k, **kwargs )

    def sinc( self, x, norm="none", **kwargs):
        """        
        """
        k = Sinc( x,norm=norm)
        return self._smooth( k, **kwargs )

    def cone( self, x, norm="area", **kwargs):
        """        
        """
        k = Cone( x, norm=norm)
        return self._smooth( k, **kwargs )

    def pyramid( self, x, norm="area", **kwargs):
        """        
        """
        k = Pyramid( x, norm=norm)
        return self._smooth( k, **kwargs )

    def sphere( self, x, norm="area", **kwargs):
        """        
        """
        k = Sphere( x,norm=norm)
        return self._smooth( k, **kwargs )

    def hanning( self, x, norm="area", **kwargs):
        """        
        """
        k = Hanning( x,norm=norm)
        return self._smooth( k, **kwargs )
 
    def racos( self, x, dampen, norm="area", **kwargs):
        """        
        """
        k = RaisedCosine( x,dampen, norm=norm)
        return self._smooth( k, **kwargs )

    def _adaptive_smooth( self, kernel, counts):
        """
        : adaptive_smooth( AdaptiveKernel(), min_counts )

            Perform an adaptive smooth.  Increase the kernel size
            until the min_counts value is reached.  This 
            uses the dmimgadapt tool.
            
            The available smoothing kernels are:

            AdaptGaussian( [minrad [,maxrad [,numrad [,radscale]]]] )
            AdaptTophat( [minrad [,maxrad [,numrad [,radscale]]]] )
            AdaptBoxcar( [minrad [,maxrad [,numrad [,radscale]]]] )
            AdaptCone( [minrad [,maxrad [,numrad [,radscale]]]] )
            AdaptPyramid( [minrad [,maxrad [,numrad [,radscale]]]] )
            AdaptSphere( [minrad [,maxrad [,numrad [,radscale]]]] )
            AdatpQuad( [minrad [,maxrad [,numrad [,radscale]]]] )
            AdaptExp( [minrad [,maxrad [,numrad [,radscale]]]] )

            The defaults are
                minrad = 1
                maxrad = 100
                numrad = 100
                radscale = "linear"

            Examples:
            
            >>> img.adaptive_smooth( AdaptCone(), 10) # uses 100 cones 
                                                    w/ radii from 1 to 100
            >>> img.adaptive_smooth( AdaptGaussian(1,20), 15 ) # uses 100 
                                              Gaussians between 1 and 20
            >>> img.adaptive_smooth( AdaptSphere(numrad=20), 0.01 ) # 20 
                                              spheres from 1 to 100

        """
        dmimgadapt = make_tool("dmimgadapt")
        return _run_cmd( dmimgadapt, self, counts=counts, function=kernel.function, 
            minrad=kernel.minrad, maxrad=kernel.maxrad, numrad=kernel.numrad,
            radscal=kernel.radscale )

    def agaus( self, counts, **kwargs):
        k = AdaptGaussian( **kwargs )
        return self._adaptive_smooth( k, counts )
        
    def atophat( self, counts, **kwargs):
        k = AdaptTophat( **kwargs )
        return self._adaptive_smooth( k, counts )
        
    def abox( self, counts, **kwargs):
        k = AdaptBoxcar( **kwargs )
        return self._adaptive_smooth( k, counts )
        
    def acone( self, counts, **kwargs):
        k = AdaptCone( **kwargs )
        return self._adaptive_smooth( k, counts )
        
    def apyramid( self, counts, **kwargs):
        k = AdaptPyramid( **kwargs )
        return self._adaptive_smooth( k, counts )
        
    def asphere( self, counts, **kwargs):
        k = AdaptSphere( **kwargs )
        return self._adaptive_smooth( k, counts )
        
    def aquad( self, counts, **kwargs):
        k = AdaptQuad( **kwargs )
        return self._adaptive_smooth( k, counts )
        
    def aexp( self, counts, **kwargs):
        k = AdaptExp( **kwargs )
        return self._adaptive_smooth( k, counts )
 
    def thresh(self, lo=0, hi=None, val=0):
        """
        Lo and hi pixel threshold
        """
        dmimgthresh = make_tool("dmimgthresh")
        if None == val:
            val = "INDEF" 
        cut = "" if None == lo else str(lo)
        if None != hi:
            cut += ":{}".format(hi)
        return _run_cmd( dmimgthresh, self, cut=cut, value=val )

    
    def filter( self, region, nullval=0, coord="sky" ):
        """
        Spatial filter
        """
        dmcopy = make_tool("dmcopy")
        if None == nullval:
            nullval = "NaN"
        
        vfspec="[{}={}][opt full,null={}]".format(coord,region, nullval )
        return _run_cmd( dmcopy, self, vfspec=vfspec )
    

    def shift( self, dx, dy, **kwargs ):
        dmregrid2 = make_tool("dmregrid2")
        return _run_cmd( dmregrid2, self, xoffset=-1.0*float(dx), yoffset=-1.0*float(dy), theta=0, rotx=0, roty=0, **kwargs )

        
    def rotate( self, angle, center=None, **kwargs ):
        dmregrid2 = make_tool("dmregrid2")
        if None == center:
            center = [ x/2.0 for x in self.get_shape() ]            
        return _run_cmd( dmregrid2, self, theta=angle, rotx=center[0], roty=center[1], **kwargs )

        
    def scale( self, sx, sy=None, **kwargs ):
        dmregrid2 = make_tool("dmregrid2")
        if None == sy:
            sy = sx
        center = [ x/2.0 for x in self.get_shape() ]            
        return _run_cmd( dmregrid2, self, xscale=sx, yscale=sy, rotx=center[0], roty=center[1], **kwargs )

    def blob( self, threshold=None, srconly=True):
        dmimgblob=make_tool("dmimgblob")
        if None == threshold:
            from sys import float_info
            threshold = float_info.epsilon
        return _run_cmd( dmimgblob, self, threshold=threshold, srconly=srconly )
                
    def adaptive_bin( self, snr, **kwargs ):
        dmnautilus = make_tool("dmnautilus")
        return _run_cmd( dmnautilus, self, snr=snr, **kwargs )

    def distance( self, edgevalue):
        dmimgdist = make_tool("dmimgdist")
        return _run_cmd( dmimgdist, self, tolerance=edgevalue )
        
    def fill(self, src, bkg=None, method="poisson"):
        dmfilth = make_tool("dmfilth")
        
        method = method.lower()
        
        if method not in ["poly", "dist", "global", "poisson", "bilint"]:
            raise ValueError("Invalid fill method: {}".format(method))
        
        if method in ["poly", "dist", "poisson", "bilint"]:
            if None == bkg:
                raise ValueError("Background region must be supplied with this method")
        
        return( _run_cmd( dmfilth, self, method=method.upper(), srclist=src.__str__(), bkglist=bkg.__str__()))

    def powerspectrum( self ):
        """
        """
        apowerspectrum = make_tool("apowerspectrum")
        return( _run_cmd( apowerspectrum, self, inp="infilereal"))

    def correlate(self, other=None):
        """
        """
        acrosscorr = make_tool("acrosscorr") 
        if None == other:
            return( _run_cmd( acrosscorr, self, inp="infile1"))

        with serialize_temp_crate( other) as tmppsf:
            oo = _run_cmd( acrosscorr, self, inp="infile1", infile2=tmppsf.name)
        return(oo) 

    def convolve( self, psf, **kwargs ):
        aconvolve = make_tool("aconvolve")
        with serialize_temp_crate( psf ) as tmppsf:
            oo = _run_cmd( aconvolve, self, kernelspec="file:{}".format(tmppsf.name), **kwargs)
        return(oo) 
    
    def deconvolve(self, psf, numiter=100, method="lucy", xc="INDEF", yc="INDEF"):
        """
        """
        if method not in [ 'lucy' ]:
            raise NotImplementedError("Unsupported deconvolution algorithm")

        arestore = make_tool("arestore")
        with serialize_temp_crate( psf) as tmppsf:
            oo = _run_cmd( arestore, self, psffile=tmppsf.name, numiter=numiter, psf_x_center=str(xc), psf_y_center=str(yc))        
        return(oo) 
    
    def csmooth( self, sigmin=3, sigmax=5, kernel="gauss", sclmin="INDEF", sclmax=20):
        csmooth = make_tool("csmooth")
        return( _run_cmd( csmooth, self, sclmap="", outsigfile="", outsclfile="", conmeth="fft", conkerneltype=kernel, sigmin=sigmin, sigmax=sigmax, sclmin=str(sclmin), sclmax=str(sclmax),  sclmode="compute") )


    def maskbin( self, maskfile):
        dmmaskbin = make_tool("dmmaskbin")
        with serialize_temp_crate( maskfile) as mask:
            oo = _run_cmd( dmmaskbin, self, maskfile=mask.name)        
        return(oo) 


    def match( self, tobe, **kwargs ):
        ri = make_tool("reproject_image")
        with serialize_temp_crate( tobe) as matchfile:
            oo = _run_cmd( ri, self, matchfile=matchfile.name, **kwargs)        
        return(oo) 
        
    # -------------- Image returns a Region ------------------------
    def contour( self, levels):
        dmcontour = make_tool("dmcontour")
        return(_get_reg( dmcontour, self, levels=levels ))
                
    def ellipse( self, fraction, **kwargs):
        dmellipse = make_tool("dmellipse")
        return(_get_reg( dmellipse, self, fraction=fraction, **kwargs ))

    def hull( self, tolerance=0):
        dmimghull = make_tool("dmimghull")
        return(_get_reg( dmimghull, self, tolerance=tolerance ))

    def lasso( self, xpos, ypos, low_value=0, high_value="INDEF", **kwargs):
        dmimglasso = make_tool("dmimglasso")
        return( _get_reg( dmimglasso, self, xpos=xpos, ypos=ypos, low_value=low_value, high_value=high_value, **kwargs))

    def get_src_region( self):
        raise NotImplementedError("This command is not currently supported")


    # ------------------- Image returns a Table -------------------------
    def histogram(self, grid="1"):
        # Todo: make grid more pythonic (list of lo/hi pars, min:max:bin, grid() etc
        dmimghist = make_tool("dmimghist")
        return( _get_table( dmimghist, self, hist=grid) )

        
    def project(self, axis):
        # Todo: make separate xproject and yproject routines 
        dmimgproject = make_tool("dmimgproject")
        if axis not in ['x','y']:
            raise ValueError("axis must be either 'x' or 'y'")
        return( _get_table( dmimgproject, self, axis=axis ))

    def radial( self, xpos, ypos, inner=0, outer=1000, step=10, coord="sky", **kwargs ):
        """
        Create a radial profile from a fixed grid inner:outer:step
        """
        dmextract = make_tool("dmextract")
        
        vfspec="[bin {}=annulus({},{},{}:{}:{})]".format(coord,xpos, ypos, inner, outer, step )
        return _get_table( dmextract, self, vfspec=vfspec, opt="generic", **kwargs )


    def extract( self, reg, coord="sky", **kwargs ):
        """
        Extract counts in region
        """

        # TODO: bkg requires infile name, more work for Crateify routine 

        dmextract = make_tool("dmextract")

        vfspec="[bin {}={}]".format(coord,reg )
        return _get_table( dmextract, self, vfspec=vfspec, opt="generic", **kwargs )


    def celldetect(self):
        raise NotImplementedError("This command is not currently supported")

    def vtpdetect(self):
        raise NotImplementedError("This command is not currently supported")
        
    def wavdetect(self):
        raise NotImplementedError("This command is not currently supported")
        
    def dmimgpick(self):
        raise NotImplementedError("This command is not currently supported")

    # ------------- Returns a par/dict ------------------
        
    @staticmethod    
    def _get_dmcoords_values(dmcoords):
        """
        """
        retval = {}
        for vv in [ 'chip_id', 'chipx', 'chipy', 'tdetx', 'tdety', 'detx', 'dety', 'x', 'y', 'logicalx', 'logicaly', 'ra', 'dec', 'theta', 'phi']:
            retval[vv] = float(getattr( dmcoords, vv )  )
        return retval
            

    def world(self, ra, dec):
        # TODO: Convert ra/dec to deg
        dmcoords = make_tool("dmcoords")
        with serialize_temp_crate( self) as tmpimg:
            dmcoords.punlearn()
            dmcoords( tmpimg.name, op="cel", celfmt="deg", ra=ra, dec=dec)
        return self._get_dmcoords_values( dmcoords )
    
    def physical(self, xx, yy):
        dmcoords = make_tool("dmcoords")
        with serialize_temp_crate( self) as tmpimg:
            dmcoords.punlearn()
            dmcoords( tmpimg.name, op="sky", celfmt="deg", x=xx, y=yy)
        return self._get_dmcoords_values( dmcoords )

    def logical(self, lx, ly):
        dmcoords = make_tool("dmcoords")
        with serialize_temp_crate( self) as tmpimg:
            dmcoords.punlearn()
            dmcoords( tmpimg.name, op="logical", celfmt="deg", logicalx=lx, logicaly=ly)
        return self._get_dmcoords_values( dmcoords )

    def det(self, detx, dety):
        dmcoords = make_tool("dmcoords")
        with serialize_temp_crate( self) as tmpimg:
            dmcoords.punlearn()
            dmcoords( tmpimg.name, op="det", celfmt="deg", detx=detx, dety=dety)
        return self._get_dmcoords_values( dmcoords )

    def msc(self, theta, phi):
        dmcoords = make_tool("dmcoords")
        with serialize_temp_crate( self) as tmpimg:
            dmcoords.punlearn()
            dmcoords( tmpimg.name, op="msc", celfmt="deg", theta=theta, phi=phi)
        return self._get_dmcoords_values( dmcoords )

    def chip(self, chipid, chipx, chipy):
        dmcoords = make_tool("dmcoords")
        with serialize_temp_crate( self) as tmpimg:
            dmcoords.punlearn()
            dmcoords( tmpimg.name, op="chip", celfmt="deg", chip_id=chipid, chipx=chipx, chipy=chipy)
        return self._get_dmcoords_values( dmcoords )

    def stats(self):
        """
        
        """
        dmstat = make_tool("dmstat")
        with serialize_temp_crate( self) as tmpimg:
            dmstat.punlearn()
            dmstat( tmpimg.name, centroid=False, sigma=True, median=True)
            retval = {
              'min' : float(dmstat.out_min),
              'max' : float(dmstat.out_max),
              'min_loc' : [float(x) for x in dmstat.out_min_loc.split(",") ],
              'max_loc' : [float(x) for x in dmstat.out_max_loc.split(",") ],
              'mean' : float(dmstat.out_mean),
              'median' : float(dmstat.out_median),
              'stddev' : float(dmstat.out_sigma),
              'sum' : float(dmstat.out_sum),
              'pixelarea' : int( dmstat.out_good ),
              'numnull' : int( dmstat.out_null )
              }
            dmstat( tmpimg.name, centroid=True, sigma=False, median=False)
            retval["centroid"] = [float(x) for x in dmstat.out_cntrd_phys.split(",") ]
            retval["lcentroid"] = [float(x) for x in dmstat.out_cntrd_log.split(",") ]

        return retval
        
    def moments(self):
        """
        """
        mom = make_tool("imgmoment")
        with serialize_temp_crate( self) as tmpimg:
            mom.punlearn()
            mom( tmpimg.name )
        retval = {}
        retval['centroid'] = [ mom.x_mu, mom.y_mu ]
        retval['angle' ] = mom.phi
        retval['major'] = mom.xsig
        retval['minor'] = mom.ysig
        retval['matrix'] = np.array ( [ [mom.m_0_0, mom.m_0_1, mom.m_0_2],
                                        [mom.m_1_0, mom.m_1_1, mom.m_1_2],
                                        [mom.m_2_0, mom.m_2_1, mom.m_2_2] ] )
        return retval

        
    def extent( self ):
        #srcextent
        raise NotImplementedError("This command is not currently supported")

    # ------------- sherpa --------------------
    def fit(self):
        #sherpa
        raise NotImplementedError("Sherpa fitting not yet available")
    
    # ----------------- ds9 -----------------------
    def ds9(self):
        try:
            from pyds9 import DS9
            d = DS9()
            
            with serialize_temp_crate( self) as img:
                d.set("file {}".format( img.name ))
            
            return d
            
        except ImportError:
            print "Please install pyds9 to enable this feature"
            return None
        except:
            raise




    
    def _test_(self):

        if True:
            for mth in ["min","max","mean","median","mode","mid","sigma","extreme",
                "locheq","kuwahara","unsharp","range","variance","nmode",
                "q25","q33","q67","q75","mcv","sum","rclip","peak","valley",
                "count","olympic","pmean","mu3","mu4","jitter","rms","nslope",
                "sig3mean","sig3median" ]:
                    getattr( self, mth)(Box(5)).write( mth+".out", clobber=True )
                
        if True:
            for mth in ["cos", "sin", "tan", "acos", "asin", "atan", "cosh", 
                "sinh", "tanh", "exp", "log", "ln", "sqrt", "fabs",
                "asinh", "acosh", "atanh",
                "byte", "short", "long", "ushort", "ulong", "float", "double", "int" ]:
                    getattr( self, mth)().write(mth+".out", clobber=True )
                
        if True:
            for k in [ Gaussian(3), Gaussian(3,5),
                Boxcar(5), Boxcar(9,5),
                Tophat(4), Tophat( 10,2),
                Mexhat(2,norm="none"),Mexhat(3,3,norm="none"),
                Exp(-1), Exp(1,-1),
                Power(-1),Power(-1,-2),
                Beta(-1),
                Sinc(2),
                Cone(5),
                Pyramid(7),
                Hanning(10),
                RaisedCosine(0.4, 5)]:
                    n=k.__repr__()
                    n=n.split(":")[1].replace(",","_").replace("(","_").replace(")","")
                    self._smooth(k).write( n+".out", clobber=True)

        if True:
            for k in [AdaptGaussian(maxrad=10), AdaptBoxcar(maxrad=10), AdaptCone(maxrad=10), 
                AdaptExp(maxrad=10), AdaptPyramid(maxrad=10), AdaptQuad(maxrad=10), 
                AdaptSphere(maxrad=10), AdaptTophat(maxrad=10)]:
                    self._adaptive_smooth( k, 10 ).write( k.function+".out", clobber=True)

        if True:
            self.thresh().write("thresh0.out", clobber=True)
            self.thresh(1).write("thresh1.out", clobber=True)
            self.thresh(None,3, val=3).write("thresh_max3.out",clobber=True)
            self._smooth(Gaussian(3)).thresh(0.1,3,val=None).write("thresh_range.out",clobber=True)

        if True:
            self.filter(Field(), coord="(#1,#2)").write("field.out", clobber=True)
            self.filter("field()", coord="(#1,#2)").write("field_str.out", clobber=True)
            self.filter(Circle(200,(300,300)), coord="(#1,#2)").write("circle.out", clobber=True)
            self.filter(Annulus(100,200,(300,300)), coord="(#1,#2)").write("annulus.out", clobber=True)
            self.filter(Box(100,center=(200,200)), coord="(#1,#2)").write("box100.out", clobber=True)
            self.filter(Box(100,200,(300,400)), coord="(#1,#2)").write("box100_200.out", clobber=True)
            self.filter(Box(100,center=(300,300),angle=45), coord="(#1,#2)").write("box100_rot45.out",clobber=True)
            self.filter(Ellipse(100,200,(300,400)), coord="(#1,#2)").write("ellipse100_200.out", clobber=True)
            self.filter(Point(center=(546,333)), coord="(#1,#2)").write("point.out",clobber=True)
            self.filter(Pie(30,100,-45,45,center=(300,300)), coord="(#1,#2)").write("pie.out",clobber=True)
            self.filter(Sector(10,70), coord="(#1,#2)").write("sector.out",clobber=True)
            self.filter(Rectangle( 100,100,300,400), coord="(#1,#2)").write("rectangle.out",clobber=True)
            self.filter(Polygon( (197,410),(370,410),(496,242),(366,280),(306,150),(197,237)), coord="(#1,#2)").write("polygon.out", clobber=True)

            two_boxes = Box(10,center=(200,200)) + Box(30,10,center=(400,400))
            self.filter(two_boxes, coord="(#1,#2)").write("2boxes.out",clobber=True)

            pacmac = Circle(200,center=(300,300)) - Circle(100, center=(450,300))
            self.filter(pacmac, coord="(#1,#2)").write("pacmac.out",clobber=True)

            el1 = Ellipse( 50, 200, center=(300,300), angle=45)
            el2 = Ellipse( 50, 200, center=(300,300), angle= 135)        
            elxor = el1-el2 + el2-el1
            eland = el1 * el2
            self.filter(elxor, coord="(#1,#2)").write("xor.out", clobber=True)
            self.filter(eland, coord="(#1,#2)").write("and.out", clobber=True)
        
        if True:
            self.fill( circle(4210,4055,37.922803), annulus(4210,4054,52.038447,74.69442)).write("fill_poisson.out", clobber=True)
            self.fill( circle(4210,4055,37.922803), annulus(4210,4054,52.038447,74.69442), method="dist").write("fill_dist.out", clobber=True)
            self.fill( circle(4210,4055,37.922803), method="global").write("fill_global.out", clobber=True)
            self.fill( circle(4210,4055,37.922803), annulus(4210,4054,52.038447,74.69442), method="bilint").write("fill_bilint.out", clobber=True)
            self.fill( circle(4210,4055,37.922803), circle(4210,4054,74.69442), method="poly").write("fill_poly.out", clobber=True)
            

        if True:
            (self+0).write("add.out", clobber=True)
            (self-0).write("sub.out", clobber=True)
            (self*1).write("mul.out", clobber=True)
            (self/1).write("div.out", clobber=True)

            (0+self+0).write("add2.out", clobber=True)
            (0-self-0).write("sub2.out", clobber=True)
            (1*self*1).write("mul2.out", clobber=True)
            (1.0/self/1).write("div2.out", clobber=True)

            (self+self).write("c_add.out", clobber=True)
            (self-self).write("c_sub.out", clobber=True)
            (self*self).write("c_mul.out", clobber=True)
            (self/self).write("c_div.out", clobber=True)
        
        
        if True:

            self.power(2).power(0.5).write("power_square_sqrt.out", clobber=True)
            self.blob().write("blob_def.out", clobber=True)
            self._smooth(Gaussian(3)).blob(0.05).write("blob_sm0p05.out", clobber=True)
            self._smooth(Gaussian(3)).distance(0.05).write("dist_sm0p05.out", clobber=True)
            self.adaptive_bin(3).write("abin_3.out",clobber=True)
            self._smooth(Boxcar(5)).adaptive_bin(4).write("abin_bsm_4.out",clobber=True)

            self.scale(2).write("scale_2.out", clobber=True)
            self.scale(0.5).write("scale_0.5.out",clobber=True)
            self.scale(2,4).write("scale_2x4.out",clobber=True)
            self.rotate(45).write("rotate_45.out",clobber=True)
            self.rotate(-45,center=(0,0)).write("rotate_about00.out",clobber=True)
            
            self.shift( 200,200).write("shift_200x200.out",clobber=True)
            self.shift(-200,200).write("shift_m200x200.out",clobber=True)


            pass
            
        if True:
            self.gaus(3).powerspectrum().write("powerspectrum.out",clobber=True)
            self.correlate().write("auto_correlate.out", clobber=True )
            self.correlate( self.tophat(3) ).write("correlate.out", clobber=True)
            
            self.gaus(3).deconvolve( self ).write("deconvolve.out", clobber=True)
            self.csmooth( 3,6,sclmax=10).write("csmooth.out", clobber=True)
            
        if True:
            self.gaus(3).contour("1").write( "contour_reg.out")
            self.box(10).ellipse( 0.5).write( "ellipse_reg.out")
            self.box(10).ellipse(0.5, shape="rotbox").write("rotbox_reg.out")
            self.filter(circle(4096,4096,20)).hull().write("hull_reg.out")
            
        if True:
            self.histogram().write("histogram_tab.out", clobber=True)
            self.histogram("0:100:2").write("histogram_0_100_2_tab.out", clobber=True)
            self.project("x").write("project_x_tab.out", clobber=True )
            self.project("y").write("project_y_tab.out", clobber=True )
        


@Crateify(CIAOImage)
def _run_cmd( mycmd, infile, outfile, vfspec="", inp="infile", **kwargs):

    mycmd.punlearn()
    setattr( mycmd, inp, infile+vfspec )
    setattr( mycmd, "outfile", outfile )
    setattr( mycmd, "clobber", True )
        
    for k in kwargs:
        setattr( mycmd, k, kwargs[k] )

    print mycmd
    x = mycmd()
    if x : print x 


@Crateify(Region)
def _get_reg( mycmd, infile, outfile, vfspec="", inp="infile", **kwargs):

    mycmd.punlearn()
    setattr( mycmd, inp, infile+vfspec )
    setattr( mycmd, "outfile", outfile )
    setattr( mycmd, "clobber", True )
        
    for k in kwargs:
        setattr( mycmd, k, kwargs[k] )

    print mycmd
    x = mycmd()
    if x : print x 


@Crateify(HistoryTABLECrate)
def _get_table( mycmd, infile, outfile, vfspec="", inp="infile", **kwargs):

    mycmd.punlearn()
    setattr( mycmd, inp, infile+vfspec )
    setattr( mycmd, "outfile", outfile )
    setattr( mycmd, "clobber", True )
        
    for k in kwargs:
        setattr( mycmd, k, kwargs[k] )

    print mycmd
    x = mycmd()
    if x : print x 



import contextlib
@contextlib.contextmanager
def serialize_temp_crate( data2save ):
    """
    
    """
    import tempfile as tempfile 
    nn = tempfile.NamedTemporaryFile( )
    try:
        data2save.write(nn.name, clobber=True )
    except:
        from crates_contrib.utils import make_image_crate
        make_image_crate( data2save ).write( nn.name, clobber=True )
    yield nn


        





#img = CIAOImage("img.fits")
#img._test_()

