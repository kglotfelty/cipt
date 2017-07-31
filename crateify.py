from __future__ import print_function

import pycrates

class Crateify( object ):
    """
    A decorator to make tools seem "crate-aware"
    
    The CIAO tools require files to work, but when in python
    and reading file into crates means that a lot of the CIAO
    functionaltiy is hard to access since it requires 
    writing the crate out, running the tool, then
    reading the file back in.
    
    This decorator can be used to do just this.
    
    Example:
    
      @Crateify()
      def bin_to_image( bincmd, infile, outfile )
            bincmd(infile=infile+"[bin sky=1]", outfile=outfile, clobber=True)
            
      from ciao_contriub.runtool import dmcopy
      tab = read_file("evt.fits")
      img = bin_to_image(dmcopy, tab)
      
    
    """
    def __init__( self, readfn=pycrates.read_file ):
        """        
        """
        self.reader = readfn
    
    def __call__( self, tool ):
        """        
        """
        def write_read_crate( mycmd, mycrate, **kwargs ):
            import tempfile 

            try:
                with tempfile.NamedTemporaryFile() as outfile:
                    with tempfile.NamedTemporaryFile() as infile:
                        mycrate.write(  infile.name, clobber=True )
                        tool( mycmd, infile.name, outfile.name, **kwargs )
                        newcrate = self.reader( outfile.name )
                        return newcrate
            except Exception as E:
                print(E)
                raise E
        return write_read_crate




