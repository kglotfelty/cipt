


__all__ = [ "HistoryIMAGECrate", "HistoryTABLECrate" ]

from pycrates import IMAGECrate, TABLECrate


class HistoryCrate(object):
    """
    File HISTORY transport mechanism
    
    This class provides a way to copy all the HISTORY records from
    an input Crate and copy them to the output.
    
    All the HISTORY records are copied from wherever they appear 
    in the header, but are always written at the end of the
    output header.
    """

    def get_history(self, filename ):
        """
        Use CXCDM to get history records
        """
        from cxcdm import dmBlockOpen, dmBlockClose, dmBlockReadComment, dmBlockMoveToKey, dmBlockGetNoKeys

        self.history = []
        tab = dmBlockOpen( filename )
        for ii in xrange( 1, dmBlockGetNoKeys(tab)+1 ):
            dmBlockMoveToKey(tab,ii)
            while True:
                hist = dmBlockReadComment(tab)
                if None == hist:
                    break
                elif 2 != len(hist):
                    break
                elif not any(hist):
                    break
                else:
                    self.history.append(hist)
        dmBlockClose(tab)
            

    def put_history(self, infile ):
        """
        Use the CXCDM to write history records
        """
        from cxcdm import dmBlockOpen, dmBlockClose, dmBlockWriteComment, dmBlockMoveToKey, dmBlockGetNoKeys
        tab = dmBlockOpen( infile, update=True)
        dmBlockMoveToKey( tab, dmBlockGetNoKeys(tab))
        for t,v in self.history:
            dmBlockWriteComment( tab, t, v )
            
        dmBlockClose(tab)


class HistoryIMAGECrate( IMAGECrate, HistoryCrate):
    """
    A HISTORY aware IMAGECrate
    
    
    """
    def __init__( self, filename, mode="r" ):        
        IMAGECrate.__init__(self,filename, mode=mode)
        self.get_history( filename)

    def write( self, outfile, *args, **kwargs ):
        """
        Overide crates write to include history
        """
        IMAGECrate.write( self, outfile, *args, **kwargs )
        self.put_history( outfile )


class HistoryTABLECrate( TABLECrate, HistoryCrate):
    """
    A HISTORY aware TABLECrate
    
    
    """
    def __init__( self, filename, mode="r" ):        
        TABLECrate.__init__(self,filename, mode=mode)
        self.get_history( filename)

    def write( self, outfile, *args, **kwargs ):
        """
        Overide crates write to include history
        """
        TABLECrate.write( self, outfile, *args, **kwargs )
        self.put_history( outfile )

