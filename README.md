
# CIAO Image Processing Toolkit

Wraps various CIAO image processin tasks together so that
they can be called within python in such a way that 
the actual file i/o is hidden from the user.

```python
>>> from ciao_contrib.cipt import *
>>> myimg = CIAOImage("foo.fits")
>>> scld = myimg.asinh()
>>> diff = myimg - scld
>>> diff.write("zoo.fits", clobber=True)
```


