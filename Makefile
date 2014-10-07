
ROOT = /data/lenin2/Scripts/MyStuff/ciao46/ciao-4.6/contrib
DEST = lib/python2.7/site-packages/ciao_contrib/cipt
DEV  = /data/da/Docs/scripts/dev

CP_F = /bin/cp -f

RDIR = $(ROOT)/$(DEST)
DDIR = $(DEV)/$(DEST)

PYCODE = cipt.py crateify.py filter_mask.py history_crate.py smooth_kernels.py __init__.py all.py


all: $(PYCODE)

install: all
	mkdir -p $(RDIR)
	$(CP_F) $(PYCODE) $(RDIR)/

install-dev: all
	mkdir -p $(DDIR)
	$(CP_F) $(PYCODE) $(DDIR)/
