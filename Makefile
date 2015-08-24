
ROOT = /data/lenin2/Scripts/MyStuff/ciao46/ciao-4.6/contrib
DEST = lib/python2.7/site-packages/ciao_contrib/cipt
DEST2 = lib/python2.7/site-packages/chips_contrib
DEV  = /data/da/Docs/scripts/dev


CP_F = /bin/cp -fv

RDIR = $(ROOT)/$(DEST)
R2DIR = $(ROOT)/$(DEST2)
DDIR = $(DEV)/$(DEST)
D2DIR = $(DEV)/$(DEST2)


PYCODE = cipt.py crateify.py enhanced_region.py history_crate.py smooth_kernels.py __init__.py all.py
PLTCODE = plot_shapes.py


all: $(PYCODE)

install: all
	@mkdir -p $(RDIR)
	@$(CP_F) $(PYCODE) $(RDIR)/
	@mkdir -p $(R2DIR)
	@$(CP_F) $(PLTCODE) $(R2DIR)/

install-dev: all
	@mkdir -p $(DDIR)
	@$(CP_F) $(PYCODE) $(DDIR)/
	@mkdir -p $(D2DIR)
	@$(CP_F) $(PLTCODE) $(D2DIR)/

