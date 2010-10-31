EXTRA_DIST = libplump_wrap.cxx libplump.py

if PYTHONBINDING
all-local: libplump_wrap.cxx libplump.py
	$(PYTHON) setup.py build
	$(PYTHON) setup.py build_ext --inplace

clean-local:
	$(PYTHON) setup.py clean --all
	rm libplump_wrap.cxx libplump.py _libplump.so

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(DESTDIR)$(prefix)

if USESWIG
libplump_wrap.cxx libplump.py: ../libplump.i
	$(SWIG) -c++ -python -I$(top_srcdir)/src -outdir . -o libplump_wrap.cxx  $(srcdir)/../libplump.i
endif

endif
