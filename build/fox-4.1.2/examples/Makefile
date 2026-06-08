#!/usr/bin/make

include ../arch.make

default: $(BUILD_TARGETS)

INCFLAGS=`../FoX-config --fcflags`

EXAMPLES=dom_example_2$(EXEEXT) dom_example_3$(EXEEXT) sax_example$(EXEEXT) sax_example_2$(EXEEXT) wcml_example$(EXEEXT) wxml_example$(EXEEXT) wkml_example$(EXEEXT)

dom_lib: sax_lib wxml_lib dom_example_2$(EXEEXT) dom_example_3$(EXEEXT) 
sax_lib: sax_example$(EXEEXT) sax_example_2$(EXEEXT)
wcml_lib: wxml_lib wcml_example$(EXEEXT)
wkml_lib: wxml_lib wkml_example$(EXEEXT)
wxml_lib: wxml_example$(EXEEXT)

# Note that we cannot use $< below since IBM make doesn't understand it.

dom_canonicalize.ns.no.exe: dom_canonicalize.ns.no.$(OBJEXT) ../objs/lib/libFoX_dom.$(LIBEXT) ../objs/lib/libFoX_sax.$(LIBEXT) ../objs/lib/libFoX_wxml.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ dom_canonicalize.ns.no.$(OBJEXT) $$(../FoX-config --libs --dom)
dom_canonicalize.ns.yes.exe: dom_canonicalize.ns.yes.$(OBJEXT) ../objs/lib/libFoX_dom.$(LIBEXT) ../objs/lib/libFoX_sax.$(LIBEXT) ../objs/lib/libFoX_wxml.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ dom_canonicalize.ns.yes.$(OBJEXT) $$(../FoX-config --libs --dom)
dom_configuration.exe: dom_configuration.$(OBJEXT) ../objs/lib/libFoX_dom.$(LIBEXT) ../objs/lib/libFoX_sax.$(LIBEXT) ../objs/lib/libFoX_wxml.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT) 
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ dom_configuration.$(OBJEXT) $$(../FoX-config --libs --dom)
dom_example_2$(EXEEXT): dom_example_2.$(OBJEXT) ../objs/lib/libFoX_dom.$(LIBEXT) ../objs/lib/libFoX_sax.$(LIBEXT) ../objs/lib/libFoX_wxml.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ dom_example_2.$(OBJEXT) $$(../FoX-config --libs --dom)
dom_example_3$(EXEEXT): dom_example_3.$(OBJEXT) ../objs/lib/libFoX_dom.$(LIBEXT) ../objs/lib/libFoX_sax.$(LIBEXT) ../objs/lib/libFoX_wxml.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ dom_example_3.$(OBJEXT) $$(../FoX-config --libs --dom)
#
sax_example$(EXEEXT): sax_example.$(OBJEXT) ../objs/lib/libFoX_sax.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ sax_example.$(OBJEXT) $$(../FoX-config --libs --sax)
sax_example_2$(EXEEXT): sax_example_2.$(OBJEXT) ../objs/lib/libFoX_sax.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ sax_example_2.$(OBJEXT) $$(../FoX-config --libs --sax)
#
wxml_example$(EXEEXT): wxml_example.$(OBJEXT) ../objs/lib/libFoX_wxml.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ wxml_example.$(OBJEXT) $$(../FoX-config --libs --wxml)
#
wcml_example$(EXEEXT): wcml_example.$(OBJEXT)  ../objs/lib/libFoX_wcml.$(LIBEXT) ../objs/lib/libFoX_wxml.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ wcml_example.$(OBJEXT) $$(../FoX-config --libs --wcml)
#

wkml_example$(EXEEXT): wkml_example.$(OBJEXT)  ../objs/lib/libFoX_wkml.$(LIBEXT) ../objs/lib/libFoX_wxml.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ wkml_example.$(OBJEXT) $$(../FoX-config --libs --wkml)


clean:
	rm -f *.$(OBJEXT) *.$(MOD_EXT) $(EXAMPLES) wkml_example.kml
