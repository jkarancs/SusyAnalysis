include Makefile.arch


#------------------------------------------------------------------------------
B2GPLOTTERO      = B2GPlotter.$(ObjSuf)
B2GPLOTTERS      = B2GPlotter.$(SrcSuf)
B2GPLOTTER       = B2GPlotter$(ExeSuf)

OBJS          = $(B2GPLOTTERO)
PROGRAMS      = $(B2GPLOTTER)

#------------------------------------------------------------------------------
B2GANALYZERO      = B2GAnalyzer.$(ObjSuf)
B2GANALYZERS      = B2GAnalyzer.$(SrcSuf)
B2GANALYZER       = B2GAnalyzer$(ExeSuf)

OBJS          += $(B2GANALYZERO)
PROGRAMS      += $(B2GANALYZER)

#------------------------------------------------------------------------------
SUSYHISTOMAKERO      = SusyHistoMaker.$(ObjSuf)
SUSYHISTOMAKERS      = SusyHistoMaker.$(SrcSuf)
SUSYHISTOMAKER       = SusyHistoMaker$(ExeSuf)

OBJS          += $(SUSYHISTOMAKERO)
PROGRAMS      += $(SUSYHISTOMAKER)

#------------------------------------------------------------------------------
DRAWEVENTO      = DrawEvent.$(ObjSuf)
DRAWEVENTS      = DrawEvent.$(SrcSuf)
DRAWEVENT       = DrawEvent$(ExeSuf)

OBJS          += $(DRAWEVENTO)
PROGRAMS      += $(DRAWEVENT)

#------------------------------------------------------------------------------


.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRAMS)

$(B2GPLOTTER):     $(B2GPLOTTERO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"
$(B2GANALYZER):     $(B2GANALYZERO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"
$(SUSYHISTOMAKER):     $(SUSYHISTOMAKERO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"
$(DRAWEVENT):     $(DRAWEVENTO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

clean:
		@rm -f $(OBJS) core


distclean:      clean
		@rm -f $(PROGRAMS) *Dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d *.log .def so_locations
		@rm -rf cxx_repository

.SUFFIXES: .$(SrcSuf)

###

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<

