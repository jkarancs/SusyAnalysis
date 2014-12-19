include Makefile.arch


##------------------------------------------------------------------------------
#PLOTTER_OLDO      = Plotter_old.$(ObjSuf)
#PLOTTER_OLDS      = Plotter_old.$(SrcSuf)
#PLOTTER_OLD       = Plotter_old$(ExeSuf)
#
#OBJS          = $(PLOTTER_OLDO)
#PROGRAMS      = $(PLOTTER_OLD)

#------------------------------------------------------------------------------
B2GPLOTTERO      = B2GPlotter.$(ObjSuf)
B2GPLOTTERS      = B2GPlotter.$(SrcSuf)
B2GPLOTTER       = B2GPlotter$(ExeSuf)

OBJS          += $(B2GPLOTTERO)
PROGRAMS      += $(B2GPLOTTER)

#------------------------------------------------------------------------------


.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRAMS)

#$(PLOTTER_OLD):     $(PLOTTER_OLDO)
#		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
#		$(MT_EXE)
#		@echo "$@ done"

$(B2GPLOTTER):     $(B2GPLOTTERO)
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

