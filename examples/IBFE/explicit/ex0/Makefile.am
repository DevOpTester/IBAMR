## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.cpp
EXTRA_DIST = input2d input3d

EXTRA_PROGRAMS =
if LIBMESH_ENABLED
if SAMRAI2D_ENABLED
EXTRA_PROGRAMS += main2d
endif
#if SAMRAI3D_ENABLED
#EXTRA_PROGRAMS += main3d
#endif
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(SOURCES)

#main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
#main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
#main3d_SOURCES = $(SOURCES)

## Google Test Application
check_PROGRAMS = test_ex0

test_ex0_SOURCES = testex0.cpp $(SOURCES)
test_ex0_CXXFLAGS = $(main2d_CXXFLAGS)
test_ex0_LDADD = $(main2d_LDADD) -lgtest

TESTS = $(check_PROGRAMS)

examples: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input2d ; \
	fi ;
check: $(check_PROGRAMS)
