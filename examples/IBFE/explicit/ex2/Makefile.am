## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
EXAMPLE_DRIVER = ex_main.cpp
GTEST_DRIVER   = test_main.cpp
CMB_EX_DRIVER  = cmb_ex_main.cpp
EXTRA_DIST     = input2d input3d input2d.test input3d.test

EXAMPLES =
GTESTS =
CMB_EXAMPLES =
EXTRA_PROGRAMS =
if LIBMESH_ENABLED
if SAMRAI2D_ENABLED
EXAMPLES += main2d 
GTESTS += test2d
CMB_EXAMPLES += cmb_main2d
EXTRA_PROGRAMS += $(EXAMPLES) $(GTESTS) $(CMB_EXAMPLES)
endif
if SAMRAI3D_ENABLED
EXAMPLES += main3d 
GTESTS += test3d
CMB_EXAMPLES += cmb_main3d
EXTRA_PROGRAMS += $(EXAMPLES) $(GTESTS) $(CMB_EXAMPLES)
endif
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(EXAMPLE_DRIVER)

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = $(EXAMPLE_DRIVER)

cmb_main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
cmb_main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
cmb_main2d_SOURCES = $(CMB_EX_DRIVER)

cmb_main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
cmb_main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
cmb_main3d_SOURCES = $(CMB_EX_DRIVER)

test2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
test2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
test2d_SOURCES = $(GTEST_DRIVER)

test3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
test3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
test3d_SOURCES = $(GTEST_DRIVER)

cmb: $(CMB_EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/cmb_input* $(PWD) ; \
	fi ;

examples: $(EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input* $(PWD) ; \
	fi ;

gtest: $(GTESTS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input* $(PWD) ; \
	fi ;
	./test2d input2d.test

gtest-long:
	make gtest

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input* ; \
	fi ;
