# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_GTEST],[
echo
echo "====================================="
echo "Configuring required package Google Test"
echo "====================================="
PACKAGE_SETUP_ENVIRONMENT
CONTRIB_SRCDIR=$1
CONTRIB_BUILDDIR=$2
AC_LANG_PUSH([C++])
AC_ARG_ENABLE([gtest],
  AS_HELP_STRING(--enable-gtest,enable support for the optional Google Test library @<:@default=no@:>@),
                 [case "$enableval" in
                    yes)  GTEST_ENABLED=yes ;;
                    no)   GTEST_ENABLED=no ;;
                    *)    AC_MSG_ERROR(--enable-gtest=$enableval is invalid; choices are "yes" and "no") ;;
                  esac],[GTEST_ENABLED=no])

AM_CONDITIONAL([GTEST_ENABLED],[test "$GTEST_ENABLED" = yes])
AC_ARG_WITH([gtest],
  AS_HELP_STRING(--with-gtest=PATH,location of required Google Test installation),
  [if test ! -d "$withval" ; then
     AC_MSG_ERROR([it is necessary to specify an existing directory when using --with-gtest=PATH])
   fi
   GTEST_DIR=$withval])
if test x$GTEST_DIR != x ; then
  if test -d "${GTEST_DIR}/include" ; then
    GTEST_CPPFLAGS="-I${GTEST_DIR}/include"
    USING_BUNDLED_GTEST=no
  fi
  if test -d "${GTEST_DIR}/lib" ; then
    GTEST_LDFLAGS="-L${GTEST_DIR}/lib"
  fi
fi
if test "$GTEST_ENABLED" = yes; then
  if test "$USING_BUNDLED_GTEST" = no ; then
    echo "first check for a system Google Test library; if not found, revert to bundled Google Test library"
    CPPFLAGS_PREPEND($GTEST_CPPFLAGS)
    $as_unset ac_cv_header_gtest_h
    AC_CHECK_HEADER([${GTEST_DIR}/googletest/include/gtest/gtest.h],[
    HAVE_GTEST=yes
    LDFLAGS_PREPEND($GTEST_LDFLAGS)
    AC_LIB_HAVE_LINKFLAGS([gtest])
    if test "$HAVE_LIBGTEST" = no ; then
    AC_LIB_HAVE_LINKFLAGS([gtest])
    if test "$HAVE_LIBGTEST" = no ; then
      AC_MSG_ERROR([could not find working Google Test])
    fi
    fi
    ])
  else
  AC_MSG_NOTICE([Optional package Google Test is enabled, using bundled])
  USING_BUNDLED_GTEST=yes
  GTEST_DIR=$CONTRIB_SRCDIR
  GTEST_BUILDDIR=$CONTRIB_BUILDDIR
  GTEST_CPPFLAGS="-I${GTEST_DIR}/googletest/include"
  GTEST_LDFLAGS=""
  HAVE_GTEST=yes
  LIBGTEST="$GTEST_BUILDDIR/googletest/lib/libgtest.la"
  AC_DEFINE(HAVE_LIBGTEST, 1, [Define if you have the libgtest library.])
  fi
else
  AC_MSG_NOTICE([Optional package Google Test is DISABLED])
  USING_BUNDLED_GTEST=no
  HAVE_GTEST=no
fi
AM_CONDITIONAL([USING_BUNDLED_GTEST],[test "$USING_BUNDLED_GTEST" = yes])
PACKAGE_CPPFLAGS_PREPEND($GTEST_CPPFLAGS)
PACKAGE_LDFLAGS_PREPEND($GTEST_LDFLAGS)
PACKAGE_CONTRIB_LIBS_PREPEND($LIBGTEST)
AC_LANG_POP([C++])
PACKAGE_RESTORE_ENVIRONMENT
])
