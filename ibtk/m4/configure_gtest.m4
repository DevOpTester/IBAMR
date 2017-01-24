# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_GTEST],[
echo
echo "==================================================="
echo "Configuring optional package Google Test Framework"
echo "==================================================="

PACKAGE_SETUP_ENVIRONMENT

AC_ARG_ENABLE([gtest],
  AS_HELP_STRING(--enable-gtest,enable support for the optional GTEST library @<:@default=no@:>@),
                 [case "$enableval" in
                    yes)  GTEST_ENABLED=yes ;;
                    no)   GTEST_ENABLED=no ;;
                    *)    AC_MSG_ERROR(--enable-gtest=$enableval is invalid; choices are "yes" and "no") ;;
                  esac],[GTEST_ENABLED=no])

AM_CONDITIONAL([GTEST_ENABLED],[test "$GTEST_ENABLED" = yes])

AC_ARG_WITH([gtest],
  AS_HELP_STRING(--with-gtest=PATH,location of optional GTEST installation),
  [if test "$GTEST_ENABLED" = no ; then
     AC_MSG_WARN([--with-gtest is specified, but support for gtest is disabled])
   else
     if test ! -d "$withval" ; then
       AC_MSG_ERROR([it is necessary to specify an existing directory when using --with-gtest=PATH])
     fi
     GTEST_DIR=$withval
   fi])

if test "$GTEST_ENABLED" = yes; then
  if test x$GTEST_DIR != x ; then
    if test -d "${GTEST_DIR}/include" ; then
      GTEST_CPPFLAGS="-I${GTEST_DIR}/include"
    fi
    if test -d "${GTEST_DIR}/lib" ; then
      GTEST_LDFLAGS="-L${GTEST_DIR}/lib"
    fi
  fi

  CPPFLAGS_PREPEND($GTEST_CPPFLAGS)
  AC_CHECK_HEADER([${GTEST_DIR}/include/gtest/gtest.h],,AC_MSG_ERROR([could not find header file gtest.h]))

  LDFLAGS_PREPEND($GTEST_LDFLAGS)
  AC_LIB_HAVE_LINKFLAGS([gtest])
  if test "$HAVE_LIBGTEST" = no ; then
    AC_MSG_ERROR([could not find working libgtest])
  fi

  PACKAGE_CPPFLAGS_PREPEND($GTEST_CPPFLAGS)
  PACKAGE_LDFLAGS_PREPEND($GTEST_LDFLAGS)
  PACKAGE_LIBS_PREPEND("$LIBGTEST")
 
else
  AC_MSG_NOTICE([Optional package Google Test Framework is DISABLED])
fi

PACKAGE_RESTORE_ENVIRONMENT

])
