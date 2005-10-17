# $Id: configure.ac,v 1.5 2005/10/17 06:34:02 afrolov Exp $
# Process this file with autoconf to produce a configure script.

# Initialization
AC_PREREQ(2.59)
AC_INIT(SCARSE, 0.4-alpha, [Andrei Frolov <afrolov@stanford.edu>])
#AC_PREFIX_DEFAULT([$HOME/local])
AC_CONFIG_LIBOBJ_DIR([lib])
AC_CONFIG_SRCDIR([src/ipb.c])
AC_CONFIG_HEADER([src/config.h])

CFLAGS_SET=${CFLAGS+yes}
CFLAGS_SAVE=$CFLAGS

# Optional features
AC_ARG_ENABLE(debug,
	[AC_HELP_STRING([--enable-debug], [Compile debugging code in])],
	AC_DEFINE(DEBUG, [], [Define if you want debugging code compiled in])
)

# Libtiff location
AC_ARG_WITH(libtiff,
	[AC_HELP_STRING([--with-libtiff=PATH], [specify location of libtiff installation tree])],
	[	test -d $withval/include && LIBTIFF_INC=$withval/include
		test -d $withval/lib && LIBTIFF_LIB=$withval/lib]
)
AC_ARG_WITH(libtiff-inc,
	[AC_HELP_STRING([--with-libtiff-inc=PATH], [specify location of libtiff headers])],
	[test -d $withval && LIBTIFF_INC=$withval]
)
AC_ARG_WITH(libtiff-lib,
	[AC_HELP_STRING([--with-libtiff-lib=PATH], [specify location of libtiff libraries])],
	[test -d $withval && LIBTIFF_LIB=$withval]
)

# Checks for programs.
AC_PROG_CC([icc gcc cc])
AC_PROG_MAKE_SET
AC_CHECK_PROG([AR], [ar], [ar -r])
AC_PROG_RANLIB
AC_PROG_AWK
AC_CHECK_PROG([PERL], [perl], [perl])
AC_CHECK_PROGS([ZIP], ["gzip -9" compress], [true])
AC_CHECK_PROGS([UXP], ["upx --best"], [true])
AC_PROG_INSTALL

# Default flags
CFLAGS='-O3'

if test "$CC" = "icc"; then
	CFLAGS="$CFLAGS -ipo -axPN"
	LDFLAGS="$LDFLAGS -i-static"
	LIBS="$LIBS -lsvml"
fi

case `uname -i` in ( i?86 ) CFLAGS="$CFLAGS -march=i386" ;; esac
case `uname -p` in ( i686 ) CFLAGS="$CFLAGS -mcpu=pentium4" ;; esac

if test "$enable_debug" = "yes"; then
	CFLAGS="$CFLAGS -g -pg -Wall"
	LDFLAGS="$LDFLAGS -g -pg"
fi

test "$CFLAGS_SET" = "yes" && CFLAGS=$CFLAGS_SAVE

test "x$LIBTIFF_INC" = "x" || CFLAGS="$CFLAGS -I$LIBTIFF_INC"
test "x$LIBTIFF_INC" = "x" || CPPFLAGS="$CPPFLAGS -I$LIBTIFF_INC"
test "x$LIBTIFF_LIB" = "x" || LDFLAGS="$LDFLAGS -L$LIBTIFF_LIB"

# Checks for libraries.
AC_SEARCH_LIBS([pow], [svml m])
AC_HAVE_LIBRARY([tiff],, [AC_MSG_ERROR([Can't find tiff libraries, giving up...])])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADER([tiffio.h],, [AC_MSG_ERROR([Can't find tiff headers, giving up...])])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_C_CONST

# Checks for library functions.
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
#AC_FUNC_STRTOD
#AC_FUNC_VPRINTF

# fnmatch replacement
AC_FUNC_FNMATCH
if test $ac_cv_func_fnmatch_works = yes; then
	rm -f lib/fnmatch.h
else
	AC_LIBOBJ(fnmatch)
	AC_CONFIG_LINKS(lib/fnmatch.h:lib/fnmatch.hin)
fi

# sscanf replacement
AC_MSG_CHECKING([for working GNU sscanf])
AC_RUN_IFELSE([AC_LANG_PROGRAM([[#include <stdio.h>]], [[
	int i, n;
	double x;
	char *s = NULL, *r = "Test31";
	char *test = " 1643 175.34e-5 Test31   ";
	
	if ((sscanf(test, "%i %lf %as %n", &i, &x, &s, &n) < 3) ||
		(i != 1643) || (x != 175.34e-5) || strcmp(s, r) || n != 25) return 1;
	else return 0;
]])], [AC_MSG_RESULT([yes])], [
AC_MSG_RESULT([no])
AC_LIBOBJ(sscanf)
], [])

# asprintf replacement
AC_CHECK_FUNC([asprintf],
	[AC_DEFINE(HAVE_ASPRINTF, [1], [Define to 1 if you have the 'asprintf' function])],
	[AC_LIBOBJ(snprintf)]
)

# getline replacement
AC_CHECK_FUNC([getline],
	[AC_DEFINE(HAVE_GETLINE, [1], [Define to 1 if you have the 'getline' function])],
	[AC_LIBOBJ(getline)]
)


#AC_CHECK_FUNCS([floor pow sqrt erf strcasecmp strchr strdup strstr])
AC_CHECK_FUNCS([erf])


AC_CONFIG_FILES([Makefile src/Makefile icclib/Makefile data/Makefile data/profiles/Makefile])
AC_OUTPUT
