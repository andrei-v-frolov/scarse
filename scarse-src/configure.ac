# $Id: configure.ac,v 1.3 2005/10/05 06:29:16 afrolov Exp $
# Process this file with autoconf to produce a configure script.

# Initialization
AC_PREREQ(2.59)
AC_INIT(SCARSE, 0.4-dev, [Andrei Frolov <afrolov@stanford.edu>])
#AC_PREFIX_DEFAULT([$HOME/local])
AC_CONFIG_SRCDIR([src/ipb.c])
AC_CONFIG_HEADER([src/config.h])

CFLAGS_SET=${CFLAGS+yes}
CFLAGS_SAVE=$CFLAGS

# Optional features
AC_ARG_ENABLE(debug,
	[AC_HELP_STRING([--enable-debug], [Compile debugging code in])],
	AC_DEFINE(DEBUG, [], [Define if you want debugging code compiled in])
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

# Checks for libraries.
AC_SEARCH_LIBS([pow], [svml m])
AC_CHECK_LIB([tiff], [TIFFOpen])

# Checks for header files.
#AC_HEADER_STDC
#AC_CHECK_HEADERS([fcntl.h float.h stddef.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
#AC_C_CONST
#AC_TYPE_SIZE_T
#AC_STRUCT_TM

# Checks for library functions.
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
#AC_FUNC_MEMCMP
#AC_FUNC_STRTOD
#AC_FUNC_VPRINTF
#AC_REPLACE_FNMATCH
#AC_CHECK_FUNCS([floor pow sqrt erf getline getdelim memset strcasecmp strchr strdup strstr])
AC_CHECK_FUNCS([erf])

AC_CONFIG_FILES([Makefile src/Makefile icclib/Makefile data/Makefile data/profiles/Makefile])
AC_OUTPUT
