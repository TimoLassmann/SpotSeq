EXTRA_DIST = autogen.sh
# Build the library in the hand subdirectory first.
LIB_TLDEVELDIR = @LIB_TLDEVELDIR@ 

# top_builddir = @top_builddir@

SUBDIRS = tldevel


MYLIBDIRS = ./tldevel/libthrpool.a \
	./tldevel/libtldevel.a


MYINCDIRS = -I${top_builddir}/${LIB_TLDEVELDIR} \
            -I${top_srcdir}/${LIB_TLDEVELDIR}

LIBS = -lm -lpthread 


AM_CPPFLAGS=${MYINCDIRS} 


noinst_LIBRARIES = libihmm.a

# AUTOMAKE_OPTIONS = subdir-objects
libihmm_a_SOURCES = \
ihmm.h \
ihmm.c \
distributions.h \
distributions.c \
randomkit.h \
randomkit.c

libihmm_a_LIBADD  =  ${MYLIBDIRS}

