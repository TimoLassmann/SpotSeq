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



bin_PROGRAMS = build_ihmm 

build_ihmm_SOURCES = \
model_nuc_seq.c \
matrix_io.h \
matrix_io.c \
ihmm_io.c \
outdir.h\
outdir.c

build_ihmm_LDADD =  ${MYLIBDIRS}  libihmm.a



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


check_PROGRAMS = libihmm_ITEST

libihmm_ITEST_SOURCES = $(libihmm_a_SOURCES)
libihmm_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITEST
libihmm_ITEST_LDADD = $(libihmm_a_LIBADD)
