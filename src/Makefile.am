# Build the library in the hand subdirectory first.
LIB_TLDEVELDIR = @LIB_TLDEVELDIR@ 

# top_builddir = @top_builddir@


MYLIBDIRS = ../tldevel/libthrpool.a \
	../tldevel/libtldevel.a


MYINCDIRS = -I${top_builddir}/${LIB_TLDEVELDIR} \
            -I${top_srcdir}/${LIB_TLDEVELDIR}

LIBS = -lm -lpthread 


AM_CPPFLAGS=${MYINCDIRS} 



bin_PROGRAMS = build_ihmm sim_seq

sim_seq_SOURCES = \
outdir.h \
outdir.c \
matrix_io.h \
matrix_io.c \
sim_seq.c

sim_seq_LDADD =  ${MYLIBDIRS} 


build_ihmm_SOURCES = \
model_nuc_seq.c \
matrix_io.h \
matrix_io.c \
ihmm_io.h \
ihmm_io.c \
outdir.h \
outdir.c \
make_dot_file.h \
make_dot_file.c

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


check_PROGRAMS = libihmm_ITEST fast_transition_ITEST

libihmm_ITEST_SOURCES = $(libihmm_a_SOURCES)
libihmm_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITEST
libihmm_ITEST_LDADD = $(libihmm_a_LIBADD)

fast_transition_ITEST_SOURCES = fast_hmm_param.c fast_hmm_param.h distributions.h \
distributions.c \
randomkit.h \
randomkit.c

fast_transition_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITEST
fast_transition_ITEST_LDADD = $(MYLIBDIRS)

