# Build the library in the hand subdirectory first.
LIB_TLDEVELDIR = @LIB_TLDEVELDIR@ 


CLEANFILES = ihmm_seq_itest_read_test.fa test_model_file.txt test.lfa
# top_builddir = @top_builddir@


MYLIBDIRS = ../tldevel/libthrpool.a \
	../tldevel/libtldevel.a


MYINCDIRS = -I${top_builddir}/${LIB_TLDEVELDIR} \
            -I${top_srcdir}/${LIB_TLDEVELDIR}

LIBS = -lm -lpthread 

RANDOMKIT_FILES = distributions.h \
distributions.c \
randomkit.h \
randomkit.c

AM_CPPFLAGS=${MYINCDIRS} 



ISEQSOURCE = ihmm_seq.h ihmm_seq.c


FASTHMMSOURCE = fast_hmm_param.c fast_hmm_param.h  \
fast_hmm_param_test_functions.h \
fast_hmm_param_test_functions.c

MODELSOURCE = model.h model.c model_io.c global.h  

BEAMSOURCE = beam_sample.h beam_sample.c 



bin_PROGRAMS = build_ihmm sim_seq build_model plot_model

plot_model_SOURCES = \
plot_model.c $(BEAMSOURCE) $(ISEQSOURCE) $(MODELSOURCE) $(FASTHMMSOURCE) $(RANDOMKIT_FILES) 
plot_model_LDADD = ${MYLIBDIRS}

build_model_SOURCES = \
build_model.c $(BEAMSOURCE) $(ISEQSOURCE) $(MODELSOURCE) $(FASTHMMSOURCE) $(RANDOMKIT_FILES) 
build_model_LDADD = ${MYLIBDIRS}

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
libihmm_a_SOURCES = $(RANDOMKIT_FILES) ihmm.h ihmm.c 

libihmm_a_LIBADD  =  ${MYLIBDIRS}


TESTS = fast_hmm_param_ITEST model_ITEST iseq_ITEST beam_ITEST dtest_ITEST

TESTS_ENVIRONMENT = \
	$(VALGRIND)

check_PROGRAMS = libihmm_ITEST fast_hmm_param_ITEST model_ITEST iseq_ITEST beam_ITEST dtest_ITEST

libihmm_ITEST_SOURCES = $(libihmm_a_SOURCES)
libihmm_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITEST
libihmm_ITEST_LDADD = $(libihmm_a_LIBADD)


iseq_ITEST_SOURCES = $(ISEQSOURCE) $(RANDOMKIT_FILES) 
iseq_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTSEQ
iseq_ITEST_LDADD = $(libihmm_a_LIBADD)

dtest_ITEST_SOURCES  = dtest.c $(RANDOMKIT_FILES) 
dtest_ITEST_CPPFLAGS = $(AM_CPPFLAGS) 
dtest_ITEST_LDADD = $(MYLIBDIRS)


fast_hmm_param_ITEST_SOURCES = $(FASTHMMSOURCE) $(RANDOMKIT_FILES) 
fast_hmm_param_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTFASTPARAM
fast_hmm_param_ITEST_LDADD = $(MYLIBDIRS)


model_ITEST_SOURCES =  $(MODELSOURCE) $(ISEQSOURCE) $(RANDOMKIT_FILES) 
model_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTMODEL
model_ITEST_LDADD = $(libihmm_a_LIBADD)


beam_ITEST_SOURCES = $(BEAMSOURCE) $(ISEQSOURCE) $(MODELSOURCE) $(FASTHMMSOURCE) $(RANDOMKIT_FILES) 
beam_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTBEAM
beam_ITEST_LDADD = $(libihmm_a_LIBADD)

