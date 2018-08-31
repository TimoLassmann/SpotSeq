# Build the library in the hand subdirectory first.
LIB_TLDEVELDIR = @LIB_TLDEVELDIR@ 


CLEANFILES = ihmm_seq_itest_read_test.fa test_model_file.txt test.lfa  libhdf5gluetest_exp.h5 libhdf5gluetest.h5 model_at_*h5

# top_builddir = @top_builddir@


MYLIBDIRS = ../tldevel/libtldevel.a


MYINCDIRS = -I${top_builddir}/${LIB_TLDEVELDIR} \
            -I${top_srcdir}/${LIB_TLDEVELDIR} 


LIBS = -lm -lpthread


#     AC_SUBST(HDF5_CFLAGS)
#     AC_SUBST(HDF5_CPPFLAGS)
#     AC_SUBST(HDF5_LDFLAGS)
#     AC_SUBST(HDF5_LIBS)

AM_CPPFLAGS=${MYINCDIRS} ${HDF5_CPPFLAGS}


RANDOMKIT_FILES = distributions.h \
distributions.c \
randomkit.h \
randomkit.c


ISEQSOURCE = ihmm_seq.h ihmm_seq.c

FASTHMMSOURCE = fast_hmm_param.c fast_hmm_param.h  \
fast_hmm_param_test_functions.h \
fast_hmm_param_test_functions.c

MODELSOURCE = model.h model.c model_io.c global.h

BEAMSOURCE = beam_sample.h beam_sample.c

THREADSOURCE = thread_data.c thread_data.h

HDF5GLUESOURCE = hdf5_glue.c hdf5_glue.h

CONVERSION = hmm_conversion.h hmm_conversion.c

FINITEHMM = finite_hmm.c finite_hmm.h

bin_PROGRAMS = build_ihmm sim_seq spotseq_model spotseq_plot spotseq_score spotseq_lcs spotseq_emit



spotseq_plot_SOURCES = \
plot_model.c  $(CONVERSION) $(ISEQSOURCE) $(FINITEHMM)  $(MODELSOURCE) $(FASTHMMSOURCE) $(RANDOMKIT_FILES)  $(HDF5GLUESOURCE)
spotseq_plot_LDADD = ${MYLIBDIRS} $(HDF5_LDFLAGS)  $(HDF5_LIBS)

spotseq_model_SOURCES = \
build_model.c $(BEAMSOURCE) $(CONVERSION) $(ISEQSOURCE) $(FINITEHMM) $(MODELSOURCE) $(FASTHMMSOURCE) $(RANDOMKIT_FILES) $(HDF5GLUESOURCE) $(THREADSOURCE)
spotseq_model_LDADD = ${MYLIBDIRS} $(HDF5_LDFLAGS)  $(HDF5_LIBS)

spotseq_score_SOURCES = \
score_sequences.c $(CONVERSION)  $(ISEQSOURCE) $(FINITEHMM) $(MODELSOURCE) $(FASTHMMSOURCE)  $(RANDOMKIT_FILES) $(HDF5GLUESOURCE) $(THREADSOURCE)
spotseq_score_LDADD = ${MYLIBDIRS} $(HDF5_LDFLAGS)  $(HDF5_LIBS)

spotseq_lcs_SOURCES = $(ISEQSOURCE) $(HDF5GLUESOURCE) $(RANDOMKIT_FILES)  $(FINITEHMM) lcs.c lcs.h
spotseq_lcs_LDADD = $(MYLIBDIRS) $(HDF5_LDFLAGS)  $(HDF5_LIBS)

spotseq_emit_SOURCES = $(ISEQSOURCE) $(HDF5GLUESOURCE) $(RANDOMKIT_FILES)  $(FINITEHMM) emit_random.c emit_random.h
spotseq_emit_LDADD = $(MYLIBDIRS) $(HDF5_LDFLAGS)  $(HDF5_LIBS)


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


#TESTS = fast_hmm_param_ITEST model_ITEST 	iseq_ITEST beam_ITEST dtest_ITEST   hdf5glue_ITEST

TESTS_ENVIRONMENT = \
	$(VALGRIND)

check_PROGRAMS = pst_ITEST libihmm_ITEST fast_hmm_param_ITEST model_ITEST iseq_ITEST beam_ITEST dtest_ITEST  hdf5glue_ITEST

libihmm_ITEST_SOURCES = $(libihmm_a_SOURCES)
libihmm_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITEST
libihmm_ITEST_LDADD = $(libihmm_a_LIBADD)


iseq_ITEST_SOURCES = $(ISEQSOURCE) $(RANDOMKIT_FILES) $(HDF5GLUESOURCE)
iseq_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTSEQ
iseq_ITEST_LDADD = $(libihmm_a_LIBADD) $(HDF5_LDFLAGS)  $(HDF5_LIBS)

dtest_ITEST_SOURCES  = dtest.c $(RANDOMKIT_FILES)
dtest_ITEST_CPPFLAGS = $(AM_CPPFLAGS)
dtest_ITEST_LDADD = $(MYLIBDIRS)


fast_hmm_param_ITEST_SOURCES = $(FASTHMMSOURCE) $(RANDOMKIT_FILES)
fast_hmm_param_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTFASTPARAM
fast_hmm_param_ITEST_LDADD = $(MYLIBDIRS)


hdf5glue_ITEST_SOURCES = $(HDF5GLUESOURCE)
hdf5glue_ITEST_CPPFLAGS  = $(AM_CPPFLAGS) $(HDF5_CPPFLAGS)  -DITEST
hdf5glue_ITEST_LDADD  = $(MYLIBDIRS) $(HDF5_LDFLAGS)  $(HDF5_LIBS)


model_ITEST_SOURCES =  $(MODELSOURCE) $(ISEQSOURCE) $(RANDOMKIT_FILES) $(HDF5GLUESOURCE)
model_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTMODEL
model_ITEST_LDADD = $(libihmm_a_LIBADD) $(HDF5_LDFLAGS)  $(HDF5_LIBS)


beam_ITEST_SOURCES = $(BEAMSOURCE) $(THREADSOURCE) $(ISEQSOURCE) $(MODELSOURCE)  $(HDF5GLUESOURCE)  $(FASTHMMSOURCE) $(RANDOMKIT_FILES)  $(CONVERSION) $(FINITEHMM)
beam_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTBEAM
beam_ITEST_LDADD = $(libihmm_a_LIBADD) $(HDF5_LDFLAGS)  $(HDF5_LIBS)

pst_ITEST_SOURCES = $(ISEQSOURCE) $(HDF5GLUESOURCE) $(RANDOMKIT_FILES)  $(FINITEHMM) lcs.c lcs.h
pst_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTLST
pst_ITEST_LDADD = $(MYLIBDIRS) $(HDF5_LDFLAGS)  $(HDF5_LIBS)

