ACLOCAL_AMFLAGS= -I m4
# Build the library in the hand subdirectory first.
LIB_TLDEVELDIR = @LIB_TLDEVELDIR@

@VALGRIND_CHECK_RULES@

CLEANFILES = ihmm_seq_itest_read_test.fa test_model_file.txt test.lfa test_dna.lfa libhdf5gluetest_exp.h5 libhdf5gluetest.h5 model_at_*h5 test1.h5 test2.h5

# top_builddir = @top_builddir@
AM_CFLAGS = @OPENMP_CFLAGS@  @TLDEVEL_CPPFLAGS@

AM_LDFLAGS = @DEAD_STRIP@
AM_LDFLAGS += -Wno-undef
AM_LDFLAGS += -static

LIBS = @TLDEVEL_LIB@ $(HDF5_LDFLAGS)  $(HDF5_LIBS) -lm

#LIBS = ${MYLIBDIRS} $(HDF5_LDFLAGS)  $(HDF5_LIBS) -lpthread


#     AC_SUBST(HDF5_CFLAGS)
#     AC_SUBST(HDF5_CPPFLAGS)
#     AC_SUBST(HDF5_LDFLAGS)
#     AC_SUBST(HDF5_LIBS)

#AM_CPPFLAGS=${MYINCDIRS} ${HDF5_CPPFLAGS}


RANDOMKIT_FILES = distributions.h \
distributions.c \
randomkit.h \
randomkit.c \
randomkit_tl_add.h \
randomkit_tl_add.c \
randomkit_io.h \
randomkit_io.c

ADJUSTEDRANDINDEXSOURCE = adjusted_rand_index.c adjusted_rand_index.h

MOTIFREFINE = motif_refinement.h motif_refinement.c

NULLMODEL = null_model_emission.h null_model_emission.c

ISEQSOURCE = ihmm_seq.h ihmm_seq.c

FASTHMMSOURCE = fast_hmm_param.c fast_hmm_param.h  \
fast_hmm_param_test_functions.h \
fast_hmm_param_test_functions.c

KALIGNSOURCE = kalign.h kalign.c

MODELSOURCE = \
model_core.h \
model_core.c \
model_io.h \
model_io.c \
model_struct.h \
model_alloc.h \
model_alloc.c \
model_help.h \
model_help.c \
global.h \
$(NULLMODEL)

BEAMSOURCE = beam_sample.h beam_sample.c

THREADSOURCE = \
thread_data.h \
thread_data.c \
thread_data_io.h \
thread_data_io.c


# HDF5GLUESOURCE = hdf5_glue.c hdf5_glue.h #save_random_state.c

SCORESOURCE = run_score.h run_score.c

EMITSOURCE = run_emit.c emit_random.h

CALIBRATESOURCE = calibrate.c calibrate.h

CONVERSION = hmm_conversion.h hmm_conversion.c

FINITEHMM = \
finite_hmm.c \
finite_hmm.h \
finite_hmm_struct.h \
finite_hmm_alloc.h \
finite_hmm_alloc.c \
finite_hmm_io.h \
finite_hmm_io.c


DIJSTRA_SOURCES = dijkstra.c dijkstra.h

LOGO_SOURCES = motif_logo.c motif_logo.h


bin_PROGRAMS = sim_seq seqer_model seqer_plot seqer_score  seqer_iclu seqer_emit cal tfbs_bench composition_bench seqer_select seqer_ari pst_test seqer_emit_motif

#seqer_lcs
if ENABLE_CAIRO
bin_PROGRAMS +=  cairo_test seqer_pplot seqer_seqmotif
AM_CFLAGS += $(CAIRO_CFLAGS)
cairo_test_LDADD = $(CAIRO_LIBS) $(MYLIBDIRS)
seqer_pplot_LDADD = $(CAIRO_LIBS) $(MYLIBDIRS)
seqer_seqmotif_LDADD = $(CAIRO_LIBS) $(MYLIBDIRS)
endif


seqer_model_SOURCES = \
build_model.c $(BEAMSOURCE) $(SCORESOURCE) $(CONVERSION) $(ISEQSOURCE) $(FINITEHMM) $(MODELSOURCE) $(FASTHMMSOURCE) $(RANDOMKIT_FILES) $(THREADSOURCE)


seqer_ari_SOURCES = model_ari_comparison.c $(MODELSOURCE) $(ISEQSOURCE) $(ADJUSTEDRANDINDEXSOURCE) $(FINITEHMM) $(RANDOMKIT_FILES)

seqer_seqmotif_SOURCES = draw_sequence_logo.c $(LOGO_SOURCES) $(SCORESOURCE) $(MODELSOURCE) $(THREADSOURCE) $(ISEQSOURCE) $(FINITEHMM) $(RANDOMKIT_FILES)

 seqer_pplot_SOURCES = pplot.c   $(LOGO_SOURCES)  $(DIJSTRA_SOURCES) $(CONVERSION) $(ISEQSOURCE) $(FINITEHMM)  $(MODELSOURCE) $(FASTHMMSOURCE) $(RANDOMKIT_FILES)

cairo_test_SOURCES = cairo_test.c

cal_SOURCES = $(CALIBRATESOURCE) $(ISEQSOURCE) $(FINITEHMM) $(EMITSOURCE) $(SCORESOURCE) $(RANDOMKIT_FILES) $(THREADSOURCE)
cal_LDADD = ${MYLIBDIRS} $(HDF5_LDFLAGS)  $(HDF5_LIBS)


seqer_plot_SOURCES = \
plot_model.c  $(CONVERSION) $(ISEQSOURCE) $(FINITEHMM)  $(MODELSOURCE) $(FASTHMMSOURCE) $(RANDOMKIT_FILES)  $(HDF5GLUESOURCE)
#seqwise_plot_LDADD = ${MYLIBDIRS} $(HDF5_LDFLAGS)  $(HDF5_LIBS)



seqer_select_SOURCES = \
select_model.c  $(SCORESOURCE) $(CONVERSION) $(ISEQSOURCE) $(FINITEHMM) $(MODELSOURCE) $(FASTHMMSOURCE) $(RANDOMKIT_FILES) $(THREADSOURCE)

#seqwise_kmer_SOURCES = \
#kmer_enrichment.c $(KALIGNSOURCE) $(MOTIFREFINE) $(SCORESOURCE) $(EMITSOURCE) $(CONVERSION) $(ISEQSOURCE) $(FINITEHMM) $(MODELSOURCE) $(RANDOMKIT_FILES) $(HDF5GLUESOURCE) $(THREADSOURCE)


seqer_score_SOURCES = \
score_sequences.c $(CONVERSION)  $(SCORESOURCE) $(ISEQSOURCE) $(FINITEHMM) $(MODELSOURCE) $(FASTHMMSOURCE)  $(RANDOMKIT_FILES) $(THREADSOURCE)
#seqer_score_LDADD = ${MYLIBDIRS} $(HDF5_LDFLAGS)  $(HDF5_LIBS)

# seqer_lcs_SOURCES = $(ISEQSOURCE) $(RANDOMKIT_FILES)  $(FINITEHMM) lcs.c lcs.h
#seqer_lcs_LDADD = $(MYLIBDIRS) $(HDF5_LDFLAGS)  $(HDF5_LIBS)


seqer_iclu_SOURCES = $(SCORESOURCE) $(MODELSOURCE) $(MOTIFREFINE) $(THREADSOURCE) $(ISEQSOURCE) $(RANDOMKIT_FILES)  $(FINITEHMM) label_suffix_array.h label_suffix_array.c  infoclust.c infoclust.h
#seqer_iclu_LDADD = $(MYLIBDIRS) $(HDF5_LDFLAGS)  $(HDF5_LIBS)

seqer_emit_motif_SOURCES = $(MODELSOURCE) $(RANDOMKIT_FILES)  $(FINITEHMM) emit_motif.c
#seqer_iclu_LDADD = $(MYLIBDIRS) $(HDF5_LDFLAGS)  $(HDF5_LIBS)


seqer_emit_SOURCES = $(ISEQSOURCE) $(MODELSOURCE)  $(RANDOMKIT_FILES)  $(FINITEHMM) $(EMITSOURCE) emit_random.c

#seqer_emit_LDADD = $(MYLIBDIRS) $(HDF5_LDFLAGS)  $(HDF5_LIBS)


sim_seq_SOURCES = $(RANDOMKIT_FILES) \
outdir.h \
outdir.c \
matrix_io.h \
matrix_io.c \
sim_seq.c \
make_dot_file.h \
make_dot_file.c

sim_seq_LDADD =  ${MYLIBDIRS}

tfbs_bench_SOURCES = $(RANDOMKIT_FILES) \
outdir.h \
outdir.c \
matrix_io.h \
matrix_io.c \
benchmark_seq.h \
benchmark_seq.c \
create_TFBS_benchmark_datasets.c

tfbs_bench_LDADD =  ${MYLIBDIRS}



composition_bench_SOURCES = $(RANDOMKIT_FILES) \
outdir.h \
outdir.c \
matrix_io.h \
matrix_io.c \
benchmark_seq.h \
benchmark_seq.c \
create_composition_benchmark_dataset.c

composition_bench_LDADD =  ${MYLIBDIRS}
#build_ihmm_SOURCES = \
model_nuc_seq.c \
matrix_io.h \
matrix_io.c \
ihmm_io.h \
ihmm_io.c \
outdir.h \
outdir.c \
make_dot_file.h \
make_dot_file.c

#build_ihmm_LDADD =  ${MYLIBDIRS}  libihmm.a



#noinst_LIBRARIES = libihmm.a

# AUTOMAKE_OPTIONS = subdir-objects
#libihmm_a_SOURCES = $(RANDOMKIT_FILES) ihmm.h ihmm.c

#libihmm_a_LIBADD  =  ${MYLIBDIRS}


TESTS =  kalign_ITEST dijkstra_ITEST fast_hmm_param_ITEST model_ITEST 	iseq_ITEST beam_ITEST dtest_ITEST   ari_ITEST

check_PROGRAMS = kalign_ITEST dijkstra_ITEST fast_hmm_param_ITEST model_ITEST iseq_ITEST beam_ITEST dtest_ITEST ari_ITEST randomkit_tl_test sequences_TEST
 # pst_ITEST

#libihmm_ITEST_SOURCES = $(libihmm_a_SOURCES)
#libihmm_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITEST
#libihmm_ITEST_LDADD = $(libihmm_a_LIBADD)


dijkstra_ITEST_SOURCES = $(DIJSTRA_SOURCES)
dijkstra_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITEST
dijkstra_ITEST_LDADD = $(libihmm_a_LIBADD)


iseq_ITEST_SOURCES = $(ISEQSOURCE) $(RANDOMKIT_FILES)
iseq_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTSEQ
iseq_ITEST_LDADD = $(libihmm_a_LIBADD) $(HDF5_LDFLAGS)  $(HDF5_LIBS)

dtest_ITEST_SOURCES  = dtest.c $(RANDOMKIT_FILES)
dtest_ITEST_CPPFLAGS = $(AM_CPPFLAGS)
dtest_ITEST_LDADD = $(MYLIBDIRS)


fast_hmm_param_ITEST_SOURCES = $(FASTHMMSOURCE) $(RANDOMKIT_FILES)
fast_hmm_param_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTFASTPARAM
fast_hmm_param_ITEST_LDADD = $(MYLIBDIRS)


# hdf5glue_ITEST_SOURCES = $(HDF5GLUESOURCE)
# hdf5glue_ITEST_CPPFLAGS  = $(AM_CPPFLAGS) $(HDF5_CPPFLAGS)  -DITEST
# hdf5glue_ITEST_LDADD  = $(MYLIBDIRS) $(HDF5_LDFLAGS)  $(HDF5_LIBS)


model_ITEST_SOURCES =  $(MODELSOURCE) $(ISEQSOURCE) $(RANDOMKIT_FILES)
model_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTMODEL
model_ITEST_LDADD = $(libihmm_a_LIBADD) $(HDF5_LDFLAGS)  $(HDF5_LIBS)


beam_ITEST_SOURCES = $(BEAMSOURCE) $(THREADSOURCE) $(ISEQSOURCE) $(MODELSOURCE) $(FASTHMMSOURCE) $(RANDOMKIT_FILES)  $(CONVERSION) $(FINITEHMM)
beam_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTBEAM
beam_ITEST_LDADD = $(libihmm_a_LIBADD) $(HDF5_LDFLAGS)  $(HDF5_LIBS)

# pst_ITEST_SOURCES = $(ISEQSOURCE)  $(RANDOMKIT_FILES)  $(FINITEHMM) lcs.c lcs.h
# pst_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTLST
# pst_ITEST_LDADD = $(MYLIBDIRS) $(HDF5_LDFLAGS)  $(HDF5_LIBS)


ari_ITEST_SOURCES = $(ADJUSTEDRANDINDEXSOURCE)
ari_ITEST_CPPFLAGS = $(AM_CPPFLAGS) -DITEST
ari_ITEST_LDADD = $(MYLIBDIRS)

kalign_ITEST_SOURCES = $(KALIGNSOURCE)
kalign_ITEST_CPPFLAGS = $(AM_CPPFLAGS) -DITEST
kalign_ITEST_LDADD = $(MYLIBDIRS)

randomkit_tl_test_SOURCES = $(RANDOMKIT_FILES)
randomkit_tl_test_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_COPY_RK
randomkit_tl_test_LDADD = $(MYLIBDIRS)


pst_test_SOURCES = \
pst.h \
pst.c \
pst_structs.h \
pst_test.c \
sim_seq_lib.h \
sim_seq_lib.c

sequences_TEST_SOURCES = \
sequence_alloc.h \
sequence_alloc.c \
sequence_io.h \
sequence_io.c \
sequence_prep.h \
sequence_prep.c \
sequence_test.c

