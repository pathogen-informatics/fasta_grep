# ship these headers in the "make dist" target
dist_data_DATA = kseq.h parse_fasta.h


# fasta_grep is our top level progra
bin_PROGRAMS = fasta_grep
fasta_grep_SOURCES = main.c
fasta_grep_LDADD=-lm libfasta_grep.la -lz

# libfasta_grep.so is our library
lib_LTLIBRARIES = libfasta_grep.la
libfasta_grep_la_SOURCES = parse_fasta.c
libfasta_grep_la_LDFLAGS= -version-info 1:0

# "make check" target
TESTS = $(check_PROGRAMS)
check_PROGRAMS = run_all_tests
run_all_tests_SOURCES = \
	../tests/check_parse_fasta.c \
	../tests/run_all_tests.c
run_all_tests_CFLAGS = -I../tests
run_all_tests_LDADD = -lm -lcheck libfasta_grep.la -lz