CFLAGS=		-g -O0 -Wall -fno-strict-aliasing
CPPFLAGS=
INCLUDES=
OBJS=
PROG=       yahs juicer agp_to_fasta ONEview
PROG_EXTRA=
LIBS=		-lm -lz

.PHONY:all extra clean depend
.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all: $(PROG)

extra: all $(PROG_EXTRA)

debug: $(PROG)
debug: CFLAGS += -DDEBUG

yahs: asset.c bamlite.c break.c graph.c kalloc.c kopen.c link.c sdict.c cov.c binomlite.c enzyme.c yahs.c ONElib.c
		$(CC) $(CFLAGS) asset.c bamlite.c break.c graph.c kalloc.c kopen.c link.c sdict.c cov.c binomlite.c enzyme.c yahs.c ONElib.c -o $@ -L. $(LIBS)

juicer: asset.c bamlite.c kalloc.c kopen.c sdict.c juicer.c ONElib.c
		$(CC) $(CFLAGS) asset.c bamlite.c kalloc.c kopen.c sdict.c juicer.c ONElib.c -o $@ -L. $(LIBS)

agp_to_fasta: asset.c kalloc.c kopen.c sdict.c agp_to_fasta.c ONElib.c
		$(CC) $(CFLAGS) asset.c kalloc.c kopen.c sdict.c agp_to_fasta.c ONElib.c -o $@ -L. $(LIBS)

ONEview: ONEview.c ONElib.c utils.c
		$(CC) $(CFLAGS) ONEview.c ONElib.c utils.c -o $@

clean:
		rm -fr *.o a.out $(PROG) $(PROG_EXTRA)

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

