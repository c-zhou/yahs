CFLAGS=		-g -O3 -Wall -fno-strict-aliasing
CPPFLAGS=
INCLUDES=
OBJS=
PROG=       yahs juicer_pre
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

yahs: asset.c bamlite.c break.c graph.c kalloc.c kopen.c link.c sdict.c yahs.c
		$(CC) $(CFLAGS) asset.c bamlite.c break.c graph.c kalloc.c kopen.c link.c sdict.c yahs.c -o $@ -L. $(LIBS)

juicer_pre: asset.c bamlite.c kalloc.c kopen.c link.c sdict.c juicer_pre.c
		$(CC) $(CFLAGS) asset.c bamlite.c kalloc.c kopen.c link.c sdict.c juicer_pre.c -o $@ -L. $(LIBS)

clean:
		rm -fr *.o a.out $(PROG) $(PROG_EXTRA)

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

