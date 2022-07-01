/*********************************************************************************
 * MIT License                                                                   *
 *                                                                               *
 * Copyright (c) 2021 Chenxi Zhou <chnx.zhou@gmail.com>                          *
 *                                                                               *
 * Permission is hereby granted, free of charge, to any person obtaining a copy  *
 * of this software and associated documentation files (the "Software"), to deal *
 * in the Software without restriction, including without limitation the rights  *
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     *
 * copies of the Software, and to permit persons to whom the Software is         *
 * furnished to do so, subject to the following conditions:                      *
 *                                                                               *
 * The above copyright notice and this permission notice shall be included in    *
 * all copies or substantial portions of the Software.                           *
 *                                                                               *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, *
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE *
 * SOFTWARE.                                                                     *
 *********************************************************************************/

/********************************** Revision History *****************************
 *                                                                               *
 * 15/04/21 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "khash.h"
#include "asset.h"
#include "sdict.h"
#include "ksort.h"
#include "kseq.h"

#undef DEBUG

KSEQ_INIT(gzFile, gzread, gzseek)

void *kopen(const char *fn, int *_fd);
int kclose(void *a);

sdict_t *sd_init(void)
{
    sdict_t *d;
    d = (sdict_t *) calloc(1, sizeof(sdict_t));
    d->n = 0;
    d->m = 16;
    d->s = (sd_seq_t *) malloc(d->m * sizeof(sd_seq_t));
    d->h = kh_init(sdict);
    return d;
}

asm_dict_t *asm_init(sdict_t *sdict)
{
    asm_dict_t *d;
    d = (asm_dict_t *) calloc(1, sizeof(asm_dict_t));
    d->n = 0;
    d->m = 16;
    d->s = (sd_aseq_t *) malloc(d->m * sizeof(sd_aseq_t));
    d->h = kh_init(sdict);
    d->u = 0;
    d->v = 16;
    d->seg = (sd_seg_t *) malloc(d->v * sizeof(sd_seg_t));
    d->a = (uint32_t *) malloc(sdict->n * sizeof(uint32_t));
    d->index = 0;
    d->sdict = sdict;
    return d;
}

void sd_destroy(sdict_t *d)
{
    uint32_t i;
    if (d == 0)
        return;
    if (d->h)
        kh_destroy(sdict, d->h);
    if(d->s) {
        for (i = 0; i < d->n; ++i) {
            free(d->s[i].name);
            if (d->s[i].seq)
                free(d->s[i].seq);
        }
        free(d->s);
    }
    free(d);
}

void asm_destroy(asm_dict_t *d)
{
    uint32_t i;
    if (d == 0)
        return;
    if (d->s) {
        for (i = 0; i < d->n; ++i)
            free(d->s[i].name);
        free(d->s);
    }
    if (d->seg)
        free(d->seg);
    if (d->a)
        free(d->a);
    if (d->index)
        free(d->index);
    if (d->h)
        kh_destroy(sdict, d->h);
    free(d);
}

uint32_t sd_put(sdict_t *d, const char *name, uint32_t len)
{
    if (!name)
        return UINT32_MAX;
    sdhash_t *h = d->h;
    khint_t k;
    int absent;
    k = kh_put(sdict, h, name, &absent);
    if (absent) {
        sd_seq_t *s;
        if (d->n == d->m) {
            d->m = d->m? d->m<<1 : 16;
            d->s = (sd_seq_t *) realloc(d->s, d->m * sizeof(sd_seq_t));
        }
        s = &d->s[d->n];
        s->len = len;
        s->seq = 0;
        kh_key(h, k) = s->name = strdup(name);
        kh_val(h, k) = d->n++;
    }
    return kh_val(h, k);
}

uint32_t sd_put1(sdict_t *d, const char *name, const char *seq, uint32_t len)
{
    uint32_t k = sd_put(d, name, len);
    d->s[k].seq = strdup(seq);
    return k;
}

uint32_t sd_get(sdict_t *d, const char *name)
{
    sdhash_t *h = d->h;
    khint_t k;
    k = kh_get(sdict, h, name);
    return k == kh_end(h)? UINT32_MAX : kh_val(h, k);
}

void sd_hash(sdict_t *d)
{
    uint32_t i;
    sdhash_t *h;
    if (d->h)
        return;
    d->h = h = kh_init(sdict);
    for (i = 0; i < d->n; ++i) {
        int absent;
        khint_t k;
        k = kh_put(sdict, h, d->s[i].name, &absent);
        kh_val(h, k) = i;
    }
}

sdict_t *make_sdict_from_fa(const char *f, uint32_t min_len)
{
    int fd, l;
    gzFile fp;
    kseq_t *ks;
    void *ko = 0;

    ko = kopen(f, &fd);
    if (ko == 0) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }
    fp = gzdopen(fd, "r");
    ks = kseq_init(fp);
    
    sdict_t *d;
    d = sd_init();
    // WARNING: might introduce bug here
    // l is int type by definition in kseq.h
    while ((l = kseq_read(ks)) >= 0)
        if (strlen(ks->seq.s) >= min_len)
            sd_put1(d, ks->name.s, ks->seq.s, strlen(ks->seq.s));

    kseq_destroy(ks);
    gzclose(fp);
    kclose(ko);

    return d;
}

sdict_t *make_sdict_from_index(const char *f, uint32_t min_len)
{
    FILE *fp;
    char *line = NULL;
    size_t ln = 0;
    ssize_t read;
    char name[4096];
    uint32_t len;

    fp = fopen(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }
    
    sdict_t *d;
    d = sd_init();
    while ((read = getline(&line, &ln, fp)) != -1) {
        sscanf(line, "%s %u", name, &len);
        if (len >= min_len)
            sd_put(d, name, len);
    }
    fclose(fp);
    return d;
}

sdict_t *make_sdict_from_gfa(const char *f, uint32_t min_len)
{
    FILE *fp;
    char *line = NULL;
    size_t ln = 0;
    ssize_t read;
    char name[4096], lens[4096];
    uint32_t len;

    fp = fopen(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    sdict_t *d;
    d = sd_init();
    while ((read = getline(&line, &ln, fp)) != -1) {
        if (line[0] == 'S') {
            sscanf(line, "%*s %s %*s %s", name, lens);
            len = strtoul(lens + 5, NULL, 10);
            if (len >= min_len)
                sd_put(d, name, len);
        }
    }
    fclose(fp);

    return d;
}

uint32_t asm_put(asm_dict_t *d, const char *name, uint64_t len, uint32_t n, uint32_t s)
{
    if (!name) 
        return UINT32_MAX;
    sdhash_t *h;
    khint_t k;
    int absent;

    h = d->h;
    k = kh_put(sdict, h, name, &absent);
    if (absent) {
        sd_aseq_t *a;
        if (d->n == d->m) {
            d->m = d->m? d->m<<1 : 16;
            d->s = (sd_aseq_t *) realloc(d->s, d->m * sizeof(sd_aseq_t));
        }
        a = &d->s[d->n];
        a->len = len;
        a->n = n;
        a->s = s;
        kh_key(h, k) = a->name = strdup(name);
        kh_val(h, k) = d->n++;
    }
    return kh_val(h, k);
}

uint32_t asm_sd_get(asm_dict_t *d, const char *name)
{
    sdhash_t *h = d->h;
    khint_t k;
    k = kh_get(sdict, h, name);
    return k == kh_end(h)? UINT32_MAX : kh_val(h, k);
}

static void seg_put(asm_dict_t *d, uint32_t s, uint32_t k, uint64_t a, uint32_t c, uint32_t x, uint32_t y)
{
    if (d->u == d->v) {
        d->v = d->v? d->v<<1 : 16;
        d->seg = (sd_seg_t *) realloc(d->seg, d->v * sizeof(sd_seg_t));
    }
    sd_seg_t seg = {s, k, a, c, x, y};
    d->seg[d->u++] = seg;
}

asm_dict_t *make_asm_dict_from_sdict(sdict_t *sdict)
{
    uint32_t i;
    asm_dict_t *d;
    d = asm_init(sdict);
    d->index = (uint64_t *) calloc(sdict->n, sizeof(uint64_t));
    uint32_t *a = d->a;
    for (i = 0; i < sdict->n; ++i) {
        asm_put(d, sdict->s[i].name, sdict->s[i].len, 1, i);
        seg_put(d, i, 0, 0, i<<1, 0, sdict->s[i].len);
        a[i] = i;
        d->index[i] = (uint64_t) sdict->s[i].len << 32 | i;
    }
    return d;
}

typedef struct {uint64_t x, y;} c_pair_t;

#define pair_key(a) ((a).x)
KRADIX_SORT_INIT(pair, c_pair_t, pair_key, 8)

void asm_index(asm_dict_t *d)
{
    int32_t i, s, c, c1;
    sd_seg_t seg;
    uint64_t x, y;
    uint32_t *a;

    s = 0;
    for (i = 0; i < d->n; ++i)
        s += d->s[i].n;
    
    c_pair_t *c_pairs;
    c_pairs = (c_pair_t *) malloc(s * sizeof(c_pair_t));
    for(i = 0; i < s; ++i) {
        seg = d->seg[i];
        x = (uint64_t) seg.c >> 1 << 32 | (seg.x + seg.y);
        y = (uint64_t) (seg.x + seg.y) << 32 | i;
        c_pair_t pair = {x, y};
        c_pairs[i] = pair;
    }

    radix_sort_pair(c_pairs, c_pairs + s);
    
    if (d->index)
        free(d->index);
    d->index = (uint64_t *) malloc(s * sizeof(uint64_t));
    a = d->a;
    c1 = INT32_MAX;
    for (i = 0; i < s; ++i) {
        d->index[i] = c_pairs[i].y;
        c = c_pairs[i].x >> 32;
        if (c != c1) {
            a[c] = i;
            c1 = c;
        }
    }

    free(c_pairs);
}

asm_dict_t *make_asm_dict_from_agp(sdict_t *sdict, const char *f)
{
    FILE *fp;
    char *line = NULL;
    size_t ln = 0;
    ssize_t read;
    char sname[256], type[4], cname[256], cstarts[16], cends[16], oris[4];
    char *name = NULL;
    uint64_t a;
    uint32_t l, cstart, cend;
    uint32_t c, s, n;

    fp = fopen(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    asm_dict_t *d;
    d = asm_init(sdict);
    a = 0;
    s = n = 0;
    while ((read = getline(&line, &ln, fp)) != -1) {
        sscanf(line, "%s %*s %*s %*s %s %s %s %s %s", sname, type, cname, cstarts, cends, oris);
        if (!strncmp(type, "N", 1))
            continue;
        if (!name)
            name = strdup(sname);
        if (strcmp(sname, name)) {
            asm_put(d, name, a, n, s - n);
            a = 0;
            n = 0;
            name = strdup(sname);
        }
        cstart = strtoul(cstarts, NULL, 10);
        cend = strtoul(cends, NULL, 10);
        c = sd_get(sdict, cname);
        if (c == UINT32_MAX) {
            fprintf(stderr, "[E::%s] sequence %s not found\n", __func__, cname);
            return 0;
        }
        c <<= 1;
        if (strncmp(oris, "+", 1)) 
            c |= 1;
        l = cend - cstart + 1;
        seg_put(d, d->n, n, a, c, cstart - 1, l);
        a += l;
        ++s;
        ++n;
    }
    asm_put(d, name, a, n, s - n);
    d->sdict = sdict;
#ifdef DEBUG
    for (int i = 0; i < s; ++i)
        fprintf(stderr, "%u %lu %u %u %u\n", d->seg[i].s, d->seg[i].a, d->seg[i].c, d->seg[i].x, d->seg[i].y);
    for (int i = 0; i < d->n; ++i)
        fprintf(stderr, "%s %lu %u %u\n", d->s[i].name, d->s[i].len, d->s[i].n, d->s[i].s);
#endif
    asm_index(d);

    fclose(fp);
    if (line)
        free(line);
    if (name)
        free(name);
    
    return d;
}

void add_unplaced_short_seqs(asm_dict_t *d, uint32_t min_len)
{
    uint32_t i;
    sdict_t *sdict;
    
    sdict = d->sdict;
    for (i = 0; i < sdict->n; ++i) {
        if (sdict->s[i].len >= min_len)
            continue;
        seg_put(d, d->n, 0, 0, i<<1, 0, sdict->s[i].len);
        asm_put(d, sdict->s[i].name, sdict->s[i].len, 1, d->u - 1);
    }
    asm_index(d);
}


/* contig coordinates to scaffold coordinates */
// one-based
int sd_coordinate_conversion(asm_dict_t *d, uint32_t id, uint32_t pos, uint32_t *s, uint64_t *p, int count_gap)
{
    if (id == UINT32_MAX) {
        *s = UINT32_MAX;
        return 1;
    }
    uint32_t i;
    uint64_t *index = d->index;
    i = d->a[id];
    while (index[i]>>32 < pos) 
        ++i;
    sd_seg_t seg = d->seg[(uint32_t) index[i]];
    *s = seg.s;
    *p = seg.c & 1? seg.a + seg.x + seg.y - pos + 1 : seg.a + pos - seg.x;
    if (count_gap)
        *p = *p + seg.k * GAP_SZ;
    return 0;
}

/* scaffold coordinates to contig coordinates */
// one-based
int sd_coordinate_rev_conversion(asm_dict_t *d, uint32_t id, uint64_t pos, uint32_t *s, uint32_t *p, int count_gap)
{
    if (id == UINT32_MAX) {
        *s = UINT32_MAX;
        return 1;    
    }

    sd_seg_t *seg = d->seg + d->s[id].s;
    uint32_t n = d->s[id].n;
    int64_t l = pos;

    uint32_t i = 0;
    int in_gap = 0;
    while (i < n && l > 0) {
        l -= seg[i].y;
        in_gap = 0;
        ++i;
        if (count_gap && i < n && l > 0) {
            l -= GAP_SZ;
            in_gap = 1;
        }
    }
    
    if (l > 0 || in_gap) {
        *s = UINT32_MAX;
        return 1;
    }

    --i;
    *s = seg[i].c >> 1;
    *p = seg[i].c & 1? seg[i].x - l + 1 : seg[i].x + seg[i].y + l;

    return 0;
}

int cmp_uint64_d (const void *a, const void *b) {
    // decreasing order
    uint64_t x, y;
    x = *(uint64_t *) a;
    y = *(uint64_t *) b;
    return x == y? 0 : (x < y? 1 : -1);
}

static void nl_stats(uint64_t *s, uint32_t n, uint64_t *n_stats, uint32_t *l_stats)
{
    uint32_t i, j;
    double b, a, bs;
    
    bs = 0;
    for (i = 0; i < n; ++i)
        bs += s[i];

    qsort(s, n, sizeof(uint64_t), cmp_uint64_d);

    b = 0;
    j = 0;
    a = .1 * (++j) * bs;
    for (i = 0; i < n; ++i) {
        b += s[i];
        while (b >= a) {
            l_stats[j - 1] = i + 1;
            n_stats[j - 1] = s[i];
            a = .1 * (++j) * bs;
        }
    }
}

void sd_stats(sdict_t *d, uint64_t *n_stats, uint32_t *l_stats)
{
    // n_stats and l_stats are of at least size 10
    uint32_t i, n;
    uint64_t *s;

    n = d->n;
    s = (uint64_t *) calloc(n, sizeof(uint64_t));
    for (i = 0; i < n; ++i)
        s[i] = d->s[i].len;
    nl_stats(s, n, n_stats, l_stats);
    
    free(s);
}

void asm_sd_stats(asm_dict_t *d, uint64_t *n_stats, uint32_t *l_stats)
{
    // n_stats and l_stats are of at least size 10
    uint32_t i, n;
    uint64_t *s;
    
    n = d->n;
    s = (uint64_t *) calloc(n, sizeof(uint64_t));
    for (i = 0; i < n; ++i)
        s[i] = d->s[i].len;
    nl_stats(s, n, n_stats, l_stats);

    free(s);
}

char comp_table[] = {
      0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
     16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
     32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
     48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
     64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
     64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
    'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

char nucl_toupper[] = {
      0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
     16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
     32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
     48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
     64, 'A', 'N', 'C', 'N', 'N', 'N', 'G', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'T', 'N', 'N', 'N', 'N', 'N', 'N',  91,  92,  93,  94,  95,
     64, 'A', 'N', 'C', 'N', 'N', 'N', 'G', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'T', 'N', 'N', 'N', 'N', 'N', 'N', 123, 124, 125, 126, 127
};

void write_fasta_file_from_agp(const char *fa, const char *agp, FILE *fo, int line_wd, int un_oris)
{
    FILE *agp_in;
    char *line = NULL;
    sdict_t *dict;
    sd_seq_t s;
    size_t ln = 0;
    ssize_t read;
    char sname[256], type[4], cname[256], cstarts[16], cends[16], oris[16];
    char *name = NULL;
    uint32_t c;
    uint64_t i, l, cstart, cend;

    agp_in = fopen(agp, "r");
    if (agp_in == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, agp);
        exit(EXIT_FAILURE);
    }

    dict = make_sdict_from_fa(fa, 0);
    l = 0;
    while ((read = getline(&line, &ln, agp_in)) != -1) {
        sscanf(line, "%s %*s %*s %*s %s %s %s %s %s", sname, type, cname, cstarts, cends, oris);
        if (!strncmp(sname, "#", 1))
            // header lines
            continue;
        if (!strncmp(type, "N", 1) || !strncmp(type, "U", 1)) {
            cend = strtoul(cname, NULL, 10);
            for (i = 0; i < cend; ++i) {
                fputc('N', fo);
                if (++l % line_wd == 0)
                    fputc('\n', fo);
            }
            continue;
        }
        if (!name) {
            name = strdup(sname);
            fprintf(fo, ">%s\n", name);
            l = 0;
        }
        if (strcmp(sname, name)) {
            free(name);
            name = strdup(sname);
            if (l % line_wd != 0)
                fputc('\n', fo);
            fprintf(fo, ">%s\n", name);
            l = 0;
        }

        cstart = strtoul(cstarts, NULL, 10);
        cend = strtoul(cends, NULL, 10);
        c = sd_get(dict, cname);
        if (c == UINT32_MAX) {
            fprintf(stderr, "[E::%s] sequence %s not found\n", __func__, cname);
            exit(EXIT_FAILURE);
        }
        s = dict->s[c];
        if (!strncmp(oris, "+", 1) || (un_oris && (!strncmp(oris, "?", 1) || !strncmp(oris, "0", 1) || !strncmp(oris, "na", 2)))) {
            // forward
            for (i = cstart - 1; i < cend; ++i) {
                fputc(s.seq[i], fo);
                if (++l % line_wd == 0)
                    fputc('\n', fo);
            }
        } else if (!strncmp(oris, "-", 1)) {
            // reverse
            for (i = cend; i >= cstart; --i) {
                fputc(comp_table[(int) s.seq[i - 1]], fo);
                if (++l % line_wd == 0)
                    fputc('\n', fo);
            }
        } else {
            fprintf(stderr, "[E::%s] unknown orientation of sequence component %s:%lu-%lu on %s: '%s'\n", __func__, cname, cstart, cend, sname, oris);
            if (un_oris) {
                fprintf(stderr, "[E::%s] valid identifiers for unorientated sequence include: '?', '0' and 'na'\n", __func__);
                fprintf(stderr, "[E::%s] see https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/#FORMAT\n", __func__);
            }
            exit(EXIT_FAILURE);
        }
    }
    if (l % line_wd != 0)
        fputc('\n', fo);
    free(name);
    free(dict);

    fclose(agp_in);
}

void write_segs_to_agp(sd_seg_t *segs, uint32_t n, sdict_t *sd, uint32_t s, FILE *fp)
{
    uint64_t len;
    sd_seg_t seg;
    uint32_t i, t;
    len = 0;
    t = 0;
    for (i = 0; i < n; ++i) {
        seg = segs[i];
        fprintf(fp, "scaffold_%u\t%lu\t%lu\t%d\tW\t%s\t%u\t%u\t%c\n", s, len + 1, len + seg.y, ++t, sd->s[seg.c >> 1].name, seg.x + 1, seg.x + seg.y, "+-"[seg.c & 1]);
        len += seg.y;
        if (i != n - 1) {
            fprintf(fp, "scaffold_%u\t%lu\t%lu\t%d\tN\t%d\tscaffold\tyes\tna\n", s, len + 1, len + GAP_SZ, ++t, GAP_SZ);
            len += GAP_SZ;
        }
    }
}

void write_sorted_agp(asm_dict_t *dict, FILE *fo)
{
    int i, j;
    uint32_t s;
    c_pair_t *c_pairs;

    c_pairs = (c_pair_t *) malloc(dict->n * sizeof(c_pair_t));
    for(i = 0; i < dict->n; ++i) {
        c_pair_t pair = {dict->s[i].len, i};
        c_pairs[i] = pair;
    }
    radix_sort_pair(c_pairs, c_pairs + dict->n);
    s = 0;
    for (i = dict->n - 1; i >= 0; --i) {
        j = c_pairs[i].y;
        write_segs_to_agp(dict->seg + dict->s[j].s, dict->s[j].n, dict->sdict, ++s, fo);
    }
    free(c_pairs);
}

void write_sdict_to_agp(sdict_t *sdict, char *out)
{
    int i;
    FILE *agp_out;
    agp_out = fopen(out, "w");
    if (agp_out == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < sdict->n; ++i)
        fprintf(agp_out, "scaffold_%u\t1\t%u\t1\tW\t%s\t1\t%u\t+\n", i + 1, sdict->s[i].len, sdict->s[i].name, sdict->s[i].len);
    fclose(agp_out);
}

void write_asm_dict_to_agp(asm_dict_t *dict, char *out)
{
    uint64_t len = 0;
    sd_seg_t seg;
    uint32_t i, j, n, s, t;
    FILE *agp_out;
    agp_out = fopen(out, "w");
    if (agp_out == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < dict->n; ++i) {
        len = 0;
        t = 0;
        s = dict->s[i].s;
        n = dict->s[i].n;
        for (j = 0; j < n; ++j) {
            seg = dict->seg[s + j];
            fprintf(agp_out, "%s\t%lu\t%lu\t%u\tW\t%s\t%u\t%u\t%c\n", dict->s[i].name, len + 1, len + seg.y, ++t, dict->sdict->s[seg.c >> 1].name, seg.x + 1, seg.x + seg.y, "+-"[seg.c & 1]);
            len += seg.y;
            if (j != n - 1) {
                fprintf(agp_out, "%s\t%lu\t%lu\t%d\tN\t%d\tscaffold\tyes\tna\n", dict->s[i].name, len + 1, len + GAP_SZ, ++t, GAP_SZ);
                len += GAP_SZ;
            }
        }
    }
        
    fclose(agp_out);
}

