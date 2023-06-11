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
#include "kvec.h"

#undef DEBUG_DICT

KSEQ_INIT(gzFile, gzread, gzseek)

void *kopen(const char *fn, int *_fd);
int kclose(void *a);

AGP_CT_t DEFAULT_AGP_SEQ_COMPONENT_TYPE = AGP_CT_W;
AGP_CT_t DEFAULT_AGP_GAP_COMPONENT_TYPE = AGP_CT_U;
AGP_GT_t DEFAULT_AGP_GAP_TYPE = AGP_GT_SCAFFOLD;
AGP_LE_t DEFAULT_AGP_LINKAGE_EVIDENCE = AGP_LE_PROXIMITY_LIGATION;
int DEFAULT_AGP_GAP_SIZE = DEFAULT_AGP_U_GAP_SIZE;

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
    d->a = (uint64_t *) malloc(sdict->n * sizeof(uint64_t));
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
        for (i = 0; i < d->n; ++i) {
            free(d->s[i].name);
            free(d->s[i].gs);
            free(d->s[i].gsize);
            free(d->s[i].gcomp);
            free(d->s[i].gtype);
            free(d->s[i].glink);
            free(d->s[i].gevid);
        }
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
    int fd;
    int64_t l;
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
    while ((l = kseq_read(ks)) >= 0) {
        if (l > UINT32_MAX) {
            fprintf(stderr, "[E::%s] >4G sequence chunks are not supported: %s [%ld]\n", __func__, ks->name.s, l);
            exit(EXIT_FAILURE);
        }
        if (strlen(ks->seq.s) >= min_len)
            sd_put1(d, ks->name.s, ks->seq.s, strlen(ks->seq.s));
    }

    kseq_destroy(ks);
    gzclose(fp);
    kclose(ko);

    return d;
}

static int is_empty_line(char *line)
{
    char *c;
    c = line;
    while (isspace((unsigned char) *c))
        c++;
    if(*c == 0)
        return 1;
    return 0;
}

sdict_t *make_sdict_from_index(const char *f, uint32_t min_len)
{
    FILE *fp;
    char *line = NULL;
    size_t ln = 0;
    ssize_t read;
    char name[4096];
    int64_t len;

    fp = fopen(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }
    
    sdict_t *d;
    d = sd_init();
    while ((read = getline(&line, &ln, fp)) != -1) {
        if (is_empty_line(line))
            continue;
        sscanf(line, "%s %ld", name, &len);
        if (len > UINT32_MAX) {
            fprintf(stderr, "[E::%s] >4G sequence chunks are not supported: %s [%ld]\n", __func__, name, len);
            exit(EXIT_FAILURE);
        }
        if (len >= min_len)
            sd_put(d, name, len);
    }
    if (line)
        free(line);
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
    uint64_t len;

    fp = fopen(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    sdict_t *d;
    d = sd_init();
    while ((read = getline(&line, &ln, fp)) != -1) {
        if (is_empty_line(line))
            continue;
        if (line[0] == 'S') {
            sscanf(line, "%*s %s %*s %s", name, lens);
            len = strtoul(lens + 5, NULL, 10);
            if (len > UINT32_MAX) {
                fprintf(stderr, "[E::%s] >4G sequence chunks are not supported: %s [%ld]\n", __func__, name, len);
                exit(EXIT_FAILURE);
            }
            if (len >= min_len)
                sd_put(d, name, len);
        }
    }
    if (line)
        free(line);
    fclose(fp);

    return d;
}

static inline void *memcpy1(const void *src, size_t count)
{
    if (count == 0) return 0;
    void *dest = malloc(count);
    memcpy(dest, src, count);
    return dest;
}

uint32_t asm_put(asm_dict_t *d, const char *name, uint64_t len, uint64_t gap, uint32_t n, uint32_t s, uint32_t gn, 
        uint32_t *gs, uint32_t *gsize, uint32_t *gcomp, uint32_t *gtype, uint32_t *glink, uint32_t *gevid)
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
        a->gap = gap;
        a->n = n;
        a->s = s;
        a->gn = gn;
        a->gs = memcpy1(gs, sizeof(uint32_t) * gn);
        a->gsize = memcpy1(gsize, sizeof(uint32_t) * gn);
        a->gcomp = memcpy1(gcomp, sizeof(uint32_t) * gn);
        a->gtype = memcpy1(gtype, sizeof(uint32_t) * gn);
        a->glink = memcpy1(glink, sizeof(uint32_t) * gn);
        a->gevid = memcpy1(gevid, sizeof(uint32_t) * gn);
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

void seg_put(asm_dict_t *d, uint32_t s, uint32_t k, uint64_t a, uint64_t b,
        uint32_t c, uint32_t x, uint32_t y, uint32_t t, uint32_t r)
{
    if (d->u == d->v) {
        d->v = d->v? d->v<<1 : 16;
        d->seg = (sd_seg_t *) realloc(d->seg, d->v * sizeof(sd_seg_t));
    }
    sd_seg_t seg = {s, k, a, b, c, x, y, t, r};
    d->seg[d->u++] = seg;
}

asm_dict_t *make_asm_dict_from_sdict(sdict_t *sdict)
{
    uint32_t i;
    asm_dict_t *d;
    d = asm_init(sdict);
    d->index = (uint64_t *) calloc(sdict->n, sizeof(uint64_t));
    uint64_t *a = d->a;
    for (i = 0; i < sdict->n; ++i) {
        asm_put(d, sdict->s[i].name, sdict->s[i].len, 0, 1, i, 0, 0, 0, 0, 0, 0, 0);
        seg_put(d, i, 0, 0, 0, i<<1, 0, sdict->s[i].len, DEFAULT_AGP_SEQ_COMPONENT_TYPE, AGP_OT_PLUS);
        a[i] = (uint64_t) i << 32 | 1;
        d->index[i] = (uint64_t) sdict->s[i].len << 32 | i;
    }
    return d;
}

typedef struct {uint64_t x, y;} c_pair_t;

#define pair_key(a) ((a).x)
KRADIX_SORT_INIT(pair, c_pair_t, pair_key, 8)

void asm_index(asm_dict_t *d)
{
    int32_t i, s, c;
    uint64_t *a, last;
    sd_seg_t seg;

    s = d->u;
    if (s == 0) return;

    c_pair_t *c_pairs;
    c_pairs = (c_pair_t *) malloc(s * sizeof(c_pair_t));
    for(i = 0; i < s; ++i) {
        seg = d->seg[i];
        c_pairs[i].x = (uint64_t) seg.c >> 1 << 32 | (seg.x + seg.y);
        c_pairs[i].y = (uint64_t) (seg.x + seg.y) << 32 | i;
    }

    radix_sort_pair(c_pairs, c_pairs + s);
    
    if (d->index)
        free(d->index);
    d->index = (uint64_t *) malloc(s * sizeof(uint64_t));
    a = d->a;
    memset(a, 0, sizeof(uint64_t) * d->sdict->n);
    d->index[0] = c_pairs[0].y;
    c = c_pairs[0].x >> 32;
    for (i = 1, last = 0; i < s; ++i) {
        d->index[i] = c_pairs[i].y;
        if (c != c_pairs[i].x >> 32) {
            a[c] = last << 32 | (i - last);
            last = i;
            c = c_pairs[i].x >> 32;
        }
    }
    a[c] = last << 32 | (i - last);

    free(c_pairs);
}

asm_dict_t *make_asm_dict_from_agp(sdict_t *sdict, const char *f, int allow_unknown_oris)
{
    FILE *fp;
    char *line = NULL;
    size_t ln = 0;
    ssize_t read;
    char sname[256], type[4], cname[256], cstarts[16], cends[16], oris[64];
    char *name = NULL;
    uint64_t a, b;
    uint32_t l, cstart, cend, ctype, coris;
    uint32_t c, s, n;
    kvec_t(uint32_t) gs, gsize, gcomp, gtype, glink, gevid;

    fp = fopen(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    asm_dict_t *d;
    sd_seq_t *sd;
    d = asm_init(sdict);
    a = b = 0;
    s = n = 0;
    kv_init(gs);
    kv_init(gsize);
    kv_init(gcomp);
    kv_init(gtype);
    kv_init(glink);
    kv_init(gevid);
    while ((read = getline(&line, &ln, fp)) != -1) {
        if (is_empty_line(line) || !strncmp(line, "#", 1))
            // header or empty lines
            continue;
        
        sscanf(line, "%s %*s %*s %*s %s %s %s %s %s", sname, type, cname, cstarts, cends, oris);
        if (!name)
            name = strdup(sname);
        if (strcmp(sname, name)) {
            asm_put(d, name, a, b - a, n, s - n, gs.n, gs.a, gsize.a, gcomp.a, gtype.a, glink.a, gevid.a);
            a = 0;
            b = 0;
            n = 0;
            gs.n = 0;
            gsize.n = 0;
            gtype.n = 0;
            glink.n = 0;
            gevid.n = 0;
            free(name);
            name = strdup(sname);
        }

        ctype = agp_component_type_key(type);
        if (ctype == AGP_CT_N || ctype == AGP_CT_U) {
            // gap
            kv_push(uint32_t, gs, n);
            kv_push(uint32_t, gcomp, ctype);
            c = strtoul(cname, NULL, 10); kv_push(uint32_t, gsize, c);
            if (ctype == AGP_CT_U && c != DEFAULT_AGP_U_GAP_SIZE)
                fprintf(stderr, "[W::%s] type 'U' gap size is not %d\n", __func__, DEFAULT_AGP_U_GAP_SIZE);
            b += c;
            c = agp_gap_type_key(cstarts); kv_push(uint32_t, gtype, c);
            c = agp_linkage_key(cends); kv_push(uint32_t, glink, c);
            c = agp_linkage_evidence_key(oris); kv_push(uint32_t, gevid, c);
            continue;
        }

        cstart = strtoul(cstarts, NULL, 10);
        cend = strtoul(cends, NULL, 10);
        c = sd_get(sdict, cname);
        if (c == UINT32_MAX) {
            fprintf(stderr, "[E::%s] sequence %s not found\n", __func__, cname);
            return 0;
        }
        sd = &sdict->s[c];
        if (cstart < 1 || cstart > cend || cend > sd->len) {
            fprintf(stderr, "[E::%s] invalid sequence component %s:%u-%u on %s\n", __func__, cname, cstart, cend, sname);
            if (cend > sd->len)
                fprintf(stderr, "[E::%s] sequence end position (%u) greater than sequence length (%u)\n", __func__, cend, sd->len);
            exit(EXIT_FAILURE);
        }
        coris = agp_orientation_key(oris);
        if (!allow_unknown_oris && coris != AGP_OT_PLUS && coris != AGP_OT_MINUS) {
            // unknown orientation is not allowed
            fprintf(stderr, "[E::%s] unknown orientation of sequence component %s:%u-%u on %s: '%s'\n", __func__, cname, cstart, cend, sname, oris);
            if (allow_unknown_oris) {
                fprintf(stderr, "[E::%s] valid identifiers for unorientated sequence include: '?', '0' and 'na'\n", __func__);
                fprintf(stderr, "[E::%s] see %s\n", __func__, AGP_SPEC_ONLINE_DOC);
            } else {
                fprintf(stderr, "[E::%s] run program with '-u' option to include sequence components\n", __func__);
            }
            exit(EXIT_FAILURE);
        }
        c <<= 1;
        if (coris == AGP_OT_MINUS) c |= 1;
        l = cend - cstart + 1;
        seg_put(d, d->n, n, a, b, c, cstart - 1, l, ctype, coris);
        a += l;
        b += l;
        ++s;
        ++n;
    }
    asm_put(d, name, a, b - a, n, s - n, gs.n, gs.a, gsize.a, gcomp.a, gtype.a, glink.a, gevid.a);
    d->sdict = sdict;
#ifdef DEBUG_DICT
    for (int i = 0; i < s; ++i)
        fprintf(stderr, "[DEBUG_DICT::%s] %u %lu %u %u %u\n", __func__, d->seg[i].s, d->seg[i].a, d->seg[i].c, d->seg[i].x, d->seg[i].y);
    for (int i = 0; i < d->n; ++i)
        fprintf(stderr, "[DEBUG_DICT::%s] %s %lu %u %u\n", __func__, d->s[i].name, d->s[i].len, d->s[i].n, d->s[i].s);
#endif
    asm_index(d);

    fclose(fp);
    if (line)
        free(line);
    if (name)
        free(name);
    kv_destroy(gs);
    kv_destroy(gsize);
    kv_destroy(gcomp);
    kv_destroy(gtype);
    kv_destroy(glink);
    kv_destroy(gevid);

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
        seg_put(d, d->n, 0, 0, 0, i<<1, 0, sdict->s[i].len, DEFAULT_AGP_SEQ_COMPONENT_TYPE, AGP_OT_PLUS);
        asm_put(d, sdict->s[i].name, sdict->s[i].len, 0, 1, d->u - 1, 0, 0, 0, 0, 0, 0, 0);
    }
    asm_index(d);
}

/* contig coordinates to scaffold coordinates */
// 0-based
int sd_coordinate_conversion(asm_dict_t *d, uint32_t id, uint32_t pos, uint32_t *s, uint64_t *p, int count_gap)
{
    *s = UINT32_MAX;
    *p = UINT64_MAX;

    if (id >= d->sdict->n)
        return SEQ_NOT_FOUND;

    uint32_t i, n;
    uint64_t *index = d->index;
    i = (uint32_t) (d->a[id] >> 32);
    n = (uint32_t) (d->a[id]) + i;

    while (i < n && index[i]>>32 <= pos)
        ++i;
    if (i == n)
        return POS_NOT_IN_RANGE;

    sd_seg_t *seg = &d->seg[(uint32_t) index[i]];
    if (pos < seg->x || pos >= seg->x + seg->y)
        return POS_NOT_IN_RANGE;

    uint64_t offset = count_gap? seg->b : seg->a;
    *s = seg->s;
    *p = seg->c & 1? offset + seg->x + seg->y - 1 - pos : offset + pos - seg->x;

    return CC_SUCCESS;
}

/* scaffold coordinates to contig coordinates */
// 0-based
int sd_coordinate_rev_conversion(asm_dict_t *d, uint32_t id, uint64_t pos, uint32_t *s, uint32_t *p, int count_gap)
{
    *s = UINT32_MAX;
    *p = UINT32_MAX;

    if (id >= d->n)
        return SEQ_NOT_FOUND;

    sd_seg_t *seg = d->seg + d->s[id].s;
    uint32_t i = 0, n = d->s[id].n;
    uint64_t offset = 0;

    if (count_gap) {
        while (i < n && pos >= seg[i].b)
            ++i;
        seg = &seg[i-1];
        offset = pos - seg->b;
    } else {
        while (i < n && pos >= seg[i].a)
            ++i;
        seg = &seg[i-1];
        offset = pos - seg->a;
    }
    
    // check if in gap or exceeds seq length
    if (offset >= seg->y)
        return POS_NOT_IN_RANGE;

    *s = seg->c >> 1;
    *p = seg->c & 1? seg->y - 1 - offset : seg->x + offset;

    return CC_SUCCESS;
}

int cmp_uint64_d (const void *a, const void *b) {
    // decreasing order
    uint64_t x, y;
    x = *(uint64_t *) a;
    y = *(uint64_t *) b;
    return (x < y) - (x > y);
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


// deprecated
// this can now be done in two steps: first read AGP file to an assembly dictionary and then write to a FASTA file
void _write_fasta_file_from_agp_(const char *fa, const char *agp, FILE *fo, int line_wd, int allow_unknown_oris)
{
    FILE *agp_in;
    char *line = NULL;
    sdict_t *dict;
    sd_seq_t s;
    size_t ln = 0;
    ssize_t read;
    char sname[256], type[4], cname[256], cstarts[16], cends[16], oris[256];
    char *name = NULL;
    uint32_t c;
    int64_t i, l, cstart, cend, ns, L;

    agp_in = fopen(agp, "r");
    if (agp_in == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, agp);
        exit(EXIT_FAILURE);
    }

    dict = make_sdict_from_fa(fa, 0);
    l = L = ns = 0;
    while ((read = getline(&line, &ln, agp_in)) != -1) {
        if (is_empty_line(line) || !strncmp(line, "#", 1))
            // header or empty lines
            continue;
        sname[0] = type[0] = cname[0] = cstarts[0] = cends[0] = oris[0] = '\0';
        sscanf(line, "%s %*s %*s %*s %s %s %s %s %s", sname, type, cname, cstarts, cends, oris);
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
            ++ns;
            L += l;
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
            fprintf(stderr, "[E::%s] sequence component '%s' not found\n", __func__, cname);
            exit(EXIT_FAILURE);
        }
        s = dict->s[c];
        if (cstart < 1 || cstart > cend || cend > s.len) {
            fprintf(stderr, "[E::%s] invalid sequence component %s:%ld-%ld on %s\n", __func__, cname, cstart, cend, sname);
            if (cend > s.len)
                fprintf(stderr, "[E::%s] sequence end position (%ld) greater than sequence length (%u)\n", __func__, cend, s.len);
            exit(EXIT_FAILURE);
        }
        if (!strncmp(oris, "+", 1) || (allow_unknown_oris && (!strncmp(oris, "?", 1) || !strncmp(oris, "0", 1) || !strncmp(oris, "na", 2)))) {
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
            fprintf(stderr, "[E::%s] unknown orientation of sequence component %s:%ld-%ld on %s: '%s'\n", __func__, cname, cstart, cend, sname, oris);
            if (allow_unknown_oris) {
                fprintf(stderr, "[E::%s] valid identifiers for unorientated sequence include: '?', '0' and 'na'\n", __func__);
                fprintf(stderr, "[E::%s] see %s\n", __func__, AGP_SPEC_ONLINE_DOC);
            }
            exit(EXIT_FAILURE);
        }
    }
    if (l % line_wd != 0)
        fputc('\n', fo);
    ++ns;
    L += l;
    fprintf(stderr, "[I::%s] Number sequences: %ld\n", __func__, ns);
    fprintf(stderr, "[I::%s] Number bases: %ld\n", __func__, L);

    if (line)
        free(line);
    free(name);
    fclose(agp_in);
    sd_destroy(dict);
}

static inline void write_agp_seq(FILE *fo, char *s_name, uint64_t s_beg, uint64_t s_end, uint32_t b, 
        uint32_t ctype, char *c_name, uint32_t c_beg, uint32_t c_end, uint32_t c_oris)
{
    fprintf(fo, "%s\t%lu\t%lu\t%d\t%s\t%s\t%u\t%u\t%s\n", s_name, s_beg, s_end, b, 
            agp_component_type_val(ctype), c_name, c_beg, c_end, agp_orientation_val(c_oris));
}

static inline void write_agp_gap(FILE *fo, char *s_name, uint64_t s_beg, uint64_t s_end, uint32_t b, 
        uint32_t gcomp, uint32_t gsize, uint32_t gtype, uint32_t glink, uint32_t gevid)
{
    fprintf(fo, "%s\t%lu\t%lu\t%d\t%s\t%d\t%s\t%s\t%s\n", s_name, s_beg, s_end, b, 
            agp_component_type_val(gcomp), gsize, agp_gap_type_val(gtype), 
            agp_linkage_val(glink), agp_linkage_evidence_val(gevid));
}

void write_segs_to_agp(sd_seg_t *segs, uint32_t n, sd_aseq_t *aseq, sdict_t *sd, char *name, FILE *fo)
{
    uint64_t len;
    sd_seg_t seg;
    uint32_t i, j, t;
    len = 0;
    t = 0;
    if (!name) name = aseq? aseq->name : 0;
    if (!name) {
        fprintf(stderr, "[E::%s] sequence name required\n", __func__);
        exit(EXIT_FAILURE);
    } 
    if (aseq) {
        // with gap information
        for (i = 0, j = 0; i < n; ++i) {
            seg = segs[i];
            while (j < aseq->gn && aseq->gs[j] <= i) {
                // print gap
                write_agp_gap(fo, name, len + 1, len + aseq->gsize[j], ++t, aseq->gcomp[j], aseq->gsize[j], 
                        aseq->gtype[j], aseq->glink[j], aseq->gevid[j]);
                len += aseq->gsize[j];
                ++j;
            }
            write_agp_seq(fo, name, seg.b + 1, seg.b + seg.y, ++t, seg.t, sd->s[seg.c >> 1].name, 
                    seg.x + 1, seg.x + seg.y, seg.r);
            len += seg.y;
        }
        while (j < aseq->gn) {
            // print trailing gap
            write_agp_gap(fo, name, len + 1, len + aseq->gsize[j], ++t, aseq->gcomp[j], aseq->gsize[j], 
                    aseq->gtype[j], aseq->glink[j], aseq->gevid[j]);
            len += aseq->gsize[j];
            ++j;
        }
        if (len != aseq->len + aseq->gap) {
            fprintf(stderr, "[E::%s] sequence length conflict: %lu %lu\n", __func__, len, aseq->len + aseq->gap);
            //exit(EXIT_FAILURE);
        }
    } else {
        // no gap information provided
        // using default settings
        for (i = 0; i < n; ++i) {
            seg = segs[i];
            write_agp_seq(fo, name, len + 1, len + seg.y, ++t, seg.t, sd->s[seg.c >> 1].name, seg.x + 1, seg.x + seg.y, seg.r);
            /***
            fprintf(fo, "%s\t%lu\t%lu\t%d\t%s\t%s\t%u\t%u\t%s\n", 
                    name, len + 1, len + seg.y, ++t,
                    agp_component_type_val(DEFAULT_AGP_SEQ_COMPONENT_TYPE),
                    sd->s[seg.c >> 1].name, seg.x + 1, seg.x + seg.y,
                    seg.c & 1? agp_orientation_val(AGP_OT_MINUS) : agp_orientation_val(AGP_OT_PLUS));
            **/
            len += seg.y;
            if (i != n - 1) {
                write_agp_gap(fo, name, len + 1, len + DEFAULT_AGP_GAP_SIZE, ++t, DEFAULT_AGP_GAP_COMPONENT_TYPE, 
                        DEFAULT_AGP_GAP_SIZE, DEFAULT_AGP_GAP_TYPE, AGP_LG_YES, DEFAULT_AGP_LINKAGE_EVIDENCE);
                /***
                fprintf(fo, "%s\t%lu\t%lu\t%d\t%s\t%d\t%s\t%s\t%s\n", 
                        name, len + 1, len + DEFAULT_AGP_GAP_SIZE, ++t,
                        agp_component_type_val(DEFAULT_AGP_GAP_COMPONENT_TYPE), 
                        DEFAULT_AGP_GAP_SIZE,
                        agp_gap_type_val(DEFAULT_AGP_GAP_TYPE),
                        agp_linkage_val(AGP_LG_YES),
                        agp_linkage_evidence_val(DEFAULT_AGP_LINKAGE_EVIDENCE));
                **/
                len += DEFAULT_AGP_GAP_SIZE;
            }
        }
    }
}

void write_asm_dict_to_agp(asm_dict_t *dict, FILE *fo)
{
    uint32_t i;
    for (i = 0; i < dict->n; ++i) {
        write_segs_to_agp(dict->seg + dict->s[i].s, dict->s[i].n, &dict->s[i], dict->sdict, 0, fo);
    }
}

void write_sorted_agp(asm_dict_t *dict, FILE *fo)
{
    int i, j;
    uint32_t s;
    char s_name[32];
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
        sprintf(s_name, "scaffold_%u", ++s);
        write_segs_to_agp(dict->seg + dict->s[j].s, dict->s[j].n, &dict->s[j], dict->sdict, s_name, fo);
    }
    free(c_pairs);
}

void write_sdict_to_agp(sdict_t *sdict, FILE *fo)
{
    int i;
    for (i = 0; i < sdict->n; ++i)
        fprintf(fo, "scaffold_%u\t1\t%u\t1\t%s\t%s\t1\t%u\t+\n", i + 1, sdict->s[i].len, 
                agp_component_type_val(DEFAULT_AGP_SEQ_COMPONENT_TYPE), sdict->s[i].name, sdict->s[i].len);
}

uint64_t write_segs_to_fasta(sd_seg_t *segs, uint32_t n, sd_aseq_t *aseq, sdict_t *sd, char *name, int line_wd, FILE *fo)
{
    uint32_t i, j;
    uint64_t l, len;
    int64_t k, m;
    sd_seg_t seg;
    char *s;
    if (!name) name = aseq? aseq->name : 0;
    if (!name) {
        fprintf(stderr, "[E::%s] sequence name required\n", __func__);
        exit(EXIT_FAILURE);
    }
    fprintf(fo, ">%s\n", name);
    l = len = 0;
    if (aseq) {
        // with gap information
        for (i = 0, j = 0; i < n; ++i) {
            while (j < aseq->gn && aseq->gs[j] <= i) {
                // print gap
                for (k = 0, m = aseq->gsize[j]; k < m; ++k) {
                    fputc('N', fo);
                    if (++l % line_wd == 0)
                        fputc('\n', fo);
                }
                len += m;
                ++j;
            }
            seg = segs[i];
            s = sd->s[seg.c >> 1].seq;
            if ((seg.c&1) == 0) {
                // forward
                for (k = seg.x, m = seg.x + seg.y; k < m; ++k) {
                    fputc(s[k], fo);
                    if (++l % line_wd == 0)
                        fputc('\n', fo);
                }
            } else {
                // reverse
                for (k = seg.x + seg.y - 1, m = seg.x; k >= m; --k) {
                    fputc(comp_table[(int) s[k]], fo);
                    if (++l % line_wd == 0)
                        fputc('\n', fo);
                }
            }
            len += seg.y;
        }
        while (j < aseq->gn) {
            // print trailing gap
            for (k = 0, m = aseq->gsize[j]; k < m; ++k) {
                fputc('N', fo);
                if (++l % line_wd == 0)
                    fputc('\n', fo);
            }
            len += m;
            ++j;
        }
        if (len != aseq->len + aseq->gap) {
            fprintf(stderr, "[E::%s] sequence length conflict\n", __func__);
            exit(EXIT_FAILURE);
        }
    } else {
        // no gap information provided
        // using default settings
        for (i = 0; i < n; ++i) {
            seg = segs[i];
            s = sd->s[seg.c >> 1].seq;
            if ((seg.c&1) == 0) {
                // forward
                for (k = seg.x, m = seg.x + seg.y; k < m; ++k) {
                    fputc(s[k], fo);
                    if (++l % line_wd == 0)
                        fputc('\n', fo);
                }
            } else {
                // reverse
                for (k = seg.x + seg.y - 1, m = seg.x; k >= m; --k) {
                    fputc(comp_table[(int) s[k]], fo);
                    if (++l % line_wd == 0)
                        fputc('\n', fo);
                }
            }
            len += seg.y;
            
            if (i != n - 1) {
                for (k = 0, m = DEFAULT_AGP_GAP_SIZE; k < m; ++k) {
                    fputc('N', fo);
                    if (++l % line_wd == 0)
                        fputc('\n', fo);
                }
                len += m;
            }
        }
    }

    if (l % line_wd != 0)
        fputc('\n', fo);
    
    return len;
}

uint64_t write_asm_dict_to_fasta(asm_dict_t *dict, int line_wd, FILE *fo)
{
    uint32_t i;
    uint64_t len = 0;
    for (i = 0; i < dict->n; ++i)
        len += write_segs_to_fasta(dict->seg + dict->s[i].s, dict->s[i].n, 
                &dict->s[i], dict->sdict, 0, line_wd, fo);
    return len;
}

void write_fasta_file_from_agp(const char *fa, const char *agp, FILE *fo, int line_wd, int allow_unknown_oris)
{
    sdict_t *sdict = make_sdict_from_fa(fa, 0);
    asm_dict_t *dict = make_asm_dict_from_agp(sdict, agp, allow_unknown_oris);
    uint64_t len = write_asm_dict_to_fasta(dict, line_wd, fo);
    fprintf(stderr, "[I::%s] Number sequences: %u\n", __func__, dict->n);
    fprintf(stderr, "[I::%s] Number bases: %lu\n", __func__, len);
    asm_destroy(dict);
    sd_destroy(sdict);
}

uint64_t write_sequence_dictionary(FILE *fo, sdict_t *dict)
{
    uint32_t i;
    uint64_t l;
    // write number of sequences
    fwrite(&dict->n, sizeof(uint32_t), 1, fo);
    // write sequence lengths
    for (i = 0; i < dict->n; ++i)
        fwrite(&dict->s[i].len, sizeof(uint32_t), 1, fo);
    l = dict->n;
    for (i = 0; i < dict->n; ++i)
        l += strlen(dict->s[i].name);
    // write total length of sequence names
    fwrite(&l, sizeof(uint64_t), 1, fo);
    // write sequence name with a appending white space
    const unsigned char space = ' ';
    for (i = 0; i < dict->n; ++i) {
        fwrite(dict->s[i].name, 1, strlen(dict->s[i].name), fo);
        fwrite(&space, 1, 1, fo);
    }

    return sizeof(uint32_t) * (dict->n + 3) + l;
}

void file_seek_skip_sdict(FILE *fp)
{
    uint32_t n;
    uint64_t l;
    fread(&n, sizeof(uint32_t), 1, fp);
    fseek(fp, sizeof(uint32_t) * n, SEEK_CUR);
    fread(&l, sizeof(uint64_t), 1, fp);
    fseek(fp, l, SEEK_CUR);
}

int file_sdict_match(char *f, sdict_t *dict)
{
    uint32_t i, m, n, s, *lens;
    uint64_t l;
    int64_t magic_number;
    char *names, *p;
    FILE *fp;

    fp = fopen(f, "r");
    if (fp == NULL)
        return 1;

    m = fread(&magic_number, sizeof(int64_t), 1, fp);
    if (!m || !is_valid_bin_header(magic_number)) 
        return 2;

    fread(&n, sizeof(uint32_t), 1, fp);
    if (n != dict->n)
        return 3;
    
    lens = (uint32_t *) calloc(n, sizeof(uint32_t));
    fread(lens, sizeof(uint32_t), n, fp);
    for (i = 0; i < n; ++i) {
        if (dict->s[i].len != lens[i]) {
            free(lens);
            return 4;
        }
    }
    free(lens);
    
    fread(&l, sizeof(uint64_t), 1, fp);
    names = (char *) calloc(l, 1);
    m = fread(names, 1, l, fp);
    for (p = names, i = 0; i < n; ++i) {
        s = strlen(dict->s[i].name);
        if (p - names >= l || strncmp(p, dict->s[i].name, s)) {
            free(names);
            return 5;
        }
        p += s + 1;
    }
    free(names);

    fclose(fp);

    return 0;
}

