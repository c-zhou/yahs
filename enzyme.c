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
 * 14/12/21 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <ctype.h>
#include <float.h>

#include "kvec.h"
#include "enzyme.h"
#include "sdict.h"
#include "asset.h"

#undef DEBUG_ENZ

static double MIN_RE_DENS = .1;
static double MAX_RE_DENS = DBL_MAX;

re_cuts_t *re_cuts_init(uint32_t n)
{
    re_cuts_t *re_cuts = (re_cuts_t *) malloc(sizeof(re_cuts_t));
    re_cuts->n = n;
    re_cuts->density = 0;
    re_cuts->re = (re_t *) malloc(n * sizeof(re_t));
    uint32_t i;
    for (i = 0; i < n; ++i)
        re_cuts->re[i].sites = 0;
    return re_cuts;
}

void re_cuts_destroy(re_cuts_t *re_cuts)
{
    uint32_t i;
    for (i = 0; i < re_cuts->n; ++i)
        if (re_cuts->re[i].sites)
            free(re_cuts->re[i].sites);
    free(re_cuts->re);
    free(re_cuts);
}

typedef struct {size_t n, m; uint32_t *a;} u32_v;
int u32_cmp(const void *p, const void *q)
{
    uint32_t a, b;
    a = *(uint32_t *) p;
    b = *(uint32_t *) q;
    return a == b? 0 : (a > b? 1 : -1);
}

re_cuts_t *find_re_from_seqs(const char *f, uint32_t ml, char **enz_cs, int enz_n)
{
    // now find all RE cutting sites
    int i, j, c, c0, c1;
    uint32_t len, n, p;
    int64_t n_re, genome_size;
    char *pch, *seq, *re;
    sdict_t *sdict;
    re_cuts_t *re_cuts;

    sdict = make_sdict_from_fa(f, ml);
    n = sdict->n;
    n_re = genome_size = 0;
    re_cuts = re_cuts_init(n);
    for (i = 0; i < n; ++i) {
        u32_v enz_cs_pos = {0, 0, 0};

        seq = sdict->s[i].seq;
        len = sdict->s[i].len;
        for (p = 0; p < len; ++p) {
            c = seq[p];
            if (!isalpha(c)) {
                fprintf(stderr, "[E::%s] non-alphabetic chacrater in FASTA file: %c\n", __func__, c);
                re_cuts_destroy(re_cuts);
                sd_destroy(sdict);
                return 0;
            }
            seq[p] = nucl_toupper[c];
        }
        
        for (j = 0; j < enz_n; ++j) {
            re = enz_cs[j];
            pch = strstr(seq, re);
            while (pch != NULL) {
                kv_push(uint32_t, enz_cs_pos, pch - seq);
                pch = strstr(pch + 1, re);
            }
        }

        // reverse complement sequence
        for (p = 0; p < len>>1; ++p) {
            c0 = comp_table[(int) seq[p]];
            c1 = comp_table[(int) seq[len - 1 - p]];
            seq[p] = c1;
            seq[len - 1 - p] = c0;
        }
        if (len & 1) // complement the remaining base
            seq[len>>1] = comp_table[(int) seq[len>>1]];
        
        for (j = 0; j < enz_n; ++j) {
            re = enz_cs[j];
            pch = strstr(seq, re);
            while (pch != NULL) {
                kv_push(uint32_t, enz_cs_pos, len - 1 - (pch - seq));
                pch = strstr(pch + 1, re);
            }
        }

        qsort(enz_cs_pos.a, enz_cs_pos.n, sizeof(uint32_t), u32_cmp);
        re_cuts->re[i].sites = enz_cs_pos.a;
        re_cuts->re[i].n = enz_cs_pos.n;
        re_cuts->re[i].l = len;
        
        n_re += enz_cs_pos.n;
        genome_size += len;
    }

    re_cuts->density = (double) n_re / genome_size;

    fprintf(stderr, "[I::%s] NO. restriction enzyme cutting sites found in sequences: %ld\n", __func__, n_re);
    fprintf(stderr, "[I::%s] restriction enzyme cutting sites density: %.6f\n", __func__, re_cuts->density);
#ifdef DEBUG_ENZ
    fprintf(stderr, "[DEBUG_ENZ::%s] restriction enzyme cutting sites for individual sequences (n = %d)\n", __func__, n);
    for (i = 0; i < n; ++i)
        fprintf(stderr, "[DEBUG_ENZ::%s] %s %u %u %.6f\n", __func__, sdict->s[i].name, sdict->s[i].len, re_cuts->re[i].n, (double) re_cuts->re[i].n / sdict->s[i].len);
#endif

    sd_destroy(sdict);

    return re_cuts;
}

double **calc_re_cuts_density(re_cuts_t *re_cuts, uint32_t resolution)
{
    if (!re_cuts)
        return 0;

    uint32_t i, j, b;
    double **dens, *ds;
    re_t re;
    dens = (double **) malloc(re_cuts->n * sizeof(double *));
    for (i = 0; i < re_cuts->n; ++i) {
        re = re_cuts->re[i];
        b = div_ceil(re.l, resolution);
        ds = (double *) calloc(b, sizeof(double));
        for (j = 0; j < re.n; ++j)
            ds[(MAX(re.sites[j], 1) - 1) / resolution] += 1.;
        for (j = 0; j < b - 1; ++j)
            ds[j] /= (double) resolution * re_cuts->density;
        ds[b - 1] /= ((double) re.l - (double) (b - 1) * resolution) * re_cuts->density;
        
        for (j = 0; j < b; ++j)
            if (ds[j] < MIN_RE_DENS || ds[j] > MAX_RE_DENS)
                ds[j] = .0;
        
        dens[i] = ds;
#ifdef DEBUG_ENZ
        fprintf(stderr, "[DEBUG_ENZ::%s] DENS [%u/%u] (%u):", __func__, i, re_cuts->n, b);
        for (j = 0; j < b; ++j)
            fprintf(stderr, " %.6f", ds[j]);
        fprintf(stderr, "\n");
#endif
    }
    return dens;
}

static uint32_t bin_search(uint32_t *a, uint32_t n, uint32_t s)
{
    uint32_t low, high, mid;
    low = 0;
    high = n;
    while (low != high) {
        mid = (low >> 1) + (high >> 1);
        if (a[mid] < s)
            low = mid + 1;
        else
            high = mid;
    }
    return low;
} 

double **calc_re_cuts_density1(re_cuts_t *re_cuts, uint32_t resolution, asm_dict_t *dict)
{
    if (!re_cuts)
        return 0;

    uint32_t i, j, b, n, e, a;
    double **dens, *ds;
    re_t re;
    sd_aseq_t seq;
    sd_seg_t seg;
    
    n = dict->n;
    dens = (double **) malloc(n * sizeof(double *));
    for (i = 0; i < n; ++i) {
        seq = dict->s[i];
        b = div_ceil(seq.len, resolution);
        ds = (double *) calloc(b, sizeof(double));
        for (j = 0; j < seq.n; ++j) {
            seg = dict->seg[seq.s + j];
            re = re_cuts->re[seg.c >> 1]; // get seq re cuts
            e = seg.x + seg.y; // seq end position, exclusive
            // binary search to get starting re cuts
            // first cutting site no smaller than seg.x
            // not many seq breaks - linear search might be faster on average?
            a = bin_search(re.sites, re.n, seg.x);
            
            if (seg.c & 1) {
                // reverse complement
                while (a < re.n && re.sites[a] < e) {
                    ds[(MAX(seg.a + seg.y - (re.sites[a] - seg.x), 1) - 1) / resolution] += 1.;
                    ++a;
                }
            } else {
                while (a < re.n && re.sites[a] < e) {
                    ds[(MAX(seg.a + re.sites[a] - seg.x, 1) - 1) / resolution] += 1.;
                    ++a;
                }
            }
        }
        for (j = 0; j < b - 1; ++j)
            ds[j] /= (double) resolution * re_cuts->density;
        ds[b - 1] /= ((double) seq.len - (double) (b - 1) * resolution) * re_cuts->density;
        
        for (j = 0; j < b; ++j)
            if (ds[j] < MIN_RE_DENS || ds[j] > MAX_RE_DENS)
                ds[j] = .0;
        
        dens[i] = ds;
#ifdef DEBUG_ENZ
        fprintf(stderr, "[DEBUG_ENZ::%s] DENS1 [%u/%u] (%u):", __func__, i, n, b);
        for (j = 0; j < b; ++j)
            fprintf(stderr, " %.6f", ds[j]);
        fprintf(stderr, "\n");
#endif
    }

    return dens;
}

double **calc_re_cuts_density2(re_cuts_t *re_cuts, uint32_t resolution, asm_dict_t *dict)
{
    if (!re_cuts)
        return 0;

    uint32_t i, j, b, n, e, a;
    uint64_t l, p;
    double **dens, *ds;
    re_t re;
    sd_aseq_t seq;
    sd_seg_t seg;
    n = dict->n;
    dens = (double **) malloc(n * sizeof(double *));
    for (i = 0; i < n; ++i) {
        seq = dict->s[i];
        l = div_ceil(seq.len, 2); // split sequence into two parts
        b = div_ceil(l, resolution);
        ds = (double *) calloc(b << 1, sizeof(double));
        for (j = 0; j < seq.n; ++j) {
            seg = dict->seg[seq.s + j];
            re = re_cuts->re[seg.c >> 1]; // get seq re cuts
            e = seg.x + seg.y; // seq end position, exclusive
            // binary search to get starting re cuts
            // first cutting site no smaller than seg.x
            // not many seq breaks - linear search might be faster on average?
            a = bin_search(re.sites, re.n, seg.x);
            if (seg.c & 1) {
                // reverse complement
                while (a < re.n && re.sites[a] < e) {
                    p = seg.a + seg.y - (re.sites[a] - seg.x); // position on seq
                    if (p < l)
                        ds[(MAX(p, 1) - 1) / resolution << 1] += 1.;
                    else
                        ds[(MAX(seq.len - p, 1) - 1) / resolution << 1 | 1] += 1.;
                    ++a;
                }
            } else {
                while (a < re.n && re.sites[a] < e) {
                    p = seg.a + re.sites[a] - seg.x; // position on seq
                    if (p < l)
                        ds[(MAX(p, 1) - 1) / resolution << 1] += 1.;
                    else
                        ds[(MAX(seq.len - p, 1) - 1) / resolution << 1 | 1] += 1.;
                    ++a;
                }
            }
        }
        for (j = 0; j < (b - 1) << 1; ++j)
            ds[j] /= (double) resolution * re_cuts->density;
        ds[(b - 1) << 1] /= ((double) l - (double) (b - 1) * resolution) * re_cuts->density;
        ds[(b - 1) << 1 | 1] /= ((double) l - (double) (b - 1) * resolution) * re_cuts->density;
        
        for (j = 0; j < b << 1; ++j)
            if (ds[j] < MIN_RE_DENS || ds[j] > MAX_RE_DENS)
                 ds[j] = .0;

        dens[i] = ds;
#ifdef DEBUG_ENZ
        fprintf(stderr, "[DEBUG_ENZ::%s] DENS2 [%u/%u] (%u):", __func__, i, n, b);
        for (j = 0; j < b << 1; ++j)
            fprintf(stderr, " %.6f", ds[j]);
        fprintf(stderr, "\n");
#endif

    }

    return dens;
}

