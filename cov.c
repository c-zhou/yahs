/*********************************************************************************
 * MIT License                                                                   *
 *                                                                               *
 * Copyright (c) 2022 Chenxi Zhou <chnx.zhou@gmail.com>                          *
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
 * 02/09/22 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>

#include "kvec.h"
#include "bamlite.h"
#include "cov.h"
#include "asset.h"

KHASH_SET_INIT_STR(str)

// read unmapped (0x4)
// not primary alignment (0x100)
// read fails platform/vendor quality checks (0x200)
// read is PCR or optical duplicate (0x400)
// supplementary alignment (0x800)
static uint16_t mask = 0xF04;
static double max_cov_scale = 1.0;

cov_t *cov_init(uint32_t n)
{
    cov_t *cov = (cov_t *) malloc(sizeof(cov_t));
    cov->n = n;
    cov->p = (pos_t *) malloc(sizeof(pos_t) * n);
    uint32_t i;
    for (i = 0; i < n; ++i)
        kv_init(cov->p[i]);

    return cov;
}

void cov_destroy(cov_t* cov)
{
    uint32_t i;
    for (i = 0; i < cov->n; ++i)
        kv_destroy(cov->p[i]);
    free(cov->p);
    free(cov);

    return;
}

void cov_norm_destroy(cov_norm_t *cov_norm)
{
    if (cov_norm) {
        free(cov_norm->norm[0]);
        free(cov_norm->norm);
        free(cov_norm);
    }

    return;
}

static int diff_target_length(uint32_t *target_len, sdict_t *sdict)
{
    uint32_t i;
    for (i = 0; i < sdict->n; ++i)
        if (target_len[i] != sdict->s[i].len)
            return 1;
    return 0;
}

static int diff_target_name(char **target_name, sdict_t *sdict)
{
    uint32_t i;
    for (i = 0; i < sdict->n; ++i)
        if (strcmp(target_name[i], sdict->s[i].name))
            return 1;
    return 0;
}

static int u64_cmpfunc (const void *a, const void *b)
{
    return (*(uint64_t *) a > *(uint64_t *) b) - (*(uint64_t *) a < *(uint64_t *) b);
}

static uint64_t pos_compression(cov_t *cov)
{
    uint32_t i, p, p1;
    int32_t c;
    size_t j, k, n;
    uint64_t m, *a;
    
    m = 0;
    for (i = 0; i < cov->n; ++i) {
        a = cov->p[i].a;
        n = cov->p[i].n;
        if (n == 0) continue;

        qsort(a, n, sizeof(uint64_t), u64_cmpfunc);
        p = (uint32_t) (a[0] >> 32);
        c = (int32_t) a[0];
        k = 0;
        for(j = 1; j < n; ++j) {
            p1 = (uint32_t) (a[j] >> 32);
            if (p != p1) {
                if (c != 0) a[k++] = (uint64_t) p << 32 | (uint32_t) c;
                p = p1;
                c = (int32_t) a[j];
            } else {
                c += (int32_t) a[j];
            }
        }
        if (c != 0) a[k++] = (uint64_t) p << 32 | (uint32_t) c;
        cov->p[i].n = k;
        m += k;
    }

    return m;
}

cov_t *bam_cstats(const char *bam, sdict_t *sdict, int match_header)
{
    uint32_t i, flag;
    uint64_t n, m, max_m, n_recs;
    bamFile fp;
    bam1_t *b1;
    bam_header_t *h;
    cov_t *covs;
    
    khash_t(str) *hmseq; // for absent sequences
    khint_t k;
    int absent;
    hmseq = kh_init(str);

    fp = bam_open(bam, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, bam);
        exit(EXIT_FAILURE);
    }

    h = bam_header_read(fp);
    if (match_header) {
        int n_targets;
        char **target_name;
        uint32_t *target_len;

        n_targets = h->n_targets;
        target_len = h->target_len;
        target_name = h->target_name;

        // check consistency between sdict and bam header
        if (sdict && (n_targets != sdict->n || 
                    diff_target_length(target_len, sdict) || 
                    diff_target_name(target_name, sdict))) {
            fprintf(stderr, "[E::%s] sequence dictionary does not match BAM header\n", __func__);
            exit(EXIT_FAILURE);
        }
    }

    covs = cov_init(sdict->n);
    b1 = bam_init1();
    max_m = 0x7FFFFFFULL; // 128MB - ~1GB mem for covs
    n_recs = 0;
    n = 0;
    while (bam_read1(fp, b1) >= 0) {
        ++n_recs;
        if (n_recs % 1000000 == 0)
            fprintf(stderr, "[I::%s] %lu million records processed\n", __func__, n_recs / 1000000);
        flag = b1->core.flag;
        if (flag & mask) continue;
        i = sd_get(sdict, h->target_name[b1->core.tid]);
        if (i == UINT32_MAX) {
            k = kh_put(str, hmseq, h->target_name[b1->core.tid], &absent);
            if (absent) {
                kh_key(hmseq, k) = strdup(h->target_name[b1->core.tid]);
                fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, h->target_name[b1->core.tid]);
            }
            continue;
        }
        kv_push(uint64_t, covs->p[i], (uint64_t) b1->core.pos << 32 | (uint32_t) 1);
        kv_push(uint64_t, covs->p[i], (uint64_t) get_target_end(b1) << 32 | (uint32_t) -1);
        n += 2;
        if (n > max_m) {
            m = pos_compression(covs); // m is the position size after compression
            fprintf(stderr, "[I::%s] position compression n = %lu, m = %lu, max_m = %lu\n", __func__, n, m, max_m);
            if (m > n>>1) {
                max_m <<= 1; // increase m limit if compression ratio smaller than 0.5
                fprintf(stderr, "[I::%s] position memory buffer expanded max_m = %lu\n", __func__, max_m);
            }
            n = m;
        }
    }
    m = pos_compression(covs);
    fprintf(stderr, "[I::%s] position compression n = %lu, m = %lu, max_m = %lu\n", __func__, n, m, max_m);
    
    fprintf(stderr, "[I::%s] total number BAM records processed %lu\n", __func__, n_recs);
    
    for (k = 0; k < kh_end(hmseq); ++k)
        if (kh_exist(hmseq, k))
            free((char *) kh_key(hmseq, k));
    kh_destroy(str, hmseq);

    bam_destroy1(b1);
    bam_header_destroy(h);
    bam_close(fp);

    return covs;
}

cov_t *bed_cstats(const char *bed, sdict_t *sdict)
{
    uint32_t i, s, e;
    uint64_t n, m, max_m, n_recs;
    char *line, cname[4096];
    size_t ln = 0;
    ssize_t read;
    FILE *fp;
    cov_t *covs;
    
    khash_t(str) *hmseq; // for absent sequences
    khint_t k;
    int absent;
    hmseq = kh_init(str);

    fp = fopen(bed, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, bed);
        exit(EXIT_FAILURE);
    }
    
    covs = cov_init(sdict->n);
    max_m = 0x7FFFFFFULL; // 128MB - ~1GB mem for covs
    n_recs = 0;
    n = 0;
    line = 0;
    while ((read = getline(&line, &ln, fp)) != -1) {
        ++n_recs;
        if (n_recs % 1000000 == 0)
            fprintf(stderr, "[I::%s] %lu million records processed\n", __func__, n_recs / 1000000);
        sscanf(line, "%s %u %u %*s %*s %*s", cname, &s, &e);
        i = sd_get(sdict, cname);
        if (i == UINT32_MAX) {
            k = kh_put(str, hmseq, cname, &absent);
            if (absent) {
                kh_key(hmseq, k) = strdup(cname);
                fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname);
            }
            continue;
        }
        kv_push(uint64_t, covs->p[i], (uint64_t) s << 32 | (uint32_t) 1);
        kv_push(uint64_t, covs->p[i], (uint64_t) e << 32 | (uint32_t) -1);
        n += 2;
        if (n > max_m) {
            m = pos_compression(covs); // m is the position size after compression
            fprintf(stderr, "[I::%s] position compression n = %lu, m = %lu, max_m = %lu\n", __func__, n, m, max_m);
            if (m > n>>1) {
                max_m <<= 1; // increase m limit if compression ratio smaller than 0.5
                fprintf(stderr, "[I::%s] position memory buffer expanded max_m = %lu\n", __func__, max_m);
            }
            n = m;
        }
    }
    m = pos_compression(covs);
    fprintf(stderr, "[I::%s] position compression n = %lu, m = %lu, max_m = %lu\n", __func__, n, m, max_m);
    
    fprintf(stderr, "[I::%s] total number BED records processed %lu\n", __func__, n_recs);

    for (k = 0; k < kh_end(hmseq); ++k)
    if (kh_exist(hmseq, k))
        free((char *) kh_key(hmseq, k));
    kh_destroy(str, hmseq);

    if (line)
        free(line);
    fclose(fp);

    return covs;
}

KHASH_MAP_INIT_INT64(ctab, int64_t)

static void kh_ctab_put1(kh_ctab_t *ctab, int32_t c, int cnt)
{
    int absent;
    khint_t k;
    k = kh_put_ctab(ctab, c, &absent);
    if (absent)
        kh_val(ctab, k) = cnt;
    else
        kh_val(ctab, k) += cnt;
}

typedef struct {int64_t x, y;} ctab_pair;

static int ctab_cmpfunc(const void *a, const void *b)
{
    int64_t x, y;
    x = ((ctab_pair *) a)->x;
    y = ((ctab_pair *) b)->x;
    return (x > y) - (x < y);
}

double calc_avg_cov(cov_t *cov, sdict_t *sdict, double q_drop)
{
    size_t i, j, s;
    int32_t c;
    uint32_t p, p1;
    uint64_t *a;
    kh_ctab_t *ctab;
    khint_t k;

    // get base coverage table
    ctab = kh_init_ctab();
    for (i = 0; i < cov->n; ++i) {
        a = cov->p[i].a;
        s = cov->p[i].n;
        
        p = 0, c = 0;
        for (j = 0; j < s; ++j) {
            p1 = (uint32_t) (a[j] >> 32);
            if (p1 > p)
                kh_ctab_put1(ctab, c, p1 - p);
            p = p1;
            c += (int32_t) a[j];
        }

        if (p < sdict->s[i].len)
            kh_ctab_put1(ctab, c, sdict->s[i].len - p);
    }

    kvec_t(ctab_pair) cnts;
    int64_t bases, covs, qbases, qcovs, b, qb;

    bases = covs = 0;
    kv_init(cnts);
    for (k = kh_begin(ctab); k != kh_end(ctab); ++k) {
        if (kh_exist(ctab, k)) {
            ctab_pair cnt = {kh_key(ctab, k), kh_val(ctab, k)};
            kv_push(ctab_pair, cnts, cnt);
            bases += cnt.y;
            covs += cnt.x * cnt.y;
        }
    }
    qsort(cnts.a, cnts.n, sizeof(ctab_pair), ctab_cmpfunc);
    
    qb = bases * q_drop;
    qbases = bases - qb * 2;
    qcovs = covs;
    for (i = 0, b = 0; i < cnts.n; ++i) {
        b += cnts.a[i].y;
        if (b < qb) {
            qcovs -= cnts.a[i].x * cnts.a[i].y;
        } else {
            qcovs -= cnts.a[i].x * (qb - b + cnts.a[i].y);
            break;
        }
    }
    for (i = cnts.n, b = 0; i != 0 ; --i) {
        b += cnts.a[i-1].y;
        if (b < qb) {
            qcovs -= cnts.a[i-1].x * cnts.a[i-1].y;
        } else {
            qcovs -= cnts.a[i-1].x * (qb - b + cnts.a[i-1].y);
            break;
        }
    }
    
    kh_destroy(ctab, ctab);
    kv_destroy(cnts);

    double avg_cov;
    avg_cov = qbases == 0? .0 : (double) qcovs / qbases;
    
    fprintf(stderr, "[I::%s] sequence coverage stats:\n", __func__);
    fprintf(stderr, "[I::%s] sequence bases: %ld\n", __func__, bases);
    fprintf(stderr, "[I::%s] read bases: %ld\n", __func__, covs);
    fprintf(stderr, "[I::%s] q drop: %.3f\n", __func__, q_drop);
    fprintf(stderr, "[I::%s] average read coverage: %.3f\n", __func__, avg_cov);

    return avg_cov;
}

cov_norm_t *calc_cov_norms(cov_t *cov, sdict_t *sdict, uint32_t window, double q_drop)
{
    size_t i, j, k, s;
    int32_t c;
    uint32_t n, p, p1, w;
    int64_t m, c1;
    uint64_t *a;
    double *norm_a, **norm, avg_cov;
    cov_norm_t *cov_norm;
    
    n = sdict->n;
    m = 0;
    for (i = 0; i < n; ++i)
        m += div_ceil(sdict->s[i].len, window);
    norm_a = (double *) calloc(m, sizeof(double));
    norm = (double **) calloc(n, sizeof(double *));
    norm[0] = norm_a;
    for (i = 1; i < n; ++i) 
        norm[i] = norm[i-1] + div_ceil(sdict->s[i-1].len, window);

    // calculate average base coverage
    avg_cov = calc_avg_cov(cov, sdict, q_drop);
    for (i = 0; i < cov->n; ++i) {
        a = cov->p[i].a;
        s = cov->p[i].n;
        
        w = p = 0, c = 0, c1 = 0, k = 0;
        for (j = 0; j < s; ++j) {
            p1 = (uint32_t) (a[j] >> 32);
            if (p1 > p) {
                if (p1 - p < window - w) {
                    c1 += (int64_t) c * (p1 - p);
                    w += p1 - p;
                    p = p1;
                } else {
                    while (p1 > p) {
                        if (p1 - p >= window - w) {
                            c1 += (int64_t) c * (window - w);
                            norm[i][k++] = MIN(max_cov_scale, c1 > 0? avg_cov * window / c1 : DBL_MAX);
                            p += window - w;
                            c1 = 0;
                            w = 0;
                        } else {
                            c1 += (int64_t) c * (p1 - p);
                            w += p1 - p;
                            p = p1;
                        }
                    }
                }
            }
            c += (int32_t) a[j];
        }

        p1 = sdict->s[i].len;
        while (p1 > p) {
            if (p1 - p >= window - w) {
                norm[i][k++] = MIN(max_cov_scale, c1 > 0? avg_cov * window / c1 : DBL_MAX);
                p += window - w;
            } else {
                norm[i][k++] = MIN(max_cov_scale, c1 > 0? avg_cov * (p1 - p + w) / c1 : DBL_MAX);
                p = p1;
            }
            c1 = 0;
            w = 0;
        }

        if (w > 0)
            norm[i][k++] = MIN(max_cov_scale, c1 > 0? avg_cov * w / c1 : DBL_MAX);
        
        assert(k == div_ceil(p1, window));
    }
    
    cov_norm = (cov_norm_t *) malloc(sizeof(cov_norm_t));
    cov_norm->n = m;
    cov_norm->w = window;
    cov_norm->norm = norm;

    return cov_norm;
}

void print_cov_in_bed(cov_t *cov, sdict_t *sdict, FILE *fo)
{
    size_t i, j, s;
    int32_t c;
    uint32_t p, p1;
    uint64_t *a;

    for (i = 0; i < cov->n; ++i) {
        a = cov->p[i].a;
        s = cov->p[i].n;

        p = 0, c = 0;
        for (j = 0; j < s; ++j) {
            p1 = (uint32_t) (a[j] >> 32);
            if (p1 > p)
                fprintf(fo, "%s\t%u\t%u\t%d\n", sdict->s[i].name, p, p1, c);
            p = p1;
            c += (int32_t) a[j];
        }

        if (p < sdict->s[i].len)
            fprintf(fo, "%s\t%u\t%u\t%d\n", sdict->s[i].name, p, sdict->s[i].len, 0);
    }
}
