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
 * 02/09/21 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#include "kdq.h"

#include "sdict.h"
#include "break.h"
#include "asset.h"

#undef DEBUG
#undef REMOVE_NOISE
#undef DEBUG_LOCAL_BREAK

link_t *link_init(uint32_t s, uint32_t n)
{
    link_t *link = (link_t *) malloc(sizeof(link_t));
    link->s = s;
    link->n = n;
    int64_t *link_c = (int64_t *) malloc(n * sizeof(int64_t));
    int64_t i;
    for (i = 0; i < n; ++i)
        link_c[i] = i << 32;
    link->link = link_c;
    return link;
}

link_mat_t *link_mat_init(asm_dict_t *dict, uint32_t b)
{
    link_mat_t *link_mat = (link_mat_t *) malloc(sizeof(link_mat_t));
    link_mat->b = b;
    link_mat->n = dict->n;
    link_mat->link = (link_t *) malloc(link_mat->n * sizeof(link_t));

    uint32_t i;
    for (i = 0; i < link_mat->n; ++i)
        link_mat->link[i] = *link_init(i, div_ceil(dict->s[i].len, b));
    
    return link_mat;
}

void link_mat_destroy(link_mat_t *link_mat)
{
    uint32_t i;
    for (i = 0; i < link_mat->n; ++i)
        free(link_mat->link[i].link);
    free(link_mat->link);
    free(link_mat);
}

uint32_t estimate_dist_thres_from_file(const char *f, asm_dict_t *dict, double min_frac, uint32_t resolution)
{
    uint32_t i;
    FILE *fp;
    long pair_c, intra_c, cum_c;
    uint64_t max_len;
    uint32_t nb, *link_c;
    uint32_t buffer[BUFF_SIZE], m, i0, i1;
    uint64_t p0, p1;

    max_len = 0;
    for (i = 0; i < dict->n; ++i)
        if (dict->s[i].len > max_len)
            max_len = dict->s[i].len;
    nb = div_ceil(max_len, resolution);
    link_c = (uint32_t *) calloc(nb, sizeof(uint32_t));

    fp = fopen(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    pair_c = 0;
    intra_c = 0;
    while (1) {
        m = fread(&buffer, sizeof(uint32_t), BUFF_SIZE, fp);
        for (i = 0; i < m; i += 4) {
            sd_coordinate_conversion(dict, buffer[i], buffer[i + 1], &i0, &p0, 0);
            sd_coordinate_conversion(dict, buffer[i + 2], buffer[i + 3], &i1, &p1, 0);
            if (i0 == i1) {
                ++link_c[labs((long) p0 - p1) / resolution];
                ++intra_c;
            }
        }
        pair_c += m / 4;

        if (m < BUFF_SIZE) {
            if (ferror(fp))
                return 0;
            break;
        }
    }
    
    fclose(fp);
    
    i = 0;
    cum_c = 0;
    while (cum_c < intra_c * min_frac)
        cum_c += link_c[i++];
    free(link_c);

#ifdef REMOVE_NOISE
      printf("[I::%s] %ld read pairs processed, intra links: %ld \n", __func__, pair_c, intra_c);
#endif
    return i * resolution;
}

static void calc_moving_average(int64_t *arr, int32_t n, int32_t a)
{
    int i, pos, a1, a2;
    if (n < a) {
        a = n;
    }
    a1 = a >> 1;
    a2 = ~a & 1;

    int64_t *buff = (int64_t *) malloc(a * sizeof(int64_t));
    int64_t sum = 0;
    for (i = 0; i < a; ++i) {
        buff[i] = arr[i];
        sum += buff[i];
        if (i >= a1)
            arr[i - a1] = sum / (i + 1);
    }
    pos = 0;

    for (i = a1; i < n - a1 + a2; ++i) {
        arr[i] = sum / a;
        if (i < n - a + a1) {
            sum -= buff[pos];
            buff[pos] = arr[i + a - a1];
            sum += buff[pos];
            if (++pos == a)
                pos = 0;
        }
    }

    for (i = n - a1 + a2; i < n; ++i) {
        sum -= buff[pos];
        arr[i] = sum / (a1 + n - i);
        if (++pos == a)
            pos = 0;
    }

    free(buff);
}

link_mat_t *link_mat_from_file(const char *f, asm_dict_t *dict, uint32_t dist_thres, uint32_t resolution, double noise, uint32_t move_avg)
{
    FILE *fp;
    uint32_t i, j, n;
    uint32_t buffer[BUFF_SIZE], m, i0, i1;
    uint64_t p0, p1;
    long pair_c, intra_c;

    fp = fopen(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    link_mat_t *link_mat = (link_mat_t *) malloc(sizeof(link_mat_t));
    link_mat->b = resolution;
    link_mat->n = dict->n;
    link_mat->link = (link_t *) malloc(link_mat->n * sizeof(link_t));
    for (i = 0; i < link_mat->n; ++i) {
        n = div_ceil(dict->s[i].len, resolution);
        link_mat->link[i].s = i;
        link_mat->link[i].n = n;
        link_mat->link[i].link = (int64_t *) calloc(n, sizeof(int64_t));
    }

    pair_c = intra_c = 0;
    while (1) {
        m = fread(&buffer, sizeof(uint32_t), BUFF_SIZE, fp);
        for (i = 0; i < m; i += 4) {
            sd_coordinate_conversion(dict, buffer[i], buffer[i + 1], &i0, &p0, 0);
            sd_coordinate_conversion(dict, buffer[i + 2], buffer[i + 3], &i1, &p1, 0);

            if (p0 > p1)
                SWAP(uint64_t, p0, p1);
            if (i0 == i1 && p1 - p0 <= dist_thres) {
                ++intra_c;
                link_mat->link[i0].link[(MAX(p0, 1) - 1) / resolution] += 1;
                link_mat->link[i1].link[(MAX(p1, 1) - 1) / resolution] -= 1;
            }
        }
        pair_c += m / 4;

        if (m < BUFF_SIZE) {
            if (ferror(fp))
                return 0;
            break;
        }
    }
    fclose(fp);

#ifdef DEBUG
    printf("[I::%s] %ld read pairs processed, intra links: %ld \n", __func__, pair_c, intra_c);
#endif
    
    int64_t *link;
    for (i = 0; i < link_mat->n; ++i) {
        link = link_mat->link[i].link;
        n = link_mat->link[i].n;
        for (j = 1; j < n; ++j)
            link[j] += link[j - 1];
    }

#ifdef REMOVE_NOISE
    double l;
    noise = noise * dist_thres * (dist_thres + resolution) / 2;
    for (i = 0; i < link_mat->n; ++i) {
        link = link_mat->link[i].link;
        n = link_mat->link[i].n;
        for (j = 0; j < n; ++j) {
            l = (double) link[j] - noise;
            link[j] = (uint32_t) MAX(l, .1);
        }
    }
#endif

    int32_t ma_k = move_avg / resolution;
    if (ma_k > 1)
        for (i = 0; i < link_mat->n; ++i)
            calc_moving_average(link_mat->link[i].link, link_mat->link[i].n, ma_k);

    for (i = 0; i < link_mat->n; ++i) {
        int64_t *link_c = link_mat->link[i].link;
        n = link_mat->link[i].n;
        for (j = 0; j < n; ++j)
            link_c[j] |= (int64_t) j << 32;
    }
    
    return link_mat;
}

// sort by position = higher 32 bits
int pos_cmp(const void *p, const void *q)
{
    int64_t d = *(int64_t *) p - *(int64_t *) q;
    return d == 0? 0 : (d > 0? 1 : -1);
}

// sort by link count = lower 32 bits
int cnt_cmp(const void *p, const void *q)
{
    int32_t d = *(int32_t *) p - *(int32_t *) q;
    return d == 0? 0 : (d > 0? 1 : -1);
}

KDQ_INIT(int64_t)

static void add_break_point(bp_t *bp, uint64_t p)
{
    if (bp->n == bp->m) {
        bp->m <<= 1;
        bp->p = (uint64_t *) realloc(bp->p, bp->m * sizeof(uint64_t));
    }
    bp->p[bp->n] = p;
    ++bp->n;
}

bp_t *detect_break_points_local_joint(link_mat_t *link_mat, uint32_t bin_size, double fold_thres, uint32_t flank_size, asm_dict_t *dict, uint32_t *bp_n)
{
    uint32_t i, j, b_n, b_m;
    double mcnt;
    int64_t *link;
    uint32_t s, e, t;
    int8_t a;
    bp_t *bp, *bp1;
    sd_seg_t *segs, seg;
    sd_aseq_t seq;
    
    segs = dict->seg;
    b_n = 0;
    b_m = 16;
    bp = (bp_t *) malloc(b_m * sizeof(bp_t));
    bp1 = 0;
    for (i = 0; i < link_mat->n; ++i) {
        seq = dict->s[i];
        link = link_mat->link[i].link;
        a = 0;
        for (j = 1; j < seq.n; j++) {
            seg = segs[seq.s + j];
            s = (MAX(seg.a - MIN(flank_size, segs[seq.s + j - 1].y), 1) - 1) / bin_size;
            e = (MAX(seg.a + MIN(flank_size, seg.y), 1) - 1) / bin_size;
            // s = (MAX(seg.a - MIN(flank_size, seg.a), 1) - 1) / bin_size;
            // e = (MAX(seg.a + MIN(flank_size, seq.len - seg.a), 1) - 1) / bin_size;
            t = e - s + 1;
            qsort(link + s, t, sizeof(int64_t), cnt_cmp);
            mcnt = t & 1? (int32_t) link[(e + s) / 2] : ((int32_t) link[(e + s) / 2] + (int32_t) link[(e + s) / 2 + 1]) / 2.;
            mcnt *= fold_thres;
            qsort(link + s, t, sizeof(int64_t), pos_cmp);

            if ((int32_t) link[(MAX(seg.a, 1) - 1) / bin_size] < mcnt) {
                if (!a) {
                    if (b_n == b_m) {
                        b_m <<= 1;
                        bp = (bp_t *) realloc(bp, b_m * sizeof(bp_t));
                    }
                    bp1 = bp + b_n;
                    bp1->s = i;
                    bp1->n = 0;
                    bp1->m = 4;
                    bp1->p = (uint64_t *) malloc(bp1->m * sizeof(uint64_t));
                    ++b_n;
                    a = 1;
                }
                add_break_point(bp1, seg.a);
#ifdef DEBUG_LOCAL_BREAK
                printf("[I::%s] break local joint: %s at %lu (link number %d < %.3f)\n", __func__, seq.name, seg.a, (int32_t) link[(MAX(seg.a, 1) - 1) / bin_size], mcnt);
#endif
            }
        }
    }

    *bp_n = b_n;

    return bp;
}

static int make_dual_break(int64_t *link, uint32_t s, uint32_t e, uint32_t d, double fold_thres, uint32_t *bp_s, uint32_t *bp_e)
{
    // make dual break if there is a bump
    // find the exact break points
    uint32_t i, l, lm, u, ls, le;

    ls = u = 0;
    lm = UINT32_MAX;
    for (i = s; i <= e; ++i) {
        l = (uint32_t) link[i];
        if (l < lm) {
            u = 0;
            ls = i;
            lm = l;
        } else if (l > MAX(10., lm / fold_thres)) {
            if (++u >= d)
                break;
        }
    }

    le = u = 0;
    lm = UINT32_MAX;
    for (i = s; i <= e; ++i) {
        l = (uint32_t) link[e - i + s];
        if (l < lm) {
            u = 0;
            le = e - i + s;
            lm = l;
        } else if (l > MAX(10., lm / fold_thres)) {
            if (++u >= d)
                break;
        }
    }
    
    if (le >= ls + d) {
        *bp_s = ls;
        *bp_e = le;
        return 1;
    }

    return 0;
}

bp_t *detect_break_points(link_mat_t *link_mat, uint32_t bin_size, uint32_t merge_size, double fold_thres, uint32_t dual_break_thres, uint32_t *bp_n)
{
    uint32_t i, j, k, n, m, d, b, b_n, b_m;
    double mcnt;
    int64_t *link;
    uint32_t s, e, t, min_c, p, bp_s, bp_e;
    kdq_t(int64_t) *q;
    bp_t *bp, *bp1;
    
    b_n = 0;
    b_m = 16;
    bp = (bp_t *) malloc(b_m * sizeof(bp_t));
    bp1 = 0;
    m = merge_size / bin_size;
    d = dual_break_thres / bin_size;
    q = kdq_init(int64_t);
    for (i = 0; i < link_mat->n; ++i) {
        link = link_mat->link[i].link;
        n = link_mat->link[i].n;
        // sort by link count
        qsort(link, n, sizeof(int64_t), cnt_cmp);
        // find median
        mcnt = n & 1? (int32_t) link[n / 2] : ((int32_t) link[n / 2] + (int32_t) link[n / 2 - 1]) / 2.;
        // find count threshold
        mcnt *= fold_thres;
        // find positions below threshold
        b = 0;
        while (b < n && (int32_t) link[b] < mcnt)
            ++b;
        if (!b || b == n)
            continue;
        // sort by position for positions below threshold
        qsort(link, b, sizeof(int64_t), pos_cmp);
        // detect blocks for break points
        kdq_clean(q);
        s = link[0] >> 32;
        e = s;
        for (j = 1; j < b; ++j) {
            t = link[j] >> 32;
            if (e + m < t) {
                // new block
                kdq_push(int64_t, q, (int64_t) s << 32 | e);
                s = t;
                e = s;
            } else {
                e = t;
            }
        }
        // add last block
        kdq_push(int64_t, q, (int64_t) s << 32 | e);
        // sort by position for all positions
        qsort(link, n, sizeof(int64_t), pos_cmp);
        // detect precise break points
        if (b_n == b_m) {
            b_m <<= 1;
            bp = (bp_t *) realloc(bp, b_m * sizeof(bp_t));
        }
        bp1 = bp + b_n;
        bp1->s = i;
        bp1->n = 0;
        bp1->m = 4;
        bp1->p = (uint64_t *) malloc(bp1->m * sizeof(uint64_t));
        ++b_n;

        for (j = 0; j < kdq_size(q); ++j) {
            s = kdq_at(q, j) >> 32;
            e = (int32_t) kdq_at(q, j);
            
            if (e - s > d && make_dual_break(link, s, e, d, fold_thres, &bp_s, &bp_e)) {
                // dual break
                if (s != 0)
                    add_break_point(bp1, bp_s);
                if (e != n - 1)
                    add_break_point(bp1, bp_e);
                continue;
            }

            // skip the first and last block
            if (s == 0 || e == n - 1)
                continue;
            
            // find the position of the minimum value in this valley
            // which will be a break point
            min_c = INT32_MAX;
            p = UINT32_MAX;
            for (k = s; k <= e; ++k) {
                if ((int32_t) link[k] < min_c) {
                    min_c = (int32_t) link[k];
                    p = k;
                }
            }

            if (p != UINT32_MAX)
                add_break_point(bp1, p);
        }

        if (bp1->n > 1) {
            // revisit to remove false positives
            kdq_clean(q);
            for (k = 0; k < bp1->n; ++k) {
                s = k > 0? bp1->p[k - 1] : 0;
                e = k < bp1->n - 1? bp1->p[k + 1] : n - 1;
                t = e - s + 1;
                qsort(link + s, t, sizeof(int64_t), cnt_cmp);
                mcnt = t & 1? (int32_t) link[(e + s) / 2] : ((int32_t) link[(e + s) / 2] + (int32_t) link[(e + s) / 2 + 1]) / 2.;
                mcnt *= fold_thres;
                qsort(link + s, t, sizeof(int64_t), pos_cmp);
                if ((int32_t) link[bp1->p[k]] > mcnt)
                    kdq_push(int64_t, q, k);
            }
            for (k = 0; k < kdq_size(q); ++k)
                bp1->p[kdq_at(q, k)] = UINT32_MAX;
        }
        
        j = 0;
        for (k = 0; k < bp1->n; ++k)
            if (bp1->p[k] != UINT32_MAX)
                bp1->p[j++] = bp1->p[k] * bin_size;
        bp1->n = j;
        
        if (!bp1->n) {
            free(bp1->p);
            --b_n;
        }
    }
    kdq_destroy(int64_t, q);
    *bp_n = b_n;
    
    return bp;
}

void print_link_mat(link_mat_t *link_mat, asm_dict_t *dict, FILE *fp)
{
    int i, j, n;
    int64_t *link;
    for (i = 0; i < link_mat->n; ++i) {
        link = link_mat->link[i].link;
        n = link_mat->link[i].n;
        fprintf(fp, "%s:", dict->s[link_mat->link[i].s].name);
        for (j = 0; j < n; ++j)
            fprintf(fp, " %d", (int32_t) link[j]);
        fprintf(fp, "\n");
    }
}

void print_break_point(bp_t *bp, asm_dict_t *dict, FILE *fp)
{
    int i;
    fprintf(fp, "%s[%u]:", dict->s[bp->s].name, bp->n);
    for (i = 0; i < bp->n; ++i)
        fprintf(fp, " %lu", bp->p[i]);
    fprintf(fp, "\n");
}

void write_break_agp(asm_dict_t *d, bp_t *breaks, uint32_t b_n, FILE *fp)
{
    uint32_t i, j, s, ns, ms;
    int64_t L, l;
    uint64_t len;
    sd_seg_t *segs, seg;
    sdict_t *sd;
    uint64_t *p;
    uint32_t p_n;

    ns = 0;
    ms = 4096;
    segs = (sd_seg_t *) malloc(ms * sizeof(sd_seg_t));
    sd = d->sdict;
    s = 0;
    for (i = 0; i < d->n; ++i) {
        if (b_n == 0 || i < breaks->s) {
            // no breaks
            write_segs_to_agp(d->seg + d->s[i].s, d->s[i].n, sd, ++s, fp);
        } else {
            // contain break points
            p = breaks->p;
            p_n = breaks->n;
            len = 0;
            for (j = 0; j < d->s[i].n; ++j) {
                seg = d->seg[d->s[i].s + j];
                if (p_n == 0 || seg.a + seg.y <= p[0]) {
                    // current seg contains no break points
                    sd_seg_t seg1 = {seg.s, ns, len, seg.c, seg.x, seg.y};
                    segs[ns] = seg1;
                    if (ns == ms) {
                        ms <<= 1;
                        segs = (sd_seg_t *) realloc(segs, ms * sizeof(sd_seg_t));
                    }
                    ++ns;
                    len += seg.y;
                } else {
                    l = (int64_t) p[0] - seg.a;
                    assert(-l < UINT32_MAX && l < UINT32_MAX);
                    if (l > 0) {
                        // split seg
                        sd_seg_t seg1 = {seg.s, ns, len, seg.c, seg.c & 1? (uint32_t) (-l + seg.x + seg.y) : seg.x, (uint32_t) l};
                        segs[ns] = seg1;
                        ++ns;
                    }

                    write_segs_to_agp(segs, ns, sd, ++s, fp);
                    ns = 0;
                    len = 0;
                    if (--p_n)
                        ++p;

                    L = l;
                    while (p_n > 0 && seg.a + seg.y > p[0]) {
                        l = (int64_t) (p[0] - p[-1]);
                        assert(l < UINT32_MAX);
                        sd_seg_t seg1 = {seg.s, ns, 0, seg.c, seg.c & 1? (uint32_t) (-L - l + seg.x + seg.y) : (uint32_t) (L + seg.x), (uint32_t) l};
                        write_segs_to_agp(&seg1, 1, sd, ++s, fp);
                        L += l;
                        if (--p_n) 
                            ++p;
                    }
                    
                    assert(-L < UINT32_MAX && L < UINT32_MAX);
                    if (L < seg.y) {
                        sd_seg_t seg1 = {seg.s, ns, 0, seg.c, seg.c & 1? seg.x : (uint32_t) (L + seg.x), (uint32_t) (-L + seg.y)};
                        len += -L + seg.y;
                        segs[ns] = seg1;
                        ++ns;
                    }
                }
            }
            if (ns > 0) {
                write_segs_to_agp(segs, ns, sd, ++s, fp);
                ns = 0;
            }
            if(--b_n) 
                ++breaks;
        }
    }

    free(segs);
}

