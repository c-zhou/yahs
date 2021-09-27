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

#include "kdq.h"

#include "sdict.h"
#include "break.h"
#include "asset.h"

#undef DEBUG
#undef REMOVE_NOISE

link_t *link_init(int32_t s, int32_t n)
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

link_mat_t *link_mat_init(asm_dict_t *dict, int32_t b)
{
    link_mat_t *link_mat = (link_mat_t *) malloc(sizeof(link_mat_t));
    link_mat->b = b;
    link_mat->n = dict->n;
    link_mat->link = (link_t *) malloc(link_mat->n * sizeof(link_t));

    int i;
    for (i = 0; i < dict->n; ++i)
        link_mat->link[i] = *link_init(i, dict->s[i].len / b + 1);
    return link_mat;
}

void link_mat_destroy(link_mat_t *link_mat)
{
    int i;
    for (i = 0; i < link_mat->n; ++i)
        free(link_mat->link[i].link);
    free(link_mat->link);
    free(link_mat);
}

int estimate_dist_thres_from_file(const char *f, sdict_t *dict, double min_frac, int32_t resolution)
{
    int i;
    FILE *fp;
    long pair_c, intra_c, cum_c;
    int32_t max_len, nb, *link_c;
    int32_t buffer[BUFF_SIZE], m;

    max_len = 0;
    for (i = 0; i < dict->n; ++i)
        if (dict->s[i].len > max_len)
            max_len = dict->s[i].len;
    nb = max_len / resolution + 1;
    link_c = (int32_t *) calloc(nb, sizeof(int32_t));

    fp = fopen(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    pair_c = 0;
    intra_c = 0;
    while (1) {
        m = fread(&buffer, sizeof(int32_t), BUFF_SIZE, fp);
        for (i = 0; i < m; i += 4) {
            if (buffer[i] == buffer[i + 2]) {
                ++link_c[(int) (abs(buffer[i + 1] - buffer[i + 3]) / resolution)];
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

    return i * resolution;
}

void calc_moving_average(int64_t *arr, int32_t n, int32_t a)
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

link_mat_t *link_mat_from_file(const char *f, asm_dict_t *dict, int32_t dist_thres, int32_t resolution, int32_t move_avg)
{
    FILE *fp;
    int32_t i, j, n, i0, i1, p0, p1;
    int32_t buffer[BUFF_SIZE], m;
    long pair_c, intra_c;
#ifdef REMOVE_NOISE
    long noise_c;
    double a, noise;
#endif

    fp = fopen(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    link_mat_t *link_mat = (link_mat_t *) malloc(sizeof(link_mat_t));
    link_mat->b = resolution;
    link_mat->n = dict->n;
    link_mat->link = (link_t *) malloc(link_mat->n * sizeof(link_t));
    for (i = 0; i < dict->n; ++i) {
        n = dict->s[i].len / resolution + 1;
        link_mat->link[i].s = i;
        link_mat->link[i].n = n;
        link_mat->link[i].link = (int64_t *) calloc(n, sizeof(int64_t));
    }

    pair_c = intra_c = 0;
#ifdef REMOVE_NOISE
    noise_c = 0;
#endif
    while (1) {
        m = fread(&buffer, sizeof(int32_t), BUFF_SIZE, fp);
        for (i = 0; i < m; i += 4) {
            sd_coordinate_conversion(dict, buffer[i], buffer[i + 1], &i0, &p0);
            sd_coordinate_conversion(dict, buffer[i + 2], buffer[i + 3], &i1, &p1);

            if (i0 == i1 && abs(p0 - p1) <= dist_thres) {
                ++intra_c;
                if (p0 > p1) 
                    SWAP(int32_t, p0, p1);
                link_mat->link[i0].link[p0 / resolution] += 1;
                link_mat->link[i0].link[p1 / resolution] -= 1;
#ifdef REMOVE_NOISE
            } else if (i0 != i1) {
                ++noise_c;
#endif
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

#ifdef REMOVE_NOISE
    a = 0.0;
    for (i = 0; i < link_mat->n; ++i)
        for (j = i + 1; j < link_mat->n; ++j)
            a += (double) dict->s[i].len / dist_thres * dict->s[j].len / dist_thres * 4;
    noise = a > 0? noise_c / a : 0;
#ifdef DEBUG
    printf("[I::%s] noise links: %ld; area: %.3f; noise estimation: %.3f\n", __func__, noise_c, a, noise);
#endif
#endif
    
    int64_t *link;
    for (i = 0; i < link_mat->n; ++i) {
        link = link_mat->link[i].link;
        n = link_mat->link[i].n;
        for (j = 1; j < n; ++j)
            link[j] += link[j - 1];
    }

#ifdef REMOVE_NOISE
    for (i = 0; i < link_mat->n; ++i) {
        link = link_mat->link[i].link;
        n = link_mat->link[i].n;
        for (j = 1; j < n; ++j) {
            link[j] -= noise;
            if (link[j] < 0)
                link[j] = 0;
        }
    }
#endif

    int32_t ma_k = move_avg / resolution;
    if (ma_k > 1)
        for (i = 0; i < link_mat->n; ++i)
            calc_moving_average(link_mat->link[i].link, link_mat->link[i].n, ma_k);

    for (i = 0; i < dict->n; ++i) {
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

static void add_break_point(bp_t *bp, int32_t p)
{
    if (bp->n == bp->m) {
        bp->m <<= 1;
        bp->p = (int32_t *) realloc(bp->p, bp->m * sizeof(int32_t));
    }
    bp->p[bp->n] = p;
    ++bp->n;
}

bp_t *detect_break_points_local_joint(link_mat_t *link_mat, int32_t bin_size, double fold_thres, int32_t flank_size, asm_dict_t *dict, int *bp_n)
{
    int i, j, b_n, b_m;
    double mcnt;
    int64_t *link;
    int32_t s, e, t;
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
            s = (seg.a - MIN(flank_size, segs[seq.s + j - 1].y)) / bin_size;
            e = (seg.a + MIN(flank_size, seg.y)) / bin_size;
            t = e - s + 1;
            qsort(link + s, t, sizeof(int64_t), cnt_cmp);
            mcnt = t & 1? (int32_t) link[(e + s) / 2] : ((int32_t) link[(e + s) / 2] + (int32_t) link[(e + s) / 2 + 1]) / 2.0;
            mcnt *= fold_thres;
            qsort(link + s, t, sizeof(int64_t), pos_cmp);

            if ((int32_t) link[seg.a / bin_size] < mcnt) {
                if (!a) {
                    if (b_n == b_m) {
                        b_m <<= 1;
                        bp = (bp_t *) realloc(bp, b_m * sizeof(bp_t));
                    }
                    bp1 = bp + b_n;
                    bp1->s = i;
                    bp1->n = 0;
                    bp1->m = 4;
                    bp1->p = (int32_t *) malloc(bp1->m * sizeof(int32_t));
                    ++b_n;
                    a = 1;
                }
                add_break_point(bp1, seg.a);
            }
        }
    }

    *bp_n = b_n;

    return bp;
}

static int make_dual_break(int64_t *link, int32_t s, int32_t e, int32_t n, double fold_thres)
{
    // make dual break if there is a bump
    // split the chunk into three parts
    // and compare the links
    int32_t i, l[3];
    float a;
    memset(l, 0, 3 * sizeof(int));
    a = (e - s + 1) / 3.0;
    for (i = s; i <= e; ++i)
        l[(int) ((i - s) / a)] += link[i];
    if (s == 0 && l[1] * fold_thres > l[2] ||
            e == n && l[1] * fold_thres > l[0] ||
            l[1] * fold_thres > l[0] && l[1] * fold_thres > l[2])
        return 1;
    return 0;
}

bp_t *detect_break_points(link_mat_t *link_mat, int32_t bin_size, int32_t merge_size, double fold_thres, int32_t dual_break_thres, int *bp_n)
{
    int i, j, k, n, m, d, b, b_n, b_m;
    double mcnt;
    int64_t *link;
    int32_t s, e, t, min_c, p;
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
        mcnt = n & 1? (int32_t) link[n / 2] : ((int32_t) link[n / 2] + (int32_t) link[n / 2 - 1]) / 2.0;
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
        bp1->p = (int32_t *) malloc(bp1->m * sizeof(int32_t));
        ++b_n;

        for (j = 0; j < kdq_size(q); ++j) {
            s = kdq_at(q, j) >> 32;
            e = (int32_t) kdq_at(q, j);
            
            if (e - s > d && make_dual_break(link, s, e, n - 1, 0.5)) {
                // dual break
                if (s != 0)
                    add_break_point(bp1, s);
                if (e != n - 1)
                    add_break_point(bp1, e);
                continue;
            }

            // skip the first and last block
            if (s == 0 || e == n - 1)
                continue;
            
            // find the position of the minimum value in this valley
            // which will be a break point
            min_c = INT32_MAX;
            p = -1;
            for (k = s; k <= e; ++k) {
                if ((int32_t) link[k] < min_c) {
                    min_c = (int32_t) link[k];
                    p = k;
                }
            }

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
                mcnt = t & 1? (int32_t) link[(e + s) / 2] : ((int32_t) link[(e + s) / 2] + (int32_t) link[(e + s) / 2 + 1]) / 2.0;
                mcnt *= fold_thres;
                qsort(link + s, t, sizeof(int64_t), pos_cmp);
                if ((int32_t) link[bp1->p[k]] > mcnt)
                    kdq_push(int64_t, q, k);
            }
            for (k = 0; k < kdq_size(q); ++k)
                bp1->p[kdq_at(q, k)] = -1;
        }
        
        j = 0;
        for (k = 0; k < bp1->n; ++k)
            if (bp1->p[k] >= 0)
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
    fprintf(fp, "%s[%d]:", dict->s[bp->s].name, bp->n);
    for (i = 0; i < bp->n; ++i)
        fprintf(fp, " %d", bp->p[i]);
    fprintf(fp, "\n");
}

static void write_segs_to_agp(sd_seg_t *segs, int n, sdict_t *sd, int s, FILE *fp)
{
    uint64_t len = 0;
    sd_seg_t seg;
    int i, t = 0;
    for (i = 0; i < n; ++i) {
        seg = segs[i];
        fprintf(fp, "scaffold_%u\t%lu\t%lu\t%u\tW\t%s\t%d\t%d\t%c\n", s, len + 1, len + seg.y, ++t, sd->s[seg.c >> 1].name, seg.x + 1, seg.x + seg.y, "+-"[seg.c & 1]);
        len += seg.y;
        if (i != n - 1) {
            fprintf(fp, "scaffold_%u\t%lu\t%lu\t%u\tN\t%d\tscaffold\tyes\tna\n", s, len + 1, len + GAP_SZ, ++t, GAP_SZ);
            len += GAP_SZ;
        }
    }
}

void write_agp(asm_dict_t *d, bp_t *breaks, int32_t b_n, FILE *fp)
{
    int i, j, L, l, s, ns, ms;
    uint64_t len;
    sd_seg_t seg, *segs;
    sdict_t *sd = d->sdict;
    int32_t *p, p_n;

    ns = 0;
    ms = 4096;
    segs = (sd_seg_t *) malloc(ms * sizeof(sd_seg_t));
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
                    sd_seg_t seg1 = {seg.s, len, seg.c, seg.x, seg.y};
                    segs[ns] = seg1;
                    if (ns == ms) {
                        ms <<= 1;
                        segs = (sd_seg_t *) realloc(segs, ms * sizeof(sd_seg_t));
                    }
                    ++ns;
                    len += seg.y;
                } else {
                    l = p[0] - seg.a;
                    if (l > 0) {
                        // split seg
                        sd_seg_t seg1 = {seg.s, len, seg.c, seg.c & 1? seg.x + seg.y - l : seg.x, l};
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
                        l = p[0] - p[-1];
                        sd_seg_t seg1 = {seg.s, 0, seg.c, seg.c & 1? seg.x + seg.y - L - l : seg.x + L, l};
                        write_segs_to_agp(&seg1, 1, sd, ++s, fp);
                        L += l;
                        if (--p_n) 
                            ++p;
                    }

                    if (L < seg.y) {
                        sd_seg_t seg1 = {seg.s, 0, seg.c, seg.c & 1? seg.x : seg.x + L, seg.y - L};
                        len += seg.y - L;
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

