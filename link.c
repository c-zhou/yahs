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
 * 23/06/21 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#include "khash.h"
#include "kvec.h"
#include "bamlite.h"
#include "sdict.h"
#include "enzyme.h"
#include "link.h"
#include "asset.h"

#undef DEBUG_NOISE
#undef DEBUG_ORIEN
#undef DEBUG_INTRA

#define USE_MEDIAN_NORM
#define MIN_RE_DENS .1
static uint32_t MAX_RADIUS = 100;

void *kopen(const char *fn, int *_fd);
int kclose(void *a);

KHASH_SET_INIT_STR(str)

void intra_link_mat_destroy(intra_link_mat_t *link_mat)
{
    uint32_t i;
    for (i = 0; i < link_mat->n; ++i)
        if (link_mat->links[i].n)
            free(link_mat->links[i].link);
    if (link_mat->links)
        free(link_mat->links);
    free(link_mat);
}

void inter_link_mat_destroy(inter_link_mat_t *link_mat)
{
    uint32_t i, j;
    for (i = 0; i < link_mat->n; ++i) {
        if (link_mat->links[i].n) {
            for (j = 0; j < 4; ++j) {
                free(link_mat->links[i].link[j]);
                free(link_mat->links[i].linkb[j]);
            }
        }
    }
    if (link_mat->links)
        free(link_mat->links);
    free(link_mat);
}

// index mapping (c0, c1) -> (i, j) -> (n * 2 - i - 3) * i / 2 + j - 1 
// 0  (0, 1)
// 1  (0, 2) 5  (1, 2)
// 2  (0, 3) 6  (1, 3) 9  (2, 3)
// 3  (0, 4) 7  (1, 4) 10 (2, 4) 12 (3, 4)
// 4  (0, 5) 8  (1, 5) 11 (2, 5) 13 (3, 5) 14 (4, 5)
// number of cells
// (b0, b1, r) -> n -> MIN(b0, r) * MIN(b1, r)
// index mapping for cells
// (i, j, r) -> p -> MAX(0, i - 1) * b1 + j
// reverse index mapping for cells
// (p, r) -> (i, j) -> (p / b1, p % b1)
//     # # # # # #
//     o x x x x x
//     o o x x x x
// b0  o o o x x x
//     o o o o x x
//     o o o o o x
//     o o o o o o
//     O o o o o o
//          b1
// TODO improve memory efficiency
// number of cells
// (b0, b1, r) -> n
// for (k = 0; k < MIN(b0, r); ++k) p += MIN(b1, r - k)
// index mapping for cells
// (i, j, r) -> p
// for (k = 0; k < i; ++k) p += MIN(b1, r - k)
// p += j
// reverse index mapping for cells
// (p, r) -> (i, j)
// ???
// (1) radius [15] >= b0 + b1 [14] >= b0 [8] >= b1 [6]
//     b0 * b1
//     o o o o o o
//     o o o o o o
//     o o o o o o
// b0  o o o o o o
//     o o o o o o
//     o o o o o o
//     o o o o o o
//     O o o o o o
//          b1
// (2) b0 + b1 [14] >= radius [10] >= b0 [8] >= b1 [6]
//     k = (b0 + b1 - radius)
//     b0 * b1 - k * (k - 1) / 2
//     o o o # # #
//     o o o o # #
//     o o o o o #
// b0  o o o o o o
//     o o o o o o
//     o o o o o o
//     o o o o o o
//     O o o o o o
//          b1
// (3) b0 [8] >= radius [7] >= b1 [6]
//     (r - b1) * 
//     # # # # # #
//     o # # # # #
//     o o # # # #
// b0  o o o # # #
//     o o o o # #
//     o o o o o #
//     o o o o o o
//     O o o o o o
//          b1
// (4) b0 [8] >= b1 [6] >= radius [4]
//     # # # # # #
//     # # # # # #
//     # # # # # #
// b0  # # # # # #
//     o # # # # #
//     o o # # # #
//     o o o # # #
//     O o o o # #
//          b1
// 
inter_link_mat_t *inter_link_mat_init(asm_dict_t *dict, uint32_t resolution, uint32_t radius)
{
    inter_link_mat_t *link_mat;
    inter_link_t *link;
    uint32_t i, j, k, b0, b1, n, m, p, r2;

    n = dict->n;
    m = (long) n * (n - 1) / 2;
    r2 = resolution * 2;
    link_mat = (inter_link_mat_t *) malloc(sizeof(inter_link_mat_t));
    link_mat->n = m;
    link_mat->r = radius;
    link_mat->links = (inter_link_t *) malloc(m * sizeof(inter_link_t));
    
    for (i = 0; i < n; ++i) {
        if (dict->s[i].len < r2) {
            for (j = i + 1; j < n; ++j)
                link_mat->links[(long) (n * 2 - i - 3) * i / 2 + j - 1].n = 0;
            continue;
        }
        
        b0 = MIN(radius, div_ceil(dict->s[i].len, r2));
        for (j = i + 1; j < n; ++j) {
            link = &link_mat->links[(long) (n * 2 - i - 3) * i / 2 + j - 1];
            if (dict->s[j].len < r2) {
                link->n = 0;
                continue;
            }
            
            b1 = MIN(radius, div_ceil(dict->s[j].len, r2));
            p = b0 * b1; 
            link->c0 = i;
            link->c1 = j;
            link->b0 = b0;
            link->b1 = b1;
            link->r = MIN(radius, (uint32_t) (b0 + b1 - 1));
            link->n = p;
            link->n0 = 0;
            link->linkt = 0;
            assert(p == (long) b0 * b1);
            for (k = 0; k < 4; ++k) {
                link->link[k] = (double *) calloc(p, sizeof(double));
                link->linkb[k] = (double *) calloc(link->r, sizeof(double));
                if (!link->link[k] || !link->linkb[k]) {
                    fprintf(stderr, "[E::%s] memory allocation failure\n", __func__);
                    exit(EXIT_FAILURE);
                }
            }

            memset(link->norms, 0, sizeof(link->norms));
        }
    }

    return link_mat;
}

long estimate_inter_link_mat_init_rss(asm_dict_t *dict, uint32_t resolution, uint32_t radius)
{
    long bytes, p, m;
    uint32_t i, j, b0, b1, n, r, r2;

    n = dict->n;
    m = (long) n * (n - 1) / 2;
    if (m > UINT32_MAX)
        return -1;

    r2 = resolution * 2;

    bytes = 0;
    bytes += sizeof(inter_link_mat_t);
    bytes += m * sizeof(inter_link_t);

    for (i = 0; i < n; ++i) {
        if (dict->s[i].len < r2)
            continue;

        b0 = MIN(radius, div_ceil(dict->s[i].len, r2));
        for (j = i + 1; j < n; ++j) {
            if (dict->s[j].len < r2)
                continue;

            b1 = MIN(radius, div_ceil(dict->s[j].len, r2));
            p = b0 * b1;
            if (p > UINT32_MAX)
                return -1;
            r = MIN(radius, (uint32_t) (b0 + b1 - 1));

            bytes += (p + r) * 4 * sizeof(double);
        }
    }

    return bytes;
}

// index mapping (i, j) -> ([n - (j - i)] + [n - 1]) / 2 * (j - i) + j -> (n * 2 - j + i - 1) * (j - i) / 2 + j
// distance k [0 ... n-1] start from (n * 2 - k - 1) * k / 2 + k
// 0  (0,0)
// 7  (0,1) 1  (1,1)
// 13 (0,2) 8  (1,2) 2  (2,2)
// 18 (0,3) 14 (1,3) 9  (2,3) 3  (3,3)
// 22 (0,4) 19 (1,4) 15 (2,4) 10 (3,4) 4  (4,4)
// 25 (0,5) 23 (1,5) 20 (2,5) 16 (3,5) 11 (4,5) 5  (5,5)
// 27 (0,6) 26 (1,6) 24 (2,6) 21 (3,6) 17 (4,6) 12 (5,6) 6  (6,6)
intra_link_mat_t *intra_link_mat_init(void *dict, uint32_t resolution, int use_gap_seq)
{
    intra_link_mat_t *link_mat;
    intra_link_t *link;
    uint32_t i, n, b, p;
    uint64_t len;

    n = use_gap_seq? ((asm_dict_t *) dict)->n : ((sdict_t *) dict)->n;
    link_mat = (intra_link_mat_t *) malloc(sizeof(intra_link_mat_t));
    link_mat->n = n;
    link_mat->links = (intra_link_t *) calloc(n, sizeof(intra_link_t));

    for (i = 0; i < n; ++i) {
        link = &link_mat->links[i];
        link->c = i;
        len = use_gap_seq? ((asm_dict_t *) dict)->s[i].len : ((sdict_t *) dict)->s[i].len;
        if (len < resolution) {
            link->n = 0;
            continue;
        }
        b = div_ceil(len, resolution);
        link->n = b;
        p = (long) b * (b + 1) / 2;
        link->link = (double *) calloc(p, sizeof(double));
        if (!link->link) {
            fprintf(stderr, "[E::%s] memory allocation failure\n", __func__);
            exit(EXIT_FAILURE);
        }

#ifdef DEBUG_INTRA
        fprintf(stderr, "[DEBUG_INTRA::%s] %s bins: %d\n", __func__, use_gap_seq? ((asm_dict_t *) dict)->s[i].name : ((sdict_t *) dict)->s[i].name, link->n);
#endif
    }

    return link_mat;
}

long estimate_intra_link_mat_init_rss(void *dict, uint32_t resolution, int use_gap_seq)
{
    long bytes, p;
    uint32_t i, n, b;
    uint64_t len;

    n = use_gap_seq? ((asm_dict_t *) dict)->n : ((sdict_t *) dict)->n;

    bytes = 0;
    bytes += sizeof(intra_link_mat_t);
    bytes += n * sizeof(intra_link_t);

    for (i = 0; i < n; ++i) {
        len = use_gap_seq? ((asm_dict_t *) dict)->s[i].len : ((sdict_t *) dict)->s[i].len;
        if (len < resolution)
            continue;
        b = div_ceil(len, resolution);
        p = (long) b * (b + 1) / 2;
        if (p > UINT32_MAX)
            return -1;

        bytes += p * sizeof(double);
    }
    return bytes;
}

static void intra_link_mat_normalise(intra_link_mat_t *link_mat, void *dict, re_cuts_t *re_cuts, uint32_t resolution, int use_gap_seq)
{
    intra_link_t *link;
    uint32_t i, j, k, n, b;
    uint64_t len;
    double a, a0, a1, s;
    double **re_dens;

    n = use_gap_seq? ((asm_dict_t *) dict)->n : ((sdict_t *) dict)->n;

    re_dens = use_gap_seq? calc_re_cuts_density1(re_cuts, resolution, dict) : calc_re_cuts_density(re_cuts, resolution);
    
    for (i = 0; i < n; ++i) {
        len = use_gap_seq? ((asm_dict_t *) dict)->s[i].len : ((sdict_t *) dict)->s[i].len;
        if (len < resolution) continue;
        
        link = &link_mat->links[i];
        b = link->n;
        // relative size of the last cell
        a = ((double) len - (double) (b - 1) * resolution) / resolution;

        for (j = 0; j < b; ++j) {
            a0 = j < b-1? 1. : a;
            for (k = j; k < b; ++k) {
                a1 = k < b-1? 1. : a;
                a1 *= a0;
                a1 = a1 < .5? .0 : a1;
                s = re_dens? re_dens[i][j] * re_dens[i][k] : 1.;
                s = s < MIN_RE_DENS? .0 : s;
                s *= a1;
                if (s < FLT_EPSILON)
                    link->link[(long) (b * 2 - k + j - 1) * (k - j) / 2 + k] = -DBL_MAX;
                else
                    link->link[(long) (b * 2 - k + j - 1) * (k - j) / 2 + k] /= s;
            }
        }
    }
    
    if (re_dens) {
        for (i = 0; i < n; ++i)
            free(re_dens[i]);
        free(re_dens);
    }
}

intra_link_mat_t *intra_link_mat_from_file(const char *f, cov_norm_t *cov_norm, asm_dict_t *dict, re_cuts_t *re_cuts, uint32_t resolution, int use_gap_seq, uint8_t mq)
{
    uint32_t i, k, m, i0, i1, b0, b1;
    uint64_t p0, p1, pair_n, pair_c, intra_c;
    int64_t magic_number;
    uint8_t buffer[BUFF_SIZE * 17];
    intra_link_mat_t *link_mat;
    intra_link_t *link;
    FILE *fp;

    fp = fopen(f, "r");
    if (fp == NULL)
        return 0;

    m = fread(&magic_number, sizeof(int64_t), 1, fp);
    if (!m || !is_valid_bin_header(magic_number)) {
        fprintf(stderr, "[E::%s] not a valid BIN file\n", __func__);
        return 0;
    }
    file_seek_skip_sdict(fp);
    m = fread(&pair_n, sizeof(uint64_t), 1, fp);

    link_mat = use_gap_seq? intra_link_mat_init(dict, resolution, use_gap_seq) : intra_link_mat_init(dict->sdict, resolution, use_gap_seq);

    pair_c = 0;
    intra_c = 0;

    while (pair_c < pair_n) {
        m = fread(buffer, sizeof(uint8_t), BUFF_SIZE * 17, fp);

        for (i = 0; i < m && pair_c < pair_n; i += 17, ++pair_c) {
            if (*(uint8_t *) (buffer + i + 16) < mq)
                continue;

            if (use_gap_seq) {
                sd_coordinate_conversion(dict, *(uint32_t *) (buffer + i),     *(uint32_t *) (buffer + i + 4),  &i0, &p0, 0);
                sd_coordinate_conversion(dict, *(uint32_t *) (buffer + i + 8), *(uint32_t *) (buffer + i + 12), &i1, &p1, 0);
                b0 = p0 / resolution;
                b1 = p1 / resolution;
            } else {
                i0 = *(uint32_t *) (buffer + i);
                i1 = *(uint32_t *) (buffer + i + 8);
                b0 = *(uint32_t *) (buffer + i + 4) / resolution;
                b1 = *(uint32_t *) (buffer + i + 12) / resolution;
            }

            if (i0 == i1 && i0 != UINT32_MAX) {
                link = &link_mat->links[i0];
                if (link->n) {
                    if (b0 > b1)
                        SWAP(uint32_t, b0, b1);
                    k = (long) (link->n * 2 - b1 + b0 - 3) * (b1 - b0) / 2 + b1;
                    // link->link[k] += 1.;
                    link->link[k] += cov_norm->norm[*(uint32_t *) (buffer + i)][*(uint32_t *) (buffer + i + 4) / cov_norm->w] *
                        cov_norm->norm[*(uint32_t *) (buffer + i + 8)][*(uint32_t *) (buffer + i + 12) / cov_norm->w];
                }

                ++intra_c;
            }
        }
    }
#ifdef DEBUG
    fprintf(stderr, "[DEBUG::%s] %lu read pairs processed, %lu intra links \n", __func__, pair_c, intra_c);
#endif
    fclose(fp);
    
    // TODO normalise links by cell size
    if (use_gap_seq)
        intra_link_mat_normalise(link_mat, dict, re_cuts, resolution, use_gap_seq);
    else
        intra_link_mat_normalise(link_mat, dict->sdict, re_cuts, resolution, use_gap_seq);

    return link_mat;
}

static void inter_link_mat_normalise(inter_link_mat_t *link_mat, asm_dict_t *dict, re_cuts_t *re_cuts, uint32_t resolution, uint32_t radius)
{
    inter_link_t *link;
    uint32_t i, j, k, l, b0, b1, n, p, r2;
    double a0, a1, a;
    double **re_dens, re;

    n = dict->n;
    r2 = resolution * 2;

    re_dens = calc_re_cuts_density2(re_cuts, resolution, dict);

    for (i = 0; i < n; ++i) {
        if (dict->s[i].len < r2) continue;
        
        b0 = MIN(radius, div_ceil(dict->s[i].len, r2));
        // relative size of the last cell
        a0 = MIN(1., (dict->s[i].len / 2. - (double) (b0 - 1) * resolution) / resolution);
        
        for (j = i + 1; j < n; ++j) {
            if (dict->s[j].len < r2) continue;
            
            b1 = MIN(radius, div_ceil(dict->s[j].len, r2));
            // relative size of the last cell
            a1 = MIN(1., (dict->s[j].len / 2. - (double) (b1 - 1) * resolution) / resolution);
            
            link = &link_mat->links[(long) (n * 2 - i - 3) * i / 2 + j - 1];
            p = link->n;
            // calculate relative areas for each cell
            for (l = 0; l < p; ++l) {
                a = 1.;
                if (l / b1 == b0 - 1)
                    a *= a0;
                if (l % b1 == b1 - 1)
                    a *= a1;
                if (a < .5)
                    a = .0;
                re = re_dens? re_dens[i][l / b1] * re_dens[j][l % b1] : 1.;
                a *= re < MIN_RE_DENS? .0 : re;
                if (a < FLT_EPSILON)
                    for (k = 0; k < 4; ++k) link->link[k][l] = -DBL_MAX;
                else
                    for (k = 0; k < 4; ++k) link->link[k][l] /= a;
            }
        }
    }

    if (re_dens) {
        for (i = 0; i < n; ++i)
            free(re_dens[i]);
        free(re_dens);
    }
}

inter_link_mat_t *inter_link_mat_from_file(const char *f, cov_norm_t *cov_norm, asm_dict_t *dict, re_cuts_t *re_cuts, uint32_t resolution, uint32_t radius, uint8_t mq)
{
    uint32_t i, j, k, m, n, i0, i1, b0, b1;
    uint64_t p0, p1, pair_n, pair_c, inter_c, radius_c, noise_c;
    int64_t magic_number;
    uint8_t buffer[BUFF_SIZE * 17];
    double l0, l1, c, a, na[4], nc[4];
    inter_link_mat_t *link_mat;
    inter_link_t *link;
    FILE *fp;

    fp = fopen(f, "r");
    if (fp == NULL)
        return 0;

    m = fread(&magic_number, sizeof(int64_t), 1, fp);
    if (!m || !is_valid_bin_header(magic_number)) {
        fprintf(stderr, "[E::%s] not a valid BIN file\n", __func__);
        return 0;
    }
    file_seek_skip_sdict(fp);
    m = fread(&pair_n, sizeof(uint64_t), 1, fp);

    n = dict->n;
    link_mat = inter_link_mat_init(dict, resolution, radius);
    pair_c = inter_c = radius_c = 0;

    while (pair_c < pair_n) {
        m = fread(buffer, sizeof(uint8_t), BUFF_SIZE * 17, fp);

        for (i = 0; i < m && pair_c < pair_n; i += 17, ++pair_c) {
            if (*(uint8_t *) (buffer + i + 16) < mq)
                continue;
            
            sd_coordinate_conversion(dict, *(uint32_t *) (buffer + i),     *(uint32_t *) (buffer + i + 4),  &i0, &p0, 0);
            sd_coordinate_conversion(dict, *(uint32_t *) (buffer + i + 8), *(uint32_t *) (buffer + i + 12), &i1, &p1, 0);

            if (i0 != i1 && i0 != UINT32_MAX) {
                c = cov_norm->norm[*(uint32_t *) (buffer + i)][*(uint32_t *) (buffer + i + 4) / cov_norm->w] *
                    cov_norm->norm[*(uint32_t *) (buffer + i + 8)][*(uint32_t *) (buffer + i + 12) / cov_norm->w];

                if (i0 > i1) {
                    SWAP(uint32_t, i0, i1);
                    SWAP(uint64_t, p0, p1);
                }

                link = &link_mat->links[(long) (n * 2 - i0 - 3) * i0 / 2 + i1 - 1];

                if (link->n == 0)
                    continue;

                l0 = dict->s[i0].len / 2.;
                l1 = dict->s[i1].len / 2.;
                
                if (p0 >= l0 && p1 < l1) {
                    // i0(-) -> i1(+)
                    b0 = (uint32_t) ((2 * l0 - p0) / resolution);
                    b1 = (uint32_t) ((double) p1 / resolution);
                    if (b0 < link->b0 && b1 < link->b1 && b0 + b1 < radius) {
                        k = (long) (MAX(1, b0) - 1) * link->b1 + b1;
                        // link->link[0][k] += 1.;
                        link->link[0][k] += c;
                        ++radius_c;
                    }
                } else if(p0 >= l0 && p1 >= l1) {
                    // i0(-) -> i1(-)
                    b0 = (uint32_t) ((2 * l0 - p0) / resolution);
                    b1 = (uint32_t) ((2 * l1 - p1) / resolution);
                    if (b0 < link->b0 && b1 < link->b1 && b0 + b1 < radius) {
                        k = (long) (MAX(1, b0) - 1) * link->b1 + b1;
                        // link->link[1][k] += 1.;
                        link->link[1][k] += c;
                        ++radius_c;
                    }
                } else if(p0 < l0 && p1 < l1) {
                    // i0(+) -> i1(+)
                    b0 = (uint32_t) ((double) p0 / resolution);
                    b1 = (uint32_t) ((double) p1 / resolution);
                    if (b0 < link->b0 && b1 < link->b1 && b0 + b1 < radius) {
                        k = (long) (MAX(1, b0) - 1) * link->b1 + b1;
                        // link->link[2][k] += 1.;
                        link->link[2][k] += c;
                        ++radius_c;
                    }
                } else if(p0 < l0 && p1 >= l1) {
                    // i0(+) -> i1(-)
                    b0 = (uint32_t) ((double) p0 / resolution);
                    b1 = (uint32_t) ((2 * l1 - p1) / resolution);
                    if (b0 < link->b0 && b1 < link->b1 && b0 + b1 < radius) {
                        k = (long) (MAX(1, b0) - 1) * link->b1 + b1;
                        // link->link[3][k] += 1.;
                        link->link[3][k] += c;
                        ++radius_c;
                    }
                }

                ++inter_c;
            }
        }
    }

#ifdef DEBUG
    fprintf(stderr, "[DEBUG::%s] %lu read pairs processed, %lu inter links \n", __func__, pair_c, inter_c);
    fprintf(stderr, "[DEBUG::%s] within radius %d: %lu\n", __func__, radius, radius_c);
#endif
    fclose(fp);

    // TODO normalise links by cell size
    inter_link_mat_normalise(link_mat, dict, re_cuts, resolution, radius);

    // calculate noise level
    noise_c = a = 0;
    for (i = 0; i < link_mat->n; ++i) {
        memset(nc, 0, sizeof(nc));
        memset(na, 0, sizeof(na));
        link = &link_mat->links[i];
        for (j = 0; j < link->n; ++j) {
            for (k = 0; k < 4; ++k) {
                if (link->link[k][j] >= 0) {
                    nc[k] += link->link[k][j];
                    na[k] += 1.;
                }
            }
        }
        
        j = 0;
        for (k = 1; k < 4; ++k)
            if (nc[k] < nc[j])
                j = k;
        noise_c += nc[j];
        a += na[j];
    }
    link_mat->noise = a > 0? noise_c / a : 0;
#ifdef DEBUG_NOISE
    fprintf(stderr, "[DEBUG_NOISE::%s] noise links: %lu; area: %.12f; noise estimation: %.12f\n", __func__, noise_c, a, link_mat->noise);
#endif

    return link_mat;
}

cov_norm_t *cov_norm_from_file(const char *f, sdict_t *dict, uint32_t window)
{
    uint32_t i;
    uint64_t m, n;
    double *norm_a, **norm;
    int64_t magic_number;
    cov_norm_t *cov_norm;
    FILE *fp;
    
    fp = fopen(f, "r");
    if (fp == NULL)
        return 0;

    m = fread(&magic_number, sizeof(int64_t), 1, fp);
    if (!m || !is_valid_bin_header(magic_number)) {
        fprintf(stderr, "[E::%s] not a valid BIN file\n", __func__);
        return 0;
    }
    file_seek_skip_sdict(fp);
    m = fread(&n, sizeof(uint64_t), 1, fp);
    fseek(fp, n * 17, SEEK_CUR);
    m = fread(&n, sizeof(uint64_t), 1, fp);
    norm_a = (double *) malloc(sizeof(double) * n);
    fread(norm_a, sizeof(double), n, fp);

    norm = (double **) malloc(sizeof(double *) * dict->n);
    norm[0] = norm_a;
    for (i = 1; i < dict->n; ++i)
        norm[i] = norm[i-1] + div_ceil(dict->s[i-1].len, window);
    
    fclose(fp);
    
    cov_norm = (cov_norm_t *) malloc(sizeof(cov_norm_t));
    cov_norm->n = n;
    cov_norm->w = window;
    cov_norm->norm = norm;

    return cov_norm;
}

void norm_destroy(norm_t *norm)
{
    //uint32_t i;
    free(norm->bs);
    //for (i = 0; i < norm->n; ++i)
    //    free(norm->link[i]);
    //free(norm->link);
    free(norm->linkc);
    free(norm->norms);
    free(norm);
}

int dcmp (const void *a, const void *b) {
    double cmp = *(double *) a - *(double *) b;
    return cmp > 0? 1 : ( cmp < 0? -1 : 0);
}

norm_t *calc_norms(intra_link_mat_t *link_mat, uint32_t d_min_cell, double d_mass_frac)
{
    uint32_t i, j, n, b, r, r0, t;
    uint32_t *bs;
    double intra_c, tmp_c;
    double *norms, *link, *linkc;
    norm_t *norm;
    
    // find maximum numebr of bands
    n = 0;
    for (i = 0; i < link_mat->n; ++i)
        n = MAX(link_mat->links[i].n, n);
    if (n <= 1) {
        fprintf(stderr, "[E::%s] no bands for norm calculation, try a higher resolution\n", __func__);
        return 0;
    }
    n -= 1; // do not use the last one as it might be an incomplete cell
    // calculate number of cells for each band
    bs = (uint32_t *) calloc(n, sizeof(uint32_t));
    for (i = 0; i < link_mat->n; ++i) {
        if (link_mat->links[i].n) {
            b = link_mat->links[i].n - 1;
            for (j = 0; j < b; ++j)
                bs[j] += b - j;
        }
    }

    r0 = 0;
    while (bs[r0] >= d_min_cell)
        ++r0;
    if (r0 < 10) {
        fprintf(stderr, "[E::%s] no enough bands (%d) for norm calculation, try a higher resolution\n", __func__, r0);
        free(bs);
        return 0;
    }

    // caluclate links in each band and radius
    // calculate norms - using median or mean?
    norms = (double *) malloc(n * sizeof(double));
    linkc = (double *) malloc(n * sizeof(double));
    intra_c = .0;
    for (i = 0; i < n; ++i) {
        link = (double *) malloc(bs[i] * sizeof(double));
        t = 0;
        for (j = 0; j < link_mat->n; ++j) {
            b = link_mat->links[j].n;
            if (b > i + 1) {
                memcpy(link + t, link_mat->links[j].link + (long) (b * 2 - i - 1) * i / 2 + i, (b - 1 - i) * sizeof(double));
                t += b - 1 - i;
            }
        }

        tmp_c = .0;
        for (j = 0; j < bs[i]; ++j)
            tmp_c += link[j];
        linkc[i] = tmp_c;
        intra_c += tmp_c;
        
#ifdef USE_MEDIAN_NORM
        qsort(link, bs[i], sizeof(double), dcmp);
        norms[i] = bs[i] & 1? link[bs[i] / 2] : (link[bs[i] / 2] + link[bs[i] / 2 - 1]) / 2;
#else
        norms[i] = MAX(linkc[i], 1.) / bs[i];
#endif

        free(link);
    }

    intra_c -= linkc[0];
    tmp_c = .0;
    for (r = 1; r < n; ++r) {
        tmp_c += linkc[r];
        if (tmp_c / intra_c >= d_mass_frac)
            break;
    }
    
    r = MIN(MIN(r, r0), MAX_RADIUS);

#ifdef USE_MEDIAN_NORM
    // adjust radius by norms
    // change radius if norm drops to zero before reaching r
    // find the first position with 10 consective zeros
    int zeros = 0;
    int p = -1;
    for (i = 0; i < n; ++i) {
        if (norms[i] < 1) {
            ++zeros;
            if (zeros >= 10) {
                p = i - 9;
                break;
            }
        } else {
            zeros = 0;
        }
    }
    if (p >= 0 && p < r)
        r = p;
    for (i = 0; i <= r; ++i) {
        if(norms[i] < 1)
            norms[i] = 1;
    }
#endif

    norm = (norm_t *) malloc(sizeof(norm_t));
    norm->n = n;
    norm->r = r;
    norm->bs = bs;
    //norm->link = link;
    norm->linkc = linkc;
    norm->norms = norms;
    return norm;
}

void print_inter_link_norms(FILE *fp, inter_link_mat_t *link_mat, asm_dict_t *dict)
{
    uint32_t i, n;
    inter_link_t *link;

    n = link_mat->n;
    for (i = 0; i < n; ++i) {
        link = &link_mat->links[i];
        if (link->n)
            fprintf(fp, "NORM: %s %s (%u) [%.3f %.3f %.3f %.3f]\n", dict->s[link->c0].name, dict->s[link->c1].name, link->n0, link->norms[0], link->norms[1], link->norms[2], link->norms[3]);
    }
}

void print_inter_link_bands(FILE *fp, inter_link_mat_t *link_mat, asm_dict_t *dict)
{
    uint32_t i, j, k, n;
    inter_link_t *link;
    
    n = link_mat->n;
    for (i = 0; i < n; ++i) {
        link = &link_mat->links[i];
        if (link->n) {
            for (j = 0; j < 4; ++j) {
                fprintf(fp, "BAND [%u %u]: %s%s %s%s (%u, %u) [%.3f]", i, j, dict->s[link->c0].name, j < 2? "-" : "+", dict->s[link->c1].name, j & 1? "-" : "+", link->n0, link->r, link->norms[j]);
                for (k = 0; k < link->r; ++k)
                    fprintf(fp, " %.0f", link->linkb[j][k]);
                fprintf(fp, "\n");
            }
        }
    }
}

double *get_max_inter_norms(inter_link_mat_t *link_mat, asm_dict_t *dict)
{
    uint32_t i, j, n, c0, c1;
    double *max_norms;
    inter_link_t *link;

    n = link_mat->n;
    max_norms = (double *) calloc(dict->n, sizeof(double));
    for (i = 0; i < n; ++i) {
        link = &link_mat->links[i];
        if (link->n) {
            c0 = link->c0;
            c1 = link->c1;
            for (j = 0; j < 4; ++j) {
                max_norms[c0] = MAX(max_norms[c0], link->norms[j]);
                max_norms[c1] = MAX(max_norms[c1], link->norms[j]);
            }
        }
    }
    return max_norms;
}

void print_norms(FILE *fp, norm_t *norm)
{
    uint32_t i/**, j, n**/;
    //uint64_t link_count;
    //double *link;

    //link_count = 0;

    fprintf(fp, "norms:");
    for (i = 0; i < norm->n; ++i)
        fprintf(fp, " %.1f", norm->norms[i]);
    fprintf(fp, "\n");
    
    fprintf(fp, "links:");
    for (i = 0; i < norm->n; ++i)
        fprintf(fp, " %.0f", norm->linkc[i]);
    fprintf(fp, "\n");

    /***
    for (i = 0; i < norm->n; ++i) {
        fprintf(fp, "band %d:", i);
        n = norm->bs[i];
        link = norm->link[i];
        for (j = 0; j < n; ++j) {
            link_count += link[j];
            fprintf(fp, " %.0f", link[j]);
        }
        fprintf(fp, "\n");
    }
    **/
}

int estimate_noise(inter_link_mat_t *link_mat)
{
    uint32_t i, j, k, n;
    uint32_t bin, *hist;
    uint64_t total;
    double l;
    inter_link_t *inter_link;

    bin = 4096;
    hist = (uint32_t *) calloc(bin, sizeof(uint32_t));
    for (i = 0; i < link_mat->n; ++i) {
        inter_link = &link_mat->links[i];
        n = inter_link->n;
        for (j = 0; j < n; ++j) {
            for (k = 0; k < 4; ++k) {
                l = inter_link->link[k][j];
                if (l < .0) 
                    continue;
                ++hist[MIN(bin - 1, (int) l)];
            }
        }
    }
   
    hist[0] = 0;

    total = 0;
    for (i = 0; i < bin; ++i)
        total += hist[i];
    total >>= 1;

#ifdef DEBUG_NOISE
    for (i = 0; i < bin; ++i)
        fprintf(stderr, "[DEBUG_NOISE::%s] noise bin %d %d\n", __func__, i, hist[i]);
#endif

    n = 0;
    while (hist[n] < total) {
        ++n;
        if (n == bin)
            break;
        hist[n] += hist[n - 1];
    }

    if (n == 4096)
        fprintf(stderr, "[W::%s] noise esimation might be incorrect %d\n", __func__, n);
    
    free(hist);

    return n;
}

int inter_link_norms(inter_link_mat_t *link_mat, norm_t *norm, int use_estimated_noise, double max_noise_ratio, double *la)
{
    uint32_t i, j, k, b, b1, r, n, n0;
    uint32_t *fn;
    int8_t bl;
    double l, *norms, *fr, noise;
    double c, c0, t;
    inter_link_t *inter_link;

    noise = .0;
    if (use_estimated_noise) {
        // noise = estimate_noise(link_mat);
        noise = link_mat->noise;
        fprintf(stderr, "[I::%s] using noise level %.3f\n", __func__, noise);
    }

    r = link_mat->r;
    norms = (double *) calloc(r + 1, sizeof(double));
    for (i = 0; i <= r; ++i) {
        if (norm->norms[i] * max_noise_ratio > noise) {
            norms[i] = norm->norms[i] - noise;
        } else {
            fprintf(stderr, "[I::%s] using new radius %u (%u) as noise ratio exceeds threshold %.3f\n", __func__, i, r, max_noise_ratio);
            r = i;
            break;
        }
    }
    if (r <= 1) return 1;

    /***
    // no weighting
    c = c0 = 0;
    for (i = 0; i < link_mat->n; ++i) {
        inter_link = &link_mat->links[i];
        n = inter_link->n;
        if (n == 0)
            continue;
        r = inter_link->r;
        b1 = inter_link->b1;
        n0 = 0;
        for (j = 0; j < n; ++j) {
            b = j / b1 + j % b1;
            if (b < r) {
                if (norms[b + 1] > 0) {
                    bl = 0;
                    for (k = 0; k < 4; ++k) {
                        l = inter_link->link[k][j];
                        if (l < .0)
                            break;
                        if (l >= norms[b + 1] + noise) {
                            inter_link->norms[k] += 1.;
                            ++c0;
                        }
                        bl = 1;
                    }
                    if (bl)
                        ++n0;
                }
            }
        }
        for (k = 0; k < 4; ++k)
            inter_link->norms[k] /= n0;
        inter_link->n0 = n0;
        c += n0;
    }
    **/
    
    c = c0 = 0;
    for (i = 0; i < link_mat->n; ++i) {
        inter_link = &link_mat->links[i];
        n = inter_link->n;
        if (n == 0)
            continue;
        r = inter_link->r;
        b1 = inter_link->b1;

        fr = (double *) calloc(4 * r, sizeof(double));
        fn = (uint32_t *) calloc(r, sizeof(uint32_t));

        for (j = 0; j < n; ++j) {
            b = j / b1 + j % b1;
            if (b < r) {
                if (norms[b + 1] > 0) {
                    bl = 0;
                    for (k = 0; k < 4; ++k) {
                        l = inter_link->link[k][j];
                        if(l < .0)
                            break;
                        l = MAX(.0, l - noise);
                        fr[k * r + b] += MIN(1., l / norms[b + 1]);
                        bl = 1;
                    }
                    if (bl)
                        ++fn[b];
                }
            }
        }
        if (fn[0] == 0)
            goto loop_i_end;

        // cumsum
        for (k = 0; k < 4; ++k)
            inter_link->norms[k] = fr[k * r];
        for (b = 1; b < r; ++b) {
            for (k = 0; k < 4; ++k) {
                inter_link->norms[k] += fr[k * r + b] * fr[k * r + b - 1] / fn[b - 1];
                fr[k * r + b] += fr[k * r + b - 1];
            }
            fn[b] += fn[b - 1];
        }

        n0 = fn[r - 1];
        t = 0;
        for (k = 0; k < 4; ++k) {
            if (inter_link->norms[k] > t)
                t = inter_link->norms[k];
            inter_link->norms[k] /= n0;
        }
        inter_link->n0 = n0;
        c += n0;
        c0 += t;
loop_i_end:
        free(fr);
        free(fn);
    }
    
    *la = c == 0? .0 : c0 / c;
    fprintf(stderr, "[I::%s] average link count: %.3f %.3f %.3f\n", __func__, c0, c, *la);
    
    free(norms);

    return 0;
}

int8_t *calc_link_directs_from_file(const char *f, asm_dict_t *dict, uint8_t mq)
{
    uint32_t i, j, k, m, n, na, i0, i1, b0, b1, b, ma, sma, n_ma;
    uint64_t p0, p1, pair_n, pair_c, inter_c;
    int64_t magic_number;
    uint8_t buffer[BUFF_SIZE * 17];
    uint32_t *link, l;
    int8_t *directs;
    FILE *fp;

    fp = fopen(f, "r");
    if (fp == NULL)
        return 0;
    
    m = fread(&magic_number, sizeof(int64_t), 1, fp);
    if (!m || !is_valid_bin_header(magic_number)) {
        fprintf(stderr, "[E::%s] not a valid BIN file\n", __func__);
        return 0;
    }
    file_seek_skip_sdict(fp);
    m = fread(&pair_n, sizeof(uint64_t), 1, fp);

    n = dict->n;
    na = (long) n * (n - 1) / 2;
    link = (uint32_t *) calloc(na << 2, sizeof(uint32_t));
    pair_c = inter_c = 0;

    while (pair_c < pair_n) {
        m = fread(buffer, sizeof(uint8_t), BUFF_SIZE * 17, fp);
        
        for (i = 0; i < m && pair_c < pair_n; i += 17, ++pair_c) {
            if (*(uint8_t *) (buffer + i + 16) < mq)
                continue;

            sd_coordinate_conversion(dict, *(uint32_t *) (buffer + i),     *(uint32_t *) (buffer + i + 4),  &i0, &p0, 0);
            sd_coordinate_conversion(dict, *(uint32_t *) (buffer + i + 8), *(uint32_t *) (buffer + i + 12), &i1, &p1, 0);

            if (i0 != i1 && i0 != UINT32_MAX) {
                if (i0 > i1) {
                    SWAP(uint32_t, i0, i1);
                    SWAP(uint64_t, p0, p1);
                }
                b0 = !((double) p0 / dict->s[i0].len > .5);
                b1 = !((double) p1 / dict->s[i1].len > .5);
                b = (b0 << 1) | (!b1);
                ++link[b * na + (long) (n * 2 - i0 - 1) * i0 / 2 + i1 - i0 - 1];
                
                ++inter_c;
            }
        }
    }

#ifdef DEBUG
    fprintf(stderr, "[DEBUG::%s] %lu read pairs processed, %lu inter links \n", __func__, pair_c, inter_c);
#endif
    fclose(fp);
    
    directs = (int8_t *) malloc(na * sizeof(int8_t));
    for (i = 0; i < na; ++i) {
        ma = sma = n_ma = j = 0;
        for (k = 0; k < 4; ++k) {
            l = link[k * na + i];
            if (l > ma) {
                sma = ma;
                ma = l;
                n_ma = 1;
                j = k;
            } else if (l == ma) {
                ++n_ma;
            } else if (l > sma) {
                sma = l;
            }
        }
        
        directs[i] = n_ma == 1 && ma * .9 > sma? j : -1;
    }
    free(link);
#ifdef DEBUG_ORIEN
    int8_t d;
    for (i0 = 0; i0 < n; ++i0) {
        for (i1 = i0 + 1; i1 < n; ++i1) {
            d = directs[(long) (n * 2 - i0 - 1) * i0 / 2 + i1 - i0 - 1];
            if (d >=0)
                fprintf(stderr, "[DEBUG_ORIEN::%s] #Oriens: %s %s %hhi \n", __func__, dict->s[i0].name, dict->s[i1].name, d);
        }
    }
#endif
    return directs;
}

void calc_link_directs(inter_link_mat_t *link_mat, double min_norm, asm_dict_t *dict, int8_t *directs)
{
    uint32_t i, j, k, n, b, b0, b1, r, n_ma;
    int8_t t;
    inter_link_t *inter_link;
    double area[4]; // area under the cumsum of linkb
    double ma, sma;
    uint32_t max_b = UINT32_MAX;

    r = link_mat->r;
    for (i = 0; i < link_mat->n; ++i) {
        inter_link = &link_mat->links[i];
        n = inter_link->n;
        if (n == 0)
            continue;
        b0 = inter_link->b0;
        b1 = inter_link->b1;
        for (j = 0; j < n; ++j) {
            b = j / b1 + j % b1;
            if (b <= b0 && b <= b1 && b <= max_b && b < r)
                for (k = 0; k < 4; ++k)
                    if (inter_link->link[k][j] > .0)
                        inter_link->linkb[k][b] += inter_link->link[k][j];
        }
    }

    for (i = 0; i < link_mat->n; ++i) {
        inter_link = &link_mat->links[i];
        if (inter_link->n == 0)
            continue;
        r = inter_link->r;
        memset(area, 0, sizeof(area));
        for (j = 0; j < r; ++j)
            for (k = 0; k < 4; ++k)
                area[k] += inter_link->linkb[k][j]; //* (r - j);
        // select the one with maximum cumsum area as linkt
        t = 0;
        ma = sma = 0;
        n_ma = 0;
        for (k = 0; k < 4; ++k) {
            if (area[k] > ma) {
                sma = ma;
                ma = area[k];
                n_ma = 1;
                t = k;
            } else if (area[k] == ma) {
                ++n_ma;
            } else if (area[k] > sma) {
                sma = area[k];
            }
        }
        if (n_ma == 1 && ma * .9 > sma) {
            inter_link->linkt = 1 << t;
        } else {
            t = 0;
            for (k = 0; k < 4; ++k)
                if (inter_link->norms[k] >= min_norm)
                    t |= (1 << k);
            inter_link->linkt = t;
        }

#ifdef DEBUG_ORIEN
        fprintf(stderr, "[DEBUG_ORIEN::%s] %d %.0f %.0f %d %d [%.3f %.3f %.3f %.3f] [%.0f  %.0f  %.0f  %.0f] %s %s\n", __func__, inter_link->linkt, ma, sma, inter_link->b0, inter_link->b1, inter_link->norms[0], inter_link->norms[1], inter_link->norms[2], inter_link->norms[3], area[0], area[1], area[2], area[3], dict->s[inter_link->c0].name, dict->s[inter_link->c1].name);
#endif
    }
}

static char *parse_bam_rec(bam1_t *b, bam_header_t *h, int32_t *s, int32_t *e, uint8_t *q, char **cname)
{
    // 0x4 0x100 0x200 0x400 0x800
    if (b->core.flag & 0xF04)
        return 0;
    *cname = h->target_name[b->core.tid];
    *s = b->core.pos;
    *e = get_target_end(b);
    *q = b->core.qual;

    return strdup(bam1_qname(b));
}

static int parse_bam_rec1(bam1_t *b, bam_header_t *h, char **cname0, int32_t *s0, int32_t *e0, char **cname1, int32_t *s1)
{
    // 0x4 0x8 [0x40] 0x100 0x200 0x400 0x800
    // if (b->core.flag & 0xF4C)
    if (b->core.flag & 0xF0C)
        return -1;

    *cname0 = h->target_name[b->core.tid];
    *s0 = b->core.pos;
    *e0 = get_target_end(b);
    *cname1 = h->target_name[b->core.mtid];
    *s1 = b->core.mpos;

    return (b->core.flag & 0x40)? 1 : 0;
}

void dump_links_from_bam_file(const char *f, const char *fai, uint32_t ml, uint8_t mq, uint32_t wd, double q_drop, const char *out)
{
    bamFile fp;
    FILE *fo;
    bam_header_t *h;
    bam1_t *b;
    char *cname0, *cname1, *rname0, *rname1;
    int32_t s0, s1, e0, e1, first;
    uint32_t i0, i1, p0, p1;
    uint8_t q, q0, q1;
    int8_t buff;
    uint64_t rec_c, pair_c, inter_c, intra_c, sd_l, n, m, max_m;
    enum bam_sort_order so;
    cov_t *cov;
    
    khash_t(str) *hmseq; // for absent sequences
    khint_t k;
    int absent;
    hmseq = kh_init(str);

    sdict_t *dict = make_sdict_from_index(fai, ml);
    
    fp = bam_open(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }
    
    fo = fopen(out, "w");
    if (fo == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }
    write_bin_header(fo);
    sd_l = write_sequence_dictionary(fo, dict);

    h = bam_header_read(fp);
    b = bam_init1();
    so = bam_hrecs_sort_order(h);
    cname0 = cname1 = rname0 = rname1 = 0;
    s0 = s1 = e0 = e1 = 0;
    i0 = i1 = p0 = p1 = 0;
    rec_c = pair_c = inter_c = intra_c = 0;
    buff = 0;
    cov = cov_init(dict->n);
    n = 0;
    max_m = 0x7FFFFFFULL; // 128MB - ~1GB mem for covs

    fwrite(&pair_c, sizeof(uint64_t), 1, fo);

    if (so == ORDER_NAME) {
        // sorted by read names
        while (bam_read1(fp, b) >= 0 ) {
            
            if (++rec_c % 1000000 == 0)
                fprintf(stderr, "[I::%s] %lu million records processed, %lu read pairs \n", __func__, rec_c / 1000000, pair_c);

            if (buff == 0) {
                q0 = 255;
                rname0 = parse_bam_rec(b, h, &s0, &e0, &q0, &cname0);
                if (!rname0)
                    continue;
                ++buff;
            } else if (buff == 1) {
                q1 = 255;
                rname1 = parse_bam_rec(b, h, &s1, &e1, &q1, &cname1);
                if (!rname1)
                    continue;
                if (strcmp(rname0, rname1) == 0) {
                    buff = 0;

                    i0 = sd_get(dict, cname0);
                    i1 = sd_get(dict, cname1);
                    q = MIN(q0, q1);

                    if (i0 == UINT32_MAX) {
                        k = kh_put(str, hmseq, cname0, &absent);
                        if (absent) {
                            kh_key(hmseq, k) = strdup(cname0);
                            fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname0);
                        }
                    } else if (i1 == UINT32_MAX) {
                        k = kh_put(str, hmseq, cname1, &absent);
                        if (absent) {
                            kh_key(hmseq, k) = strdup(cname1);
                            fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname1);
                        }
                    } else if (q >= mq) {
                        // update read coverage
                        // fix bad alignments
                        s0 = MIN_MAX(s0, 0, dict->s[i0].len - 1);
                        e0 = MIN_MAX(e0, 0, dict->s[i0].len - 1);
                        s1 = MIN_MAX(s1, 0, dict->s[i1].len - 1);
                        e1 = MIN_MAX(e1, 0, dict->s[i1].len - 1);
                        kv_push(uint64_t, cov->p[i0], (uint64_t) s0 << 32 | (uint32_t) 1);
                        kv_push(uint64_t, cov->p[i0], (uint64_t) e0 << 32 | (uint32_t) -1);
                        kv_push(uint64_t, cov->p[i1], (uint64_t) s1 << 32 | (uint32_t) 1);
                        kv_push(uint64_t, cov->p[i1], (uint64_t) e1 << 32 | (uint32_t) -1);
                        n += 4;
                        if (n > max_m) {
                            m = pos_compression(cov); // m is the position size after compression
                            fprintf(stderr, "[I::%s] position compression n = %lu, m = %lu, max_m = %lu\n", __func__, n, m, max_m);
                            if (m > n>>1) {
                                max_m <<= 1; // increase m limit if compression ratio smaller than 0.5
                                fprintf(stderr, "[I::%s] position memory buffer expanded max_m = %lu\n", __func__, max_m);
                            }
                            n = m;
                        }

                        // from zero-based to one-based
                        p0 = s0 / 2 + e0 / 2 + (s0 & 1 && e0 & 1);
                        p1 = s1 / 2 + e1 / 2 + (s1 & 1 && e1 & 1);
                        if (i0 > i1) {
                            SWAP(uint32_t, i0, i1);
                            SWAP(uint32_t, p0, p1);
                        }
                        fwrite(&i0, sizeof(uint32_t), 1, fo);
                        fwrite(&p0, sizeof(uint32_t), 1, fo);
                        fwrite(&i1, sizeof(uint32_t), 1, fo);
                        fwrite(&p1, sizeof(uint32_t), 1, fo);
                        fwrite(&q, sizeof(uint8_t), 1, fo);

                        if (i0 == i1)
                            ++intra_c;
                        else
                            ++inter_c;

                        ++pair_c;
                    }
                    
                    free(rname0);
                    free(rname1);
                    rname0 = 0;
                    rname1 = 0;
                } else {
                    buff = 1;

                    cname0 = cname1;
                    s0 = s1;
                    e0 = e1;
                    q0 = q1;
                    free(rname0);
                    rname0 = rname1;
                    rname1 = 0;
                }
            }
        }
    } else {
        // sorted by coordinates or others
        if (mq > 0)
            fprintf(stderr, "[W::%s] BAM file is not sorted by read name. Filtering by mapping quality %hhu suppressed \n", __func__, mq);
        q = 255;

        while (bam_read1(fp, b) >= 0) {
            
            if (++rec_c % 1000000 == 0)
                fprintf(stderr, "[I::%s] %lu million records processed, %lu read pairs \n", __func__, rec_c / 1000000, pair_c);

            first = parse_bam_rec1(b, h, &cname0, &s0, &e0, &cname1, &s1);

            if (first < 0) continue; // read not paired or not primary

            i0 = sd_get(dict, cname0);
            if (i0 != UINT32_MAX) {
                // update read coverage
                // fix bad alignments
                s0 = MIN_MAX(s0, 0, dict->s[i0].len - 1);
                e0 = MIN_MAX(e0, 0, dict->s[i0].len - 1);
                kv_push(uint64_t, cov->p[i0], (uint64_t) s0 << 32 | (uint32_t) 1);
                kv_push(uint64_t, cov->p[i0], (uint64_t) e0 << 32 | (uint32_t) -1);
                n += 2;
                if (n > max_m) {
                    m = pos_compression(cov); // m is the position size after compression
                    fprintf(stderr, "[I::%s] position compression n = %lu, m = %lu, max_m = %lu\n", __func__, n, m, max_m);
                    if (m > n>>1) {
                        max_m <<= 1; // increase m limit if compression ratio smaller than 0.5
                        fprintf(stderr, "[I::%s] position memory buffer expanded max_m = %lu\n", __func__, max_m);
                    }
                    n = m;
                }
            }

            if (first == 0) { // second read only for each read pair
                // i0 = sd_get(dict, cname0);
                i1 = sd_get(dict, cname1);

                if (i0 == UINT32_MAX) {
                    k = kh_put(str, hmseq, cname0, &absent);
                    if (absent) {
                        kh_key(hmseq, k) = strdup(cname0);
                        fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname0);
                    }
                } else if (i1 == UINT32_MAX) {
                    k = kh_put(str, hmseq, cname1, &absent);
                    if (absent) {
                        kh_key(hmseq, k) = strdup(cname1);
                        fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname1);
                    }
                } else {
                    // fix bad alignments
                    s1 = MIN_MAX(s1, 0, dict->s[i1].len - 1);
                    if (i0 > i1) {
                        SWAP(uint32_t, i0, i1);
                        SWAP(uint32_t, s0, s1);
                    }
                    fwrite(&i0, sizeof(uint32_t), 1, fo);
                    fwrite(&s0, sizeof(uint32_t), 1, fo);
                    fwrite(&i1, sizeof(uint32_t), 1, fo);
                    fwrite(&s1, sizeof(uint32_t), 1, fo);
                    fwrite(&q, sizeof(uint8_t), 1, fo);

                    if (i0 == i1)
                        ++intra_c;
                    else
                        ++inter_c;

                    ++pair_c;
                }
            }
        }
    }

    m = pos_compression(cov);
    fprintf(stderr, "[I::%s] position compression n = %lu, m = %lu, max_m = %lu\n", __func__, n, m, max_m);

    for (k = 0; k < kh_end(hmseq); ++k)
        if (kh_exist(hmseq, k))
            free((char *) kh_key(hmseq, k));
    kh_destroy(str, hmseq);

    if (rname0)
        free(rname0);
    if (rname1)
        free(rname1);
    bam_destroy1(b);
    bam_header_destroy(h);
    bam_close(fp);
    // write pair number
    fseek(fo, sizeof(uint64_t) + sd_l, SEEK_SET);
    fwrite(&pair_c, sizeof(uint64_t), 1, fo);
    fseek(fo, pair_c * 17, SEEK_CUR);

    fprintf(stderr, "[I::%s] dumped %lu read pairs from %lu records: %lu intra links + %lu inter links \n", __func__, pair_c, rec_c, intra_c, inter_c);

    // cov_t *cov = bam_cstats(f, dict, 0);
    cov_norm_t *cov_norm = calc_cov_norms(cov, dict, wd, q_drop);
    fwrite(&cov_norm->n, sizeof(uint64_t), 1, fo);
    fwrite(cov_norm->norm[0], sizeof(double), cov_norm->n, fo);
    cov_destroy(cov);
    cov_norm_destroy(cov_norm);

    sd_destroy(dict);
    fclose(fo);
}

void dump_links_from_bed_file(const char *f, const char *fai, uint32_t ml, uint8_t mq, uint32_t wd, double q_drop, const char *out)
{
    FILE *fp, *fo;
    int fd;
    void *fh;
    char *line = NULL;
    size_t ln = 0;
    ssize_t read;
    char cname0[4096], cname1[4096], rname0[4096], rname1[4096];
    uint32_t s0, s1, e0, e1, i0, i1, p0, p1;
    uint8_t q, q0, q1;
    int8_t buff;
    uint64_t rec_c, pair_c, inter_c, intra_c, sd_l, n, m, max_m;
    cov_t *cov;

    khash_t(str) *hmseq; // for absent sequences
    khint_t k;
    int absent;
    hmseq = kh_init(str);

    sdict_t *dict = make_sdict_from_index(fai, ml);

    fh = kopen(f, &fd);
    fp = fdopen(fd, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    fo = fopen(out, "w");
    if (fo == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }
    write_bin_header(fo);
    sd_l = write_sequence_dictionary(fo, dict);

    rec_c = pair_c = inter_c = intra_c = 0;
    buff = 0;
    cov = cov_init(dict->n);
    n = 0;
    max_m = 0x7FFFFFFULL; // 128MB - ~1GB mem for covs

    fwrite(&pair_c, sizeof(uint64_t), 1, fo);

    while ((read = getline(&line, &ln, fp)) != -1) {
        
        if (++rec_c % 1000000 == 0)
            fprintf(stderr, "[I::%s] %lu million records processed, %lu read pairs \n", __func__, rec_c / 1000000, pair_c);
    
        if (buff == 0) {
            cname0[0] = rname0[0] = '\0';
            s0 = e0 = UINT32_MAX;
            q0 = 255;
            sscanf(line, "%s %u %u %s %hhu", cname0, &s0, &e0, rname0, &q0);
            ++buff;
        } else if (buff == 1) {
            cname1[0] = rname1[0] = '\0';
            s1 = e1 = UINT32_MAX;
            q1 = 255;
            sscanf(line, "%s %u %u %s %hhu", cname1, &s1, &e1, rname1, &q1);
            if (is_read_pair(rname0, rname1)) {
                buff = 0;

                i0 = sd_get(dict, cname0);
                i1 = sd_get(dict, cname1);
                q = MIN(q0, q1);

                if (i0 == UINT32_MAX) {
                    k = kh_put(str, hmseq, cname0, &absent);
                    if (absent) {
                        kh_key(hmseq, k) = strdup(cname0);
                        fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname0);
                    }
                } else if (i1 == UINT32_MAX) {
                    k = kh_put(str, hmseq, cname1, &absent);
                    if (absent) {
                        kh_key(hmseq, k) = strdup(cname1);
                        fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname1);
                    }
                } else if (q >= mq) {
                   // update read coverage
                   // fix bad alignments
                   s0 = MIN_MAX(s0, 0, dict->s[i0].len - 1);
                   e0 = MIN_MAX(e0, 0, dict->s[i0].len - 1);
                   s1 = MIN_MAX(s1, 0, dict->s[i1].len - 1);
                   e1 = MIN_MAX(e1, 0, dict->s[i1].len - 1);
                   kv_push(uint64_t, cov->p[i0], (uint64_t) s0 << 32 | (uint32_t) 1);
                   kv_push(uint64_t, cov->p[i0], (uint64_t) e0 << 32 | (uint32_t) -1);
                   kv_push(uint64_t, cov->p[i1], (uint64_t) s1 << 32 | (uint32_t) 1);
                   kv_push(uint64_t, cov->p[i1], (uint64_t) e1 << 32 | (uint32_t) -1);
                   n += 4;
                   if (n > max_m) {
                       m = pos_compression(cov); // m is the position size after compression
                       fprintf(stderr, "[I::%s] position compression n = %lu, m = %lu, max_m = %lu\n", __func__, n, m, max_m);
                       if (m > n>>1) {
                           max_m <<= 1; // increase m limit if compression ratio smaller than 0.5
                           fprintf(stderr, "[I::%s] position memory buffer expanded max_m = %lu\n", __func__, max_m);
                       }
                       n = m;
                   }

                    p0 = s0 / 2 + e0 / 2 + (s0 & 1 && e0 & 1);
                    p1 = s1 / 2 + e1 / 2 + (s1 & 1 && e1 & 1);
                    if (i0 > i1) {
                        SWAP(uint32_t, i0, i1);
                        SWAP(uint32_t, p0, p1);
                    }
                    fwrite(&i0, sizeof(uint32_t), 1, fo);
                    fwrite(&p0, sizeof(uint32_t), 1, fo);
                    fwrite(&i1, sizeof(uint32_t), 1, fo);
                    fwrite(&p1, sizeof(uint32_t), 1, fo);
                    fwrite(&q, sizeof(uint8_t), 1, fo);
                
                    if (i0 == i1)
                        ++intra_c;
                    else
                        ++inter_c;

                    ++pair_c;
                }
            } else {
                buff = 1;

                strcpy(cname0, cname1);
                s0 = s1;
                e0 = e1;
                q0 = q1;
                strcpy(rname0, rname1);
            }
        }
    }

    m = pos_compression(cov);
    fprintf(stderr, "[I::%s] position compression n = %lu, m = %lu, max_m = %lu\n", __func__, n, m, max_m);

    for (k = 0; k < kh_end(hmseq); ++k)
        if (kh_exist(hmseq, k))
            free((char *) kh_key(hmseq, k));
    kh_destroy(str, hmseq);

    if (line)
        free(line);
    fclose(fp);
    kclose(fh);
    // write pair number
    fseek(fo, sizeof(uint64_t) + sd_l, SEEK_SET);
    fwrite(&pair_c, sizeof(uint64_t), 1, fo);
    fseek(fo, pair_c * 17, SEEK_CUR);

    fprintf(stderr, "[I::%s] dumped %lu read pairs from %lu records: %lu intra links + %lu inter links \n", __func__, pair_c, rec_c, intra_c, inter_c);

    // cov_t *cov = bed_cstats(f, dict);
    cov_norm_t *cov_norm = calc_cov_norms(cov, dict, wd, q_drop);
    fwrite(&cov_norm->n, sizeof(uint64_t), 1, fo);
    fwrite(cov_norm->norm[0], sizeof(double), cov_norm->n, fo);
    cov_destroy(cov);
    cov_norm_destroy(cov_norm);

    sd_destroy(dict);
    fclose(fo);
}

void dump_links_from_pa5_file(const char *f, const char *fai, uint32_t ml, uint8_t mq, uint32_t rl, uint32_t wd, double q_drop, const char *out)
{
    FILE *fp, *fo;
    int fd;
    void *fh;
    char *line = NULL;
    size_t ln = 0;
    ssize_t read;
    char cname0[4096], cname1[4096];
    uint32_t s0, s1, e0, e1, i0, i1, p0, p1;
    uint8_t q, q0, q1;
    uint64_t rec_c, pair_c, inter_c, intra_c, sd_l, n, m, max_m;
    cov_t *cov;

    khash_t(str) *hmseq; // for absent sequences
    khint_t k;
    int absent;
    hmseq = kh_init(str);

    sdict_t *dict = make_sdict_from_index(fai, ml);

    fh = kopen(f, &fd);
    fp = fdopen(fd, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    fo = fopen(out, "w");
    if (fo == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }
    write_bin_header(fo);
    sd_l = write_sequence_dictionary(fo, dict);

    rl >>= 1;
    rec_c = pair_c = inter_c = intra_c = 0;
    cov = cov_init(dict->n);
    n = 0;
    max_m = 0x7FFFFFFULL; // 128MB - ~1GB mem for covs

    fwrite(&pair_c, sizeof(uint64_t), 1, fo);

    while ((read = getline(&line, &ln, fp)) != -1) {

        if (++rec_c % 1000000 == 0)
            fprintf(stderr, "[I::%s] %lu million records processed, %lu read pairs \n", __func__, rec_c / 1000000, pair_c);

        cname0[0] = cname1[0] = '\0';
        p0 = p1 = UINT32_MAX;
        q0 = q1 = 255;
        sscanf(line, "%*s %s %u %s %u %hhu %hhu", cname0, &p0, cname1, &p1, &q0, &q1);

        i0 = sd_get(dict, cname0);
        i1 = sd_get(dict, cname1);
        q = MIN(q0, q1);

        if (i0 == UINT32_MAX) {
            k = kh_put(str, hmseq, cname0, &absent);
            if (absent) {
                kh_key(hmseq, k) = strdup(cname0);
                fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname0);
            }
        } else if (i1 == UINT32_MAX) {
            k = kh_put(str, hmseq, cname1, &absent);
            if (absent) {
                kh_key(hmseq, k) = strdup(cname1);
                fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname1);
            }
        } else if (q >= mq) {
            // update read coverage
            // fix bad alignments
            s0 = MIN_MAX(p0 - MIN(p0, rl), 0, dict->s[i0].len - 1);
            e0 = MIN_MAX(p0 + rl, 0, dict->s[i0].len - 1);
            s1 = MIN_MAX(p1 - MIN(p1, rl), 0, dict->s[i1].len - 1);
            e1 = MIN_MAX(p1 + rl, 0, dict->s[i1].len - 1);
            kv_push(uint64_t, cov->p[i0], (uint64_t) s0 << 32 | (uint32_t) 1);
            kv_push(uint64_t, cov->p[i0], (uint64_t) e0 << 32 | (uint32_t) -1);
            kv_push(uint64_t, cov->p[i1], (uint64_t) s1 << 32 | (uint32_t) 1);
            kv_push(uint64_t, cov->p[i1], (uint64_t) e1 << 32 | (uint32_t) -1);
            n += 4;
            if (n > max_m) {
                m = pos_compression(cov); // m is the position size after compression
                fprintf(stderr, "[I::%s] position compression n = %lu, m = %lu, max_m = %lu\n", __func__, n, m, max_m);
                if (m > n>>1) {
                    max_m <<= 1; // increase m limit if compression ratio smaller than 0.5
                    fprintf(stderr, "[I::%s] position memory buffer expanded max_m = %lu\n", __func__, max_m);
                }
                n = m;
            }

            p0 = MIN(p0, dict->s[i0].len - 1);
            p1 = MIN(p1, dict->s[i1].len - 1);
            if (i0 > i1) {
                SWAP(uint32_t, i0, i1);
                SWAP(uint32_t, p0, p1);
            }
            fwrite(&i0, sizeof(uint32_t), 1, fo);
            fwrite(&p0, sizeof(uint32_t), 1, fo);
            fwrite(&i1, sizeof(uint32_t), 1, fo);
            fwrite(&p1, sizeof(uint32_t), 1, fo);
            fwrite(&q, sizeof(uint8_t), 1, fo);

            if (i0 == i1)
                ++intra_c;
            else
                ++inter_c;

            ++pair_c;
        }
    }

    m = pos_compression(cov);
    fprintf(stderr, "[I::%s] position compression n = %lu, m = %lu, max_m = %lu\n", __func__, n, m, max_m);

    for (k = 0; k < kh_end(hmseq); ++k)
        if (kh_exist(hmseq, k))
            free((char *) kh_key(hmseq, k));
    kh_destroy(str, hmseq);

    if (line)
        free(line);
    fclose(fp);
    kclose(fh);
    // write pair number
    fseek(fo, sizeof(uint64_t) + sd_l, SEEK_SET);
    fwrite(&pair_c, sizeof(uint64_t), 1, fo);
    fseek(fo, pair_c * 17, SEEK_CUR);

    fprintf(stderr, "[I::%s] dumped %lu read pairs from %lu records: %lu intra links + %lu inter links \n", __func__, pair_c, rec_c, intra_c, inter_c);

    // cov_t *cov = bed_cstats(f, dict);
    cov_norm_t *cov_norm = calc_cov_norms(cov, dict, wd, q_drop);
    fwrite(&cov_norm->n, sizeof(uint64_t), 1, fo);
    fwrite(cov_norm->norm[0], sizeof(double), cov_norm->n, fo);
    cov_destroy(cov);
    cov_norm_destroy(cov_norm);

    sd_destroy(dict);
    fclose(fo);
}

