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

#include "bamlite.h"
#include "sdict.h"
#include "link.h"
#include "asset.h"

#undef DEBUG
#undef DEBUG_NOISE
#undef DEBUG_ORIEN

static int USE_MEDIAN_NORM = 1;
static uint32_t MAX_RADIUS = 100;

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

static inline uint32_t div_ceil(uint32_t x, uint32_t y)
{
    return 1 + ((MAX(x, 1) - 1) / y);
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
    uint32_t i, j, k, l, b0, b1, n, m, p, r2;
    double a0, a1, a;

    n = dict->n;
    m = n * (n - 1) / 2;
    r2 = resolution * 2;
    link_mat = (inter_link_mat_t *) malloc(sizeof(inter_link_mat_t));
    link_mat->n = m;
    link_mat->r = radius;
    link_mat->links = (inter_link_t *) malloc(m * sizeof(inter_link_t));

    for (i = 0; i < n; ++i) {
        if (dict->s[i].len < r2) {
            for (j = i + 1; j < n; ++j)
                link_mat->links[(n * 2 - i - 3) * i / 2 + j - 1].n = 0;
            continue;
        }
        
        b0 = MIN(radius, div_ceil(dict->s[i].len, r2));
        // relative edge length of the last cell
        a0 = MIN(1.0, (dict->s[i].len / 2.0 - (b0 - 1) * resolution) / resolution);
        for (j = i + 1; j < n; ++j) {
            link = &link_mat->links[(n * 2 - i - 3) * i / 2 + j - 1];
            if (dict->s[j].len < r2) {
                link->n = 0;
                continue;
            }
            
            b1 = MIN(radius, div_ceil(dict->s[j].len, r2));
            // relative edge length of the last cell
            a1 = MIN(1.0, (dict->s[j].len / 2.0 - (b1 - 1) * resolution) / resolution);
            p = b0 * b1; 
            link->c0 = i;
            link->c1 = j;
            link->b0 = b0;
            link->b1 = b1;
            link->r = MIN(radius, (uint32_t) (b0 + b1 - 1));
            link->n = p;
            link->n0 = 0;
            link->linkt = 0;
            
            for (k = 0; k < 4; ++k) {
                link->link[k] = (double *) calloc(p, sizeof(double));
                link->linkb[k] = (double *) calloc(link->r, sizeof(double));
            }

            // calculate relative areas for each cell 
            for (l = 0; l < p; ++l) {
                a = 1.0;
                if (l / b1 == b0 - 1)
                    a *= a0;
                if (l % b1 == b1 - 1)
                    a *= a1;
                if (a < 0.5)
                    a = -DBL_MAX;
                for (k = 0; k < 4; ++k)
                    link->link[k][l] = a;
            }
            memset(link->norms, 0, sizeof(link->norms));
        }
    }
    
    return link_mat;
}

long estimate_inter_link_mat_init_rss(asm_dict_t *dict, uint32_t resolution, uint32_t radius)
{
    long bytes;
    uint32_t i, j, b0, b1, n, m, p, r, r2;

    n = dict->n;
    m = n * (n - 1) / 2;
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
            r = MIN(radius, (uint32_t) (b0 + b1 - 1));

            bytes += 4 * (p + r) * sizeof(double);
        }
    }

    return bytes;
}

intra_link_mat_t *intra_link_mat_init(asm_dict_t *dict, uint32_t resolution)
{
    intra_link_mat_t *link_mat;
    intra_link_t *link;
    uint32_t i, j, n, b, p;
    double a;

    n = dict->n;
    link_mat = (intra_link_mat_t *) malloc(sizeof(intra_link_mat_t));
    link_mat->n = n;
    link_mat->links = (intra_link_t *) calloc(n, sizeof(intra_link_t));

    for (i = 0; i < n; ++i) {
        link = &link_mat->links[i];
        link->c = i;
        if (dict->s[i].len < resolution) {
            link->n = 0;
            continue;
        }
        b = div_ceil(dict->s[i].len, resolution);
        link->n = b;
        // relative edge length of the last cell
        a = ((double) dict->s[i].len - (b - 1) * resolution) / resolution;
        p = b * (b + 1) / 2;
        link->link = (double *) calloc(p, sizeof(double));
        for (j = 0; j < p; ++j)
            link->link[j] = 1.0;
        // calculate relative area for each cell
        for (j = 0; j < b; ++j)
            link->link[(b + j) * (b - j - 1) / 2 + b - 1] *= a;
        // the last triangle in the diagonal
        link->link[b - 1] *= a;
        // need a least half size to count
        for (j = 0; j < p; ++j)
            if (link->link[j] < 0.5)
                link->link[j] = -DBL_MAX;
#ifdef DEBUG
        printf("[I::%s] %s bins: %d\n", __func__, dict->s[i].name, link->n);
#endif
    }

    return link_mat;
}

long estimate_intra_link_mat_init_rss(asm_dict_t *dict, uint32_t resolution)
{
    long bytes;
    uint32_t i, n, b, p;

    n = dict->n;
    
    bytes = 0;
    bytes += sizeof(intra_link_mat_t);
    bytes += n * sizeof(intra_link_t);

    for (i = 0; i < n; ++i) {
        if (dict->s[i].len < resolution)
            continue;
        b = div_ceil(dict->s[i].len, resolution);
        p = b * (b + 1) / 2;
        bytes += p * sizeof(double);
    }
    return bytes;
}

intra_link_mat_t *intra_link_mat_init_sdict(sdict_t *dict, uint32_t resolution)
{
    intra_link_mat_t *link_mat;
    intra_link_t *link;
    uint32_t i, j, n, b, p;
    double a;

    n = dict->n;
    link_mat = (intra_link_mat_t *) malloc(sizeof(intra_link_mat_t));
    link_mat->n = n;
    link_mat->links = (intra_link_t *) malloc(n * sizeof(intra_link_t));

    for (i = 0; i < n; ++i) {
        link = &link_mat->links[i];
        link->c = i;
        if (dict->s[i].len < resolution) {
            link->n = 0;
            continue;
        }
        b = div_ceil(dict->s[i].len, resolution);
        link->n = b;
        // relative edge length of the last cell
        a = ((double) dict->s[i].len - (b - 1) * resolution) / resolution;
        p = b * (b + 1) / 2;
        link->link = (double *) calloc(p, sizeof(double));
        for (j = 0; j < p; ++j)
            link->link[j] = 1.0;
        // calculate relative area for each cell
        for (j = 0; j < b; ++j)
            link->link[(b + j) * (b - j - 1) / 2 + b - 1] *= a;
        // the last triangle in the diagonal
        link->link[b - 1] *= a;
        // need a least half size to count
        for (j = 0; j < p; ++j)
            if (link->link[j] < 0.5)
                link->link[j] = -DBL_MAX;
#ifdef DEBUG
        printf("[I::%s] %s bins: %d\n", __func__, dict->s[i].name, link->n);
#endif
    }

    return link_mat;
}

long estimate_intra_link_mat_init_sdict_rss(sdict_t *dict, uint32_t resolution)
{
    long bytes;
    uint32_t i, n, b, p;

    n = dict->n;
    
    bytes = 0;
    bytes += sizeof(intra_link_mat_t);
    bytes += n * sizeof(intra_link_t);

    for (i = 0; i < n; ++i) {
        if (dict->s[i].len < resolution)
            continue;
        b = div_ceil(dict->s[i].len, resolution);
        p = b * (b + 1) / 2;
        bytes += p * sizeof(double);
    }
    return bytes;
}


static inline void normalise_by_size(double *l)
{
    if (*l < 0) 
        return;
    double f, i;
    // f is the float part, i.e., area
    // i is the integer part, i.e., links
    f = modf(*l, &i);
    if (f == 0.0) {
        f = 1.0;
        i -= 1.0;
    }
    assert(i >= 0);
    // no normalisation for small cells to avoid extreme cases
    if (f >= 0.1) 
        i /= f;
    *l = i;
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
intra_link_mat_t *intra_link_mat_from_file(const char *f, asm_dict_t *dict, uint32_t resolution, int use_gap_seq)
{
    uint32_t i, j, n, b0, b1;
    uint32_t buffer[BUFF_SIZE], m, i0, i1, p0, p1;
    long pair_c, intra_c;
    intra_link_mat_t *link_mat;
    intra_link_t *link;
    FILE *fp;

    fp = fopen(f, "r");
    if (fp == NULL)
        return 0;

    link_mat = use_gap_seq? intra_link_mat_init(dict, resolution) : intra_link_mat_init_sdict(dict->sdict, resolution);

    pair_c = 0;
    intra_c = 0;

    while (1) {
        m = fread(&buffer, sizeof(uint32_t), BUFF_SIZE, fp);
        for (i = 0; i < m; i += 4) {
            if (use_gap_seq) {
                sd_coordinate_conversion(dict, buffer[i], buffer[i + 1], &i0, &p0);
                sd_coordinate_conversion(dict, buffer[i + 2], buffer[i + 3], &i1, &p1);
                b0 = (MAX(p0, 1) - 1) / resolution;
                b1 = (MAX(p1, 1) - 1) / resolution;
            } else {
                i0 = buffer[i];
                i1 = buffer[i + 2];
                b0 = (MAX(buffer[i + 1], 1) - 1) / resolution;
                b1 = (MAX(buffer[i + 3], 1) - 1) / resolution;
            }

            if (i0 == i1) {
                ++intra_c;
                link = &link_mat->links[i0];
                if (link->n) {
                    if (b0 > b1)
                        SWAP(uint32_t, b0, b1);
                    link->link[(link->n * 2 - b1 + b0 - 3) * (b1 - b0) / 2 + b1] += 1.0;
                }
            }
        }
        pair_c += m / 4;

        if (m < BUFF_SIZE) {
            if (ferror(fp))
                return 0;
            break;
        }
    }
#ifdef DEBUG
    printf("[I::%s] %ld read pairs processed, %ld intra links \n", __func__, pair_c, intra_c);
#endif
    fclose(fp);

    // normalise links by cell size
    for (i = 0; i < link_mat->n; ++i) {
        link = &link_mat->links[i];
        n = link->n;
        if (n == 0)
            continue;
        n = n * (n + 1) / 2;
        for (j = 0; j < n; ++j)
            normalise_by_size(&link->link[j]);
    }

    return link_mat;
}

inter_link_mat_t *inter_link_mat_from_file(const char *f, asm_dict_t *dict, uint32_t resolution, uint32_t radius)
{
    uint32_t i, j, n, k, b0, b1;
    uint32_t buffer[BUFF_SIZE], m, i0, i1, p0, p1;
    double l0, l1, a, na[4], nc[4];
    long pair_c, inter_c, radius_c, noise_c;
    inter_link_mat_t *link_mat;
    inter_link_t *link;
    FILE *fp;

    fp = fopen(f, "r");
    if (fp == NULL)
        return 0;

    n = dict->n;
    link_mat = inter_link_mat_init(dict, resolution, radius);
    pair_c = inter_c = radius_c = 0;

    while (1) {
        m = fread(&buffer, sizeof(uint32_t), BUFF_SIZE, fp);
        for (i = 0; i < m; i += 4) {
            sd_coordinate_conversion(dict, buffer[i], buffer[i + 1], &i0, &p0);
            sd_coordinate_conversion(dict, buffer[i + 2], buffer[i + 3], &i1, &p1);

            if (i0 != i1) {
                ++inter_c;
                
                if (i0 > i1) {
                    SWAP(uint32_t, i0, i1);
                    SWAP(uint32_t, p0, p1);
                }

                link = &link_mat->links[(n * 2 - i0 - 3) * i0 / 2 + i1 - 1];

                if (link->n == 0)
                    continue;

                l0 = dict->s[i0].len / 2.0;
                l1 = dict->s[i1].len / 2.0;
                
                if (p0 >= l0 && p1 < l1) {
                    // i0(-) -> i1(+)
                    b0 = (uint32_t) ((2 * l0 - p0) / resolution);
                    b1 = (uint32_t) ((double) p1 / resolution);
                    if (b0 < link->b0 && b1 < link->b1 && b0 + b1 < radius) {
                        link->link[0][(MAX(1, b0) - 1) * link->b1 + b1] += 1.0;
                        ++radius_c;
                    }
                } else if(p0 >= l0 && p1 >= l1) {
                    // i0(-) -> i1(-)
                    b0 = (uint32_t) ((2 * l0 - p0) / resolution);
                    b1 = (uint32_t) ((2 * l1 - p1) / resolution);
                    if (b0 < link->b0 && b1 < link->b1 && b0 + b1 < radius) {
                        link->link[1][(MAX(1, b0) - 1) * link->b1 + b1] += 1.0;
                        ++radius_c;
                    }
                } else if(p0 < l0 && p1 < l1) {
                    // i0(+) -> i1(+)
                    b0 = (uint32_t) ((double) p0 / resolution);
                    b1 = (uint32_t) ((double) p1 / resolution);
                    if (b0 < link->b0 && b1 < link->b1 && b0 + b1 < radius) {
                        link->link[2][(MAX(1, b0) - 1) * link->b1 + b1] += 1.0;
                        ++radius_c;
                    }
                } else if(p0 < l0 && p1 >= l1) {
                    // i0(+) -> i1(-)
                    b0 = (uint32_t) ((double) p0 / resolution);
                    b1 = (uint32_t) ((2 * l1 - p1) / resolution);
                    if (b0 < link->b0 && b1 < link->b1 && b0 + b1 < radius) {
                        link->link[3][(MAX(1, b0) - 1) * link->b1 + b1] += 1.0;
                        ++radius_c;
                    }
                }
            }
        }
        pair_c += m / 4;

        if (m < BUFF_SIZE) {
            if (ferror(fp))
                return 0;
            break;
        }
    }

#ifdef DEBUG
    printf("[I::%s] %ld read pairs processed, %ld inter links \n", __func__, pair_c, inter_c);
    printf("[I::%s] within radius %d: %ld\n", __func__, radius, radius_c);
#endif
    fclose(fp);

    // normalise links by cell size
    for (i = 0; i < link_mat->n; ++i) {
        link = &link_mat->links[i];
        for (j = 0; j < link->n; ++j)
            for (k = 0; k < 4; ++k)
                normalise_by_size(&link->link[k][j]);
    }

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
                    na[k] += 1.0;
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
#ifdef DEBUG
    printf("[I::%s] noise links: %ld; area: %.3f; noise estimation: %.3f\n", __func__, noise_c, a, link_mat->noise);
#endif

    return link_mat;
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

int dcmp (const void * a, const void * b) {
   double cmp = *(double *) a - *(double *) b;
   return cmp > 0? 1 : ( cmp < 0? -1 : 0);
}

norm_t *calc_norms(intra_link_mat_t *link_mat)
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
    while (bs[r0] >= 30) 
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
    intra_c = 0.0;
    for (i = 0; i < n; ++i) {
        link = (double *) malloc(bs[i] * sizeof(double));
        t = 0;
        for (j = 0; j < link_mat->n; ++j) {
            b = link_mat->links[j].n;
            if (b > i + 1) {
                memcpy(link + t, link_mat->links[j].link + (b * 2 - i - 1) * i / 2 + i, (b - 1 - i) * sizeof(double));
                t += b - 1 - i;
            }
        }

        tmp_c = 0.0;
        for (j = 0; j < bs[i]; ++j)
            tmp_c += link[j];
        linkc[i] = tmp_c;
        intra_c += tmp_c;
        
        if (USE_MEDIAN_NORM) {
            qsort(link, bs[i], sizeof(double), dcmp);
            norms[i] = bs[i] & 1? link[bs[i] / 2] : (link[bs[i] / 2] + link[bs[i] / 2 - 1]) / 2;
        } else {
            norms[i] = MAX(linkc[i], 1.0) / bs[i];
        }

        free(link);
    }

    intra_c -= linkc[0];
    tmp_c = 0.0;
    for (r = 1; r < n; ++r) {
        tmp_c += linkc[r];
        if (tmp_c / intra_c >= 0.99)
            break;
    }
    
    r = MIN(MIN(r, r0), MAX_RADIUS);

    if (USE_MEDIAN_NORM) {
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
    }

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
                if (l < 0.0) 
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
        printf("[I::%s] noise bin %d %d\n", __func__, i, hist[i]);
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

void inter_link_norms(inter_link_mat_t *link_mat, norm_t *norm, int use_estimated_noise, double *la)
{
    uint32_t i, j, k, b, b1, r, n, n0;
    uint32_t *fn;
    int8_t bl;
    double l, *norms, *fr, noise;
    double c, c0, t;
    inter_link_t *inter_link;

    noise = 0.0;
    if (use_estimated_noise) {
        // noise = estimate_noise(link_mat);
        noise = link_mat->noise;
        fprintf(stderr, "[I::%s] using noise level %.3f\n", __func__, noise);
    }

    r = link_mat->r;
    norms = (double *) malloc((r + 1) * sizeof(double));
    for (i = 0; i <= r; ++i)
        norms[i] = norm->norms[i] - noise;

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
                        if (l < 0.0)
                            break;
                        if (l >= norms[b + 1] + noise) {
                            inter_link->norms[k] += 1.0;
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
    
    double c1 = 0;
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
                        if(l < 0.0)
                            break;
                        l = MAX(0.0, l - noise);
                        fr[k * r + b] += MIN(1.0, l / norms[b + 1]);
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
        c1 += t;
loop_i_end:
        free(fr);
        free(fn);
    }
    
    *la = c1 / c;
    fprintf(stderr, "[I::%s] average link count: %.3f %.3f %.3f\n", __func__, c1, c, *la);
    
    free(norms);
}

int8_t *calc_link_directs_from_file(const char *f, asm_dict_t *dict)
{
    uint32_t i, j, n, na, k, b0, b1, b, ma, sma, n_ma;
    uint32_t buffer[BUFF_SIZE], m, i0, i1, p0, p1;
    uint32_t *link, l;
    int8_t *directs;
    long pair_c, inter_c;
    FILE *fp;

    fp = fopen(f, "r");
    if (fp == NULL)
        return 0;

    n = dict->n;
    na = n * (n - 1) / 2;
    link = (uint32_t *) calloc(na << 2, sizeof(uint32_t));
    pair_c = inter_c = 0;

    while (1) {
        m = fread(&buffer, sizeof(uint32_t), BUFF_SIZE, fp);
        for (i = 0; i < m; i += 4) {
            sd_coordinate_conversion(dict, buffer[i], buffer[i + 1], &i0, &p0);
            sd_coordinate_conversion(dict, buffer[i + 2], buffer[i + 3], &i1, &p1);

            if (i0 != i1) {
                ++inter_c;

                if (i0 > i1) {
                    SWAP(uint32_t, i0, i1);
                    SWAP(uint32_t, p0, p1);
                }
                b0 = !((double) p0 / dict->s[i0].len > 0.5);
                b1 = !((double) p1 / dict->s[i1].len > 0.5);
                b = (b0 << 1) | (!b1);
                ++link[b * na + (n * 2 - i0 - 1) * i0 / 2 + i1 - i0 - 1];
            }
        }
        pair_c += m / 4;
        if (m < BUFF_SIZE) {
            if (ferror(fp))
                return 0;
            break;
        }
    }

#ifdef DEBUG
    printf("[I::%s] %ld read pairs processed, %ld inter links \n", __func__, pair_c, inter_c);
#endif
    fclose(fp);
    
    i0 = asm_sd_get(dict, "scaffold_11");
    i1 = asm_sd_get(dict, "scaffold_111");
    if (i0 > i1)
        SWAP(uint32_t, i0, i1);
    j = (n * 2 - i0 - 1) * i0 / 2 + i1 - i0 - 1;
    printf("XXXXX: %u %u %u %u %u %u\n", dict->s[i0].len, dict->s[i1].len, link[j], link[na + j], link[na * 2 + j], link[na * 3 + j]);

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
        
        directs[i] = n_ma == 1 && ma * 0.9 > sma? j : -1;
    }
    free(link);
#ifdef DEBUG_ORIEN
    int8_t d;
    for (i0 = 0; i0 < n; ++i0) {
        for (i1 = i0 + 1; i1 < n; ++i1) {
            d = directs[(n * 2 - i0 - 1) * i0 / 2 + i1 - i0 - 1];
            if (d >=0)
                printf("[I::%s] #Oriens: %s %s %hhi \n", __func__, dict->s[i0].name, dict->s[i1].name, d);
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
                    if (inter_link->link[k][j] > 0.0)
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
        if (n_ma == 1 && ma * 0.9 > sma) {
            inter_link->linkt = 1 << t;
        } else {
            t = 0;
            for (k = 0; k < 4; ++k)
                if (inter_link->norms[k] >= min_norm)
                    t |= (1 << k);
            inter_link->linkt = t;
        }

#ifdef DEBUG_ORIEN
        printf("[I::%s] %d %.0f %.0f %d %d [%.3f %.3f %.3f %.3f] [%.0f  %.0f  %.0f  %.0f] %s %s\n", __func__, inter_link->linkt, ma, sma, inter_link->b0, inter_link->b1, inter_link->norms[0], inter_link->norms[1], inter_link->norms[2], inter_link->norms[3], area[0], area[1], area[2], area[3], dict->s[inter_link->c0].name, dict->s[inter_link->c1].name);
#endif
    }
}

static int32_t get_target_end(uint32_t *cigar, int n_cigar, int32_t s)
{
    int i;
    uint8_t c;
    for (i = 0; i < n_cigar; ++i) {
        c = cigar[i] & BAM_CIGAR_MASK;
        if (c == BAM_CMATCH || c == BAM_CDEL)
            s += cigar[i] >> BAM_CIGAR_SHIFT;
    }
    return s;
}

static char *parse_bam_rec(bam1_t *b, bam_header_t *h, uint8_t q, int32_t *s, int32_t *e, char **cname)
{
    // 0x4 0x100 0x400 0x800
    if (b->core.flag & 0xD04)
        return 0;
    *cname = h->target_name[b->core.tid];
    if (b->core.qual < q) {
        *s = -1;
        *e = -1;
    } else {
        *s = b->core.pos + 1;
        *e = get_target_end(bam1_cigar(b), b->core.n_cigar, b->core.pos) + 1;
    }
    
    return strdup(bam1_qname(b));
}

void dump_links_from_bam_file(const char *f, const char *fai, uint8_t mq, const char *out)
{
    bamFile fp;
    FILE *fo;
    bam_header_t *h;
    bam1_t *b;
    char *cname0, *cname1, *rname0, *rname1;
    int32_t s0, s1, e0, e1;
    uint32_t i0, i1, p0, p1;
    int8_t buff;
    long rec_c, pair_c, inter_c, intra_c;
    sdict_t *dict = make_sdict_from_index(fai);

    fp = bam_open(f, "r"); // sorted by read name
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open fail %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }
    
    fo = fopen(out, "w");
    if (fo == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }

    h = bam_header_read(fp);
    b = bam_init1();
    cname0 = cname1 = rname0 = rname1 = 0;
    s0 = s1 = e0 = e1 = 0;
    i0 = i1 = p0 = p1 = 0;
    rec_c = pair_c = inter_c = intra_c = 0;
    buff = 0;
    while (bam_read1(fp, b) >= 0 ) {
        if (buff == 0) {
            rname0 = parse_bam_rec(b, h, mq, &s0, &e0, &cname0);
            if (!rname0)
                continue;
            ++buff;
        } else if (buff == 1) {
            rname1 = parse_bam_rec(b, h, mq, &s1, &e1, &cname1);
            if (!rname1)
                continue;
            if (strcmp(rname0, rname1) == 0) {
                ++pair_c;

                if (s0 > 0 && s1 > 0) {
                    i0 = sd_get(dict, cname0);
                    i1 = sd_get(dict, cname1);

                    if (i0 == UINT32_MAX || i1 == UINT32_MAX) {
                        fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, i0 < 0? cname0 : cname1);
                    } else {
                        if (i0 == i1)
                            ++intra_c;
                        else
                            ++inter_c;
                        p0 = (s0 + e0) / 2;
                        p1 = (s1 + e1) / 2;
                        if (i0 > i1) {
                            SWAP(uint32_t, i0, i1);
                            SWAP(uint32_t, p0, p1);
                        }
                        fwrite(&i0, sizeof(uint32_t), 1, fo);
                        fwrite(&p0, sizeof(uint32_t), 1, fo);
                        fwrite(&i1, sizeof(uint32_t), 1, fo);
                        fwrite(&p1, sizeof(uint32_t), 1, fo);
                    }
                }
                free(rname0);
                free(rname1);
                rname0 = 0;
                rname1 = 0;
                buff = 0;
            } else {
                cname0 = cname1;
                s0 = s1;
                e0 = e1;
                free(rname0);
                rname0 = rname1;
                rname1 = 0;
                buff = 1;
            }
        }

        if (++rec_c % 1000000 == 0)
            fprintf(stderr, "[I::%s] %ld million records processed, %ld read pairs \n", __func__, rec_c / 1000000, pair_c);
    }
    if (rname0)
        free(rname0);
    if (rname1)
        free(rname1);
    bam_destroy1(b);
    bam_header_destroy(h);
    sd_destroy(dict);
    bam_close(fp);
    fclose(fo);

    fprintf(stderr, "[I::%s] dumped %ld read pairs: %ld intra links + %ld inter links \n", __func__, pair_c, intra_c, inter_c);
}

void dump_links_from_bed_file(const char *f, const char *fai, uint8_t mq, const char *out)
{
    FILE *fp, *fo;
    char *line = NULL;
    size_t ln = 0;
    ssize_t read;
    char cname0[4096], cname1[4096], rname0[4096], rname1[4096];
    uint32_t s0, s1, e0, e1, i0, i1, p0, p1;
    int8_t buff;
    long rec_c, pair_c, inter_c, intra_c;
    sdict_t *dict = make_sdict_from_index(fai);

    fp = fopen(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }
    fo = fopen(out, "w");
    if (fo == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }
    
    s0 = s1 = e0 = e1 = 0;
    i0 = i1 = p0 = p1 = 0;
    rec_c = pair_c = inter_c = intra_c = 0;
    buff = 0;
    while ((read = getline(&line, &ln, fp)) != -1) {
        if (buff == 0) {
            sscanf(line, "%s %u %u %s", cname0, &s0, &e0, rname0);
            ++buff;
        } else if (buff == 1) {
            sscanf(line, "%s %u %u %s", cname1, &s1, &e1, rname1);
            if (is_read_pair(rname0, rname1)) {
                ++pair_c;

                i0 = sd_get(dict, cname0);
                i1 = sd_get(dict, cname1);

                if (i0 == UINT32_MAX || i1 == UINT32_MAX) {
                    fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, i0 < 0? cname0 : cname1);
                } else {
                    if (i0 == i1)
                        ++intra_c;
                    else
                        ++inter_c;
                    p0 = (s0 + e0) / 2;
                    p1 = (s1 + e1) / 2;
                    if (i0 > i1) {
                        SWAP(uint32_t, i0, i1);
                        SWAP(uint32_t, p0, p1);
                    }
                    fwrite(&i0, sizeof(uint32_t), 1, fo);
                    fwrite(&p0, sizeof(uint32_t), 1, fo);
                    fwrite(&i1, sizeof(uint32_t), 1, fo);
                    fwrite(&p1, sizeof(uint32_t), 1, fo);
                }
                buff = 0;
            } else {
                strcpy(cname0, cname1);
                s0 = s1;
                e0 = e1;
                strcpy(rname0, rname1);
                buff = 1;
            }
        }

        if (++rec_c % 1000000 == 0)
            fprintf(stderr, "[I::%s] %ld million records processed, %ld read pairs \n", __func__, rec_c / 1000000, pair_c);
    }
    
    if (line)
        free(line);
    sd_destroy(dict);
    fclose(fp);
    fclose(fo);

    fprintf(stderr, "[I::%s] dumped %ld read pairs: %ld intra links + %ld inter links \n", __func__, pair_c, intra_c, inter_c);
}

