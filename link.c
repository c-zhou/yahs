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
static int MAX_RADIUS = 100;

void intra_link_mat_destroy(intra_link_mat_t *link_mat)
{
    int i;
    for (i = 0; i < link_mat->n; ++i)
        if (link_mat->links[i].n)
            free(link_mat->links[i].link);
    if (link_mat->links)
        free(link_mat->links);
    free(link_mat);
}

void inter_link_mat_destroy(inter_link_mat_t *link_mat)
{
    int i, j;
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
inter_link_mat_t *inter_link_mat_init(asm_dict_t *dict, int resolution, int radius)
{
    inter_link_mat_t *link_mat;
    inter_link_t *link;
    int32_t i, j, k, l, b0, b1, n, m, p;
    double a0, a1, a;

    n = dict->n;
    m = n * (n - 1) / 2;
    link_mat = (inter_link_mat_t *) malloc(sizeof(inter_link_mat_t));
    link_mat->n = m;
    link_mat->r = radius;
    link_mat->links = (inter_link_t *) malloc(m * sizeof(inter_link_t));

    for (i = 0; i < n; ++i) {
        if (dict->s[i].len < 2 * resolution) {
            for (j = i + 1; j < n; ++j)
                link_mat->links[(n * 2 - i - 3) * i / 2 + j - 1].n = 0;
            continue;
        }
        
        b0 = MIN(radius, ceil((double) dict->s[i].len / 2.0 / resolution));
        // relative edge length of the last cell
        a0 = MIN(1.0, ((double) dict->s[i].len / 2.0 - (b0 - 1) * resolution) / resolution);
        for (j = i + 1; j < n; ++j) {
            link = &link_mat->links[(n * 2 - i - 3) * i / 2 + j - 1];
            if (dict->s[j].len < 2 * resolution) {
                link->n = 0;
                continue;
            }
            
            b1 = MIN(radius, ceil((double) dict->s[j].len / 2.0 / resolution));
            // relative edge length of the last cell
            a1 = MIN(1.0, ((double) dict->s[j].len / 2.0 - (b1 - 1) * resolution) / resolution);
            p = b0 * b1; 
            link->c0 = i;
            link->c1 = j;
            link->b0 = b0;
            link->b1 = b1;
            link->r = MIN(radius, b0 + b1 - 1);
            link->n = p;
            link->n0 = -1;
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

intra_link_mat_t *intra_link_mat_init(asm_dict_t *dict, int resolution)
{
    intra_link_mat_t *link_mat;
    intra_link_t *link;
    int32_t i, j, n, b, p;
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
        b = ceil((double) dict->s[i].len / resolution);
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

intra_link_mat_t *intra_link_mat_init_sdict(sdict_t *dict, int resolution)
{
    intra_link_mat_t *link_mat;
    intra_link_t *link;
    int32_t i, j, n, b, p;
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
        b = ceil((double) dict->s[i].len / resolution);
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
intra_link_mat_t *intra_link_mat_from_file(const char *f, asm_dict_t *dict, int resolution, int use_gap_seq)
{
    int32_t i, j, n, i0, i1, b0, b1, p0, p1;
    int32_t buffer[BUFF_SIZE], m;
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
        m = fread(&buffer, sizeof(int32_t), BUFF_SIZE, fp);
        for (i = 0; i < m; i += 4) {
            if (use_gap_seq) {
                sd_coordinate_conversion(dict, buffer[i], buffer[i + 1], &i0, &p0);
                sd_coordinate_conversion(dict, buffer[i + 2], buffer[i + 3], &i1, &p1);
                b0 = p0 / resolution;
                b1 = p1 / resolution;
            } else {
                i0 = buffer[i];
                i1 = buffer[i + 2];
                b0 = buffer[i + 1] / resolution;
                b1 = buffer[i + 3] / resolution;
            }

            if (i0 == i1) {
                ++intra_c;
                link = &link_mat->links[i0];
                n = link->n - 1;
                if (n < 0)
                    continue;
                if (b0 > b1)
                    SWAP(int, b0, b1);
                link->link[(n * 2 - b1 + b0 -1) * (b1 - b0) / 2 + b1] += 1.0;
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

inter_link_mat_t *inter_link_mat_from_file(const char *f, asm_dict_t *dict, int resolution, int radius)
{
    int32_t i, j, n, k, i0, i1, b0, b1, p0, p1;
    int32_t buffer[BUFF_SIZE], m;
    double l0, l1, a;
    long pair_c, inter_c, radius_c, noise_c;
    inter_link_mat_t *link_mat;
    inter_link_t *link;
    FILE *fp;

    fp = fopen(f, "r");
    if (fp == NULL)
        return 0;

    n = dict->n;
    link_mat = inter_link_mat_init(dict, resolution, radius);
    pair_c = inter_c = radius_c = noise_c = 0;

    while (1) {
        m = fread(&buffer, sizeof(int32_t), BUFF_SIZE, fp);
        for (i = 0; i < m; i += 4) {
            sd_coordinate_conversion(dict, buffer[i], buffer[i + 1], &i0, &p0);
            sd_coordinate_conversion(dict, buffer[i + 2], buffer[i + 3], &i1, &p1);

            if (i0 != i1) {
                ++inter_c;
                
                if (i0 > i1) {
                    SWAP(int32_t, i0, i1);
                    SWAP(int32_t, p0, p1);
                }

                link = &link_mat->links[(n * 2 - i0 - 3) * i0 / 2 + i1 - 1];

                if (link->n == 0)
                    continue;

                ++noise_c;

                l0 = dict->s[i0].len / 2.0;
                l1 = dict->s[i1].len / 2.0;
                
                if (p0 >= l0 && p1 < l1) {
                    // i0(-) -> i1(+)
                    b0 = (int) ((2 * l0 - p0) / resolution);
                    b1 = (int) ((double) p1 / resolution);
                    if (b0 < link->b0 && b1 < link->b1 && b0 + b1 < radius) {
                        link->link[0][MAX(0, b0 - 1) * link->b1 + b1] += 1.0;
                        ++radius_c;
                    }
                } else if(p0 >= l0 && p1 >= l1) {
                    // i0(-) -> i1(-)
                    b0 = (int) ((2 * l0 - p0) / resolution);
                    b1 = (int) ((2 * l1 - p1) / resolution);
                    if (b0 < link->b0 && b1 < link->b1 && b0 + b1 < radius) {
                        link->link[1][MAX(0, b0 - 1) * link->b1 + b1] += 1.0;
                        ++radius_c;
                    }
                } else if(p0 < l0 && p1 < l1) {
                    // i0(+) -> i1(+)
                    b0 = (int) ((double) p0 / resolution);
                    b1 = (int) ((double) p1 / resolution);
                    if (b0 < link->b0 && b1 < link->b1 && b0 + b1 < radius) {
                        link->link[2][MAX(0, b0 - 1) * link->b1 + b1] += 1.0;
                        ++radius_c;
                    }
                } else if(p0 < l0 && p1 >= l1) {
                    // i0(+) -> i1(-)
                    b0 = (int) ((double) p0 / resolution);
                    b1 = (int) ((2 * l1 - p1) / resolution);
                    if (b0 < link->b0 && b1 < link->b1 && b0 + b1 < radius) {
                        link->link[3][MAX(0, b0 - 1) * link->b1 + b1] += 1.0;
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
    a = 0.0;
    for (i = 0; i < link_mat->n; ++i) {
        link = &link_mat->links[i];
        if (link->n == 0)
            continue;
        l0 = (double) dict->s[link->c0].len / resolution;
        l1 = (double) dict->s[link->c1].len / resolution;
        a += l0 * l1;
    }
    link_mat->noise = a > 0? noise_c / a : 0;
#ifdef DEBUG
    printf("[I::%s] noise links: %ld; area: %.3f; noise estimation: %.3f\n", __func__, noise_c, a, link_mat->noise);
#endif

    return link_mat;
}

void norm_destroy(norm_t *norm)
{
    int i;
    free(norm->bs);
    for (i = 0; i < norm->n; ++i)
        free(norm->link[i]);
    free(norm->link);
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
    int32_t i, j, n, b, r, r0;
    int32_t *bs, *t;
    double **link;
    double intra_c, tmp_c;
    double *norms, *intra_links, *linkc;
    norm_t *norm;
    // find maximum numebr of bands
    n = 0;
    for (i = 0; i < link_mat->n; ++i)
        n = MAX(link_mat->links[i].n, n);
    n -= 1; // do not use the last one as it might be an incomplete cell
    if (n <= 0) {
        fprintf(stderr, "[E::%s] no bands (%d) for norm calculation, try a higher resolution\n", __func__, n);
        return 0;
    }
    // calculate number of cells for each band
    bs = (int32_t *) calloc(n, sizeof(int32_t));
    for (i = 0; i < link_mat->n; ++i) {
        b = link_mat->links[i].n - 1;
        for (j = 0; j < b; ++j)
            bs[j] += b - j;
    }

    r0 = 0;
    while (bs[r0] >= 30) 
        ++r0;
    if (r0 < 2) {
        fprintf(stderr, "[E::%s] no enough bands (%d) for norm calculation, try a higher resolution\n", __func__, r0);
        free(bs);
        return 0;
    }

    // initialise cells
    link = (double **) malloc(n * sizeof(double *));
    for (i = 0; i < n; ++i)
        link[i] = (double *) malloc(bs[i] * sizeof(double));
    // fill cells
    t = (int32_t *) calloc(n, sizeof(int32_t));
    for (i = 0; i < link_mat->n; ++i) {
        b = link_mat->links[i].n;
        intra_links = link_mat->links[i].link;
        for (j = 0; j < b - 1; ++j) {
            memcpy(link[j] + t[j], intra_links + (b * 2 - j - 1) * j / 2 + j, (b - 1 - j) * sizeof(double));
            t[j] += b - 1 - j;
        }
    }
    free(t);
    // caluclate links in each band and radius
    linkc = (double *) malloc(n * sizeof(double));
    intra_c = 0.0;
    for (i = 0; i < n; ++i) {
        intra_links = link[i];
        tmp_c = 0.0;
        for (j = 0; j < bs[i]; ++j)
            tmp_c += intra_links[j];
        linkc[i] = tmp_c;
        intra_c += tmp_c;
    }
    intra_c -= linkc[0];
    tmp_c = 0.0;
    for (r = 1; r < n; ++r) {
        tmp_c += linkc[r];
        if (tmp_c / intra_c >= 0.95)
            break;
    }
    
    r = MIN(MIN(r, r0), MAX_RADIUS);
    // calculate norms
    // using median or mean?
    norms = (double *) malloc(n * sizeof(double));
    for (i = 0; i < n; ++i) {
        if (USE_MEDIAN_NORM) {
            qsort(link[i], bs[i], sizeof(double), dcmp);
            norms[i] = bs[i] & 1? link[i][bs[i] / 2] : (link[i][bs[i] / 2] + link[i][bs[i] / 2 - 1]) / 2;
        } else {
            norms[i] = MAX(linkc[i], 1.0) / bs[i];
        }
    }

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
    norm->link = link;
    norm->linkc = linkc;
    norm->norms = norms;
    return norm;
}

void print_inter_link_norms(FILE *fp, inter_link_mat_t *link_mat, asm_dict_t *dict)
{
    int i, n;
    inter_link_t *link;
    char *cname0, *cname1;
    n = link_mat->n;
    for (i = 0; i < n; ++i) {
        link = &link_mat->links[i];
        if (link->n == 0)
            continue;
        cname0 = dict->s[link->c0].name;
        cname1 = dict->s[link->c1].name;
        fprintf(fp, "NORM: %s %s (%d) [%.3f %.3f %.3f %.3f]\n", cname0, cname1, link->n0, link->norms[0], link->norms[1], link->norms[2], link->norms[3]);
    }
}

void print_inter_link_bands(FILE *fp, inter_link_mat_t *link_mat, asm_dict_t *dict)
{
    int i, j, k, n;
    inter_link_t *link;
    char *cname0, *cname1;
    n = link_mat->n;
    for (i = 0; i < n; ++i) {
        link = &link_mat->links[i];
        if (link->n == 0)
            continue;
        cname0 = dict->s[link->c0].name;
        cname1 = dict->s[link->c1].name;
        for (j = 0; j < 4; ++j) {
            fprintf(fp, "BAND: %s%s %s%s (%d, %d) [%.3f]", cname0, j < 2? "-" : "+", cname1, j & 1? "-" : "+", link->n0, link->r, link->norms[j]);
            for (k = 0; k < link->r; ++k)
                fprintf(fp, " %.0f", link->linkb[j][k]);
            fprintf(fp, "\n");
        }
    }
}

double *get_max_inter_norms(inter_link_mat_t *link_mat, asm_dict_t *dict)
{
    double *max_norms;
    inter_link_t *link;
    int i, j, n;
    int32_t c0, c1;
    n = link_mat->n;
    max_norms = (double *) calloc(dict->n, sizeof(double));
    for (i = 0; i < n; ++i) {
        link = &link_mat->links[i];
        if (link->n == 0)
            continue;
        c0 = link->c0;
        c1 = link->c1;
        for (j = 0; j < 4; ++j) {
            max_norms[c0] = MAX(max_norms[c0], link->norms[j]);
            max_norms[c1] = MAX(max_norms[c1], link->norms[j]);
        }
    }
    return max_norms;
}

void print_norms(FILE *fp, norm_t *norm)
{
    int i, j, n;
    int64_t link_count;
    double *link;
    link_count = 0;

    fprintf(fp, "norms:");
    for (i = 0; i < norm->n; ++i)
        fprintf(fp, " %.1f", norm->norms[i]);
    fprintf(fp, "\n");
    
    fprintf(fp, "links:");
    for (i = 0; i < norm->n; ++i)
        fprintf(fp, " %.0f", norm->linkc[i]);
    fprintf(fp, "\n");

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
}

static int estimate_noise(inter_link_mat_t *link_mat)
{
    int32_t i, j, k, n, l;
    uint32_t bin, *hist;
    uint64_t total;
    inter_link_t *inter_link;

    bin = 4096;
    hist = calloc(bin, sizeof(uint32_t));
    for (i = 0; i < link_mat->n; ++i) {
        inter_link = &link_mat->links[i];
        n = inter_link->n;
        for (j = 0; j < n; ++j) {
            for (k = 0; k < 4; ++k) {
                l = inter_link->link[k][j];
                if (l < 0.0) 
                    continue;
                l = (uint64_t) l;
                if (l >= bin)
                    l = bin - 1;
                ++hist[l];
            }
        }
    }
   
    hist[0] = 0;

    total = 0;
    for (i = 0; i < bin; ++i)
        total += hist[i];
    total >>= 1;

#ifdef DEBUG_NOISE
    for (i = 0; i < bin; ++i) {
        printf("%d %d\n", i, hist[i]);
    }
#endif

    n = 0;
    while (hist[n] < total) {
        ++n;
        if (n == bin)
            break;
        hist[n] += hist[n-1];
    }

    if (n == 4096)
        fprintf(stderr, "[W::%s] noise esimation might be incorrect %d\n", __func__, n);
    
    free(hist);

    return n;
}

void inter_link_norms(inter_link_mat_t *link_mat, norm_t *norm, int use_estimated_noise)
{
    int32_t i, j, k, b, b1, r, n, n0;
    double l, *norms, *fr, noise;
    int32_t *fn;
    inter_link_t *inter_link;

    noise = 0.0;
    if (use_estimated_noise)
        // noise = estimate_noise(link_mat);
        noise = link_mat->noise;

    r = link_mat->r;
    norms = (double *) malloc((r + 1) * sizeof(double));
    for (i = 0; i <= r; ++i)
        norms[i] = norm->norms[i] - noise;
    
    fr = (double *) malloc(4 * r * sizeof(double));
    fn = (int32_t *) malloc(r * sizeof(int32_t));

    for (i = 0; i < link_mat->n; ++i) {
        inter_link = &link_mat->links[i];
        n = inter_link->n;
        if (n == 0)
            continue;
        b1 = inter_link->b1;
        n0 = 0;
        memset(fr, 0, 4 * r * sizeof(double));
        memset(fn, 0, r * sizeof(int32_t));

        for (j = 0; j < n; ++j) {
            b = j / b1 + j % b1;
            if (b < r && norms[b + 1] > 0) {
                for (k = 0; k < 4; ++k) {
                    l = inter_link->link[k][j];
                    if(l < 0.0) 
                        break;
                    l = MAX(0.0, l - noise);
                    fr[k * r + b] += MIN(1.0, l / norms[b + 1]);
                }
                if (l >= 0.0)
                    ++fn[b];
            }
        }

        // cumsum
        for (b = 1; b < r; ++b) {
            fn[b] += fn[b - 1];
            for (k = 0; k < 4; ++k)
                fr[k * r + b] += fr[k * r + b - 1];
        }
        // calculate filling rate
        for (b = 1; b < r; ++b)
            if (fn[b])
                for (k = 0; k < 4; ++k)
                    fr[k * r + b] /= fn[b];
        for (k = 0; k < 4; ++k)
            fr[k * r] = 1.0;

        for (j = 0; j < n; ++j) {
            b = j / b1 + j % b1;
            if (b < r && norms[b + 1] > 0) {
                for (k = 0; k < 4; ++k) {
                    l = inter_link->link[k][j];
                    if(l < 0.0) 
                        break;
                    l = MAX(0.0, l - noise);
                    inter_link->norms[k] += fr[k * r + b] * MIN(1.0, l / norms[b + 1]);
                }
                if (l >= 0.0)
                    ++n0;
            }
        }

        for (k = 0; k < 4; ++k)
            if (n0)
                inter_link->norms[k] /= n0;
        inter_link->n0 = n0;
    }

    free(fr);
    free(fn);
    free(norms);
}

void calc_link_directs(inter_link_mat_t *link_mat, double min_norm, asm_dict_t *sdict)
{
    int32_t i, j, k, n, b, b1, r, t, n_ma;
    inter_link_t *inter_link;
    double area[4]; // area under the cumsum of linkb
    double ma, sma;

    r = link_mat->r;
    for (i = 0; i < link_mat->n; ++i) {
        inter_link = &link_mat->links[i];
        n = inter_link->n;
        if (n == 0)
            continue;
        b1 = inter_link->b1;
        for (j = 0; j < n; ++j) {
            b = j / b1 + j % b1;
            if (b < r)
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
                area[k] += inter_link->linkb[k][j] * (r - j);
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
        if (n_ma == 1 && ma * 0.8 > sma) {
            inter_link->linkt = 1 << t;
        } else {
            t = 0;
            for (k = 0; k < 4; ++k)
                if (inter_link->norms[k] >= min_norm)
                    t |= (1 << k);
            inter_link->linkt = t;
        }

#ifdef DEBUG_ORIEN
        printf("[I::%s] %d %.0f %.0f %d %d [%.3f %.3f %.3f %.3f] [%.0f  %.0f  %.0f  %.0f] %s %s\n", __func__, inter_link->linkt, ma, sma, inter_link->b0, inter_link->b1, inter_link->norms[0], inter_link->norms[1], inter_link->norms[2], inter_link->norms[3], area[0], area[1], area[2], area[3], sdict->s[inter_link->c0].name, sdict->s[inter_link->c1].name);
#endif
    }
}

static uint32_t get_target_end(uint32_t *cigar, int n_cigar, uint32_t s)
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

static char *parse_bam_rec(bam1_t *b, bam_header_t *h, int32_t *s, int32_t *e, char **cname)
{
    *cname = h->target_name[b->core.tid];
    if (b->core.flag & 0x4 || b->core.flag & 0x400) {
        *s = -1;
        *e = -1;
    } else {
        *s = b->core.pos + 1;
        *e = get_target_end(bam1_cigar(b), b->core.n_cigar, b->core.pos) + 1;
    }
    
    return strdup(bam1_qname(b));
}

void dump_links_from_bam_file(const char *f, const char *fai, const char *out)
{
    bamFile fp;
    FILE *fo;
    bam_header_t *h;
    bam1_t *b;
    char *cname0, *cname1, *rname0, *rname1;
    int32_t s0, s1, e0, e1, i0, i1, p0, p1;
    int8_t buff;
    long pair_c, inter_c, intra_c;
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

    pair_c = inter_c = intra_c = buff = 0;
    while (bam_read1(fp, b) >= 0 ) {
        if (buff == 0) {
            rname0 = parse_bam_rec(b, h, &s0, &e0, &cname0);
            ++buff;
        } else if (buff == 1) {
            rname1 = parse_bam_rec(b, h, &s1, &e1, &cname1);
            if (strcmp(rname0, rname1) == 0) {
                if (++pair_c % 1000000 == 0)
                    fprintf(stderr, "[I::%s] %ld million read pairs processed \n", __func__, pair_c / 1000000);
                
                if (s0 > 0 && s1 > 0) {
                    i0 = sd_get(dict, cname0);
                    i1 = sd_get(dict, cname1);

                    if (i0 < 0 || i1 < 0) {
                        fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, i0 < 0? cname0 : cname1);
                    } else {
                        if (i0 == i1)
                            ++intra_c;
                        else
                            ++inter_c;
                        p0 = (s0 + e0) / 2;
                        p1 = (s1 + e1) / 2;
                        if (i0 > i1) {
                            SWAP(int32_t, i0, i1);
                            SWAP(int32_t, p0, p1);
                        }
                        fwrite(&i0, sizeof(int32_t), 1, fo);
                        fwrite(&p0, sizeof(int32_t), 1, fo);
                        fwrite(&i1, sizeof(int32_t), 1, fo);
                        fwrite(&p1, sizeof(int32_t), 1, fo);
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
                buff = 1;
            }
        }
    }

    if (rname0)
        free(rname0);
    if (rname1)
        free(rname1);
    bam_destroy1(b);
    bam_header_destroy(h);
    bam_close(fp);
    fclose(fo);

    fprintf(stderr, "[I::%s] dumped %ld read pairs: %ld intra links + %ld inter links \n", __func__, pair_c, intra_c, inter_c);
}

void dump_links_from_bed_file(const char *f, const char *fai, const char *out)
{
    FILE *fp, *fo;
    char *line = NULL;
    size_t ln = 0;
    ssize_t read;
    char cname0[4096], cname1[4096], rname0[4096], rname1[4096];
    int32_t s0, s1, e0, e1, i0, i1, p0, p1;
    int8_t buff;
    long pair_c, inter_c, intra_c;
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

    pair_c = inter_c = intra_c = buff = 0;
    while ((read = getline(&line, &ln, fp)) != -1) {
        if (buff == 0) {
            sscanf(line, "%s %d %d %s", cname0, &s0, &e0, rname0);
            ++buff;
        } else if (buff == 1) {
            sscanf(line, "%s %d %d %s", cname1, &s1, &e1, rname1);
            if (is_read_pair(rname0, rname1)) {
                if (++pair_c % 1000000 == 0)
                    fprintf(stderr, "[I::%s] %ld million read pairs processed \n", __func__, pair_c / 1000000);
                i0 = sd_get(dict, cname0);
                i1 = sd_get(dict, cname1);

                if (i0 < 0 || i1 < 0) {
                    fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, i0 < 0? cname0 : cname1);
                } else {
                    if (i0 == i1)
                        ++intra_c;
                    else
                        ++inter_c;
                    p0 = (s0 + e0) / 2;
                    p1 = (s1 + e1) / 2;
                    if (i0 > i1) {
                        SWAP(int32_t, i0, i1);
                        SWAP(int32_t, p0, p1);
                    }
                    fwrite(&i0, sizeof(int32_t), 1, fo);
                    fwrite(&p0, sizeof(int32_t), 1, fo);
                    fwrite(&i1, sizeof(int32_t), 1, fo);
                    fwrite(&p1, sizeof(int32_t), 1, fo);
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
    }
    
    if (line)
        free(line);
    fclose(fp);
    fclose(fo);

    fprintf(stderr, "[I::%s] dumped %ld read pairs: %ld intra links + %ld inter links \n", __func__, pair_c, intra_c, inter_c);
}

