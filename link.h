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
#ifndef LINK_H_
#define LINK_H_

#include <stdlib.h>
#include <stdint.h>

#include "sdict.h"
#include "enzyme.h"
#include "cov.h"

#define SQRT2 1.41421356237
#define SQRT2_2 .70710678118

typedef struct {
    uint32_t c;
    uint32_t n;
    double *link;
} intra_link_t;

typedef struct {
    uint32_t c0, c1; // sequence id
    uint32_t b0, b1; // number of bands
    uint32_t r; // real radius MIN(R, b0 + b1 - 1)
    uint32_t n, n0; // number of cells, n = b0 * b1, n0 is the real number of cells
    int8_t linkt; // 0...3 determine the join direction
    // link[0]: i0(-) -> i1(+)
    // link[1]: i0(-) -> i1(-)
    // link[2]: i0(+) -> i1(+)
    // link[3]: i0(+) -> i1(-)
    double *link[4]; // size 4 x n
    double *linkb[4]; // total number of links in each band, of size 4 x r
    double norms[4];
} inter_link_t;

typedef struct {
    uint32_t n;
    intra_link_t *links;
} intra_link_mat_t;

typedef struct {
    uint32_t n;
    uint32_t r; // radius
    double noise; // noise
    inter_link_t *links;
} inter_link_mat_t;

typedef struct {
    uint32_t n; // number of bands
    uint32_t *bs; // number of cells in each band [1 x n]
    // double **link; // number of links in each cell [n x b]
    double *norms; // norms [1 x n]
    double *linkc; // link count in each band [1 x n]
    uint32_t r; // number of first r bands contains at least 90% links adjusted by norms
} norm_t;

#ifdef __cplusplus 
extern "C" {
#endif

intra_link_mat_t *intra_link_mat_init(void *dict, uint32_t resolution, int use_gap_seq);
inter_link_mat_t *inter_link_mat_init(asm_dict_t *dict, uint32_t resolution, uint32_t radius);
intra_link_mat_t *intra_link_mat_from_file(const char *f, cov_norm_t *cov_norm, asm_dict_t *dict, re_cuts_t *re_cuts, uint32_t resolution, int use_gap_seq, uint8_t mq);
inter_link_mat_t *inter_link_mat_from_file(const char *f, cov_norm_t *cov_norm, asm_dict_t *dict, re_cuts_t *re_cuts, uint32_t resolution, uint32_t radius, uint8_t mq);
cov_norm_t *cov_norm_from_file(const char *f, sdict_t *dict, uint32_t window);
intra_link_t *get_intra_link(intra_link_mat_t *link_mat, uint32_t i, uint32_t j);
inter_link_t *get_inter_link(inter_link_mat_t *link_mat, uint32_t i, uint32_t j);
norm_t *calc_norms(intra_link_mat_t *link_mat);
int inter_link_norms(inter_link_mat_t *link_mat, norm_t *norm, int use_estimated_noise, double max_noise_ratio, double *la);
void inter_link_weighted_norms(inter_link_mat_t *link_mat, norm_t *norm);
void print_norms(FILE *fp, norm_t *norm);
void print_intra_links(FILE *fp, intra_link_mat_t *link_mat, sdict_t *dict);
void print_inter_links(FILE *fp, inter_link_mat_t *link_mat, asm_dict_t *dict);
void print_inter_link_norms(FILE *fp, inter_link_mat_t *link_mat, asm_dict_t *dict);
void print_inter_link_bands(FILE *fp, inter_link_mat_t *link_mat, asm_dict_t *dict);
void intra_link_mat_destroy(intra_link_mat_t *link_mat);
void inter_link_mat_destroy(inter_link_mat_t *link_mat);
void norm_destroy(norm_t *norm);
double *get_max_inter_norms(inter_link_mat_t *link_mat, asm_dict_t *dict);
int8_t *calc_link_directs_from_file(const char *f, asm_dict_t *dict, uint8_t mq);
void calc_link_directs(inter_link_mat_t *link_mat, double min_norm, asm_dict_t *dict, int8_t *directs);
void dump_links_from_bam_file(const char *f, const char *fai, uint32_t ml, uint8_t mq, uint32_t wd, double q_drop, const char *out);
void dump_links_from_bed_file(const char *f, const char *fai, uint32_t ml, uint8_t mq, uint32_t wd, double q_drop, const char *out);
long estimate_intra_link_mat_init_rss(void *dict, uint32_t resolution, int use_gap_seq);
long estimate_inter_link_mat_init_rss(asm_dict_t *dict, uint32_t resolution, uint32_t radius);
#ifdef __cplusplus
}
#endif

#endif /* LINK_H_ */

