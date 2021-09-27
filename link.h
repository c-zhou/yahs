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

#include "sdict.h"

#define SQRT2 1.41421356237
#define SQRT2_2 0.70710678118

typedef struct {
    int32_t c;
    int32_t n;
    double *link;
} intra_link_t;

typedef struct {
    int32_t c0, c1; // sequence id
    int32_t b0, b1; // number of bands
    int32_t r; // real radius MIN(R, b0 + b1 - 1)
    int32_t n, n0; // number of cells, n = b0 * b1, n0 is the real number of cells
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
    int32_t n;
    intra_link_t *links;
} intra_link_mat_t;

typedef struct {
    int32_t n;
    int32_t r; // radius
    double noise; // noise
    inter_link_t *links;
} inter_link_mat_t;

typedef struct {
    int32_t n; // number of bands
    int32_t *bs; // number of cells in each band [1 x n]
    double **link; // number of links in each cell [n x b]
    double *norms; // norms [1 x n]
    double *linkc; // link count in each band [1 x n]
    int32_t r; // number of first r bands contains at least 90% links adjusted by norms
} norm_t;

#ifdef __cplusplus 
extern "C" {
#endif

intra_link_mat_t *intra_link_mat_init(asm_dict_t *dict, int resolution);
intra_link_mat_t *intra_link_mat_init_sdict(sdict_t *dict, int resolution);
inter_link_mat_t *inter_link_mat_init(asm_dict_t *dict, int resolution, int radius);
intra_link_mat_t *intra_link_mat_from_file(const char *f, asm_dict_t *dict, int resolution, int use_gap_seq);
inter_link_mat_t *inter_link_mat_from_file(const char *f, asm_dict_t *dict, int resolution, int radius);
intra_link_t *get_intra_link(intra_link_mat_t *link_mat, int32_t i, int32_t j);
inter_link_t *get_inter_link(inter_link_mat_t *link_mat, int32_t i, int32_t j);
norm_t *calc_norms(intra_link_mat_t *link_mat);
void inter_link_norms(inter_link_mat_t *link_mat, norm_t *norm, int use_estimated_noise);
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
void calc_link_directs(inter_link_mat_t *link_mat, double min_norm, asm_dict_t *dict);
void dump_links_from_bam_file(const char *f, const char *fai, const char *out);
void dump_links_from_bed_file(const char *f, const char *fai, const char *out);

#ifdef __cplusplus
}
#endif

