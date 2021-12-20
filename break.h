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
#ifndef BREAK_H_
#define BREAK_H_

#include <stdlib.h>
#include <stdint.h>
#include "sdict.h"

#define SQRT2 1.41421356237
#define SQRT2_2 0.70710678118

typedef struct {
    uint32_t s; // seq id
    uint32_t n; // link bin number
    int64_t *link; // bin pos << 32 | link count 
} link_t;

typedef struct {
    uint32_t b; // bin size
    uint32_t n; // number seqs
    link_t *link;
} link_mat_t;

typedef struct {
    uint32_t n, m, s; //s: seq id
    uint64_t *p;
} bp_t;

#ifdef __cplusplus
extern "C" {
#endif

link_t *link_init(uint32_t s, uint32_t n);
link_mat_t *link_mat_init(asm_dict_t *dict, uint32_t b);
link_mat_t *link_mat_from_file(const char *f, asm_dict_t *dict, uint32_t dist_thres, uint32_t resolution, double noise, uint32_t move_avg);
uint32_t estimate_dist_thres_from_file(const char *f, asm_dict_t *dict, double min_frac, uint32_t resolution);
void link_mat_destroy(link_mat_t *link_mat);
void print_link_mat(link_mat_t *link_mat, asm_dict_t *dict, FILE *fp);
bp_t *detect_break_points(link_mat_t *link_mat, uint32_t bin_size, uint32_t merge_size, double fold_thres, uint32_t dual_break_thres, uint32_t *bp_n);
void print_break_point(bp_t *bp, asm_dict_t *dict, FILE *fp);
bp_t *detect_break_points_local_joint(link_mat_t *link_mat, uint32_t bin_size, double fold_thres, uint32_t flank_size, asm_dict_t *dict, uint32_t *bp_n);
void write_break_agp(asm_dict_t *d, bp_t *breaks, uint32_t b_n, FILE *fp);

#ifdef __cplusplus
}
#endif

#endif /* BREAK_H_ */

