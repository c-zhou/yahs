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
#ifndef COV_H_
#define COV_H_

#include <stdint.h>
#include <stdio.h>

#include "sdict.h"

typedef struct {
    size_t n, m;
    // (uint64_t) pos << 32 | (uint32_t) count
    // count can be negative
    uint64_t *a;
} pos_t;

typedef struct {
    uint32_t n; // seq number
    pos_t *p;
} cov_t;

typedef struct {
    uint64_t n; // total numbers
    uint32_t w; // window size
    double **norm;
} cov_norm_t;

#ifdef __cplusplus
extern "C" {
#endif

cov_t *cov_init(uint32_t n);
void cov_destroy(cov_t *cov);
void cov_norm_destroy(cov_norm_t *cov_norm);
cov_t *bam_cstats(const char *bam, sdict_t *sdict, int match_header);
cov_t *bed_cstats(const char *bed, sdict_t *sdict);
double calc_avg_cov(cov_t *cov, sdict_t *sdict, double q_drop);
cov_norm_t *calc_cov_norms(cov_t *cov, sdict_t *sdict, uint32_t window, double q_drop);
void print_cov_in_bed(cov_t *cov, sdict_t *sdict, FILE *fo);

#ifdef __cplusplus
}
#endif

#endif /* COV_H_ */

