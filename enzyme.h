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
 * 14/12/21 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#ifndef ENZYME_H_
#define ENZYME_H_

#include <stdint.h>

#include "sdict.h"

typedef struct {
    uint32_t l, n; // seq len, number cuts
    uint32_t *sites; // cutting sites
} re_t;

typedef struct {
    uint32_t n; // sequence number
    double density; // number cuts per base
    re_t *re; // cutting sites
} re_cuts_t;

#ifdef __cplusplus
extern "C" {
#endif

re_cuts_t *re_cuts_init(uint32_t n);
void re_cuts_destroy(re_cuts_t *re_cuts);
re_cuts_t *find_re_from_seqs(const char *f, uint32_t ml, char **enz_cs, int enz_n);
double **calc_re_cuts_density(re_cuts_t *re_cuts, uint32_t resolution);
double **calc_re_cuts_density1(re_cuts_t *re_cuts, uint32_t resolution, asm_dict_t *dict);
double **calc_re_cuts_density2(re_cuts_t *re_cuts, uint32_t resolution, asm_dict_t *dict);

#ifdef __cplusplus
}
#endif

#endif /* ENZYME_H_ */

