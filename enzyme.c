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
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <ctype.h>

#include "kvec.h"
#include "enzyme.h"
#include "sdict.h"

#undef DEBUG_ENZ

re_cuts_t *re_cuts_init(uint32_t n)
{
    re_cuts_t *re_cuts = (re_cuts_t *) malloc(sizeof(re_cuts_t));
    re_cuts->n = n;
    re_cuts->density = 0;
    re_cuts->re = (re_t *) malloc(n * sizeof(re_t));
    uint32_t i;
    for (i = 0; i < n; ++i) {
        re_cuts->re[i].sites = 0;
        re_cuts->re[i].dens = 0;
    }
    return re_cuts;
}

void re_cuts_destroy(re_cuts_t *re_cuts)
{
    uint32_t i;
    for (i = 0; i < re_cuts->n; ++i) {
        if (re_cuts->re[i].sites)
            free(re_cuts->re[i].sites);
        if (re_cuts->re[i].dens)
            free(re_cuts->re[i].dens);
    }
    free(re_cuts->re);
    free(re_cuts);
}

typedef struct {size_t n, m; uint32_t *a;} u32_v;
int u32_cmp(const void *p, const void *q)
{
    int64_t d = *(int64_t *) p - *(int64_t *) q;
    return d == 0? 0 : (d > 0? 1 : -1);
}

re_cuts_t *find_re_from_seqs(const char *f, char **enz_cs, int enz_n)
{
    // now find all RE cutting sites
    int i, j, c, c0, c1;
    uint32_t len, n, p, n_re;
    int64_t genome_size;
    char *pch, *seq, *re;
    sdict_t *sdict;
    re_cuts_t *re_cuts;

    sdict = make_sdict_from_fa(f);
    n = sdict->n;
    n_re = 0;
    genome_size = 0;
    re_cuts = re_cuts_init(n);
    for (i = 0; i < n; ++i) {
        u32_v enz_cs_pos = {0, 0, 0};

        seq = sdict->s[i].seq;
        len = sdict->s[i].len;
        for (p = 0; p < len; ++p) {
            c = seq[p];
            if (!isalpha(c)) {
                fprintf(stderr, "[E::%s] non-alphabetic chacrater in FASTA file: %c\n", __func__, c);
                re_cuts_destroy(re_cuts);
                sd_destroy(sdict);
                return 0;
            }
            seq[p] = nucl_toupper[c];
        }
        
        for (j = 0; j < enz_n; ++j) {
            re = enz_cs[j];
            pch = strstr(seq, re);
            while (pch != NULL) {
                kv_push(uint32_t, enz_cs_pos, pch - seq);
                pch = strstr(pch + 1, re);
            }
        }

        // reverse complement sequence
        for (p = 0; p < len>>1; ++p) {
            c0 = comp_table[(int) seq[p]];
            c1 = comp_table[(int) seq[len - 1 - p]];
            seq[p] = c1;
            seq[len - 1 - p] = c0;
        }
        if (len & 1) // complement the remaining base
            seq[len>>1] = comp_table[(int) seq[len>>1]];
        
        for (j = 0; j < enz_n; ++j) {
            re = enz_cs[j];
            pch = strstr(seq, re);
            while (pch != NULL) {
                kv_push(uint32_t, enz_cs_pos, len - 1 - (pch - seq));
                pch = strstr(pch + 1, re);
            }
        }

        qsort(enz_cs_pos.a, enz_cs_pos.n, sizeof(uint32_t), u32_cmp);
        re_cuts->re[i].sites = enz_cs_pos.a;
        re_cuts->re[i].n = enz_cs_pos.n;
        re_cuts->re[i].l = len;
        
        n_re += enz_cs_pos.n;
        genome_size += len;
    }

    re_cuts->density = (double) n_re / genome_size;

    fprintf(stderr, "[I::%s] NO. restriction enzyme cutting sites found in sequences: %u\n", __func__, n_re);
    fprintf(stderr, "[I::%s] restriction enzyme cutting sites density: %.6f\n", __func__, re_cuts->density);
#ifdef DEBUG_ENZ
    printf("[I::%s] restriction enzyme cutting sites for individual sequences (n = %d)\n", __func__, n);
    for (i = 0; i < n; ++i)
        printf("[I::%s] %s %u %u %.6f\n", __func__, sdict->s[i].name, sdict->s[i].len, re_cuts->re[i].n, (double) re_cuts->re[i].n / sdict->s[i].len);
#endif

    sd_destroy(sdict);

    return re_cuts;
}

