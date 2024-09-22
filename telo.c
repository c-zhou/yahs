/*********************************************************************************
 * MIT License                                                                   *
 *                                                                               *
 * Copyright (c) 2024 Chenxi Zhou <chnx.zhou@gmail.com>                          *
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
 * 11/09/24 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/

/** 
 *  A telomere sequence finder adopted from seqtk https://github.com/lh3/seqtk 
 ***/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "khash.h"
#include "sdict.h"

KHASH_SET_INIT_INT64(64)

int telo_penalty  = 1;
int telo_maxDrop  = 2000;
int telo_minScore = 300;
char *telo_motif  = NULL;

char *telo_motif_db[] = {
    "AAAATTGTCCGTCC",
    "AAACCACCCT",
    "AAACCC",
    "AAACCCC",
    "AAACCCT",
    "AAAGAACCT",
    "AAATGTGGAGG",
    "AACAGACCCG",
    "AACCATCCCT",
    "AACCC",
    "AACCCAGACCC",
    "AACCCAGACCT",
    "AACCCAGACGC",
    "AACCCCAACCT",
    "AACCCGAACCT",
    "AACCCT",
    "AACCCTG",
    "AACCCTGACGC",
    "AACCT",
    "AAGGAC",
    "ACCCAG",
    "ACCTG",
    "ACGGCAGCG"
};

unsigned char seq_nt6_table[256] = {
    0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,

    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

int check_motif(char *motif)
{
	int i, c, n;

	if (motif == NULL)
		return 0;

	n = strlen(motif);
	for (i = 0; i < n; i++) {
		c = seq_nt6_table[(uint8_t) motif[i]];
		if (c < 1 || c > 4)
			return 1;
	}

	return 0;
}

static uint64_t telo_finder_core(char *sequence, int slen, kh_64_t *mtab, int mlen, int penalty, int max_drop, int min_score)
{
	int c, hit;
	uint64_t x, mask, sum_telo;
	ssize_t i, l, max_i = -1, st = 0;
	int64_t score, max;
	
	mask = (1ULL<<2*mlen) - 1;
	sum_telo = 0;
	score = max = 0, max_i = -1;
	for (i = 0, l = 0, x = 0; i < slen; ++i) { // 5'-end
		hit = 0, c = seq_nt6_table[(uint8_t) sequence[i]];
		if (c >= 1 && c <= 4) { // not N
			x = (x<<2 | (c-1)) & mask;
			if (++l >= mlen && kh_get(64, mtab, x) != kh_end(mtab)) // x is at least mlen long and is a telomere motif
				hit = 1;
		} else l = 0, x = 0; // N, ambiguous base
		if (i >= mlen) score += hit? 1 : -penalty;
		if (score > max) max = score, max_i = i;
		else if (max - score > max_drop) break;
	}
	if (max >= min_score) {
		sum_telo += (uint64_t) (max_i + 1) << 32;
		st = max_i + 1;
	}
	score = max = 0, max_i = -1;
	for (i = slen - 1, l = 0, x = 0; i >= st; --i) { // 3'-end; similar to the for loop above
		hit = 0, c = seq_nt6_table[(uint8_t) sequence[i]];
		if (c >= 1 && c <= 4) {
			x = (x<<2 | (4-c)) & mask;
			if (++l >= mlen && kh_get(64, mtab, x) != kh_end(mtab))
				hit = 1;
		} else l = 0, x = 0;
		if (slen - i >= mlen) score += hit? 1 : -penalty;
		if (score > max) max = score, max_i = i;
		else if (max - score > max_drop) break;
	}
	if (max >= min_score)
		sum_telo += slen - max_i;
	
	return sum_telo;
}

int8_t *telo_finder(const char *f, uint32_t ml)
{
	int i, j, t, c, len, absent;
	uint32_t nseq, ntelo;
	uint64_t x, telo;
	int8_t *telo_ends;
	sdict_t *sdict;
	char **telo_db, *motif;
	khash_t(64) *h;

	sdict = make_sdict_from_fa(f, ml);
    nseq  = sdict->n;
	telo_ends = (int8_t *) calloc(nseq * 2, sizeof(int8_t));

	ntelo = 0;
	telo_db = NULL;
	if (telo_motif != NULL) {
		telo_db = (char **) malloc(sizeof(char *));
		telo_db[0] = telo_motif;
		ntelo = 1;
	} else {
		telo_db = telo_motif_db;
		ntelo = sizeof(telo_motif_db) / sizeof(char *);
	}

	h = kh_init(64); // hash table for all roations of the telomere motif
	for (t = 0; t < ntelo; t++) {
		// make motif table
		motif = telo_db[t];
		len = strlen(motif);
		kh_resize(64, h, len * 2);
		for (i = 0; i < len; ++i) {
			for (j = 0, x = 0; j < len; ++j) {
				c = seq_nt6_table[(uint8_t) motif[(i + j) % len]];
				assert(c >= 1 && c <= 4);
				x = x<<2 | (c-1);
			}
			kh_put(64, h, x, &absent);
		}

		// do telo_finder for each sequence in the dictionary
		for (i = 0; i < nseq; i++) {
			telo = telo_finder_core(sdict->s[i].seq, sdict->s[i].len, h, len, 
				telo_penalty, telo_maxDrop, telo_minScore);
			if (telo == 0)
				continue;
			if ((telo>>32) > 0) {
				telo_ends[i<<1] = 1;
				fprintf(stderr, "[I::%s] found telo motif %s in sequence %s 5'-end up to position %u\n", 
					__func__, telo_db[t], sdict->s[i].name, (uint32_t) (telo>>32));
			}
			if ((uint32_t) telo > 0) {
				telo_ends[i<<1|1] = 1;
				fprintf(stderr, "[I::%s] found telo motif %s in sequence %s 3'-end from position %u\n", 
					__func__, telo_db[t], sdict->s[i].name, sdict->s[i].len - (uint32_t) telo);
			}
		}

		kh_clear(64, h);
	}

	kh_destroy(64, h);
	if (telo_motif != NULL)
		free(telo_db);
	sd_destroy(sdict);

	return telo_ends;
}

#undef TELO_FINDER_MAIN

#ifdef TELO_FINDER_MAIN
int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "[I::%s] provide at least a sequence input\n", __func__);
		exit (1);
	}

	if (argc > 2)
		telo_motif = argv[2];

	int8_t *telo_ends = telo_finder(argv[1], 0);
	free(telo_ends);
}
#endif