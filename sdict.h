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
 * 15/04/21 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#ifndef SDICT_H_
#define SDICT_H_

#include <stdint.h>

#include "khash.h"

extern char comp_table[128];
extern char nucl_toupper[128];

typedef struct {
    char *name; // seq id
    char *seq; // sequence
    uint32_t len; // seq length
} sd_seq_t;

KHASH_MAP_INIT_STR(sdict, uint32_t)
typedef khash_t(sdict) sdhash_t;

typedef struct {
    uint32_t n, m; // n: seq number, m: memory allocated
    sd_seq_t *s; // sequence dictionary
    sdhash_t *h; // sequence hash map: name -> index
} sdict_t;

typedef struct {
    uint32_t s; // seq id
    uint32_t k; // seg index on seq
    uint64_t a; // seq start
    uint32_t c, x, y; // subseq c: id << 1 | ori, x: start, y: length
} sd_seg_t;

typedef struct {
    char *name; // seq id
    uint64_t len; // seq length
    uint32_t n; // seg number
    uint32_t s; // seg start position
} sd_aseq_t;

// used to store scaffolds
// quickly convert a contig position to scaffold position
// consider contig breaks
typedef struct {
    uint32_t n, m; // n: seq number, m: seq memory allocated
    sd_aseq_t *s; // sequence dictionary
    sdhash_t *h; // sequence hash map: name -> index
    uint32_t u, v; // u: seg number, v: seg memory allocated
    sd_seg_t *seg; // segments
    uint32_t *a; // sub sequence index map: id -> start pos, need this to deal with sub seq breaks
    uint64_t *index; // sub seq end << 32 | seg index, used to find the seg index given a sub seq position
    sdict_t *sdict; // sub sequence dictionary
} asm_dict_t;

#ifdef __cplusplus
extern "C" {
#endif

sdict_t *sd_init(void);
asm_dict_t *asm_init(sdict_t *dict);
void sd_destroy(sdict_t *d);
void asm_destroy(asm_dict_t *d);
uint32_t sd_put(sdict_t *d, const char *name, uint32_t len);
uint32_t sd_put1(sdict_t *d, const char *name, const char *seq, uint32_t len);
uint32_t sd_get(sdict_t *d, const char *name);
sdict_t *make_sdict_from_fa(const char *f, uint32_t min_len);
sdict_t *make_sdict_from_index(const char *f, uint32_t min_len);
sdict_t *make_sdict_from_gfa(const char *f, uint32_t min_len);
asm_dict_t *make_asm_dict_from_sdict(sdict_t *sdict);
uint32_t asm_put(asm_dict_t *d, const char *name, uint64_t len, uint32_t n, uint32_t s);
asm_dict_t *make_asm_dict_from_agp(sdict_t *sdict, const char *f);
void add_unplaced_short_seqs(asm_dict_t *d, uint32_t min_len);
char *get_asm_seq(asm_dict_t *d, char *name);
uint32_t asm_sd_get(asm_dict_t *d, const char *name);
int sd_coordinate_conversion(asm_dict_t *d, uint32_t id, uint32_t pos, uint32_t *s, uint64_t *p, int count_gap);
void sd_stats(sdict_t *d, uint64_t *n_stats, uint32_t *l_stats);
void asm_sd_stats(asm_dict_t *d, uint64_t *n_stats, uint32_t *l_stats);
void write_fasta_file_from_agp(const char *fa, const char *agp, FILE *fo, int line_wd);
void write_segs_to_agp(sd_seg_t *segs, uint32_t n, sdict_t *sd, uint32_t s, FILE *fp);
void write_sorted_agp(asm_dict_t *dict, FILE *fo);
void write_sdict_to_agp(sdict_t *sdict, char *out);
void write_asm_dict_to_agp(asm_dict_t *dict, char *out);
#ifdef __cplusplus
}
#endif

#endif /* SDICT_H_ */

