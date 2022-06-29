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
 * 22/06/21 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdlib.h>
#include <stdio.h>

#include "khash.h"
#include "bamlite.h"
#include "ketopt.h"
#include "sdict.h"
#include "asset.h"

#define JUICER_VERSION "1.1"

static double jc_realtime0;

KHASH_SET_INIT_STR(str)

static int make_juicer_pre_file_from_bin(char *f, char *agp, char *fai, int scale, int count_gap, FILE *fo)
{
    FILE *fp;
    uint32_t i, i0, i1;
    uint64_t p0, p1;
    uint32_t buffer[BUFF_SIZE], m;
    long pair_c;

    sdict_t *sdict = make_sdict_from_index(fai, 0);
    asm_dict_t *dict = agp? make_asm_dict_from_agp(sdict, agp) : make_asm_dict_from_sdict(sdict);

    fp = fopen(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    pair_c = 0;
    while (1) {
        m = fread(&buffer, sizeof(uint32_t), BUFF_SIZE, fp);
        for (i = 0; i < m; i += 4) {
            sd_coordinate_conversion(dict, buffer[i], buffer[i + 1], &i0, &p0, count_gap);
            sd_coordinate_conversion(dict, buffer[i + 2], buffer[i + 3], &i1, &p1, count_gap);
            if (i0 == UINT32_MAX || i1 == UINT32_MAX) {
                fprintf(stderr, "[W::%s] sequence not found \n", __func__);
            } else {
                if (strcmp(dict->s[i0].name, dict->s[i1].name) <= 0)
                    fprintf(fo, "0\t%s\t%lu\t0\t1\t%s\t%lu\t1\n", dict->s[i0].name, p0 >> scale, dict->s[i1].name, p1 >> scale);
                else
                    fprintf(fo, "0\t%s\t%lu\t1\t1\t%s\t%lu\t0\n", dict->s[i1].name, p1 >> scale, dict->s[i0].name, p0 >> scale);
            }
        }
        pair_c += m / 4;

        if (m < BUFF_SIZE) {
            if (ferror(fp))
                return 1;
            break;
        }
    }

    fprintf(stderr, "[I::%s] %ld read pairs processed\n", __func__, pair_c);
    fclose(fp);
    asm_destroy(dict);
    sd_destroy(sdict);

    return 0;
}

static int make_juicer_pre_file_from_bed(char *f, char *agp, char *fai, uint8_t mq, int scale, int count_gap, FILE *fo)
{
    FILE *fp;
    char *line = NULL;
    size_t ln = 0;
    ssize_t read;
    char cname0[4096], cname1[4096], rname0[4096], rname1[4096];
    uint32_t s0, s1, e0, e1, i0, i1;
    uint64_t p0, p1;
    int8_t buff;
    long rec_c, pair_c;
    
    khash_t(str) *hmseq; // for absent sequences
    khint_t k;
    int absent;
    hmseq = kh_init(str);

    sdict_t *sdict = make_sdict_from_index(fai, 0);
    asm_dict_t *dict = agp? make_asm_dict_from_agp(sdict, agp) : make_asm_dict_from_sdict(sdict);

    fp = fopen(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    s0 = s1 = e0 = e1 = 0;
    i0 = i1 = 0;
    p0 = p1 = 0;
    rec_c = pair_c = 0;
    buff = 0;
    while ((read = getline(&line, &ln, fp)) != -1) {
        if (buff == 0) {
            sscanf(line, "%s %u %u %s", cname0, &s0, &e0, rname0);
            ++buff;
        } else if (buff == 1) {
            sscanf(line, "%s %u %u %s", cname1, &s1, &e1, rname1);
            if (is_read_pair(rname0, rname1)) {
                ++pair_c;
                
                sd_coordinate_conversion(dict, sd_get(sdict, cname0), s0 / 2 + e0 / 2 + (s0 & 1 && e0 & 1), &i0, &p0, count_gap);
                sd_coordinate_conversion(dict, sd_get(sdict, cname1), s1 / 2 + e1 / 2 + (s1 & 1 && e1 & 1), &i1, &p1, count_gap);
                if (i0 == UINT32_MAX) {
                    k = kh_put(str, hmseq, cname0, &absent);
                    if (absent) {
                        kh_key(hmseq, k) = strdup(cname0);
                        fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname0);
                    }
                } else if(i1 == UINT32_MAX) {
                    k = kh_put(str, hmseq, cname1, &absent);
                    if (absent) {
                        kh_key(hmseq, k) = strdup(cname1);
                        fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname0);
                    }
                } else {
                    if (strcmp(dict->s[i0].name, dict->s[i1].name) <= 0)
                        fprintf(fo, "0\t%s\t%lu\t0\t1\t%s\t%lu\t1\n", dict->s[i0].name, p0 >> scale, dict->s[i1].name, p1 >> scale);
                    else
                        fprintf(fo, "0\t%s\t%lu\t1\t1\t%s\t%lu\t0\n", dict->s[i1].name, p1 >> scale, dict->s[i0].name, p0 >> scale);
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

    fprintf(stderr, "[I::%s] %ld read pairs processed\n", __func__, pair_c);
    
    for (k = 0; k < kh_end(hmseq); ++k)
        if (kh_exist(hmseq, k))
            free((char *) kh_key(hmseq, k));
    kh_destroy(str, hmseq);

    if (line)
        free(line);
    fclose(fp);
    asm_destroy(dict);
    sd_destroy(sdict);

    return 0;
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

static char *parse_bam_rec1(bam1_t *b, bam_header_t *h, char **cname0, int32_t *s0, char **cname1, int32_t *s1)
{
    // 0x4 0x8 0x40 0x100 0x400 0x800
    if (b->core.flag & 0xD4C)
        return 0;
    *cname0 = h->target_name[b->core.tid];
    *s0 = b->core.pos + 1;
    *cname1 = h->target_name[b->core.mtid];
    *s1 = b->core.mpos + 1;

    return bam1_qname(b);
}

static int make_juicer_pre_file_from_bam(char *f, char *agp, char *fai, uint8_t mq, int scale, int count_gap, FILE *fo)
{
    bamFile fp;
    bam_header_t *h;
    bam1_t *b;
    char *cname0, *cname1, *rname0, *rname1;
    int32_t s0, s1, e0, e1;
    uint32_t i0, i1;
    uint64_t p0, p1;
    int8_t buff;
    long rec_c, pair_c;
    enum bam_sort_order so;

    khash_t(str) *hmseq; // for absent sequences
    khint_t k;
    int absent;
    hmseq = kh_init(str);

    sdict_t *sdict = make_sdict_from_index(fai, 0);
    asm_dict_t *dict = agp? make_asm_dict_from_agp(sdict, agp) : make_asm_dict_from_sdict(sdict);

    fp = bam_open(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }
    
    h = bam_header_read(fp);
    b = bam_init1();
    so = bam_hrecs_sort_order(h);
    cname0 = cname1 = rname0 = rname1 = 0;
    s0 = s1 = e0 = e1 = 0;
    i0 = i1 = 0;
    p0 = p1 = 0;
    rec_c = pair_c = 0;
    buff = 0;    
    
    if (so == ORDER_NAME) {
        // sorted by read names
        while (bam_read1(fp, b) >= 0 ) {
            if (++rec_c % 1000000 == 0)
                fprintf(stderr, "[I::%s] %ld million records processed, %ld read pairs \n", __func__, rec_c / 1000000, pair_c);

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

                    if (s0 > 0 && s1 >0) {
                        sd_coordinate_conversion(dict, sd_get(sdict, cname0), s0 / 2 + e0 / 2 + (s0 & 1 && e0 & 1), &i0, &p0, count_gap);
                        sd_coordinate_conversion(dict, sd_get(sdict, cname1), s1 / 2 + e1 / 2 + (s1 & 1 && e1 & 1), &i1, &p1, count_gap);

                        if (i0 == UINT32_MAX) {
                            k = kh_put(str, hmseq, cname0, &absent);
                            if (absent) {
                                kh_key(hmseq, k) = strdup(cname0);
                                fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname0);
                            }
                        } else if (i1 == UINT32_MAX) {
                            k = kh_put(str, hmseq, cname1, &absent);
                            if (absent) {
                                kh_key(hmseq, k) = strdup(cname1);
                                fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname1);
                            }
                        } else {
                            if (strcmp(dict->s[i0].name, dict->s[i1].name) <= 0)
                                fprintf(fo, "0\t%s\t%lu\t0\t1\t%s\t%lu\t1\n", dict->s[i0].name, p0 >> scale, dict->s[i1].name, p1 >> scale);
                            else
                                fprintf(fo, "0\t%s\t%lu\t1\t1\t%s\t%lu\t0\n", dict->s[i1].name, p1 >> scale, dict->s[i0].name, p0 >> scale);
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
        }
    } else {
        // sorted by coordinates or others
        if (mq > 0)
            fprintf(stderr, "[W::%s] BAM file is not sorted by read name. Filtering by mapping quality %hu suppressed \n", __func__, mq);

        while (bam_read1(fp, b) >= 0 ) {
            if (++rec_c % 1000000 == 0)
                fprintf(stderr, "[I::%s] %ld million records processed, %ld read pairs \n", __func__, rec_c / 1000000, pair_c); 
            
            if(!parse_bam_rec1(b, h, &cname0, &s0, &cname1, &s1))
                continue;

            ++pair_c;
            sd_coordinate_conversion(dict, sd_get(sdict, cname0), s0, &i0, &p0, count_gap);
            sd_coordinate_conversion(dict, sd_get(sdict, cname1), s1, &i1, &p1, count_gap);

            if (i0 == UINT32_MAX) {
                k = kh_put(str, hmseq, cname0, &absent);
                if (absent) {
                    kh_key(hmseq, k) = strdup(cname0);
                    fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname0);
                }
            } else if (i1 == UINT32_MAX) {
                k = kh_put(str, hmseq, cname1, &absent);
                if (absent) {
                    kh_key(hmseq, k) = strdup(cname1);
                    fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname1);
                }
            } else {
                if (strcmp(dict->s[i0].name, dict->s[i1].name) <= 0)
                    fprintf(fo, "0\t%s\t%lu\t0\t1\t%s\t%lu\t1\n", dict->s[i0].name, p0 >> scale, dict->s[i1].name, p1 >> scale);
                else
                    fprintf(fo, "0\t%s\t%lu\t1\t1\t%s\t%lu\t0\n", dict->s[i1].name, p1 >> scale, dict->s[i0].name, p0 >> scale);
            }
        }
    }

    fprintf(stderr, "[I::%s] %ld read pairs processed\n", __func__, pair_c);

    for (k = 0; k < kh_end(hmseq); ++k)
        if (kh_exist(hmseq, k))
            free((char *) kh_key(hmseq, k));
    kh_destroy(str, hmseq);

    if (rname0)
        free(rname0);
    if (rname1)
        free(rname1);
    bam_destroy1(b);
    bam_header_destroy(h);
    bam_close(fp);
    asm_destroy(dict);
    sd_destroy(sdict);

    return 0;
}

static uint64_t assembly_annotation(const char *f, const char *out_agp, const char *out_annot, const char *out_lift, int *scale, uint64_t max_s, uint64_t *g)
{
    FILE *fp, *fo_agp, *fo_annot, *fo_lift;
    char *line;
    size_t ln = 0;
    ssize_t read;
    char sname[256], type[4], cname[256], cstarts[16], cends[16], oris[4];
    char *name;
    uint32_t c, s, l, cstart, cend;
    uint64_t genome_size, scaled_gs;

    fp = fopen(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }
    fo_agp = fopen(out_agp, "w");
    fo_annot = fopen(out_annot, "w");
    fo_lift = fopen(out_lift, "w");
    
    line = name = NULL;
    genome_size = 0;
    c = s = 0;
    while ((read = getline(&line, &ln, fp)) != -1) {
        sscanf(line, "%s %*s %*s %*s %s %s %s %s %s", sname, type, cname, cstarts, cends, oris);
        if (!strncmp(type, "N", 1))
            continue;
        if (!name)
            name = strdup(sname);
        if (strcmp(sname, name)) {
            name = strdup(sname);
            ++s;
        }
        l = strtoul(cends, NULL, 10) - strtoul(cstarts, NULL, 10) + 1;
        fprintf(fo_agp, "assembly\t%lu\t%lu\t%u\tW\t%s\t%s\t%s\t%s\n", genome_size + 1, genome_size + l, ++c, cname, cstarts, cends, oris);
        genome_size += l;
    }
    fclose(fp);

    scaled_gs = linear_scale(genome_size, scale, max_s);

    int *seqs;
    fp = fopen(f, "r");
    line = name = NULL;
    seqs = (int *) calloc(c + s, sizeof(int));
    c = s = 0;
    while ((read = getline(&line, &ln, fp)) != -1) {
        sscanf(line, "%s %*s %*s %*s %s %s %s %s %s", sname, type, cname, cstarts, cends, oris);
        if (!strncmp(type, "N", 1))
            continue;
        if (!name)
             name = strdup(sname);
        if (strcmp(sname, name)) {
            name = strdup(sname);
            ++s;
        }
        cstart = strtoul(cstarts, NULL, 10);
        cend = strtoul(cends, NULL, 10);
        l = cend - cstart + 1;
        // l = l >> *scale;
        ++c;
        fprintf(fo_annot, ">ctg%08u.1 %u %u\n", c, c, l);
        // fprintf(fo_lift, "ctg%08u.1 %s %u %u\n", c, cname, cstart, cend);
        fprintf(fo_lift, "ctg%08u.1\t%u\t%u\t%u\tW\t%s\t%u\t%u\t+\n", c, 1, l, 1, cname, cstart, cend);
        seqs[c + s - 1] = c * (strncmp(oris, "+", 1)? -1 : 1);
    }
    uint32_t i;
    for (i = 0; i < c + s - 1; ++i) {
        if (seqs[i] == 0) {
            fprintf(fo_annot, "\n");
        } else {
            fprintf(fo_annot, "%d", seqs[i]);
            if (seqs[i + 1] != 0)
                fprintf(fo_annot, " ");
        }
    }
    fprintf(fo_annot, "%d\n", seqs[c + s - 1]);
    fclose(fp);
    
    fclose(fo_agp);
    fclose(fo_annot);
    fclose(fo_lift);

    if (line)
        free(line);
    if (name)
        free(name);
    free(seqs);
    
    *g = genome_size;

    return scaled_gs;
}

uint64_t assembly_scale_max_seq(asm_dict_t *dict, int *scale, uint64_t max_s, uint64_t *g)
{
    uint32_t i;
    uint64_t s;
    sd_aseq_t seq;
 
    s = 0;
    for (i = 0; i < dict->n; ++i) {
        seq = dict->s[i];
        s = MAX(s, seq.len + (seq.n - 1) * GAP_SZ);
    }
    
    *g = s;

    return linear_scale(s, scale, max_s);
}

static void print_help_pre(FILE *fp_help)
{
    fprintf(fp_help, "Usage: juicer pre [options] <hic.bed>|<hic.bam>|<hic.bin> <scaffolds.agp> <contigs.fa.fai>\n");
    fprintf(fp_help, "Options:\n");
    fprintf(fp_help, "    -a                preprocess for assembly mode\n");
    fprintf(fp_help, "    -q INT            minimum mapping quality [10]\n");
    fprintf(fp_help, "    -o STR            output file prefix (required for '-a' mode) [stdout]\n");
    fprintf(fp_help, "    --version         show version number\n");
}

static ko_longopt_t long_options[] = {
    { "help",           ko_no_argument, 'h' },
    { "version",        ko_no_argument, 'V' },
    { 0, 0, 0 }
};

static int main_pre(int argc, char *argv[])
{
    FILE *fo;
    char *fai, *agp, *agp1, *link_file, *out, *out1, *annot, *lift, *ext;
    int mq, asm_mode;;
    
    liftrlimit();
    jc_realtime0 = realtime();

    const char *opt_str = "q:ao:Vh";
    ketopt_t opt = KETOPT_INIT;
    int c, ret;
    FILE *fp_help = stderr;
    fai = agp = agp1 = link_file = out = out1 = annot = lift = 0;
    mq = 10;
    asm_mode = 0;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >= 0) {
        if (c == 'o') {
            out = opt.arg;
        } else if (c == 'q') {
            mq = atoi(opt.arg);
        } else if (c == 'a') {
            asm_mode = 1;
        } else if (c == 'h') {
            fp_help = stdout;   
        } else if (c == 'V') {
            puts(JUICER_VERSION);
            return 0;
        } else if (c == '?') {
            fprintf(stderr, "[E::%s] unknown option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        } else if (c == ':') {
            fprintf(stderr, "[E::%s] missing option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        }
    }

    if (fp_help == stdout) {
        print_help_pre(stdout);
        return 0;
    }

    if (asm_mode && !out) {
        fprintf(stderr, "[E::%s] missing input: -o option is required for assembly mode (-a)\n", __func__);
        return 1;
    }

    if (argc - opt.ind < 3) {
        fprintf(stderr, "[E::%s] missing input: three positional options required\n", __func__);
        print_help_pre(stderr);
        return 1;
    }
    
    if (mq < 0 || mq > 255) {
        fprintf(stderr, "[E::%s] invalid mapping quality threshold: %d\n", __func__, mq);
        return 1;
    }

    uint8_t mq8;
    mq8 = (uint8_t) mq;

    link_file = argv[opt.ind];
    agp = argv[opt.ind + 1];
    fai = argv[opt.ind + 2];

    if (out) {
        out1 = (char *) malloc(strlen(out) + 35);
        sprintf(out1, "%s.txt", out);
    }

    fo = out1 == 0? stdout : fopen(out1, "w");
    if (fo == 0) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }
    ret = 0;
    
    sdict_t *sdict;
    asm_dict_t *dict;
    int scale;
    uint64_t max_s, scaled_s;
    
    sdict = make_sdict_from_index(fai, 0);
    scale = 0;
    max_s = scaled_s = 0;
    agp1 = (char *) malloc(MAX(strlen(agp), out? strlen(out) : 0) + 35);
    if (asm_mode) {
        annot = (char *) malloc(strlen(out) + 35);
        lift = (char *) malloc(strlen(out) + 35);
        sprintf(agp1, "%s.assembly.agp", out);
        sprintf(annot, "%s.assembly", out);
        sprintf(lift, "%s.liftover.agp", out);
        scaled_s = assembly_annotation(agp, agp1, annot, lift, &scale, (uint64_t) INT_MAX, &max_s);
        dict = make_asm_dict_from_agp(sdict, agp1);
    } else {
        sprintf(agp1, "%s", agp);
        dict = make_asm_dict_from_agp(sdict, agp1);
        scaled_s = assembly_scale_max_seq(dict, &scale, (uint64_t) INT_MAX, &max_s);
    }

    ext = link_file + strlen(link_file) - 4;
    if (strcmp(ext, ".bam") == 0) {
        fprintf(stderr, "[I::%s] make juicer pre input from BAM file %s\n", __func__, link_file);
        ret = make_juicer_pre_file_from_bam(link_file, agp1, fai, mq8, scale, !asm_mode, fo);
    } else if (strcmp(ext, ".bed") == 0) {
        fprintf(stderr, "[I::%s] make juicer pre input from BED file %s\n", __func__, link_file);
        ret = make_juicer_pre_file_from_bed(link_file, agp1, fai, mq8, scale, !asm_mode, fo);
    } else if (strcmp(ext, ".bin") == 0) {
        fprintf(stderr, "[I::%s] make juicer pre input from BIN file %s\n", __func__, link_file);
        ret = make_juicer_pre_file_from_bin(link_file, agp1, fai, scale, !asm_mode, fo);
    } else {
        fprintf(stderr, "[E::%s] unknown link file format. File extension .bam, .bed or .bin is expected\n", __func__);
        exit(EXIT_FAILURE);
    }

    if (asm_mode) {
        fprintf(stderr, "[I::%s] genome size: %lu\n", __func__, max_s);
        fprintf(stderr, "[I::%s] scale factor: %d\n", __func__, 1 << scale);
        fprintf(stderr, "[I::%s] chromosome sizes for juicer_tools pre -\n", __func__);
        fprintf(stderr, "PRE_C_SIZE: assembly %lu\n", scaled_s);
        fprintf(stderr, "[I::%s] JUICER_PRE CMD: java -Xmx36G -jar ${juicer_tools} pre %s %s.hic <(echo \"assembly %lu\")\n", __func__, out1, out, scaled_s);
    } else {
        if (scale) {
            fprintf(stderr, "[W::%s] maximum scaffold length exceeds %d (=%lu)\n", __func__, INT_MAX, max_s);
            fprintf(stderr, "[W::%s] using scale factor: %d\n", __func__, scale);
        }
        fprintf(stderr, "[I::%s] chromosome sizes for juicer_tools pre -\n", __func__);
        uint32_t i;
        sd_aseq_t seq;
        for (i = 0; i < dict->n; ++i) {
            seq = dict->s[i];
            fprintf(stderr, "PRE_C_SIZE: %s %lu\n", seq.name, (seq.len + (seq.n - 1) * GAP_SZ) >> scale);
        }
    }

    asm_destroy(dict);
    sd_destroy(sdict);

    if (out != 0)
        fclose(fo);
    
    if (out1)
        free(out1);

    if (agp1)
        free(agp1);

    if (annot)
        free(annot);

    if (lift)
        free(lift);

    fprintf(stderr, "[I::%s] Version: %s\n", __func__, JUICER_VERSION);
    fprintf(stderr, "[I::%s] CMD: juicer", __func__);
    int i;
    for (i = 0; i < argc; ++i)
        fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n[I::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - jc_realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);

    return ret;
}

static int assembly_to_agp(char *assembly, char *lift, sdict_t *sdict, FILE *fo)
{
    asm_dict_t *dict;
    FILE *fp;

    fp = fopen(assembly, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, assembly);
        exit(EXIT_FAILURE);
    }

    dict = make_asm_dict_from_agp(sdict, lift);
    
    char *line = NULL;
    size_t ln = 0;
    ssize_t read;
    char cname[1024], c0[1024], c1[1024], *cstr;
    int32_t cid, sid, fid, sign;
    uint32_t clen, mlen, s, p, *coords;
    int64_t slen;
    size_t m;

    m = 4;
    coords = (uint32_t *) malloc(sizeof(uint32_t) * m * 3);
    c0[0] = c1[0] = '\0';
    mlen = 0;
    sid = 0;

    while ((read = getline(&line, &ln, fp)) != -1) {
        if (line[0] == '>') {
            sscanf(line, "%s %d %u", cname, &cid, &clen);
            cstr = strstr(cname, ":::");
            if (cstr != NULL) {
                cname[cstr - cname] = '\0';
            }
            strcpy(c1, cname + 1);
            
            if (!strcmp(c0, c1)) {
                mlen += clen;
            } else {
                mlen = clen;
                strcpy(c0, c1);
            }

            if (sd_coordinate_rev_conversion(dict, asm_sd_get(dict, c1), mlen - clen + 1, &s, &p, 0)) {
                fprintf(stderr, "[E::%s] coordinates conversion error %s %u\n", __func__, c1, mlen - clen + 1);
                exit(EXIT_FAILURE);
            }
            
            if (cid > m) {
                m <<= 1;
                coords = (uint32_t *) realloc(coords, sizeof(uint32_t) * m * 3);
            }

            cid -= 1;
            coords[cid * 3] = s;
            coords[cid * 3 + 1] = p;
            coords[cid * 3 + 2] = clen;
        } else {
            ++sid;
            fid = 0;
            slen = 0;

            char  *eptr, *fptr;
            
            cid = strtol(line, &eptr, 10);
            sign = cid > 0;
            cid = abs(cid) - 1;
            fprintf(fo, "scaffold_%d\t%ld\t%ld\t%d\tW\t%s\t%d\t%d\t%c\n", sid, slen + 1, slen + coords[cid * 3 + 2], ++fid, sdict->s[coords[cid * 3]].name, coords[cid * 3 + 1], coords[cid * 3 + 1] + coords[cid * 3 + 2] - 1, "-+"[sign]);
            slen += coords[cid * 3 + 2];
            while (*eptr != '\n') {
                fprintf(fo, "scaffold_%d\t%ld\t%ld\t%d\tN\t%d\tscaffold\tyes\tna\n", sid, slen + 1, slen + GAP_SZ, ++fid, GAP_SZ);
                slen += GAP_SZ;
                
                cid = strtol(eptr + 1, &fptr, 10);
                sign = cid > 0;
                cid = abs(cid) - 1;
                fprintf(fo, "scaffold_%d\t%ld\t%ld\t%d\tW\t%s\t%d\t%d\t%c\n", sid, slen + 1, slen + coords[cid * 3 + 2], ++fid, sdict->s[coords[cid * 3]].name, coords[cid * 3 + 1], coords[cid * 3 + 1] + coords[cid * 3 + 2] - 1, "-+"[sign]);
                slen += coords[cid * 3 + 2];

                eptr = fptr;
            }
        }
    }
    
    fclose(fp);

    asm_destroy(dict);

    free(coords);

    return 0;
}

static void print_help_post(FILE *fp_help)
{
    fprintf(fp_help, "Usage: juicer post [options] <review.assembly> <liftover.agp> <contigs.fa[.fai]>\n");
    fprintf(fp_help, "Options:\n");
    fprintf(fp_help, "    -o STR            output file prefix (required for scaffolds FASTA output) [stdout]\n");
    fprintf(fp_help, "    --version         show version number\n");
}

static int main_post(int argc, char *argv[])
{
    FILE *fo;
    char *fa, *fa1, *out, *out1;

    liftrlimit();
    jc_realtime0 = realtime();

    const char *opt_str = "o:Vh";
    ketopt_t opt = KETOPT_INIT;
    int c, ret, is_fai;
    FILE *fp_help = stderr;
    sdict_t *sdict;
    fa = fa1 = out = out1 = 0;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >= 0) {
        if (c == 'o') {
            out = opt.arg;
        } else if (c == 'h') {
            fp_help = stdout;
        } else if (c == 'V') {
            puts(JUICER_VERSION);
            return 0;
        } else if (c == '?') {
            fprintf(stderr, "[E::%s] unknown option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        } else if (c == ':') {
            fprintf(stderr, "[E::%s] missing option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        }
    }

    if (fp_help == stdout) {
        print_help_post(stdout);
        return 0;
    }

    if (argc - opt.ind < 3) {
        fprintf(stderr, "[E::%s] missing input: three positional options required\n", __func__);
        print_help_post(stderr);
        return 1;
    }

    if (out) {
        out1 = (char *) malloc(strlen(out) + 35);
        sprintf(out1, "%s.FINAL.agp", out);
    }

    fo = out1 == 0? stdout : fopen(out1, "w");
    if (fo == 0) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }

    fa = argv[opt.ind + 2];
    is_fai = strlen(fa) > 4 && !strcmp(fa + strlen(fa) - 4, ".fai");

    ret = 0;
    sdict = is_fai? make_sdict_from_index(fa, 0) : make_sdict_from_fa(fa, 0);
    ret = assembly_to_agp(argv[opt.ind], argv[opt.ind + 1], sdict, fo);
    fflush(fo);
    if (out != 0)
        fclose(fo);

    if (!ret && !is_fai && out1) {
        fa1 = (char *) malloc(strlen(out) + 35);
        sprintf(fa1, "%s.FINAL.fa", out);
        fprintf(stderr, "[I::%s] writing FASTA file for scaffolds\n", __func__);
        FILE *fo1;
        fo1 = fopen(fa1, "w");
        if (fo1 == NULL) {
            fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, fa1);
            exit(EXIT_FAILURE);
        }
        write_fasta_file_from_agp(fa, out1, fo1, 60);
        fclose(fo1);
    }

    sd_destroy(sdict);

    if (out1)
        free(out1);

    if (fa1)
        free(fa1);

    fprintf(stderr, "[I::%s] Version: %s\n", __func__, JUICER_VERSION);
    fprintf(stderr, "[I::%s] CMD: juicer", __func__);
    int i;
    for (i = 0; i < argc; ++i)
        fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n[I::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - jc_realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);

    return ret;
}

static int usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   juicer <command> <arguments>\n");
    fprintf(stderr, "Version: %s\n\n", JUICER_VERSION);
    fprintf(stderr, "Command: pre       generate files compatible with Juicebox toolset\n");
    fprintf(stderr, "         post      generate assembly files after Juicebox curation\n");
    fprintf(stderr, "\n");
    return 1;
}

int main(int argc, char *argv[])
{
    if (argc == 1)
        return usage();
    if (strcmp(argv[1], "pre") == 0)
        return main_pre(argc - 1, argv + 1);
    if (strcmp(argv[1], "post") == 0)
        return main_post(argc - 1, argv + 1);
    else {
        fprintf(stderr, "[E::%s] unrecognized command '%s'. Abort!\n", __func__, argv[1]);
        return 1;
    }
}
