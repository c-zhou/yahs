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

#include "bamlite.h"
#include "ketopt.h"
#include "sdict.h"
#include "link.h"
#include "asset.h"

int make_juicer_pre_file_from_bin(char *f, char *agp, char *fai, FILE *fo)
{
    sdict_t *sdict = make_sdict_from_index(fai);
    asm_dict_t *dict = agp? make_asm_dict_from_agp(sdict, agp) : make_asm_dict_from_sdict(sdict);

    FILE *fp;
    uint32_t i, i0, i1, p0, p1;
    uint32_t buffer[BUFF_SIZE], m;
    long pair_c;

    fp = fopen(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    pair_c = 0;
    while (1) {
        m = fread(&buffer, sizeof(uint32_t), BUFF_SIZE, fp);
        for (i = 0; i < m; i += 4) {
            sd_coordinate_conversion(dict, buffer[i], buffer[i + 1], &i0, &p0);
            sd_coordinate_conversion(dict, buffer[i + 2], buffer[i + 3], &i1, &p1);
            if (i0 == UINT32_MAX || i1 == UINT32_MAX) {
                fprintf(stderr, "[W::%s] sequence not found \n", __func__);
            } else {
                if (strcmp(dict->s[i0].name, dict->s[i1].name) <= 0)
                    fprintf(fo, "0\t%s\t%u\t0\t1\t%s\t%u\t1\n", dict->s[i0].name, p0, dict->s[i1].name, p1);
                else
                    fprintf(fo, "0\t%s\t%u\t1\t1\t%s\t%u\t0\n", dict->s[i1].name, p1, dict->s[i0].name, p0);
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

    return 0;
}

int make_juicer_pre_file_from_bed(char *f, char *agp, char *fai, uint8_t mq, FILE *fo)
{
    sdict_t *sdict = make_sdict_from_index(fai);
    asm_dict_t *dict = agp? make_asm_dict_from_agp(sdict, agp) : make_asm_dict_from_sdict(sdict);

    FILE *fp;
    char *line = NULL;
    size_t ln = 0;
    ssize_t read;
    char cname0[4096], cname1[4096], rname0[4096], rname1[4096];
    uint32_t s0, s1, e0, e1, i0, i1, p0, p1;
    int8_t buff;
    long pair_c;

    fp = fopen(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    pair_c = 0;
    buff = 0;
    while ((read = getline(&line, &ln, fp)) != -1) {
        if (buff == 0) {
            sscanf(line, "%s %u %u %s", cname0, &s0, &e0, rname0);
            ++buff;
        } else if (buff == 1) {
            sscanf(line, "%s %u %u %s", cname1, &s1, &e1, rname1);
            if (is_read_pair(rname0, rname1)) {
                if (++pair_c % 1000000 == 0)
                    fprintf(stderr, "[I::%s] %ld million read pairs processed \n", __func__, pair_c / 1000000);

                sd_coordinate_conversion(dict, sd_get(sdict, cname0), s0 / 2 + e0 / 2, &i0, &p0);
                sd_coordinate_conversion(dict, sd_get(sdict, cname1), s1 / 2 + e1 / 2, &i1, &p1);
                if (i0 == UINT32_MAX || i1 == UINT32_MAX) {
                    fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, i0 < 0? cname0 : cname1);
                } else {
                    if (strcmp(dict->s[i0].name, dict->s[i1].name) <= 0)
                        fprintf(fo, "0\t%s\t%u\t0\t1\t%s\t%u\t1\n", dict->s[i0].name, p0, dict->s[i1].name, p1);
                    else
                        fprintf(fo, "0\t%s\t%u\t1\t1\t%s\t%u\t0\n", dict->s[i1].name, p1, dict->s[i0].name, p0);
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

    fprintf(stderr, "[I::%s] %ld read pairs processed\n", __func__, pair_c);
    
    if (line)
        free(line);
    fclose(fp);

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
    *cname = h->target_name[b->core.tid];
    // 0x4 0x100 0x400 0x800
    if (b->core.flag & 0xD04 || b->core.qual < q) {
        *s = -1;
        *e = -1;
    } else {
        *s = b->core.pos + 1;
        *e = get_target_end(bam1_cigar(b), b->core.n_cigar, b->core.pos) + 1;
    }

    return strdup(bam1_qname(b));
}

int make_juicer_pre_file_from_bam(char *f, char *agp, char *fai, uint8_t mq, FILE *fo)
{
    sdict_t *sdict = make_sdict_from_index(fai);
    asm_dict_t *dict = agp? make_asm_dict_from_agp(sdict, agp) : make_asm_dict_from_sdict(sdict);

    bamFile fp;
    bam_header_t *h;
    bam1_t *b;
    char *cname0, *cname1, *rname0, *rname1;
    int32_t s0, s1, e0, e1;
    uint32_t i0, i1, p0, p1;
    int8_t buff;
    long pair_c;

    fp = bam_open(f, "r"); // sorted by read name
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open fail %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    h = bam_header_read(fp);
    b = bam_init1();
    cname0 = cname1 = rname0 = rname1 = 0;

    pair_c = 0;
    buff = 0;
    while (bam_read1(fp, b) >= 0 ) {
        if (buff == 0) {
            rname0 = parse_bam_rec(b, h, mq, &s0, &e0, &cname0);
            ++buff;
        } else if (buff == 1) {
            rname1 = parse_bam_rec(b, h, mq, &s1, &e1, &cname1);
            if (strcmp(rname0, rname1) == 0) {
                if (++pair_c % 1000000 == 0)
                    fprintf(stderr, "[I::%s] %ld million read pairs processed \n", __func__, pair_c / 1000000);

                if (s0 >= 0 && s1 >0) {
                    sd_coordinate_conversion(dict, sd_get(sdict, cname0), (s0 + e0) / 2, &i0, &p0);
                    sd_coordinate_conversion(dict, sd_get(sdict, cname1), (s1 + e1) / 2, &i1, &p1);

                    if (i0 == UINT32_MAX || i1 == UINT32_MAX) {
                        fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, i0 < 0? cname0 : cname1);
                    } else {
                        if (strcmp(dict->s[i0].name, dict->s[i1].name) <= 0)
                            fprintf(fo, "0\t%s\t%u\t0\t1\t%s\t%u\t1\n", dict->s[i0].name, p0, dict->s[i1].name, p1);
                        else
                            fprintf(fo, "0\t%s\t%u\t1\t1\t%s\t%u\t0\n", dict->s[i1].name, p1, dict->s[i0].name, p0);
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

    fprintf(stderr, "[I::%s] %ld read pairs processed\n", __func__, pair_c);

    if (rname0)
        free(rname0);
    if (rname1)
        free(rname1);
    bam_destroy1(b);
    bam_header_destroy(h);
    bam_close(fp);

    return 0;
}

static void print_help(FILE *fp_help)
{
    fprintf(fp_help, "Usage: juicer_pre [options] <hic.bed>|<hic.bam>|<hic.bin> <scaffolds.agp> <contigs.fa.fai>\n");
    fprintf(fp_help, "Options:\n");
    fprintf(fp_help, "    -q INT            minimum mapping quality [10]\n");
    fprintf(fp_help, "    -o STR            output to file [stdout]\n");
}

static ko_longopt_t long_options[] = {
    { "help",           ko_no_argument, 'h' },
    { 0, 0, 0 }
};

int main(int argc, char *argv[])
{
    if (argc < 2) {
        print_help(stderr);
        return 1;
    }

    FILE *fo;
    char *fai, *agp, *link_file, *out, *ext;
    int mq;

    const char *opt_str = "q:o:h";
    ketopt_t opt = KETOPT_INIT;
    int c, ret;
    FILE *fp_help = stderr;
    fai = agp = link_file = out = 0;
    mq = 10;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >= 0) {
        if (c == 'o') {
            out = opt.arg;
        } else if (c == 'q') {
            mq = atoi(opt.arg);
        } else if (c == 'h') {
            fp_help = stdout;
        } else if (c == '?') {
            fprintf(stderr, "[E::%s] unknown option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        } else if (c == ':') {
            fprintf(stderr, "[E::%s] missing option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        }
    }

    if (fp_help == stdout) {
        print_help(stdout);
        return 0;
    }

    if (argc - opt.ind < 3) {
        fprintf(stderr, "[E::%s] missing input: three positional options required\n", __func__);
        print_help(stderr);
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

    fo = out == 0? stdout : fopen(out, "w");
    if (fo == 0) {
        fprintf(stderr, "[E::%s] cannot open fail %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }
    ret = 0;
    ext = link_file + strlen(link_file) - 4;
    if (strcmp(ext, ".bam") == 0) {
        fprintf(stderr, "[I::%s] make juicer pre input from bam file %s\n", __func__, link_file);
        ret = make_juicer_pre_file_from_bam(link_file, agp, fai, mq8, fo);
    } else if (strcmp(ext, ".bed") == 0) {
        fprintf(stderr, "[I::%s] make juicer pre input from bed file %s\n", __func__, link_file);
        ret = make_juicer_pre_file_from_bed(link_file, agp, fai, mq8, fo);
    } else if (strcmp(ext, ".bin") == 0) {
        fprintf(stderr, "[I::%s] make juicer pre input from bin file %s\n", __func__, link_file);
        ret = make_juicer_pre_file_from_bin(link_file, agp, fai, fo);
    } else {
        fprintf(stderr, "[E::%s] unknown link file format. File extension .bam, .bed or .bin is expected\n", __func__);
        exit(EXIT_FAILURE);
    }

    if (out != 0)
        fclose(fo);
    
    return ret;
}
