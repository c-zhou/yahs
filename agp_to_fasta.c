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

#include "ketopt.h"
#include "sdict.h"

static void print_help(FILE *fp_help)
{
    fprintf(fp_help, "Usage: agp_to_fasta [options] <scaffolds.agp> <contigs.fa>\n");
    fprintf(fp_help, "Options:\n");
    fprintf(fp_help, "    -l INT            line width [60]\n");
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
    char *fa, *agp, *out;
    int line_wd;

    const char *opt_str = "o:l:h";
    ketopt_t opt = KETOPT_INIT;
    int c;
    FILE *fp_help = stderr;
    fa = agp = out = 0;
    line_wd = 60;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >= 0) {
        if (c == 'l') {
            line_wd = atoi(opt.arg);
        } else if (c == 'o') {
            out = opt.arg;
        } else if (c == 'h') {
            fp_help = stdout;
        }else if (c == '?') {
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

    if (argc - opt.ind < 2) {
        fprintf(stderr, "[E::%s] missing input: three positional options required\n", __func__);
        print_help(stderr);
        return 1;
    }

    agp = argv[opt.ind];
    fa = argv[opt.ind + 1];

    fo = out == 0? stdout : fopen(out, "w");
    if (fo == 0) {
        fprintf(stderr, "[E::%s] cannot open fail %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }
    
    write_fasta_file_from_agp(fa, agp, fo, line_wd);

    if (out != 0)
        fclose(fo);
    
    return 0;
}
