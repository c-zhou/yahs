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
#ifndef ASSET_H_
#define ASSET_H_

#include <stdint.h>
#include <stdio.h>
#include <zlib.h>

#define SWAP(T, x, y) {T tmp = x; x = y; y = tmp;}
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MIN_MAX(x, min, max) MIN((MAX((x), (min))), (max))
#define BUFF_SIZE 4096
#define BIN_H 0x5941485342494E56
#define BIN_V 0x2
#define strcasecmp(s1, s2) strcmp_case_insensitive(s1, s2)

typedef struct {
    int fd;
    gzFile fp;
    void *stream;
    void *buffer;
    void *koaux;
    int64_t nline;
} iostream_t;

#ifdef __cplusplus
extern "C" {
#endif

double realtime(void);
double cputime(void);
void liftrlimit();
long peakrss(void);
void ram_limit(long *total, long *avail);
int file_copy(char *fin, char *fout);
int8_t is_read_pair(const char *rname0, const char *rname1);
uint32_t div_ceil(uint64_t x, uint32_t y);
uint64_t linear_scale(uint64_t g, int *scale, uint64_t max_g);
void write_bin_header(FILE *fo);
int is_valid_bin_header(int64_t magic_number);
int strcmp_case_insensitive(const char *s1, const char *s2);
void positive_or_die(int num);
int is_empty_line(char *line);
iostream_t *iostream_open(const char *spath);
void iostream_close(iostream_t *iostream);
char *iostream_getline(iostream_t *iostream);
#ifdef __cplusplus
}
#endif

#endif /* ASSET_H_ */

