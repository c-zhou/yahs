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
 * 11/07/21 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#ifndef GRAPH_H_
#define GRAPH_H_

#include <stdio.h>
#include <stdint.h>

#include "sdict.h"

typedef struct {
    uint32_t v, w; // vetex_id | ori
    int32_t rank;
    double wt;
    uint64_t link_id:61, strong:1, del:1, comp:1; // link_id: a pair of dual arcs are supposed to have the same link_id
} graph_arc_t;

#define graph_arc_head(a) ((a).v)
#define graph_arc_tail(a) ((a).w)
#define graph_arc_n(g, v) ((uint32_t)(g)->idx[(v)])
#define graph_arc_a(g, v) (&(g)->arc[(g)->idx[(v)]>>32])
#define graph_n_vtx(g) ((g)->sdict->n)

typedef struct {
    // segments
    uint32_t max_rank;
    asm_dict_t *sdict;
    // links
    uint64_t m_arc, n_arc;
    graph_arc_t *arc;
    uint64_t *idx;
} graph_t;


#ifdef __cplusplus
extern "C" {
#endif

graph_t *graph_init(void);
void graph_destroy(graph_t *g);
void graph_print(const graph_t *g, FILE *fp, int no_seq);
graph_arc_t *graph_add_arc(graph_t *g, uint32_t v, uint32_t w, int64_t link_id, int comp, double wt);
int graph_arc_is_sorted(const graph_t *g);
void graph_arc_sort(graph_t *g);
void graph_arc_index(graph_t *g);
graph_t *read_graph_from_gfa(char *gfa);
void graph_clean(graph_t *g, int shear);
void graph_print_gv(const graph_t *g, FILE *fp);
graph_t *graph_print_gv_around_node(const graph_t *g, FILE *fp, uint32_t *v, int radius);
void graph_print_all_clusters(graph_t *g, FILE *fp);
int graph_remove_multi_arcs(graph_t *g);
int graph_add_symm_arcs(graph_t *g);
int trim_graph_simple_filter(graph_t *g, double min_wt, double min_diff_h, double min_diff_l, int min_len);
int trim_graph_tips(graph_t *g);
int trim_graph_blunts(graph_t *g);
int trim_graph_repeats(graph_t *g);
int trim_graph_self_loops(graph_t *g);
int trim_graph_transitive_edges(graph_t *g);
int trim_graph_pop_bubbles(graph_t *g);
int trim_graph_pop_undirected(graph_t *g);
int trim_graph_weak_edges(graph_t *g);
int trim_graph_ambiguous_edges(graph_t *g);
asm_dict_t *make_asm_dict_from_graph(graph_t *g, asm_dict_t *dict);

#ifdef __cplusplus
}
#endif

static inline void graph_arc_del(graph_t *g, uint32_t v, uint32_t w, int del)
{
    uint32_t i, nv = graph_arc_n(g, v);
    graph_arc_t *av = graph_arc_a(g, v);
    for (i = 0; i < nv; ++i)
        if (av[i].w == w) 
            av[i].del = !!del;
}

static inline graph_arc_t *graph_arc(graph_t *g, uint32_t v, uint32_t w)
{
    uint32_t i, nv = graph_arc_n(g, v);
    graph_arc_t *av = graph_arc_a(g, v);
    for (i = 0; i < nv; ++i)
        if (av[i].w == w) 
            return av + i;
    return 0;
}

static inline void graph_vtx_del(graph_t *g, uint32_t s)
{
    uint32_t k;
    for (k = 0; k < 2; ++k) {
        uint32_t i, v = s<<1 | k;
        uint32_t nv = graph_arc_n(g, v);
        graph_arc_t *av = graph_arc_a(g, v);
        for (i = 0; i < nv; ++i) {
            av[i].del = 1;
            graph_arc_del(g, av[i].w^1, v^1, 1);
        }
    }
}

static inline int graph_arc_del_existed(graph_t *g, graph_arc_t *a)
{
    if (a->del) 
        return 0;
    a->del = 1;
    graph_arc_del(g, a->w^1, a->v^1, 1);
    return 1;
}

#endif /* GRAPH_H_ */

