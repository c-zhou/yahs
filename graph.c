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
#include <stdlib.h>
#include <assert.h>

#include "ksort.h"
#include "kdq.h"
#include "graph.h"
#include "asset.h"

#undef DEBUG
#undef DEBUG_GRAPH_PRUNE

#define graph_arc_key(a) ((a).v)
KRADIX_SORT_INIT(arc, graph_arc_t, graph_arc_key, 8)
KDQ_INIT(uint32_t)

graph_t *graph_init(void)
{
    graph_t *g;
    g = (graph_t *) calloc(1, sizeof(graph_t));
    return g;
}

void graph_destroy(graph_t *g)
{
    if (g == 0)
        return;
    if (g->arc)
        free(g->arc);
    if (g->idx)
        free(g->idx);
    free(g);
}

void graph_print(const graph_t *g, FILE *fp, int no_seq)
{
    uint32_t i, n;
    uint64_t k;
    sd_aseq_t *s;

    n = g->sdict->n;
    for (i = 0; i < n; ++i) {
        s = &g->sdict->s[i];
        fprintf(fp, "S\t%s\t", s->name);
        if (!no_seq) {
            char *seq = 0;//get_asm_seq(g->sdict, s->name);
            if (seq) {
                fputs(seq, fp);
                free(seq);
            } else {
                fputc('*', fp);
            }
        } else {
            fputc('*', fp);
        }
        fprintf(fp, "\tLN:i:%lu", s->len);
        fputc('\n', fp);
    }
    
    for (k = 0; k < g->n_arc; ++k) {
        const graph_arc_t *a = &g->arc[k];
        if (a->del || a->comp) 
            continue;
        fprintf(fp, "L\t%s\t%c\t%s\t%c\t0M\tWT:f:%.6f", g->sdict->s[a->v>>1].name, "+-"[a->v&1], g->sdict->s[a->w>>1].name, "+-"[a->w&1], a->wt);
        if (a->rank >= 0) 
            fprintf(fp, "\tSR:i:%d", a->rank);
        fputc('\n', fp);
    }
}

graph_arc_t *graph_add_arc(graph_t *g, uint32_t v, uint32_t w, int64_t link_id, int comp, double wt)
{
    graph_arc_t *a;
    if (g->m_arc == g->n_arc) {
        uint64_t old_m = g->m_arc;
        g->m_arc = g->m_arc? g->m_arc<<1 : 16;
        g->arc = (graph_arc_t *) realloc(g->arc, g->m_arc * sizeof(graph_arc_t));
        memset(&g->arc[old_m], 0, (g->m_arc - old_m) * sizeof(graph_arc_t));
    }
    a = &g->arc[g->n_arc++];
    a->v = v;
    a->w = w;
    a->wt = wt;
    a->rank = -1;
    a->link_id = link_id >= 0? link_id : g->n_arc - 1;
    if (link_id >= 0) 
        a->rank = g->arc[link_id].rank; // TODO: this is not always correct!
    a->del = a->strong = 0;
    a->comp = comp;
    return a;
}

int graph_arc_is_sorted(const graph_t *g)
{
    uint64_t e;
    for (e = 1; e < g->n_arc; ++e)
        if (g->arc[e-1].v > g->arc[e].v)
            break;
    return (e == g->n_arc);
}

void graph_arc_sort(graph_t *g)
{
    radix_sort_arc(g->arc, g->arc + g->n_arc);
}

uint64_t *graph_arc_index_core(size_t max_seq, size_t n, const graph_arc_t *a)
{
    size_t i, last;
    uint64_t *idx;
    idx = (uint64_t *) calloc(max_seq * 2, 8);
    for (i = 1, last = 0; i <= n; ++i)
        if (i == n || graph_arc_head(a[i - 1]) != graph_arc_head(a[i]))
            idx[graph_arc_head(a[i - 1])] = (uint64_t) last<<32 | (i - last), last = i;
    return idx;
}

void graph_arc_index(graph_t *g)
{
    if (g->idx) 
        free(g->idx);
    g->idx = graph_arc_index_core(g->sdict->n, g->n_arc, g->arc);
}

graph_t *read_graph_from_gfa(char *gfa)
{
    FILE *fp;
    char *line = NULL;
    size_t ln = 0;
    ssize_t read;
    char c0[1024], c1[1024], s0[4], s1[4], wts[1024];
    double wt;

    graph_t *g;
    g = graph_init();
    g->sdict = make_asm_dict_from_sdict(make_sdict_from_gfa(gfa, 0));
    
    fp = fopen(gfa, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, gfa);
        exit(EXIT_FAILURE);
    }

    while ((read = getline(&line, &ln, fp)) != -1) {
        if (line[0] == 'L') {
            sscanf(line, "%*s %s %s %s %s %*s %s", c0, s0, c1, s1, wts);
            wt = strtof(wts + 5, NULL);
            graph_add_arc(g, (asm_sd_get(g->sdict, c0)<<1) | (s0[0]=='-'), (asm_sd_get(g->sdict, c1)<<1) | (s1[0]=='-'), -1, 0, wt);
        }
    }
    fclose(fp);

    graph_arc_sort(g);
    graph_arc_index(g);

    return g;
}

void graph_clean(graph_t *g, int shear)
{
    uint64_t n_arc, n;
    graph_arc_t *arc, *a;
    arc = g->arc;
    n_arc = g->n_arc;
    n = 0;
    for(a = arc; a < arc + n_arc; ++a) {
        if (!a->del) {
            memmove(arc + n, a, sizeof(graph_arc_t));
            ++n;
        }
    }
    g->n_arc = n;
#ifdef DEBUG
    printf("[I::%s] graph cleaned: #arcs %lu -> %ld\n", __func__, n_arc, n);
#endif
    if (shear) {
        uint64_t m = 16;
#ifdef DEBUG
        uint64_t m_arc = g->m_arc;
#endif
        while (m < n) 
            m <<= 1;
        g->arc = (graph_arc_t *) realloc(g->arc, m * sizeof(graph_arc_t));
        g->m_arc = m;
#ifdef DEBUG
        printf("[I::%s] memory sheared: #arcs %lu -> %ld\n", __func__, m_arc, m);
#endif
    }

    graph_arc_sort(g);
    graph_arc_index(g);
}

void graph_print_gv(const graph_t *g, FILE *fp)
{
    uint64_t k;
    char *seq;

    fprintf(fp, "digraph {\n");
    fprintf(fp, "\t{\n");
    for (k = 0; k < g->sdict->n; ++k) {
        seq = g->sdict->s[k].name;
        fprintf(fp, "\t\t\"%s%c\" [style=filled fillcolor = salmon];\n", seq, '-');
        fprintf(fp, "\t\t\"%s%c\";\n", seq, '+');
    }
    fprintf(fp, "\t}\n");
    for (k = 0; k < g->n_arc; ++k) {
        const graph_arc_t *a = &g->arc[k];
        if (a->del) 
            continue;
        fprintf(fp, "\t\"%s%c\" -> \"%s%c\"[label=\"%.3f\"];\n", g->sdict->s[a->v>>1].name, "+-"[a->v&1], g->sdict->s[a->w>>1].name, "+-"[a->w&1], a->wt);
    }
    fprintf(fp, "}\n");
}

int graph_remove_multi_arcs(graph_t *g)
{
    uint32_t nv, na, v, n_del;
    graph_arc_t *av, *a, *b;

    nv = graph_n_vtx(g);
    nv <<= 1;
    n_del = 0;
    for (v = 0; v < nv; ++v) {
        na = graph_arc_n(g, v);
        if (na > 0) {
            av = graph_arc_a(g, v);
            for (a = av; a < av + na; ++a) {
                if (a->del) 
                    continue;
                for(b = a + 1; b < av + na; ++b) {
                    if (b->del) 
                        continue;
                    if (a->w == b->w) {
                        b->del = 1;
                        ++n_del;
                    }
                }
            }
        }
    }
#ifdef DEBUG
    printf("number multi-arcs deleted: %d\n", n_del);
#endif
    graph_clean(g, 1);
    return n_del;
}

int graph_add_symm_arcs(graph_t *g)
{
    uint32_t na, v, w, n_add;
    graph_arc_t *av, *a;

    // need to do it this as graph_add_arc could reallocate g->arc
    na = g->n_arc;
    av = g->arc;
    n_add = 0;
    for (a = av; a < av + na; ++a) {
        v = a->v;
        w = a->w;
        if (!graph_arc(g, w^1, v^1))
            ++n_add;
    }
    uint64_t old_m = g->m_arc;
    while (g->m_arc < na + n_add)
        g->m_arc <<= 1;
    g->arc = (graph_arc_t *) realloc(g->arc, g->m_arc * sizeof(graph_arc_t));
    memset(&g->arc[old_m], 0, (g->m_arc - old_m) * sizeof(graph_arc_t));
    
    av = g->arc;
    for (a = av; a < av + na; ++a) {
        v = a->v;
        w = a->w;
        if (graph_arc(g, w^1, v^1))
            continue;
        graph_add_arc(g, w^1, v^1, -1, 0, .0);
    }
#ifdef DEBUG
    printf("number symm-arcs added: %d\n", n_add);
#endif
    graph_arc_sort(g);
    graph_arc_index(g);

    return n_add;
}

void graph_print_all_clusters(graph_t *g, FILE *fp)
{
    uint32_t nv, na, v, i, c;
    uint64_t lens;
    graph_arc_t *av, *a;
    kdq_t(uint32_t) *q, *s;
    
    graph_remove_multi_arcs(g);
    graph_add_symm_arcs(g);

    nv = graph_n_vtx(g);
    q = kdq_init(uint32_t);
    s = kdq_init(uint32_t);
    uint8_t *visited = calloc(nv, sizeof(uint8_t));

    c = 0;
    for (i = 0; i < nv; ++i) {
        if (visited[i]) 
            continue;
        kdq_push(uint32_t, q, i << 1);
        kdq_push(uint32_t, s, i);
        visited[i] = 1;
        while (kdq_size(q) > 0) {
            v = *kdq_shift(uint32_t, q);
            na = graph_arc_n(g, v);
            if (na) {
                av = graph_arc_a(g, v);
                for (a = av; a < av + na; ++a) {
                    if (!visited[(a->w) >> 1]) {
                        kdq_push(uint32_t, q, a->w);
                        kdq_push(uint32_t, s, a->w >> 1);
                        visited[a->w >> 1] = 1;
                    }
                }
            }
            na = graph_arc_n(g, v^1);
            if (na) {
                av = graph_arc_a(g, v^1);
                for (a = av; a < av + na; ++a) {
                    if (!visited[(a->w^1) >> 1]) {
                        kdq_push(uint32_t, q, a->w^1);
                        kdq_push(uint32_t, s, (a->w^1) >> 1);
                        visited[(a->w^1) >> 1] = 1;
                    }
                }    
            }
        }

        lens = 0;
        for (i = 0; i < kdq_size(s); ++i)
            lens += g->sdict->s[kdq_at(s, i)].len;
        fprintf(fp, "Cluster %d [%ld seqs, %ld bp]:", ++c, kdq_size(s), lens);
        while (kdq_size(s) > 0)
            fprintf(fp, " %s", g->sdict->s[*kdq_shift(uint32_t, s)].name);
        fprintf(fp, "\n");
    }
    
    kdq_destroy(uint32_t, q);
    kdq_destroy(uint32_t, s);
    free(visited);
}

int trim_graph_simple_filter(graph_t *g, double min_wt, double min_diff_h, double min_diff_l, int min_len)
{
    uint32_t v, nv, na, n_del, n_ma;
    double mwt, smwt;
    graph_arc_t *av, *a;
    nv = graph_n_vtx(g);
    nv <<= 1;
    n_del = 0;

    for (v = 0; v < nv; ++v) {
        if (g->sdict->s[v>>1].len >= min_len)
            continue;
        na = graph_arc_n(g, v);
        av = graph_arc_a(g, v);
        for (a = av; a < av + na; ++a)
            n_del += graph_arc_del_existed(g, a);
    }
    graph_clean(g, 1);

    for (v = 0; v < nv; ++v) {
        na = graph_arc_n(g, v);
        av = graph_arc_a(g, v);
        
        mwt = smwt = 0;
        n_ma = 0;
        for (a = av; a < av + na; ++a) {
            if (a->wt > mwt) {
                smwt = mwt;
                mwt = a->wt;
                n_ma = 1;
            } else if (a->wt == mwt) {
                ++n_ma;
            } else if (a->wt > smwt) {
                smwt = a->wt;
            }
        }

        for (a = av; a < av + na; ++a) {
            if ((a->wt >= min_wt && a->wt >= mwt * min_diff_h) || 
                    (a->wt < min_wt && a->wt == mwt && mwt * min_diff_l >= smwt && n_ma == 1))
                continue;
            n_del += graph_arc_del_existed(g, a);
        }
    }

    graph_clean(g, 1);

#ifdef DEBUG_GRAPH_PRUNE
    printf("#del simple filter: %u\n", n_del);
#endif

    return n_del;
}

static int is_mwt(graph_t *g, uint32_t v, uint32_t w)
{
    int b;
    uint32_t na;
    double mwt;
    graph_arc_t *a, *av;

    na = graph_arc_n(g, v);
    if (na == 0)
        return 0;
    av = graph_arc_a(g, v);
    b = 0;
    mwt = 0;
    for (a = av; a < av + na; ++a) {
        if (a->wt > mwt) {
            mwt = a->wt;
            if (a->w == w)
                b = 1;
        }
    }
    return b;
}

int trim_graph_tips(graph_t *g)
{
    uint32_t v, nv, n_del;
    graph_arc_t *av;
    nv = graph_n_vtx(g);
    nv <<= 1;
    n_del = 0;

    for (v = 0; v < nv; ++v) {
        if (graph_arc_n(g, v^1) != 1 || graph_arc_n(g, v) > 0)
            continue;
        av = graph_arc_a(g, v^1);
        if (graph_arc_n(g, av->w^1) > 1 && !is_mwt(g, av->w^1, v))
            n_del += graph_arc_del_existed(g, av);
    }

    graph_clean(g, 1);

#ifdef DEBUG_GRAPH_PRUNE
    printf("#del tips: %u\n", n_del);
#endif

    return n_del;
}

int trim_graph_blunts(graph_t *g)
{
    uint32_t v, nv, na, n_del;
    uint8_t del;
    graph_arc_t *a, *av;
    nv = graph_n_vtx(g);
    nv <<= 1;
    n_del = 0;

    for (v = 0; v < nv; ++v) {
        if (graph_arc_n(g, v^1) > 0)
            continue;
        na = graph_arc_n(g, v);
        if (na > 1) {
            av = graph_arc_a(g, v);
            del = 1;
            for (a = av; a < av + na; ++a) {
                if (graph_arc_n(g, a->w) < 1) {
                    del = 0;
                    break;
                }
            }
            if (del)
                for (a = av; a < av + na; ++a)
                    n_del += graph_arc_del_existed(g, a);
        }
    }

    graph_clean(g, 1);

#ifdef DEBUG_GRAPH_PRUNE
    printf("#del blunts: %u\n", n_del);
#endif

    return n_del;
}

int trim_graph_repeats(graph_t *g)
{
    // TODO a more sophisticated implemetation for complex repeats
    uint32_t v, v1, v2, w1, w2, nv, n_del;
    graph_arc_t *av, *aw;
    nv = graph_n_vtx(g);
    nv <<= 1;
    n_del = 0;

    for (v = 0; v < nv; ++v) {
        if (graph_arc_n(g, v) != 2 || graph_arc_n(g, v^1) != 2)
            continue;
        av = graph_arc_a(g, v^1);
        v1 = av->w;
        v2 = (av + 1)->w;
        aw = graph_arc_a(g, v);
        w1 = aw->w;
        w2 = (aw + 1)->w;

        if ((graph_arc(g, v1, w1) && graph_arc(g, v2, w2)) ||
                (graph_arc(g, v1, w2) && graph_arc(g, v2, w1))) {
            n_del += graph_arc_del_existed(g, av);
            n_del += graph_arc_del_existed(g, av + 1);
            n_del += graph_arc_del_existed(g, aw);
            n_del += graph_arc_del_existed(g, aw + 1);
        }
    }

    graph_clean(g, 1);

#ifdef DEBUG_GRAPH_PRUNE
    printf("#del repeats: %u\n", n_del);
#endif
    
    return n_del;
}

int trim_graph_self_loops(graph_t *g)
{
    uint32_t n_del;
    graph_arc_t *a, *aw;
    n_del = 0;
    for (a = g->arc; a < g->arc + g->n_arc; ++a) {
        // might have bugs here
        if (graph_arc_n(g, a->v) > 1) 
            continue;
        aw = graph_arc(g, a->w, a->v);
        if (aw) 
            n_del += graph_arc_del_existed(g, aw);
    }

    graph_clean(g, 1);

#ifdef DEBUG_GRAPH_PRUNE
    printf("#del self loops: %u\n", n_del);
#endif

    return n_del;
}

int trim_graph_transitive_edges(graph_t *g)
{
    // TODO a more sophisticated implemetation for transitive edges > 2
    uint32_t v, nv, na, n_del;
    graph_arc_t *a1, *a2, *av;
    nv = graph_n_vtx(g);
    nv <<= 1;
    n_del = 0;

    for (v = 0; v < nv; ++v) {
        na = graph_arc_n(g, v);
        if (na < 2) 
            continue;
        av = graph_arc_a(g, v);
        for (a1 = av; a1 < av + na; ++a1) {
            for (a2 = av; a2 < av + na; ++a2) {
                if (a1 == a2) 
                    continue;
                if (graph_arc(g, a1->w, a2->w))
                    n_del += graph_arc_del_existed(g, a2);
            }
        }
    }

    graph_clean(g, 1);

#ifdef DEBUG_GRAPH_PRUNE
    printf("#del transitive edges: %u\n", n_del);
#endif

    return n_del;
}

int trim_graph_pop_bubbles(graph_t *g)
{
    // TODO a more sophisticated implemetation for complex bubbles
    uint32_t v, v1, v2, nv, na, n_del;
    graph_arc_t *av, *av1, *av2;
    nv = graph_n_vtx(g);
    nv <<= 1;
    n_del = 0;

    for (v = 0; v < nv; ++v) {
        na = graph_arc_n(g, v);
        if (na != 2)
            continue;
        av = graph_arc_a(g, v);
        v1 = av->w;
        v2 = (av + 1)->w;
        if (graph_arc_n(g, v1) == 1 && graph_arc_n(g, v2) == 1) {
            av1 = graph_arc_a(g, v1);
            av2 = graph_arc_a(g, v2);
            if (av1->w == av2->w) {
                n_del += graph_arc_del_existed(g, av);
                n_del += graph_arc_del_existed(g, av + 1);
                n_del += graph_arc_del_existed(g, av1);
                n_del += graph_arc_del_existed(g, av2);
            }
        }
    }

    graph_clean(g, 1);

#ifdef DEBUG_GRAPH_PRUNE
    printf("#del bubbles: %u\n", n_del);
#endif

    return n_del;
}

int trim_graph_pop_undirected(graph_t *g)
{
    uint32_t v, w, n_del;
    graph_arc_t *a;
    n_del = 0;
    for (a = g->arc; a < g->arc + g->n_arc; ++a) {
        if (graph_arc_n(g, a->v) != 1)
            continue;
        v = a->v;
        w = a->w;
        if (graph_arc(g, v^1, w) &&
                graph_arc_n(g, v^1) > 1 &&
                graph_arc_n(g, w) > 0) {
            ++n_del;
            graph_arc_del(g, w^1, v, 1);
            graph_arc_del(g, v^1, w, 1);
        }
    }

    graph_clean(g, 1);

#ifdef DEBUG_GRAPH_PRUNE
    printf("#solved undirected: %u\n", n_del);
#endif

    return n_del;
}

static inline int exist_strong_edge(graph_t *g, uint32_t v, double wt)
{
    uint32_t na;
    graph_arc_t *a, *av;
    na = graph_arc_n(g, v);
    if (na < 2) 
        return 0;
    av = graph_arc_a(g, v);
    for (a = av; a < av + na; ++a)
        if (a->wt > wt)
            return 1;
    return 0;
}

int trim_graph_weak_edges(graph_t *g)
{
    uint32_t n_del;
    graph_arc_t *a;
    n_del = 0;
    for (a = g->arc; a < g->arc + g->n_arc; ++a)
        if (exist_strong_edge(g, a->v, a->wt) && exist_strong_edge(g, a->w^1, a->wt))
            n_del += graph_arc_del_existed(g, a);

    graph_clean(g, 1);

#ifdef DEBUG_GRAPH_PRUNE
    printf("#del weak edges: %u\n", n_del);
#endif

    return n_del;
}

int trim_graph_ambiguous_edges(graph_t *g)
{
    uint32_t n_del;
    graph_arc_t *a;
    n_del = 0;
    for (a = g->arc; a < g->arc + g->n_arc; ++a)
        if (graph_arc_n(g, a->v) > 1 || graph_arc_n(g, a->w^1) > 1)
            n_del += graph_arc_del_existed(g, a);

    graph_clean(g, 1);

#ifdef DEBUG_GRAPH_PRUNE
    printf("#del ambiguous edges: %u\n", n_del);
#endif

    return n_del;
}

int search_graph_path(graph_t *g, asm_dict_t *dict, char *out)
{
    uint32_t i, j, r, qs, v, nv, na, s;
    uint64_t len;
    graph_arc_t *av;
    kdq_t(uint32_t) *q;

    nv = graph_n_vtx(g);
    nv <<= 1;
    q = kdq_init(uint32_t);
    uint8_t *visited = calloc(nv, sizeof(uint8_t));
    sdict_t *sd = dict->sdict;
    sd_seg_t cseg;
    uint32_t nseg;
    int pst, step, k, t;
    uint32_t ori;
    FILE *agp_out;

    if (out) {
        char *agp_out_name = (char *) malloc(strlen(out) + 5);
        sprintf(agp_out_name, "%s.agp", out);
        agp_out = fopen(agp_out_name, "w");
        free(agp_out_name);
    } else {
        agp_out = fopen("scaffolds_FINAL.agp", "w");
    }

    if (agp_out == NULL) {
        fprintf(stderr, "[E::%s] fail to open file to write\n", __func__);
        return 0;
    }


    s = 0;
    for (r = 0; r < 2; ++r) {
        // second visit for circles
        for (i = 0; i < nv; ++i) {
            if (visited[i]) 
                continue;

            na = graph_arc_n(g, i);
            if (!na || r) {
                v = i^1; // v is either a path source or a singleton
                kdq_clean(q);
                while (1) {
                    if (visited[v]) 
                        break; //circle
                    kdq_push(uint32_t, q, v);
                    visited[v] = visited[v^1] = 1;
                    na = graph_arc_n(g, v);
                    av = graph_arc_a(g, v);
                    assert(na < 2);
                    if (!na) 
                        break;
                    v = av->w;
                }

                // write files
                qs = kdq_size(q);
                if (qs) {
                    ++s;
                    len = 0;
                    t = 0;
                    for (j = 0; j < qs; ++j) {
                        v = kdq_at(q, j);
                        nseg = dict->s[v>>1].n;
                        ori = v&1;
                        if (ori) {
                            pst = dict->s[v>>1].s + nseg - 1;
                            step = -1;
                        } else {
                            pst = dict->s[v>>1].s;
                            step = 1;
                        }
                        for (k = 0; k < nseg; ++k) {
                            cseg = dict->seg[pst + step * k];
                            fprintf(agp_out, "scaffold_%u\t%lu\t%lu\t%u\tW\t%s\t%u\t%u\t%c\n", s, len + 1, len + cseg.y, ++t, sd->s[cseg.c >> 1].name, cseg.x + 1, cseg.x + cseg.y, "+-"[(cseg.c & 1) ^ ori]);
                            len += cseg.y;
                            if (k != nseg - 1 || j != qs - 1) {
                                fprintf(agp_out, "scaffold_%u\t%lu\t%lu\t%u\tN\t%d\tscaffold\tyes\tna\n", s, len + 1, len + GAP_SZ, ++t, GAP_SZ);
                                len += GAP_SZ;
                            }
                        }
                    }
                }
            }
        }
    }
    kdq_destroy(uint32_t, q);
    
    fclose(agp_out);

    return 0;
}


