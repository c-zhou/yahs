/*********************************************************************************
 * MIT License                                                                   *
 *                                                                               *
 * Copyright (c) 2023 Chenxi Zhou <chnx.zhou@gmail.com>                          *
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
 * 06/06/23 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/

#ifndef __AGP_SPEC_H__
#define __AGP_SPEC_H__

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>
#include <ctype.h>

/************************ AGP specification v2.1 ***********************
 * https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/#FORMAT *
 ***********************************************************************/

#define AGP_SPEC_VERSION 2.1
#define AGPSpec_hash(str) DJB_hash(str)
#define DEFAULT_AGP_U_GAP_SIZE 100
#define DEFAULT_AGP_N_GAP_SIZE 200

static const char* const AGP_SPEC_ONLINE_DOC = "https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/#FORMAT";

// DJB hash function
// ignore case
static inline uint32_t DJB_hash(const char *str)
{
    uint32_t c, hash = 5381;
    while ((c = tolower(*str++)))
        hash = ((hash << 5) + hash) + c;
    return hash;
}

typedef enum agp_component_type {
    AGP_CT_A = 177670, // Active Finishing
    AGP_CT_D = 177673, // Draft HTG (often phase1 and phase2 are called Draft, whether or not they have the draft keyword)
    AGP_CT_F = 177675, // Finished HTG (phase3)
    AGP_CT_G = 177676, // Whole Genome Finishing
    AGP_CT_O = 177684, // Other sequence (typically means no HTG keyword)
    AGP_CT_P = 177685, // Pre Draft
    AGP_CT_W = 177692, // WGS contig
    AGP_CT_N = 177683, // gap with specified size
    AGP_CT_U = 177690  // gap of unknown size, defaulting to 100 bases
} AGP_CT_t;

typedef enum agp_gap_type {
    AGP_GT_SCAFFOLD        = 1270835015, // a gap between two sequence contigs in a scaffold (superscaffold or ultra-scaffold)
    AGP_GT_CONTIG          = 4142180361, // an unspanned gap between two sequence contigs
    AGP_GT_CENTROMERE      = 848023097,  // a gap inserted for the centromere
    AGP_GT_SHORT_ARM       = 3013023508, // a gap inserted at the start of an acrocentric chromosome
    AGP_GT_HETEROCHROMATIN = 311936241,  // a gap inserted for an especially large region of heterochromatic sequence (may also include the centromere)
    AGP_GT_TELOMERE        = 3963944450, // a gap inserted for the telomere
    AGP_GT_REPEAT          = 422440038,  // an unresolvable repeat
    AGP_GT_CONTAMINATION   = 3496389561  // a gap inserted in place of foreign sequence to maintain the coordinates
} AGP_GT_t;

typedef enum agp_linkage {
    AGP_LG_YES = 193512214, // there is evidence of linkage between the adjacent lines
    AGP_LG_NO  = 5863650    // no evidence of linkage between the adjacent lines
} AGP_LG_t;

typedef enum agp_orientation {
    AGP_OT_PLUS    = 177616,     // +
    AGP_OT_MINUS   = 177618,     // -
    AGP_OT_UNKNOWN = 177636,     // ?
    AGP_OT_ZERO    = 177621,     // 0    unknown (deprecated)
    AGP_OT_ZERO_S  = 2090971781, // zero unknown (deprecated)
    AGP_OT_NA      = 5863636     // na
} AGP_OT_t;

typedef enum agp_linkage_evidence {
    AGP_LE_NA                 = 5863636   , // used when no linkage is being asserted (column 8b is ‘no’) 
    AGP_LE_PAIRED_ENDS        = 330997937 , // paired sequences from the two ends of a DNA fragment, mate-pairs and molecular-barcoding
    AGP_LE_ALIGN_GENUS        = 4129816497, // alignment to a reference genome within the same genus
    AGP_LE_ALIGN_XGENUS       = 3807300873, // alignment to a reference genome within another genus
    AGP_LE_ALIGN_TRNSCPT      = 657295357 , // alignment to a transcript from the same species
    AGP_LE_WITHIN_CLONE       = 3548268584, // sequence on both sides of the gap is derived from the same clone, but the gap is not spanned
                                            // by paired-ends. The adjacent sequence contigs have unknown order and orientation
    AGP_LE_CLONE_CONTIG       = 231014361 , // linkage is provided by a clone contig in the tiling path (TPF). For example, a gap where there
                                            // is a known clone, but there is not yet sequence for that clone
    AGP_LE_MAP                = 193499011 , // linkage asserted using a non-sequence based map such as RH, linkage, fingerprint or optical
    AGP_LE_PCR                = 193502346 , // PCR using primers on both sides of the gap
    AGP_LE_PROXIMITY_LIGATION = 3044749424, // ligation of segments of DNA that were brought into proximity in chromatin (Hi-C and related technologies
    AGP_LE_STROBE             = 479447028 , // strobe sequencing
    AGP_LE_UNSPECIFIED        = 2572506804  // used only for gaps of type contamination and when converting old AGPs that lack a field for linkage
                                            // evidence into the new format
} AGP_LE_t;

static inline int check_collision(const char *t_val, const char *p_val)
{
    // case insensitive string comparison
    const unsigned char *p1 = (const unsigned char *) t_val;
    const unsigned char *p2 = (const unsigned char *) p_val;
    int result;
    if (p1 == p2)
        return 0;
    while ((result = tolower (*p1) - tolower(*p2++)) == 0)
        if (*p1++ == '\0')
            break;
    return result;
}

static inline void AGPSpec_die(char *format, ...)
{
    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    exit(EXIT_FAILURE);
}

static inline const char *agp_component_type_val(uint32_t key)
{
    switch ((AGP_CT_t) key) {
        case AGP_CT_A: return "A";
        case AGP_CT_D: return "D";
        case AGP_CT_F: return "F";
        case AGP_CT_G: return "G";
        case AGP_CT_O: return "O";
        case AGP_CT_P: return "P";
        case AGP_CT_W: return "W";
        case AGP_CT_N: return "N";
        case AGP_CT_U: return "U";
        default:       AGPSpec_die("[E::%s] AGP component type KEY error: %u\n", __func__, key);
    }
    return 0;
}

static inline const char *agp_gap_type_val(uint32_t key)
{
    switch ((AGP_GT_t) key) {
        case AGP_GT_SCAFFOLD:        return "scaffold";
        case AGP_GT_CONTIG:          return "contig";
        case AGP_GT_CENTROMERE:      return "centromere";
        case AGP_GT_SHORT_ARM:       return "short_arm";
        case AGP_GT_HETEROCHROMATIN: return "heterochromatin";
        case AGP_GT_TELOMERE:        return "telomere";
        case AGP_GT_REPEAT:          return "repeat";
        case AGP_GT_CONTAMINATION:   return "contamination";
        default:                     AGPSpec_die("[E::%s] AGP gap type KEY error: %u\n", __func__, key);
    }
    return 0;
}

static inline const char *agp_linkage_val(uint32_t key)
{
    switch ((AGP_LG_t) key) {
        case AGP_LG_YES: return "yes";
        case AGP_LG_NO:  return "no";
        default:         AGPSpec_die("[E::%s] AGP linkage KEY error: %u\n", __func__, key);
    }
    return 0;
}

static inline const char *agp_orientation_val(uint32_t key)
{
    switch ((AGP_OT_t) key) {
        case AGP_OT_PLUS:    return "+";
        case AGP_OT_MINUS:   return "-";
        case AGP_OT_UNKNOWN:
        case AGP_OT_ZERO:
        case AGP_OT_ZERO_S:  return "?";
        case AGP_OT_NA:      return "na";
        default:             AGPSpec_die("[E::%s] AGP orientation KEY error: %u\n", __func__, key);
    }
    return 0;
}

static inline const char *agp_linkage_evidence_val(uint32_t key)
{
    switch ((AGP_LE_t) key) {
        case AGP_LE_NA:                 return "na";
        case AGP_LE_PAIRED_ENDS:        return "paired-ends";
        case AGP_LE_ALIGN_GENUS:        return "align_genus";
        case AGP_LE_ALIGN_XGENUS:       return "align_xgenus";
        case AGP_LE_ALIGN_TRNSCPT:      return "align_trnscpt";
        case AGP_LE_WITHIN_CLONE:       return "within_clone";
        case AGP_LE_CLONE_CONTIG:       return "clone_contig";
        case AGP_LE_MAP:                return "map";
        case AGP_LE_PCR:                return "pcr";
        case AGP_LE_PROXIMITY_LIGATION: return "proximity_ligation";
        case AGP_LE_STROBE:             return "strobe";
        case AGP_LE_UNSPECIFIED:        return "unspecified";
        default:                        AGPSpec_die("[E::%s] AGP linkage evidence KEY error: %u\n", __func__, key);
    }
    return 0;
}

static inline const uint32_t agp_component_type_key(const char *val)
{
    uint32_t key = AGPSpec_hash(val);
    int c = 0;
    switch ((AGP_CT_t) key) {
        case AGP_CT_A: c = check_collision(val, "A"); break;
        case AGP_CT_D: c = check_collision(val, "D"); break;
        case AGP_CT_F: c = check_collision(val, "F"); break;
        case AGP_CT_G: c = check_collision(val, "G"); break;
        case AGP_CT_O: c = check_collision(val, "O"); break;
        case AGP_CT_P: c = check_collision(val, "P"); break;
        case AGP_CT_W: c = check_collision(val, "W"); break;
        case AGP_CT_N: c = check_collision(val, "N"); break;
        case AGP_CT_U: c = check_collision(val, "U"); break;
        default:       AGPSpec_die("[E::%s] AGP component type VAL error: %s\n", __func__, val);
    }
    if (c) AGPSpec_die("[E::%s] AGP component type VAL error (hash collision): %s\n", __func__, val);
    return key;
}

static inline const uint32_t agp_gap_type_key(const char *val)
{
    uint32_t key = AGPSpec_hash(val);
    int c = 0;
    switch ((AGP_GT_t) key) {
        case AGP_GT_SCAFFOLD:        c = check_collision(val, "scaffold"); break;
        case AGP_GT_CONTIG:          c = check_collision(val, "contig"); break;
        case AGP_GT_CENTROMERE:      c = check_collision(val, "centromere"); break;
        case AGP_GT_SHORT_ARM:       c = check_collision(val, "short_arm"); break;
        case AGP_GT_HETEROCHROMATIN: c = check_collision(val, "heterochromatin"); break;
        case AGP_GT_TELOMERE:        c = check_collision(val, "telomere"); break;
        case AGP_GT_REPEAT:          c = check_collision(val, "repeat"); break;
        case AGP_GT_CONTAMINATION:   c = check_collision(val, "contamination"); break;
        default:                     AGPSpec_die("[E::%s] AGP gap type VAL error: %s\n", __func__, val);
    }
    if (c) AGPSpec_die("[E::%s] AGP gap type VAL error (hash collision): %s\n", __func__, val);
    return key;
}

static inline const uint32_t agp_linkage_key(const char *val)
{
    uint32_t key = AGPSpec_hash(val);
    int c = 0;
    switch ((AGP_LG_t) key) {
        case AGP_LG_YES: c = check_collision(val, "yes"); break;
        case AGP_LG_NO:  c = check_collision(val, "no"); break;
        default:         AGPSpec_die("[E::%s] AGP linkage VAL error: %s\n", __func__, val);
    }
    if (c) AGPSpec_die("[E::%s] AGP linkage VAL error (hash collision): %s\n", __func__, val);
    return key;
}

static inline const uint32_t agp_orientation_key(const char *val)
{
    uint32_t key = AGPSpec_hash(val);
    int c = 0;
    switch ((AGP_OT_t) key) {
        case AGP_OT_PLUS:    c = check_collision(val, "+"); break;
        case AGP_OT_MINUS:   c = check_collision(val, "-"); break;
        case AGP_OT_UNKNOWN: c = check_collision(val, "?"); break;
        case AGP_OT_ZERO:    c = check_collision(val, "0"); break;
        case AGP_OT_ZERO_S:  c = check_collision(val, "zero"); break;
        case AGP_OT_NA:      c = check_collision(val, "na"); break;
        default:             AGPSpec_die("[E::%s] AGP orientation VAL error: %s\n", __func__, val);
    }
    if (c) AGPSpec_die("[E::%s] AGP orientation VAL error (hash collision): %s\n", __func__, val);
    return key;
}

static inline const uint32_t agp_linkage_evidence_key(const char *val)
{
    uint32_t key = AGPSpec_hash(val);
    int c = 0;
    switch ((AGP_LE_t) key) {
        case AGP_LE_NA:                 c = check_collision(val, "na"); break;
        case AGP_LE_PAIRED_ENDS:        c = check_collision(val, "paired-ends"); break;
        case AGP_LE_ALIGN_GENUS:        c = check_collision(val, "align_genus"); break;
        case AGP_LE_ALIGN_XGENUS:       c = check_collision(val, "align_xgenus"); break;
        case AGP_LE_ALIGN_TRNSCPT:      c = check_collision(val, "align_trnscpt"); break;
        case AGP_LE_WITHIN_CLONE:       c = check_collision(val, "within_clone"); break;
        case AGP_LE_CLONE_CONTIG:       c = check_collision(val, "clone_contig"); break;
        case AGP_LE_MAP:                c = check_collision(val, "map"); break;
        case AGP_LE_PCR:                c = check_collision(val, "pcr"); break;
        case AGP_LE_PROXIMITY_LIGATION: c = check_collision(val, "proximity_ligation"); break;
        case AGP_LE_STROBE:             c = check_collision(val, "strobe"); break;
        case AGP_LE_UNSPECIFIED:        c = check_collision(val, "unspecified"); break;
        default:                        AGPSpec_die("[E::%s] AGP linkage evidence VAL error: %s\n", __func__, val);
    }
    if (c) AGPSpec_die("[E::%s] AGP linkage evidence VAL error (hash collision): %s\n", __func__, val);
    return key;
}

static inline AGP_OT_t agp_orientation_rev(uint32_t key)
{
    switch ((AGP_OT_t) key) {
        case AGP_OT_PLUS:    return AGP_OT_MINUS;
        case AGP_OT_MINUS:   return AGP_OT_PLUS;
        case AGP_OT_UNKNOWN:
        case AGP_OT_ZERO:
        case AGP_OT_ZERO_S:
        case AGP_OT_NA:      return key;
        default:             AGPSpec_die("[E::%s] AGP orientation KEY error: %u\n", __func__, key);
    }
    return 0;
}

#undef AGP_SPEC_TEST

#ifdef AGP_SPEC_TEST
int main(int argc, char *argv[])
{
    char *agp_component_type_vals[] = {"A", "D", "F", "G", "O", "P", "W", "N", "U"};
    char *agp_gap_type_vals[] = {"scaffold", "contig", "centromere", "short_arm", "heterochromatin", "telomere", "repeat", "contamination"};
    char *agp_linkage_vals[] = {"yes", "no"};
    char *agp_orientation_vals[] = {"+", "-", "?", "0", "zero", "na"};
    char *agp_linkage_evidence_vals[] = {"na", "paired-ends", "align_genus", "align_xgenus", "align_trnscpt", "within_clone",
                                         "clone_contig", "map", "pcr", "proximity_ligation", "strobe", "unspecified"};
    
    int i;
    uint32_t key;
    char *val;

    for (i = 0; i < sizeof(agp_component_type_vals) / sizeof(char *); ++i)
        printf("%-18s = %-10u, //\n", agp_component_type_vals[i], AGPSpec_hash(agp_component_type_vals[i]));
    for (i = 0; i < sizeof(agp_gap_type_vals) / sizeof(char *); ++i)
        printf("%-18s = %-10u, //\n", agp_gap_type_vals[i], AGPSpec_hash(agp_gap_type_vals[i]));
    for (i = 0; i < sizeof(agp_linkage_vals) / sizeof(char *); ++i)
        printf("%-18s = %-10u, //\n", agp_linkage_vals[i], AGPSpec_hash(agp_linkage_vals[i]));
    for (i = 0; i < sizeof(agp_orientation_vals) / sizeof(char *); ++i)
        printf("%-18s = %-10u, //\n", agp_orientation_vals[i], AGPSpec_hash(agp_orientation_vals[i]));
    for (i = 0; i < sizeof(agp_linkage_evidence_vals) / sizeof(char *); ++i)
        printf("%-18s = %-10u, //\n", agp_linkage_evidence_vals[i], AGPSpec_hash(agp_linkage_evidence_vals[i]));

    for (i = 0; i < sizeof(agp_component_type_vals) / sizeof(char *); ++i) {
        val = agp_component_type_vals[i];
        key = agp_component_type_key(val);
        printf("OLD_VAL: %s NEW_VAL: %s\n", val, agp_component_type_val(key));
    }
    for (i = 0; i < sizeof(agp_gap_type_vals) / sizeof(char *); ++i) {
        val = agp_gap_type_vals[i];
        key = agp_gap_type_key(val);
        printf("OLD_VAL: %s NEW_VAL: %s\n", val, agp_gap_type_val(key));
    }
    for (i = 0; i < sizeof(agp_linkage_vals) / sizeof(char *); ++i) {
        val = agp_linkage_vals[i];
        key = agp_linkage_key(val);
        printf("OLD_VAL: %s NEW_VAL: %s\n", val, agp_linkage_val(key));
    }
    for (i = 0; i < sizeof(agp_orientation_vals) / sizeof(char *); ++i) {
        val = agp_orientation_vals[i];
        key = agp_orientation_key(val);
        printf("OLD_VAL: %s NEW_VAL: %s\n", val, agp_orientation_val(key));
    }
    for (i = 0; i < sizeof(agp_linkage_evidence_vals) / sizeof(char *); ++i) {
        val = agp_linkage_evidence_vals[i];
        key = agp_linkage_evidence_key(val);
        printf("OLD_VAL: %s NEW_VAL: %s\n", val, agp_linkage_evidence_val(key));
    }

    printf("AGP_CT_A %c= 0\n", AGP_CT_A == (AGP_CT_t) 0? '=' : '!');
    printf("AGP_CT_A %c= 1\n", AGP_CT_A == (AGP_CT_t) 1? '=' : '!');

    char *upper_case_yes = "yEs";
    key = agp_linkage_key(upper_case_yes);
    printf("OLD_VAL: %s NEW_VAL: %s\n", upper_case_yes, agp_linkage_val(key));

    char *bad_identifier = "yes_";
    key = agp_linkage_key(bad_identifier);
    printf("OLD_VAL: %s NEW_VAL: %s\n", bad_identifier, agp_linkage_val(key));
    return 0;
}
#endif // AGP_SPEC_TEST

#endif // __AGP_SPEC_H__

