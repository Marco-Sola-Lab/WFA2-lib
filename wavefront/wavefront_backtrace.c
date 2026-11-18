/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignment Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: WaveFront-Alignment module for backtracing alignments
 */

#include "utils/commons.h"
#include "wavefront_backtrace.h"

/*
 * Wavefront type
 */
#define BACKTRACE_TYPE_BITS                   4 // 4-bits for piggyback
#define BACKTRACE_TYPE_MASK 0x000000000000000Fl // Extract mask

#define BACKTRACE_PIGGYBACK_SET(offset,backtrace_type) \
  (( ((int64_t)(offset)) << BACKTRACE_TYPE_BITS) | backtrace_type)

#define BACKTRACE_PIGGYBACK_GET_TYPE(offset) \
  ((offset) & BACKTRACE_TYPE_MASK)
#define BACKTRACE_PIGGYBACK_GET_OFFSET(offset) \
  ((offset) >> BACKTRACE_TYPE_BITS)

typedef enum {
  backtrace_M       = 9,
  backtrace_D2_ext  = 8,
  backtrace_D2_open = 7,
  backtrace_D1_ext  = 6,
  backtrace_D1_open = 5,
  backtrace_I2_ext  = 4,
  backtrace_I2_open = 3,
  backtrace_I1_ext  = 2,
  backtrace_I1_open = 1,
} backtrace_type;

static const uint64_t matches_lut = 0x4D4D4D4D4D4D4D4Dul; // Matches LUT = "MMMMMMMM"
static const uint64_t insertions_lut = 0x4949494949494949ul; // Insertions LUT = "IIIIIIII"
static const uint64_t deletions_lut = 0x4444444444444444ul; // Deletions LUT = "DDDDDDDD"

/**
 * Retrieve the offset stored in the M wavefront of @p score at diagonal
 * @p k. Scores are used directly as wavefront indices (no modular memory).
 * The initialized range [wf_elements_init_min, wf_elements_init_max] is
 * used instead of [lo, hi], because heuristics trim [lo, hi] but the cells
 * they cut can still be read by later wavefront computations.
 *
 * @param wf_aligner The wavefront aligner.
 * @param score The score of the M wavefront.
 * @param k The diagonal.
 * @param offset Output: the offset at diagonal @p k.
 * @return True if the offset exists, false if @p score is out of range, the
 * wavefront is NULL, @p k is out of bounds, or the offset is null.
 */
static bool wavefront_backtrace_m_offset(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k,
    wf_offset_t* const offset) {
  if (score < 0) return false;
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  if (score >= wf_components->num_wavefronts) return false;
  wavefront_t* const mwavefront = wf_components->mwavefronts[score];
  if (mwavefront == NULL) return false;
  if (k < mwavefront->wf_elements_init_min ||
      mwavefront->wf_elements_init_max < k) return false;
  const wf_offset_t m_offset = mwavefront->offsets[k];
  if (m_offset == WAVEFRONT_OFFSET_NULL) return false;
  *offset = m_offset;
  return true;
}

/**
 * Add @p num operations of type @p lut_op (8 repetitions of a character) to the
 * CIGAR. If @p num is not multiple of 8, then @p num + (8 - (num % 8))
 * operations will be added to the CIGAR, but the begin_offset will be updated
 * using @p num. In order to do this, the cigar needs a padding of at least 7
 * elements.
 *
 * @param cigar The CIGAR.
 * @param lut_op The 64-bit representation of an operation, i.e., the same
 * operation repeated 8 times (e.g., 0x4D4D4D4D4D4D4D4Dul for 'M').
 * @param num The number of instances of the operation to add to the CIGAR.
 */
static void wavefront_backtrace_add_lut_to_cigar(
    cigar_t* const cigar,
    const uint64_t lut_op,
    int num) {
  char* operations = cigar->operations + cigar->begin_offset + 1;
  // Update offset first
  cigar->begin_offset -= num;
  // Blocks of 8-operations. We may add more than num operations if num is
  // not a multiple of 8. We need enough space for this.
  while (num > 0) {
    operations -= 8;
    memcpy(operations,&lut_op,sizeof(lut_op));
    num -= 8;
  }
}

/*
 * Backtrace Trace Patch Match/Mismsmatch
 */
int64_t wavefront_backtrace_misms(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->wf_elements_init_min <= k &&
      k <= mwavefront->wf_elements_init_max) {
    return BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k]+1,backtrace_M);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
void wavefront_backtrace_matches(
    wavefront_aligner_t* const wf_aligner,
    const int k,
    wf_offset_t offset,
    int num_matches,
    cigar_t* const cigar) {
  // Parameters
  char* operations = cigar->operations + cigar->begin_offset;
  // Update offset first
  cigar->begin_offset -= num_matches;
  // Blocks of 8-matches
  while (num_matches >= 8) {
    operations -= 8;
    *((uint64_t*)(operations+1)) = matches_lut;
    num_matches -= 8;
  }
  // Remaining matches
  int i;
  for (i=0;i<num_matches;++i) {
    *operations = 'M';
    --operations;
  }
}
/*
 * Backtrace Trace Patch Deletion
 */
int64_t wavefront_backtrace_del1_open(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->wf_elements_init_min <= k+1 &&
      k+1 <= mwavefront->wf_elements_init_max) {
    return BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k+1],backtrace_D1_open);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
int64_t wavefront_backtrace_del2_open(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->wf_elements_init_min <= k+1 &&
      k+1 <= mwavefront->wf_elements_init_max) {
    return BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k+1],backtrace_D2_open);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
int64_t wavefront_backtrace_del1_ext(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const d1wavefront = wf_aligner->wf_components.d1wavefronts[score];
  if (d1wavefront != NULL &&
      d1wavefront->wf_elements_init_min <= k+1 &&
      k+1 <= d1wavefront->wf_elements_init_max) {
    return BACKTRACE_PIGGYBACK_SET(d1wavefront->offsets[k+1],backtrace_D1_ext);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
int64_t wavefront_backtrace_del2_ext(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const d2wavefront = wf_aligner->wf_components.d2wavefronts[score];
  if (d2wavefront != NULL &&
      d2wavefront->wf_elements_init_min <= k+1 &&
      k+1 <= d2wavefront->wf_elements_init_max) {
    return BACKTRACE_PIGGYBACK_SET(d2wavefront->offsets[k+1],backtrace_D2_ext);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
/*
 * Backtrace Trace Patch Insertion
 */
int64_t wavefront_backtrace_ins1_open(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->wf_elements_init_min <= k-1 &&
      k-1 <= mwavefront->wf_elements_init_max) {
    return BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k-1]+1,backtrace_I1_open);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
int64_t wavefront_backtrace_ins2_open(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->wf_elements_init_min <= k-1 &&
      k-1 <= mwavefront->wf_elements_init_max) {
    return BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k-1]+1,backtrace_I2_open);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
int64_t wavefront_backtrace_ins1_ext(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const i1wavefront = wf_aligner->wf_components.i1wavefronts[score];
  if (i1wavefront != NULL &&
      i1wavefront->wf_elements_init_min <= k-1 &&
      k-1 <= i1wavefront->wf_elements_init_max) {
    return BACKTRACE_PIGGYBACK_SET(i1wavefront->offsets[k-1]+1,backtrace_I1_ext);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
int64_t wavefront_backtrace_ins2_ext(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const i2wavefront = wf_aligner->wf_components.i2wavefronts[score];
  if (i2wavefront != NULL &&
      i2wavefront->wf_elements_init_min <= k-1 &&
      k-1 <= i2wavefront->wf_elements_init_max) {
    return BACKTRACE_PIGGYBACK_SET(i2wavefront->offsets[k-1]+1,backtrace_I2_ext);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
/*
 * Backtrace wavefronts
 */
void wavefront_backtrace_linear(
    wavefront_aligner_t* const wf_aligner,
    const int alignment_score,
    const int alignment_k,
    const wf_offset_t alignment_offset) {
  // Parameters
  wavefront_sequences_t* const sequences = &wf_aligner->sequences;
  const int pattern_length = sequences->pattern_length;
  const int text_length = sequences->text_length;
  const wavefront_penalties_t* const penalties = &wf_aligner->penalties;
  const distance_metric_t distance_metric = penalties->distance_metric;
  // Prepare cigar
  cigar_t* const cigar = wf_aligner->cigar;
  cigar_clear(cigar);
  cigar->end_offset = cigar->max_operations - 1;
  cigar->begin_offset = cigar->max_operations - 2;
  cigar->operations[cigar->end_offset] = '\0';
  // Compute starting location
  int score = alignment_score;
  int k = alignment_k;
  int h = WAVEFRONT_H(alignment_k,alignment_offset);
  int v = WAVEFRONT_V(alignment_k,alignment_offset);
  wf_offset_t offset = alignment_offset;
  // Account for ending insertions/deletions
  if (v < pattern_length) {
    int i = pattern_length - v;
    while (i > 0) {cigar->operations[(cigar->begin_offset)--] = 'D'; --i;};
  }
  if (h < text_length) {
    int i = text_length - h;
    while (i > 0) {cigar->operations[(cigar->begin_offset)--] = 'I'; --i;};
  }
  // Trace the alignment back
  while (v > 0 && h > 0 && score > 0) {
    // Compute scores
    const int mismatch = score - penalties->mismatch;
    const int gap_open1 = score - penalties->gap_opening1;
    // Compute source offsets
    const int64_t misms = (distance_metric != indel) ?
        wavefront_backtrace_misms(wf_aligner,mismatch,k) :
        WAVEFRONT_OFFSET_NULL;
    const int64_t ins = wavefront_backtrace_ins1_open(wf_aligner,gap_open1,k);
    const int64_t del = wavefront_backtrace_del1_open(wf_aligner,gap_open1,k);
    const int64_t max_all = MAX(misms,MAX(ins,del));
    // Check source score
    if (max_all < 0) break; // No source
    // Traceback Matches
    const int max_offset = BACKTRACE_PIGGYBACK_GET_OFFSET(max_all);
    const int num_matches = offset - max_offset;
    wavefront_backtrace_matches(wf_aligner,k,offset,num_matches,cigar);
    offset = max_offset;
    // Update coordinates
    v = WAVEFRONT_V(k,offset);
    h = WAVEFRONT_H(k,offset);
    if (v <= 0 || h <= 0) break;
    // Traceback Operation
    const backtrace_type backtrace_type = BACKTRACE_PIGGYBACK_GET_TYPE(max_all);
    switch (backtrace_type) {
      case backtrace_M:
        score = mismatch;
        cigar->operations[(cigar->begin_offset)--] = 'X';
        --offset;
        break;
      case backtrace_I1_open:
        score = gap_open1;
        cigar->operations[(cigar->begin_offset)--] = 'I';
        --k; --offset;
        break;
      case backtrace_D1_open:
        score = gap_open1;
        cigar->operations[(cigar->begin_offset)--] = 'D';
        ++k;
        break;
      default:
        fprintf(stderr,"[WFA::Backtrace] Wrong type trace.4\n");
        exit(1);
        break;
    }
    // Update coordinates
    v = WAVEFRONT_V(k,offset);
    h = WAVEFRONT_H(k,offset);
  }
  // Account for last operations
  if (v > 0 && h > 0) {
    // Account for beginning series of matches
    const int num_matches = MIN(v,h);
    wavefront_backtrace_matches(wf_aligner,k,offset,num_matches,cigar);
    v -= num_matches;
    h -= num_matches;
  }
  // Account for beginning insertions/deletions
  while (v > 0) {cigar->operations[(cigar->begin_offset)--] = 'D'; --v;};
  while (h > 0) {cigar->operations[(cigar->begin_offset)--] = 'I'; --h;};
  // Set CIGAR
  ++(cigar->begin_offset);
  cigar->score = alignment_score;
}
void wavefront_backtrace_affine(
    wavefront_aligner_t* const wf_aligner,
    const affine2p_matrix_type component_begin,
    const affine2p_matrix_type component_end,
    const int alignment_score,
    const int alignment_k,
    const wf_offset_t alignment_offset) {
  // Parameters
  wavefront_sequences_t* const sequences = &wf_aligner->sequences;
  const int pattern_length = sequences->pattern_length;
  const int text_length = sequences->text_length;
  const wavefront_penalties_t* const penalties = &wf_aligner->penalties;
  const distance_metric_t distance_metric = penalties->distance_metric;
  // Prepare cigar
  cigar_t* const cigar = wf_aligner->cigar;
  cigar_clear(cigar);
  cigar->end_offset = cigar->max_operations - 1;
  cigar->begin_offset = cigar->max_operations - 2;
  cigar->operations[cigar->end_offset] = '\0';
  // Compute starting location
  affine2p_matrix_type matrix_type = component_end;
  int score = alignment_score;
  int k = alignment_k;
  int h = WAVEFRONT_H(alignment_k,alignment_offset);
  int v = WAVEFRONT_V(alignment_k,alignment_offset);
  wf_offset_t offset = alignment_offset;
  // Account for ending insertions/deletions
  if (component_end == affine2p_matrix_M) { // ends-free
    if (v < pattern_length) {
      int i = pattern_length - v;
      while (i > 0) {cigar->operations[(cigar->begin_offset)--] = 'D'; --i;};
    }
    if (h < text_length) {
      int i = text_length - h;
      while (i > 0) {cigar->operations[(cigar->begin_offset)--] = 'I'; --i;};
    }
  }
  // Trace the alignment back
  while (v > 0 && h > 0 && score > 0) {
    // Compute scores
    const int mismatch = score - penalties->mismatch;
    const int gap_open1 = score - penalties->gap_opening1 - penalties->gap_extension1;
    const int gap_open2 = score - penalties->gap_opening2 - penalties->gap_extension2;
    const int gap_extend1 = score - penalties->gap_extension1;
    const int gap_extend2 = score - penalties->gap_extension2;
    // Compute source offsets
    int64_t max_all;
    switch (matrix_type) {
      case affine2p_matrix_M: {
        const int64_t misms = wavefront_backtrace_misms(wf_aligner,mismatch,k);
        const int64_t ins1_open = wavefront_backtrace_ins1_open(wf_aligner,gap_open1,k);
        const int64_t ins1_ext  = wavefront_backtrace_ins1_ext(wf_aligner,gap_extend1,k);
        const int64_t max_ins1 = MAX(ins1_open,ins1_ext);
        const int64_t del1_open = wavefront_backtrace_del1_open(wf_aligner,gap_open1,k);
        const int64_t del1_ext  = wavefront_backtrace_del1_ext(wf_aligner,gap_extend1,k);
        const int64_t max_del1 = MAX(del1_open,del1_ext);
        if (distance_metric == gap_affine) {
          max_all = MAX(misms,MAX(max_ins1,max_del1));
          break;
        }
        const int64_t ins2_open = wavefront_backtrace_ins2_open(wf_aligner,gap_open2,k);
        const int64_t ins2_ext  = wavefront_backtrace_ins2_ext(wf_aligner,gap_extend2,k);
        const int64_t max_ins2 = MAX(ins2_open,ins2_ext);
        const int64_t del2_open = wavefront_backtrace_del2_open(wf_aligner,gap_open2,k);
        const int64_t del2_ext  = wavefront_backtrace_del2_ext(wf_aligner,gap_extend2,k);
        const int64_t max_del2 = MAX(del2_open,del2_ext);
        const int64_t max_ins = MAX(max_ins1,max_ins2);
        const int64_t max_del = MAX(max_del1,max_del2);
        max_all = MAX(misms,MAX(max_ins,max_del));
        break;
      }
      case affine2p_matrix_I1: {
        const int64_t ins1_open = wavefront_backtrace_ins1_open(wf_aligner,gap_open1,k);
        const int64_t ins1_ext  = wavefront_backtrace_ins1_ext(wf_aligner,gap_extend1,k);
        max_all = MAX(ins1_open,ins1_ext);
        break;
      }
      case affine2p_matrix_I2: {
        const int64_t ins2_open = wavefront_backtrace_ins2_open(wf_aligner,gap_open2,k);
        const int64_t ins2_ext  = wavefront_backtrace_ins2_ext(wf_aligner,gap_extend2,k);
        max_all = MAX(ins2_open,ins2_ext);
        break;
      }
      case affine2p_matrix_D1: {
        const int64_t del1_open = wavefront_backtrace_del1_open(wf_aligner,gap_open1,k);
        const int64_t del1_ext  = wavefront_backtrace_del1_ext(wf_aligner,gap_extend1,k);
        max_all = MAX(del1_open,del1_ext);
        break;
      }
      case affine2p_matrix_D2: {
        const int64_t del2_open = wavefront_backtrace_del2_open(wf_aligner,gap_open2,k);
        const int64_t del2_ext  = wavefront_backtrace_del2_ext(wf_aligner,gap_extend2,k);
        max_all = MAX(del2_open,del2_ext);
        break;
      }
      default:
        fprintf(stderr,"[WFA::Backtrace] Wrong type trace.1\n");
        exit(1);
        break;
    }
    // Check source score
    if (max_all < 0) break; // No source
    // Traceback matches
    if (matrix_type == affine2p_matrix_M) {
      const int max_offset = BACKTRACE_PIGGYBACK_GET_OFFSET(max_all);
      const int num_matches = offset - max_offset;
      wavefront_backtrace_matches(wf_aligner,k,offset,num_matches,cigar);
      offset = max_offset;
      // Update coordinates
      v = WAVEFRONT_V(k,offset);
      h = WAVEFRONT_H(k,offset);
      if (v <= 0 || h <= 0) break;
    }
    // Traceback operation
    const backtrace_type backtrace_type = BACKTRACE_PIGGYBACK_GET_TYPE(max_all);
    switch (backtrace_type) {
      case backtrace_M:
        score = mismatch;
        matrix_type = affine2p_matrix_M;
        break;
      case backtrace_I1_open:
        score = gap_open1;
        matrix_type = affine2p_matrix_M;
        break;
      case backtrace_I1_ext:
        score = gap_extend1;
        matrix_type = affine2p_matrix_I1;
        break;
      case backtrace_I2_open:
        score = gap_open2;
        matrix_type = affine2p_matrix_M;
        break;
      case backtrace_I2_ext:
        score = gap_extend2;
        matrix_type = affine2p_matrix_I2;
        break;
      case backtrace_D1_open:
        score = gap_open1;
        matrix_type = affine2p_matrix_M;
        break;
      case backtrace_D1_ext:
        score = gap_extend1;
        matrix_type = affine2p_matrix_D1;
        break;
      case backtrace_D2_open:
        score = gap_open2;
        matrix_type = affine2p_matrix_M;
        break;
      case backtrace_D2_ext:
        score = gap_extend2;
        matrix_type = affine2p_matrix_D2;
        break;
      default:
        fprintf(stderr,"[WFA::Backtrace] Wrong type trace.2\n");
        exit(1);
        break;
    }
    switch (backtrace_type) {
      case backtrace_M:
        cigar->operations[(cigar->begin_offset)--] = 'X';
        --offset;
        break;
      case backtrace_I1_open:
      case backtrace_I1_ext:
      case backtrace_I2_open:
      case backtrace_I2_ext:
        cigar->operations[(cigar->begin_offset)--] = 'I';
        --k; --offset;
        break;
      case backtrace_D1_open:
      case backtrace_D1_ext:
      case backtrace_D2_open:
      case backtrace_D2_ext:
        cigar->operations[(cigar->begin_offset)--] = 'D';
        ++k;
        break;
      default:
        fprintf(stderr,"[WFA::Backtrace] Wrong type trace.3\n");
        exit(1);
        break;
    }
    // Update coordinates
    v = WAVEFRONT_V(k,offset);
    h = WAVEFRONT_H(k,offset);
  }
  // Account for last operations
  if (matrix_type == affine2p_matrix_M) {
    if (v > 0 && h > 0) {
      // Account for beginning series of matches
      const int num_matches = MIN(v,h);
      wavefront_backtrace_matches(wf_aligner,k,offset,num_matches,cigar);
      v -= num_matches;
      h -= num_matches;
    }
    // Account for beginning insertions/deletions
    while (v > 0) {cigar->operations[(cigar->begin_offset)--] = 'D'; --v;};
    while (h > 0) {cigar->operations[(cigar->begin_offset)--] = 'I'; --h;};
  } else {
    // DEBUG
    if (v != 0 || h != 0 || (score != 0 && penalties->match == 0)) {
      fprintf(stderr,"[WFA::Backtrace] I?/D?-Beginning backtrace error\n");
      fprintf(stderr,">%.*s\n",pattern_length,sequences->pattern);
      fprintf(stderr,"<%.*s\n",text_length,sequences->text);
      exit(-1);
    }
  }
  // Set CIGAR
  ++(cigar->begin_offset);
  cigar->score = alignment_score;
}
/*
 * Backtrace (Singletrack: M-wavefronts only)
 */
/**
 * Count the consecutive matching characters that end at cell (@p v, @p h),
 * walking backwards along the diagonal. The count is bounded by
 * MIN(@p v, @p h) so it never reads before the start of either sequence.
 * Lambda inputs are compared one position at a time with the custom match
 * function.
 *
 * @param wf_aligner The wavefront aligner.
 * @param v The pattern position (vertical coordinate) of the cell.
 * @param h The text position (horizontal coordinate) of the cell.
 * @return The number of consecutive matches ending at (@p v, @p h).
 */
static int wavefront_backtrace_singletrack_count_backward_matches(
    wavefront_aligner_t* const wf_aligner,
    const int v,
    const int h) {
  wavefront_sequences_t* const sequences = &wf_aligner->sequences;
  const int max_matches = MIN(v,h);
  int nmatches = 0;
  // Custom match function (no sequence buffers)
  if (sequences->mode == wf_sequences_lambda) {
    while (nmatches < max_matches &&
           wavefront_sequences_cmp(sequences,v-nmatches-1,h-nmatches-1)) {
      ++nmatches;
    }
    return nmatches;
  }
  const char* const pattern = sequences->pattern;
  const char* const text = sequences->text;
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
  // Blocked backwards extend, bounded to avoid reading before free-start edges.
  while (nmatches + 8 <= max_matches) {
    uint64_t pattern_block, text_block;
    memcpy(&pattern_block,pattern + v - nmatches - 8,sizeof(pattern_block));
    memcpy(&text_block,text + h - nmatches - 8,sizeof(text_block));
    const uint64_t cmp = pattern_block ^ text_block;
    if (cmp != 0) {
      const int equal_left_bits = __builtin_clzll(cmp);
      nmatches += DIV_FLOOR(equal_left_bits,8);
      return nmatches;
    }
    nmatches += 8;
  }
#endif
  while (nmatches < max_matches &&
         pattern[v-nmatches-1] == text[h-nmatches-1]) {
    ++nmatches;
  }
  return nmatches;
}

/**
 * Look for an M cell at score @p gap_score from which an insertion
 * (diagonal *@p k - @p l) or a deletion (diagonal *@p k + @p l) of length
 * @p l, followed by matches, reaches an offset in [@p offset_orig, *@p offset].
 * On success, the matches and the gap are prepended to the CIGAR, and @p k
 * and @p offset are moved to that M cell.
 *
 * @param wf_aligner The wavefront aligner.
 * @param cigar The CIGAR.
 * @param gap_score The score of the M cell where the gap opens.
 * @param l The length of the gap.
 * @param offset_orig The lowest offset reachable from *@p offset by
 * backward matches.
 * @param k Input/output: the diagonal of the current cell.
 * @param offset Input/output: the offset of the current cell.
 * @return True if the gap was found.
 */
static bool wavefront_backtrace_singletrack_gap(
    wavefront_aligner_t* const wf_aligner,
    cigar_t* const cigar,
    const int gap_score,
    const int l,
    const wf_offset_t offset_orig,
    int* const k,
    wf_offset_t* const offset) {
  wf_offset_t m_offset;
  // Insertion
  if (wavefront_backtrace_m_offset(wf_aligner,gap_score,*k-l,&m_offset) &&
      m_offset + l >= offset_orig && m_offset + l <= *offset) {
    wavefront_backtrace_add_lut_to_cigar(cigar,matches_lut,*offset-(m_offset+l));
    wavefront_backtrace_add_lut_to_cigar(cigar,insertions_lut,l);
    *k -= l;
    *offset = m_offset;
    return true;
  }
  // Deletion
  if (wavefront_backtrace_m_offset(wf_aligner,gap_score,*k+l,&m_offset) &&
      m_offset >= offset_orig && m_offset <= *offset) {
    wavefront_backtrace_add_lut_to_cigar(cigar,matches_lut,*offset-m_offset);
    wavefront_backtrace_add_lut_to_cigar(cigar,deletions_lut,l);
    *k += l;
    *offset = m_offset;
    return true;
  }
  return false;
}

/**
 * Retrieve the cigar of the alignment for gap-affine and dual gap-affine
 * using exclusively the information in the M matrix (M wavefronts). This is
 * slightly slower than the general traceback, but it enables storing only
 * the M matrix when performing the alignment (I1, I2, D1 and D2 only need a
 * small scope).
 *
 * @param wf_aligner The wavefront aligner.
 * @param component_begin unused.
 * @param component_end The component (matrix) where the traceback starts.
 * @param alignment_score The score of the alignment.
 * @param alignment_k The diagonal that contains the cell (N, M), where the
 * traceback starts.
 * @param alignment_offset The offset of the cell (N, M), where the traceback
 * starts.
 */
void wavefront_backtrace_affine_singletrack(
    wavefront_aligner_t* const wf_aligner,
    const affine2p_matrix_type component_begin,
    const affine2p_matrix_type component_end,
    const int alignment_score,
    const int alignment_k,
    const wf_offset_t alignment_offset) {

  // Parameters
  wavefront_sequences_t* const sequences = &wf_aligner->sequences;
  const int pattern_length = sequences->pattern_length;
  const int text_length = sequences->text_length;
  const wavefront_penalties_t* const penalties = &wf_aligner->penalties;
  const distance_metric_t distance_metric = penalties->distance_metric;

  // Prepare cigar
  // WARNING: We want padding in the cigar so we can always use LUTs.
  cigar_t* const cigar = wf_aligner->cigar;
  cigar_clear(cigar);
  cigar->end_offset = cigar->max_operations - 1;
  cigar->begin_offset = cigar->max_operations - 2;
  cigar->operations[cigar->end_offset] = '\0';

  // Compute starting location
  int score = alignment_score;
  int k = alignment_k;
  int h = WAVEFRONT_H(alignment_k,alignment_offset);
  int v = WAVEFRONT_V(alignment_k,alignment_offset);
  wf_offset_t offset = alignment_offset;
  wf_offset_t offset_orig = offset;

  // Account for ending insertions/deletions
  if (component_end == affine2p_matrix_M) { // ends-free
    if (v < pattern_length) {
      int i = pattern_length - v;
      while (i > 0) {cigar->operations[(cigar->begin_offset)--] = 'D'; --i;};
    }
    if (h < text_length) {
      int i = text_length - h;
      while (i > 0) {cigar->operations[(cigar->begin_offset)--] = 'I'; --i;};
    }
  }

  bool in_mmatrix = component_end == affine2p_matrix_M;
  int l = 0; // Length of the current chain of insertions or deletions.

  // Trace the alignment back
  while (score != 0) {
    if (in_mmatrix) {
      v = WAVEFRONT_V(k,offset);
      h = WAVEFRONT_H(k,offset);
      const int nmatches =
          wavefront_backtrace_singletrack_count_backward_matches(wf_aligner,v,h);
      offset_orig = offset - nmatches;

      // If we come from a mismatch, then the backwards extend is correct.
      const int mismatch = score - penalties->mismatch;
      wf_offset_t m_offset;
      if (wavefront_backtrace_m_offset(wf_aligner,mismatch,k,&m_offset) &&
          m_offset + 1 == offset_orig) {
        wavefront_backtrace_add_lut_to_cigar(cigar,matches_lut,nmatches);
        cigar->operations[(cigar->begin_offset)--] = 'X';
        offset = offset_orig - 1;
        score = mismatch;
      } else {
        // Otherwise, we come from a gap. M only stores the furthest offset
        // after extension, so the gap may end anywhere in
        // [offset_orig, offset].
        in_mmatrix = false;
        l = 0;
      }
    } else {
      // Search gaps of increasing length until one opens from an M cell.
      ++l;
      const int indel1 = score - penalties->gap_opening1 -
                         l * penalties->gap_extension1;
      const int indel2 = (distance_metric == gap_affine_2p) ?
          score - penalties->gap_opening2 - l * penalties->gap_extension2 : -1;
      if (wavefront_backtrace_singletrack_gap(
            wf_aligner,cigar,indel1,l,offset_orig,&k,&offset)) {
        score = indel1;
        in_mmatrix = true;
      } else if (wavefront_backtrace_singletrack_gap(
                  wf_aligner,cigar,indel2,l,offset_orig,&k,&offset)) {
        score = indel2;
        in_mmatrix = true;
      } else if (indel1 < 0 && indel2 < 0) {
        // No source. Only valid at an ends-free start on the first row/column
        in_mmatrix = (WAVEFRONT_V(k,offset_orig) == 0 ||
                      WAVEFRONT_H(k,offset_orig) == 0);
        break;
      }
    }
  }

  // Account for last operations
  v = WAVEFRONT_V(k,offset);
  h = WAVEFRONT_H(k,offset);
  const int num_matches = MIN(v,h);
  v -= num_matches;
  h -= num_matches;
  // Leading gaps are only allowed within the begin-free limits (0 if end-to-end)
  if (!in_mmatrix ||
      v > wf_aligner->alignment_form.pattern_begin_free ||
      h > wf_aligner->alignment_form.text_begin_free) {
    fprintf(stderr,"[WFA::Backtrace] Singletrack backtrace error "
                   "(could not trace back to a valid alignment start)\n");
    if (sequences->mode != wf_sequences_lambda) {
      fprintf(stderr,">%.*s\n",pattern_length,sequences->pattern);
      fprintf(stderr,"<%.*s\n",text_length,sequences->text);
    }
    exit(-1);
  }
  // Account for beginning series of matches
  wavefront_backtrace_add_lut_to_cigar(cigar,matches_lut,num_matches);
  // Account for beginning insertions/deletions
  wavefront_backtrace_add_lut_to_cigar(cigar,deletions_lut,v);
  wavefront_backtrace_add_lut_to_cigar(cigar,insertions_lut,h);

  // Set CIGAR
  ++(cigar->begin_offset);
  cigar->score = alignment_score;
}

/*
 * Backtrace from BT-Buffer (pcigar)
 */
void wavefront_backtrace_pcigar(
    wavefront_aligner_t* const wf_aligner,
    const int alignment_k,
    const int alignment_offset,
    const pcigar_t pcigar_last,
    const bt_block_idx_t prev_idx_last) {
  // Parameters
  wf_backtrace_buffer_t* const bt_buffer =  wf_aligner->wf_components.bt_buffer;
  // Traceback pcigar-blocks
  bt_block_t bt_block_last = {
      .pcigar = pcigar_last,
      .prev_idx = prev_idx_last
  };
  bt_block_t* const init_block = wf_backtrace_buffer_traceback_pcigar(bt_buffer,&bt_block_last);
  // Fetch initial coordinate
  const int init_position_offset = init_block->pcigar;
  wf_backtrace_init_pos_t* const backtrace_init_pos =
      vector_get_elm(bt_buffer->alignment_init_pos,init_position_offset,wf_backtrace_init_pos_t);
  // Unpack pcigar blocks (packed alignment)
  const int begin_v = backtrace_init_pos->v;
  const int begin_h = backtrace_init_pos->h;
  const int end_v = WAVEFRONT_V(alignment_k,alignment_offset);
  const int end_h = WAVEFRONT_H(alignment_k,alignment_offset);
  if (wf_aligner->penalties.distance_metric <= gap_linear) {
    wf_backtrace_buffer_unpack_cigar_linear(
        bt_buffer,&wf_aligner->sequences,
        begin_v,begin_h,end_v,end_h,wf_aligner->cigar);
  } else {
    wf_backtrace_buffer_unpack_cigar_affine(
        bt_buffer,&wf_aligner->sequences,
        begin_v,begin_h,end_v,end_h,wf_aligner->cigar);
  }
}
