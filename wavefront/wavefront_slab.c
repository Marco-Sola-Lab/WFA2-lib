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
 * DESCRIPTION: WaveFront Slab for fast pre-allocated wavefronts' memory handling
 */

#include "utils/commons.h"
#include "system/mm_allocator.h"
#include "wavefront_slab.h"

/*
 * Constants
 */
#define WF_SLAB_ARENA_CHUNK_SIZE (256ULL * 1024 * 1024)

/*
 * Setup
 */
wavefront_slab_t* wavefront_slab_new(
    const int init_wf_length,
    const bool allocate_backtrace,
    const wf_slab_mode_t slab_mode,
    mm_allocator_t* const mm_allocator) {
  // Allocate
  wavefront_slab_t* const wavefront_slab =
      mm_allocator_alloc(mm_allocator,wavefront_slab_t);
  // Attributes
  wavefront_slab->allocate_backtrace = allocate_backtrace;
  wavefront_slab->slab_mode = slab_mode;
  // Wavefront Slabs
  wavefront_slab->init_wf_length = init_wf_length;
  wavefront_slab->current_wf_length = init_wf_length;
  wavefront_slab->wavefronts = NULL;
  wavefront_slab->wavefronts_free = NULL;
  // Stats
  wavefront_slab->memory_used = 0;
  // MM
  wavefront_slab->mm_allocator = mm_allocator;
  wavefront_slab->arena = mm_virtual_arena_new(WF_SLAB_ARENA_CHUNK_SIZE);
  // Return
  return wavefront_slab;
}
void wavefront_slab_reap_free(
    wavefront_slab_t* const wavefront_slab) {
  mm_virtual_arena_clear(wavefront_slab->arena);
}
void wavefront_slab_reap_repurpose(
    wavefront_slab_t* const wavefront_slab) {
  mm_virtual_arena_clear(wavefront_slab->arena);
}
void wavefront_slab_reap(
    wavefront_slab_t* const wavefront_slab) {
  wavefront_slab->current_wf_length = wavefront_slab->init_wf_length;
  mm_virtual_arena_clear(wavefront_slab->arena);
}
void wavefront_slab_clear(
    wavefront_slab_t* const wavefront_slab) {
  mm_virtual_arena_clear(wavefront_slab->arena);
}
void wavefront_slab_delete(
    wavefront_slab_t* const wavefront_slab) {
  mm_virtual_arena_delete(wavefront_slab->arena);
  mm_allocator_free(wavefront_slab->mm_allocator,wavefront_slab);
}
/*
 * Accessors
 */
void wavefront_slab_set_mode(
    wavefront_slab_t* const wavefront_slab,
    const wf_slab_mode_t slab_mode) {
  wavefront_slab->slab_mode = slab_mode;
}
/*
 * Slab Allocator
 */
wavefront_t* wavefront_slab_allocate(
    wavefront_slab_t* const wavefront_slab,
    const int min_lo,
    const int max_hi) {
  const int wf_elements_allocated = WAVEFRONT_LENGTH(min_lo,max_hi);
  mm_virtual_arena_t* const arena = wavefront_slab->arena;

  // Allocate wavefront handler
  wavefront_t* const wavefront = (wavefront_t*) mm_virtual_arena_allocate(
      arena,sizeof(wavefront_t),8);

  // Allocate diagonal offsets (64-byte aligned for AVX2)
  wavefront->wf_elements_allocated = wf_elements_allocated;
  wavefront->offsets_mem = (wf_offset_t*) mm_virtual_arena_allocate(
      arena,(size_t)wf_elements_allocated*sizeof(wf_offset_t),64);

  // Allocate backtrace memory
  if (wavefront_slab->allocate_backtrace) {
    wavefront->bt_pcigar_mem = (pcigar_t*) mm_virtual_arena_allocate(
        arena,(size_t)wf_elements_allocated*sizeof(pcigar_t),64);
    wavefront->bt_prev_mem = (bt_block_idx_t*) mm_virtual_arena_allocate(
        arena,(size_t)wf_elements_allocated*sizeof(bt_block_idx_t),64);
  } else {
    wavefront->bt_pcigar_mem = NULL;
    wavefront->bt_prev_mem = NULL;
  }

  // Init wavefront
  wavefront->status = wavefront_status_busy;
  wavefront_init(wavefront,min_lo,max_hi);
  // Return
  return wavefront;
}
void wavefront_slab_free(
    wavefront_slab_t* const wavefront_slab,
    wavefront_t* const wavefront) {
  (void)wavefront_slab;
  wavefront->status = wavefront_status_free;
}
/*
 * Utils
 */
uint64_t wavefront_slab_get_size(
    wavefront_slab_t* const wavefront_slab) {
  return mm_virtual_arena_get_used(wavefront_slab->arena);
}