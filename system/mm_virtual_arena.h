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
 * DESCRIPTION: Simple managed virtual memory arena that reduces allocation
 *   overhead by reserving demand-paged virtual memory and resetting in O(1).
 */

#ifndef MM_VIRTUAL_ARENA_H_
#define MM_VIRTUAL_ARENA_H_

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Structure
 */
typedef struct mm_virtual_chunk {
  void* base;
  uint64_t offset;
  uint64_t capacity;
  struct mm_virtual_chunk* next;
} mm_virtual_chunk_t;

typedef struct {
  mm_virtual_chunk_t* head;
  mm_virtual_chunk_t* current;
  uint64_t chunk_capacity;
  uint64_t chunk_count;
} mm_virtual_arena_t;

/*
 * Setup
 */
mm_virtual_arena_t* mm_virtual_arena_new(
    const uint64_t chunk_capacity);
void mm_virtual_arena_clear(
    mm_virtual_arena_t* const arena);
void mm_virtual_arena_delete(
    mm_virtual_arena_t* const arena);

/*
 * Allocate
 */
void* mm_virtual_arena_allocate(
    mm_virtual_arena_t* const arena,
    const uint64_t num_bytes,
    const uint64_t align_bytes);

/*
 * Accessors
 */
uint64_t mm_virtual_arena_get_used(
    const mm_virtual_arena_t* const arena);
uint64_t mm_virtual_arena_get_chunk_count(
    const mm_virtual_arena_t* const arena);

#ifdef __cplusplus
}
#endif

#endif /* MM_VIRTUAL_ARENA_H_ */