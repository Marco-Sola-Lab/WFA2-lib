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

#include "system/mm_virtual_arena.h"

#include <stdlib.h>

#if defined(__linux__) && (defined(__x86_64__) || defined(__aarch64__))
#include <sys/mman.h>
#include <unistd.h>

/*
 * Platform Primitives
 */
static inline uintptr_t mm_virtual_align_address(
    const uintptr_t address,
    const uint64_t alignment) {
  const uintptr_t mask = (uintptr_t)(alignment - 1);
  return (address + mask) & ~mask;
}

static inline uint64_t mm_virtual_page_size(void) {
  const long sz = sysconf(_SC_PAGESIZE);
  return (sz > 0) ? (uint64_t)sz : 4096;
}

static inline void* mm_virtual_os_reserve(
    const uint64_t bytes) {
  if (bytes == 0) return NULL;
  void* ptr = mmap(NULL, bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
  return (ptr == MAP_FAILED) ? NULL : ptr;
}

static inline void mm_virtual_os_release(
    void* const base,
    const uint64_t bytes) {
  if (base == NULL || bytes == 0) return;
  munmap(base, bytes);
}

/*
 * Chunks
 */
static mm_virtual_chunk_t* mm_virtual_chunk_new(
    const uint64_t capacity) {
  // Compute page-aligned capacity
  const uint64_t page_sz = mm_virtual_page_size();
  const uint64_t aligned_capacity = (uint64_t)mm_virtual_align_address((uintptr_t)capacity, page_sz);

  // Reserve virtual address space
  void* const mapped_base = mm_virtual_os_reserve(aligned_capacity);
  if (mapped_base == NULL) return NULL;

  // Allocate handler
  mm_virtual_chunk_t* const chunk = (mm_virtual_chunk_t*) malloc(sizeof(mm_virtual_chunk_t));
  if (chunk == NULL) {
    mm_virtual_os_release(mapped_base, aligned_capacity);
    return NULL;
  }

  // Setup
  chunk->base = mapped_base;
  chunk->offset = 0;
  chunk->capacity = aligned_capacity;
  chunk->next = NULL;
  return chunk;
}

static void mm_virtual_chunk_delete(
    mm_virtual_chunk_t* const chunk) {
  if (chunk == NULL) return;
  if (chunk->base != NULL) {
    mm_virtual_os_release(chunk->base, chunk->capacity);
  }
  free(chunk);
}

/*
 * Setup & Lifecycle
 */
mm_virtual_arena_t* mm_virtual_arena_new(
    const uint64_t chunk_capacity) {
  if (chunk_capacity == 0) return NULL;

  // Allocate initial chunk
  mm_virtual_chunk_t* const initial_chunk = mm_virtual_chunk_new(chunk_capacity);
  if (initial_chunk == NULL) return NULL;

  // Allocate handler
  mm_virtual_arena_t* const arena = (mm_virtual_arena_t*) malloc(sizeof(mm_virtual_arena_t));
  if (arena == NULL) {
    mm_virtual_chunk_delete(initial_chunk);
    return NULL;
  }

  // Setup
  arena->head = initial_chunk;
  arena->current = initial_chunk;
  arena->chunk_capacity = initial_chunk->capacity;
  arena->chunk_count = 1;
  return arena;
}

void mm_virtual_arena_clear(
    mm_virtual_arena_t* const arena) {
  if (arena == NULL) return;
  // Reset offsets across all chained chunks
  mm_virtual_chunk_t* curr = arena->head;
  while (curr != NULL) {
    curr->offset = 0;
    curr = curr->next;
  }
  arena->current = arena->head;
}

void mm_virtual_arena_reap(
    mm_virtual_arena_t* const arena) {
  if (arena == NULL) return;
  // Retain primary chunk and release secondary chained chunks
  mm_virtual_chunk_t* curr = arena->head->next;
  while (curr != NULL) {
    mm_virtual_chunk_t* const next = curr->next;
    mm_virtual_chunk_delete(curr);
    curr = next;
  }
  arena->head->offset = 0;
  arena->head->next = NULL;
  arena->current = arena->head;
  arena->chunk_count = 1;
}

void mm_virtual_arena_delete(
    mm_virtual_arena_t* const arena) {
  if (arena == NULL) return;
  // Free all chunks
  mm_virtual_chunk_t* curr = arena->head;
  while (curr != NULL) {
    mm_virtual_chunk_t* const next = curr->next;
    mm_virtual_chunk_delete(curr);
    curr = next;
  }
  free(arena);
}

/*
 * Allocate
 */
void* mm_virtual_arena_allocate(
    mm_virtual_arena_t* const arena,
    const uint64_t num_bytes,
    const uint64_t align_bytes) {
  if (arena == NULL || num_bytes == 0) return NULL;
  const uint64_t effective_align = (align_bytes == 0) ? sizeof(void*) : align_bytes;
  mm_virtual_chunk_t* chunk = arena->current;

  // Compute aligned address
  const uintptr_t current_addr = (uintptr_t)chunk->base + chunk->offset;
  const uintptr_t aligned_addr = mm_virtual_align_address(current_addr, effective_align);
  const uint64_t padding = (uint64_t)(aligned_addr - current_addr);
  const uint64_t total_required = padding + num_bytes;

  // Check available chunk size
  if (chunk->offset + total_required <= chunk->capacity) {
    chunk->offset += total_required;
    return (void*)aligned_addr;
  }

  // Reuse chained chunk if available
  if (chunk->next != NULL) {
    arena->current = chunk->next;
    return mm_virtual_arena_allocate(arena, num_bytes, align_bytes);
  }

  // Allocate new chained chunk
  const uint64_t new_capacity = (num_bytes > arena->chunk_capacity) ? num_bytes : arena->chunk_capacity;
  mm_virtual_chunk_t* const new_chunk = mm_virtual_chunk_new(new_capacity);
  if (new_chunk == NULL) return NULL;

  chunk->next = new_chunk;
  arena->current = new_chunk;
  ++(arena->chunk_count);
  return mm_virtual_arena_allocate(arena, num_bytes, align_bytes);
}

/*
 * Accessors
 */
uint64_t mm_virtual_arena_get_used(
    const mm_virtual_arena_t* const arena) {
  if (arena == NULL) return 0;
  uint64_t total = 0;
  mm_virtual_chunk_t* curr = arena->head;
  while (curr != NULL) {
    total += curr->offset;
    curr = curr->next;
  }
  return total;
}

uint64_t mm_virtual_arena_get_chunk_count(
    const mm_virtual_arena_t* const arena) {
  return (arena != NULL) ? arena->chunk_count : 0;
}

#else /* !defined(__linux__) || (!defined(__x86_64__) && !defined(__aarch64__)) */

/*
 * Fallback stubs for unsupported platforms (macOS, Windows, 32-bit).
 * Ensures libwfa.a links cleanly without undefined symbol errors while
 * forcing wavefront_slab to use the classic slab allocator.
 */
mm_virtual_arena_t* mm_virtual_arena_new(
    const uint64_t chunk_capacity) {
  (void)chunk_capacity;
  return NULL;
}

void mm_virtual_arena_clear(
    mm_virtual_arena_t* const arena) {
  (void)arena;
}

void mm_virtual_arena_reap(
    mm_virtual_arena_t* const arena) {
  (void)arena;
}

void mm_virtual_arena_delete(
    mm_virtual_arena_t* const arena) {
  (void)arena;
}

void* mm_virtual_arena_allocate(
    mm_virtual_arena_t* const arena,
    const uint64_t num_bytes,
    const uint64_t align_bytes) {
  (void)arena;
  (void)num_bytes;
  (void)align_bytes;
  return NULL;
}

uint64_t mm_virtual_arena_get_used(
    const mm_virtual_arena_t* const arena) {
  (void)arena;
  return 0;
}

uint64_t mm_virtual_arena_get_chunk_count(
    const mm_virtual_arena_t* const arena) {
  (void)arena;
  return 0;
}

#endif