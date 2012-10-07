/** * Copyright (c) 2012 MIT License by 6.172 Staff
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 **/

// Implements the ADT specified in bitarray.h as a packed array of bits; a bit
// array containing bit_sz bits will consume roughly bit_sz/8 bytes of
// memory.


#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>

#include <sys/types.h>

#include "./bitarray.h"

typedef u_int64_t uint64_t;

// ********************************* Types **********************************

// Concrete data type representing an array of bits.
struct bitarray {
  // The number of bits represented by this bit array.
  // Need not be divisible by 8.
  size_t bit_sz;

  // The underlying memory buffer that stores the bits in
  // packed form (8 per byte).
  char *buf;
};


// ******************** Prototypes for static functions *********************

// Rotates a subarray left by an arbitrary number of bits.
//
// bit_offset is the index of the start of the subarray
// bit_length is the length of the subarray, in bits
// bit_left_amount is the number of places to rotate the
//                    subarray left
//
// The subarray spans the half-open interval
// [bit_offset, bit_offset + bit_length)
// That is, the start is inclusive, but the end is exclusive.
static void bitarray_rotate_left(bitarray_t *const bitarray,
                                 const size_t bit_offset,
                                 const size_t bit_length,
                                 const size_t bit_left_amount);

// Rotates a subarray left by one bit.
//
// bit_offset is the index of the start of the subarray
// bit_length is the length of the subarray, in bits
//
// The subarray spans the half-open interval
// [bit_offset, bit_offset + bit_length)
// That is, the start is inclusive, but the end is exclusive.
static void bitarray_rotate_left_one(bitarray_t *const bitarray,
                                     const size_t bit_offset,
                                     const size_t bit_length);

// Portable modulo operation that supports negative dividends.
//
// Many programming languages define modulo in a manner incompatible with its
// widely-accepted mathematical definition.
// http://stackoverflow.com/questions/1907565/c-python-different-behaviour-of-the-modulo-operation
// provides details; in particular, C's modulo
// operator (which the standard calls a "remainder" operator) yields a result
// signed identically to the dividend e.g., -1 % 10 yields -1.
// This is obviously unacceptable for a function which returns size_t, so we
// define our own.
//
// n is the dividend and m is the divisor
//
// Returns a positive integer r = n (mod m), in the range
// 0 <= r < m.
static size_t modulo(const ssize_t n, const size_t m);

// Produces a mask which, when ANDed with a byte, retains only the
// bit_index th byte.
//
// Example: bitmask(5) produces the byte 0b00100000.
//
// (Note that here the index is counted from right
// to left, which is different from how we represent bitarrays in the
// tests.  This function is only used by bitarray_get and bitarray_set,
// however, so as long as you always use bitarray_get and bitarray_set
// to access bits in your bitarray, this reverse representation should
// not matter.
static char bitmask(const size_t bit_index);
#define BITMASK_L(idx) (1 << ((idx) & 7))


// ******************************* Functions ********************************

//headers
void bitarray_swap_32(uint64_t* buf64,
        const size_t idx1,
        const size_t idx2);

bitarray_t *bitarray_new(const size_t bit_sz) {
  // Allocate an underlying buffer of ceil(bit_sz/8) bytes.
  char *const buf = calloc(1, bit_sz / 8 + ((bit_sz % 8 == 0) ? 0 : 1) + 8); //FIXME: be smarter about adding 8
  if (buf == NULL) {
    return NULL;
  }

  // Allocate space for the struct.
  bitarray_t *const bitarray = malloc(sizeof(struct bitarray));
  if (bitarray == NULL) {
    free(buf);
    return NULL;
  }

  bitarray->buf = buf;
  bitarray->bit_sz = bit_sz;
  return bitarray;
}

void bitarray_free(bitarray_t *const bitarray) {
  if (bitarray == NULL) {
    return;
  }
  free(bitarray->buf);
  bitarray->buf = NULL;
  free(bitarray);
}

size_t bitarray_get_bit_sz(const bitarray_t *const bitarray) {
  return bitarray->bit_sz;
}

bool bitarray_get(const bitarray_t *const bitarray, const size_t bit_index) {
  assert(bit_index < bitarray->bit_sz);

  // We're storing bits in packed form, 8 per byte.  So to get the nth
  // bit, we want to look at the (n mod 8)th bit of the (floor(n/8)th)
  // byte.
  //
  // In C, integer division is floored explicitly, so we can just do it to
  // get the byte; we then bitwise-and the byte with an appropriate mask
  // to produce either a zero byte (if the bit was 0) or a nonzero byte
  // (if it wasn't).  Finally, we convert that to a boolean.
  return (bitarray->buf[bit_index >> 3] & bitmask(bit_index)) ?
             true : false;
}

void bitarray_set(bitarray_t *const bitarray,
                  const size_t bit_index,
                  const bool value) {
  assert(bit_index < bitarray->bit_sz);

  // We're storing bits in packed form, 8 per byte.  So to set the nth
  // bit, we want to set the (n mod 8)th bit of the (floor(n/8)th) byte.
  //
  // In C, integer division is floored explicitly, so we can just do it to
  // get the byte; we then bitwise-and the byte with an appropriate mask
  // to clear out the bit we're about to set.  We bitwise-or the result
  // with a byte that has either a 1 or a 0 in the correct place.
  bitarray->buf[bit_index >> 3] =
      (bitarray->buf[bit_index >> 3] & ~bitmask(bit_index)) |
           (value ? bitmask(bit_index) : 0);
}

void bitarray_rotate(bitarray_t *const bitarray,
                     const size_t bit_offset,
                     const size_t bit_length,
                     const ssize_t bit_right_amount) {
  assert(bit_offset + bit_length <= bitarray->bit_sz);

  if (bit_length == 0) {
    return;
  }

  // Convert a rotate left or right to a left rotate only, and eliminate
  // multiple full rotations.
  bitarray_rotate_left(bitarray, bit_offset, bit_length,
           modulo(-bit_right_amount, bit_length));
}

void bitarray_reverse(bitarray_t *const bitarray,
                     const size_t bit_offset,
                     const size_t bit_length) {
	size_t idx1 = bit_offset;
	size_t idx2 = bit_offset + bit_length - 1;
	const size_t limit = bit_offset + (bit_length >> 1) - 8;
	
	size_t idx1b;
	size_t idx2b;
	
        while (idx1 < limit + 8) {
		idx1b = idx1 >> 3;
		idx2b = idx2 >> 3;
		const char bm1 = BITMASK_L(idx1);
		const char bm2 = BITMASK_L(idx2);
		char t1 = bitarray->buf[idx1b] & bm1;
		char t2 = bitarray->buf[idx2b] & bm2;
		bitarray->buf[idx1b] = (bitarray->buf[idx1b] ^ t1 ) | (t2 ? bm1 : 0);
		bitarray->buf[idx2b] = (bitarray->buf[idx2b] & ~bm2) | (t1 ? bm2 : 0);
		idx1++;
		idx2--;
	}
}

void bitarray_reverse_fast(bitarray_t *const bitarray,
                     const size_t bit_offset,
                     const size_t bit_length) {
    size_t idx1 = bit_offset;
    size_t idx2 = bit_offset + bit_length - 1;
    uint64_t* buf64 = (uint64_t*)bitarray->buf;
    if (bit_length == 0)
        return;
    while (idx2 - idx1 > 64) {
        bitarray_swap_32(buf64, idx1, idx2-31);
        idx1 += 32;
        idx2 -= 32;
    }
    bitarray_reverse(bitarray, idx1, idx2 - idx1 + 1);
}

inline uint64_t bm_32_64(size_t b_offset) {
    uint64_t start = 0xFFFFFFFF00000000;
    return start >> b_offset;
}

uint64_t reverse_64_word(uint64_t w) {
    // swap odd and even bits
    w = ((w >> 1) & 0x5555555555555555) | ((w & 0x5555555555555555) << 1);
    // swap consecutiwe pairs
    w = ((w >> 2) & 0x3333333333333333) | ((w & 0x3333333333333333) << 2);
    // swap nibbles ... 
    w = ((w >> 4) & 0x0F0F0F0F0F0F0F0F) | ((w & 0x0F0F0F0F0F0F0F0F) << 4);
    // swap bytes
    w = ((w >> 8) & 0x00FF00FF00FF00FF) | ((w & 0x00FF00FF00FF00FF) << 8);
    // swap 2-byte pairs
    w = ((w >> 16) & 0x0000FFFF0000FFFF) | ((w & 0x0000FFFF0000FFFF) << 16);
    return w;
}

void bitarray_swap_32(uint64_t* buf64,
        const size_t idx1,
        const size_t idx2) {
    size_t idx_word1 = idx1 / sizeof(uint64_t) / 8;
    size_t idx_word2 = idx2 / sizeof(uint64_t) / 8;
    size_t idx_word_offset1 = idx1 % (sizeof(uint64_t) * 8);
    size_t idx_word_offset2 = idx2 % (sizeof(uint64_t) * 8);
    //printf("%x %x\n", &buf64[idx_word1], &buf64[idx_word2]);
    uint64_t w1 = buf64[idx_word1];
    uint64_t w2 = buf64[idx_word2];
    uint64_t bm1 = bm_32_64(idx_word_offset1);
    uint64_t bm2 = bm_32_64(idx_word_offset2);
    uint64_t reversed_bits1 = reverse_64_word(((w1 & bm1) << idx_word_offset1) | ((w2 & bm2) << idx_word_offset2 >> 32));
    uint64_t extra_bits1 = w1 & ~bm1;
    uint64_t extra_bits2 = w2 & ~bm2;
    uint64_t bitsforidx1 = (reversed_bits1 << 32) >> idx_word_offset1  | extra_bits1;
    uint64_t bitsforidx2 = (reversed_bits1 & 0xFFFFFFFF00000000) >> idx_word_offset2  | extra_bits2;
    buf64[idx_word2] = bitsforidx2;
    buf64[idx_word1] = bitsforidx1;
}

inline void bitarray_swap_block(bitarray_t *const bitarray,
        size_t idx1,
        size_t idx2,
        size_t length) {
    uint64_t* buf64 = (uint64_t*)bitarray->buf;
    size_t idx_word_offset1 = idx1 % (sizeof(uint64_t) * 8);
    size_t idx_word_offset2 = idx2 % (sizeof(uint64_t) * 8);
    size_t idx_word_offset1b = (idx_word_offset1 + 32) % 64;
    size_t idx_word_offset2b = (idx_word_offset2 + 32) % 64;
    uint64_t bm1 = bm_32_64(idx_word_offset1);
    uint64_t bm2 = bm_32_64(idx_word_offset2);
    uint64_t bm1b = bm_32_64(idx_word_offset1b);
    uint64_t bm2b = bm_32_64(idx_word_offset2b);
    while (length >= 64) {
        size_t idx_word1 = idx1 / sizeof(uint64_t) / 8;
        size_t idx_word2 = idx2 / sizeof(uint64_t) / 8;
        uint64_t w1 = buf64[idx_word1];
        uint64_t w2 = buf64[idx_word2];
        uint64_t extra_bits1 = w1 & ~bm1;
        uint64_t extra_bits2 = w2 & ~bm2;
        uint64_t bitsforidx1 = w2 << idx_word_offset2 >> idx_word_offset1 | extra_bits1;
        uint64_t bitsforidx2 = w1 << idx_word_offset1 >> idx_word_offset2  | extra_bits2;
        buf64[idx_word2] = bitsforidx2;
        buf64[idx_word1] = bitsforidx1;

        length -= 32;
        idx1 += 32;
        idx2 += 32;

        idx_word1 = idx1 / sizeof(uint64_t) / 8;
        idx_word2 = idx2 / sizeof(uint64_t) / 8;
        w1 = buf64[idx_word1];
        w2 = buf64[idx_word2];
        extra_bits1 = w1 & ~bm1b;
        extra_bits2 = w2 & ~bm2b;
        bitsforidx1 = w2 << idx_word_offset2b >> idx_word_offset1b | extra_bits1;
        bitsforidx2 = w1 << idx_word_offset1b >> idx_word_offset2b  | extra_bits2;
        buf64[idx_word2] = bitsforidx2;
        buf64[idx_word1] = bitsforidx1;
        length -= 32;
        idx1 += 32;
        idx2 += 32;
    }
    while (length >= 32) {
        size_t idx_word1 = idx1 / sizeof(uint64_t) / 8;
        size_t idx_word2 = idx2 / sizeof(uint64_t) / 8;
        size_t idx_word_offset1 = idx1 % (sizeof(uint64_t) * 8);
        size_t idx_word_offset2 = idx2 % (sizeof(uint64_t) * 8);
        uint64_t w1 = buf64[idx_word1];
        uint64_t w2 = buf64[idx_word2];
        uint64_t bm1 = bm_32_64(idx_word_offset1);
        uint64_t bm2 = bm_32_64(idx_word_offset2);
        uint64_t extra_bits1 = w1 & ~bm1;
        uint64_t extra_bits2 = w2 & ~bm2;
        uint64_t bitsforidx1 = w2 << idx_word_offset2 >> idx_word_offset1 | extra_bits1;
        uint64_t bitsforidx2 = w1 << idx_word_offset1 >> idx_word_offset2  | extra_bits2;
        buf64[idx_word2] = bitsforidx2;
        buf64[idx_word1] = bitsforidx1;
        length -= 32;
        idx1 += 32;
        idx2 += 32;
    }
    //TODO: special case
    while (length > 0) {
        bool tmp = bitarray_get(bitarray, idx1);
        bitarray_set(bitarray, idx1, bitarray_get(bitarray, idx2));
        bitarray_set(bitarray, idx2, tmp);
        idx1++;
        idx2++;
        length--;
    }
}

static void bitarray_rotate_left_reverse(bitarray_t *const bitarray,
                                 const size_t bit_offset,
                                 const size_t bit_length,
                                 const size_t bit_left_amount) {
	bitarray_reverse_fast(bitarray, bit_offset, bit_left_amount);
	bitarray_reverse_fast(bitarray, bit_offset + bit_left_amount, bit_length - bit_left_amount);
	bitarray_reverse_fast(bitarray, bit_offset, bit_length);
}

static void bitarray_rotate_left(bitarray_t *const bitarray,
                                 const size_t bit_offset,
                                 const size_t bit_length,
                                 const size_t bit_left_amount) {
    //printf("rotate\n");
    if (-bit_length + 2* bit_left_amount < 32) {
        bitarray_rotate_left_reverse(bitarray, bit_offset, bit_length, bit_left_amount);
        return;
  }
    if (bit_left_amount == 0 || bit_left_amount == bit_length)
        return;
    size_t i = bit_left_amount;
    size_t j = bit_length - bit_left_amount;
    int c = 0;
    while (i != j && i > 64 && j > 64) {
        //printf("i:%d\tj:%d\n",i,j);
        if (i < j) {
            //i is shorter
            bitarray_swap_block(bitarray,
                    bit_left_amount - i + bit_offset,
                    bit_left_amount + j - i + bit_offset, i);
            j -= i;
        } else {
            //j is shorter
            bitarray_swap_block(bitarray,
                    bit_left_amount - i + bit_offset,
                    bit_left_amount + bit_offset, j);
            i -= j;

        }
        c++;
    }
    if (i < j) {
        //i is shorter
        bitarray_swap_block(bitarray,
                bit_left_amount - i + bit_offset,
                bit_left_amount + j - i + bit_offset, i);
        bitarray_rotate_left_reverse(bitarray, bit_offset, j, i);
    } else {
        //j is shorter
        bitarray_swap_block(bitarray,
                bit_left_amount - i + bit_offset,
                bit_left_amount + bit_offset, j);
        bitarray_rotate_left_reverse(bitarray, bit_offset + j, j+i, i-j);

    }
    //swap i and j
    //bitarray_swap_block((uint64_t*)bitarray->buf,
    //        bit_left_amount - i + bit_offset,
    //        bit_left_amount + bit_offset, i);

}

static void bitarray_rotate_left_one(bitarray_t *const bitarray,
                                     const size_t bit_offset,
                                     const size_t bit_length) {
  // Grab the first bit in the range, shift everything left by one, and
  // then stick the first bit at the end.
  const bool first_bit = bitarray_get(bitarray, bit_offset);
  size_t i;
  for (i = bit_offset; i + 1 < bit_offset + bit_length; i++) {
    bitarray_set(bitarray, i, bitarray_get(bitarray, i + 1));
  }
  bitarray_set(bitarray, i, first_bit);
}

static size_t modulo(const ssize_t n, const size_t m) {
  const ssize_t signed_m = (ssize_t)m;
  assert(signed_m > 0);
  const ssize_t result = ((n % signed_m) + signed_m) % signed_m;
  assert(result >= 0);
  return (size_t)result;
}

static char bitmask(const size_t bit_index) {
  return 1 << (bit_index & 7); // % 8
}

