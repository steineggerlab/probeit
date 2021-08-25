#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <assert.h>
#include <x86intrin.h>

#define SIMD_SIZE 256
#define ALIGNED_SIZE \
  (SIMD_SIZE / CHAR_SIZE)  // For 128-bit SIMD, ALIGNED_SIZE is 16 bytes
#define CHAR_SIZE 8

// Note: We provide a constructor that can pass char* bits and we just copy this
// pointer to bits_, there is no deep copy because we donot wanna too much
// memory copy. If bitset involves copy between each other, you should control
// it yourself. Due to the same reason, the AND OR operator also return char*
// rather than Bitset.
class Bitset {
 public:
  // Support a given char*, but must specify the size
  Bitset(char* bits, int num_bit, int num_byte)
      : bits_(bits), num_bits_(num_bit), num_bytes_(num_byte) {}

  // Bitset is using char* to store data
  Bitset(int size) { Resize(size); }

  // Support null bitset.
  Bitset() : bits_(nullptr), num_bits_(0), num_bytes_(0) {}

  // TODO: free is needed, but since we did not wanna deep copy to involve
  // performance issue in constructor and assigner, just comment free here.
  ~Bitset() {
    if (bits_ != nullptr) {
      // free(bits_);
    }
  }

  // Use resize to allocate bitset
  void Resize(int size) {

    num_bits_ = size;

    // Note: aligned_alloc requires an integral multiple of alignment for
    // allocated size
    int integral = num_bits_ / SIMD_SIZE;
    int remainder = num_bits_ % SIMD_SIZE;

    // If remainder = 0, that means the current size of bit is the integral
    // of alignment. The allocated size is just the integral multiple
    if (remainder == 0) {
      num_bytes_ = integral * ALIGNED_SIZE;
    }
    // If remainder is not 0, that means the size exceed 1 ALIGNED_SIZE
    else {
      num_bytes_ = (integral + 1) * ALIGNED_SIZE;
    }

    // bits must be aligned. For 128-bit SIMD, it is 16
    posix_memalign((void **) &bits_, ALIGNED_SIZE, num_bytes_);

    // Init the memory with zero
    memset(bits_, 0, num_bytes_);
  }

  // For exampel, bits are (7, 10, 13, 45, 100)
  // and bitset are (000000000000000). This func set bitset with
  // (00000100100100001)
  void Set(std::vector<int>& bits) {
    // Ensure the bitset has already inited
    assert(num_bits_ > 0);
    // Iterate the input region
    for (auto& bit : bits) {
      int idx = bit / CHAR_SIZE;
      int offbit = bit % CHAR_SIZE;
      char offset = 0x01;
      offset = offset << offbit;
      bits_[idx] |= offset;
    }
  }

  // Note: the bit starts from 0
  void Set(int bit) {
    // Ensure the bitset has already inited
    assert(num_bits_ > 0);
    int idx = bit / CHAR_SIZE;
    int offbit = bit % CHAR_SIZE;
    char offset = 0x01;
    offset = offset << offbit;
    bits_[idx] |= offset;
  }
  
  bool isSet(int bit){
    assert(num_bits_ > 0);
    int idx = bit / CHAR_SIZE;
    int offbit = bit % CHAR_SIZE;
    return ((bits_[idx]) & (1<<(offbit)));
  } 
  // Set all bits with 0
  void Clear() { memset(bits_, 0, num_bytes_); }

  // TODO: Return Bitset where its value to value copying. But for char* it only
  // copies address, not deep copy. So in destructor, we do not free it. We
  // should fix this in the future
  // Note: If use this API pls change SIMD_SIZE to 128
  void AND(Bitset& rh_bitset) {
    // The SIMD vector
    __m256i A, B, C;

    // The char* of right hand bitset
    char* rbits = rh_bitset.Get();

    // How many chars a SIMD vector contains
    // For example, 256-bit SIMD, and char is 8, so it processes 16 bytes each
    // time (span = 32)
    int span = SIMD_SIZE / CHAR_SIZE;

    for (int i = 0; i < num_bytes_; i += span) {
      A = _mm256_load_si256((__m256i*)&bits_[i]);
      B = _mm256_load_si256((__m256i*)&rbits[i]);
      C = _mm256_and_si256(A, B);
      _mm256_store_si256((__m256i*)&bits_[i], C);
    }

  }

  // TODO: Return Bitset where its value to value copying. But for char* it only
  // copies address, not deep copy. So in destructor, we do not free it. We
  // should fix this in the future
  // Note: If use this API pls change SIMD_SIZE to 128
  void Or(Bitset& rh_bitset) {
    // The SIMD vector
    __m256i A, B, C;

    // The char* of right hand bitset
    char* rbits = rh_bitset.Get();

    // How many chars a SIMD vector contains
    // For example, 256-bit SIMD, and char is 8, so it processes 16 bytes each
    // time (span = 32)
    int span = SIMD_SIZE / CHAR_SIZE;

    for (int i = 0; i < num_bytes_; i += span) {
      A = _mm256_load_si256((__m256i*)&bits_[i]);
      B = _mm256_load_si256((__m256i*)&rbits[i]);
      C = _mm256_or_si256(A, B);
      _mm256_store_si256((__m256i*)&bits_[i], C);
    }
  }

  // Return the number of bits (which are equal to 1)
  int Count() {
    // Return count
    int count = 0;

    int count_test = 0;

    // The SIMD vector
    __m256i A;

    // How many chars a SIMD vector contains, 256/8=32
    int span = SIMD_SIZE / CHAR_SIZE;

    for (int i = 0; i < num_bytes_; i += span) {
      A = _mm256_load_si256((__m256i*)&bits_[i]);
      count += popcnt256(A);
    }

    return count;
  }


  int AndNotCountDontStore(Bitset& rh_bitset) {
    // Return count
    int count = 0;

    // The SIMD vector
    __m256i A, B, C;

    // The char* of right hand bitset
    char* rbits = rh_bitset.Get();

    // How many chars a SIMD vector contains
    // For example, 128-bit SIMD, and char is 8, so it processes 16 bytes each
    // time (span = 16)
    int span = SIMD_SIZE / CHAR_SIZE;

    for (int i = 0; i < num_bytes_; i += span) {
      A = _mm256_load_si256((__m256i*)&bits_[i]);
      B = _mm256_load_si256((__m256i*)&rbits[i]);
      C = _mm256_andnot_si256(B, A);
      count += popcnt256(C);
    }

    return count;
  }

  int AndNotCount(Bitset& rh_bitset) {
    // Return count
    int count = 0;

    // The SIMD vector
    __m256i A, B, C;

    // The char* of right hand bitset
    char* rbits = rh_bitset.Get();

    // How many chars a SIMD vector contains
    // For example, 128-bit SIMD, and char is 8, so it processes 16 bytes each
    // time (span = 16)
    int span = SIMD_SIZE / CHAR_SIZE;

    for (int i = 0; i < num_bytes_; i += span) {
      A = _mm256_load_si256((__m256i*)&bits_[i]);
      B = _mm256_load_si256((__m256i*)&rbits[i]);
      C = _mm256_andnot_si256(B, A);
      _mm256_store_si256((__m256i*)&bits_[i], C);
      count += popcnt256(C);
    }

    return count;
  }

  // Return the count for AND operation
  int AndCount(Bitset& rh_bitset) {
    // Return count
    int count = 0;

    // The SIMD vector
    __m256i A, B, C;

    // The char* of right hand bitset
    char* rbits = rh_bitset.Get();

    // How many chars a SIMD vector contains
    // For example, 128-bit SIMD, and char is 8, so it processes 16 bytes each
    // time (span = 16)
    int span = SIMD_SIZE / CHAR_SIZE;

    for (int i = 0; i < num_bytes_; i += span) {
      A = _mm256_load_si256((__m256i*)&bits_[i]);
      B = _mm256_load_si256((__m256i*)&rbits[i]);
      C = _mm256_and_si256(A, B);
      _mm256_store_si256((__m256i*)&bits_[i], C);
      count += popcnt256(C);
    }

    return count;
  }

  char* Get() { return bits_; }
  int Size() { return num_bits_; }

 private:
  inline int popcnt128(__m128i n) {
    const __m128i n_hi = _mm_unpackhi_epi64(n, n);
#ifdef _MSC_VER
    return __popcnt64(_mm_cvtsi128_si64(n)) +
           __popcnt64(_mm_cvtsi128_si64(n_hi));
#else
    return __popcntq(_mm_cvtsi128_si64(n)) + __popcntq(_mm_cvtsi128_si64(n_hi));
#endif
  }

  // Support 256
  inline int popcnt256(__m256i n) {
    uint64_t* u = (uint64_t*)&n;
    return _mm_popcnt_u64(u[0]) + _mm_popcnt_u64(u[1]) + _mm_popcnt_u64(u[2]) +
           _mm_popcnt_u64(u[3]);
  }

 private:
  // Bitset value
  char* bits_;

  // number of bits of this bitset
  int num_bits_;

  // number of bytes of this bitset
  // size_ = num_bits_/CHAR_SIZE + 1;
  int num_bytes_;
};

