// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_DATA_STRUCTURES_SHORT_ALLOC_H
#define IRL_DATA_STRUCTURES_SHORT_ALLOC_H

// short_alloc implementation from
// https://howardhinnant.github.io/short_alloc.h

// The MIT License (MIT)
//
// Copyright (c) 2015 Howard Hinnant
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// Example usage from Howard Hinnant website.
// Create a vector<T> template with a small buffer of 200 bytes.
//   Note for vector it is possible to reduce the alignment requirements
//   down to alignof(T) because vector doesn't allocate anything but T's.
//   And if we're wrong about that guess, it is a compile-time error, not
//   a run time error.
// template <class T, std::size_t BufSize = 200>
// using SmallVector = std::vector<T, short_alloc<T, BufSize, alignof(T)>>;
//
// int main()
//{
//    // Create the stack-based arena from which to allocate
//    SmallVector<int>::allocator_type::arena_type a;
//    // Create the vector which uses that arena.
//    SmallVector<int> v{a};
//    // Exercise the vector and note that new/delete are not getting called.
//    v.push_back(1);
//    memuse();
//    v.push_back(2);
//    memuse();
//    v.push_back(3);
//    memuse();
//    v.push_back(4);
//    memuse();
//    // Yes, the correct values are actually in the vector
//    for (auto i : v)
//        std::cout << i << ' ';
//    std::cout << '\n';
//}
//
// memory = 0
// alloc = 0
// memory = 0
// alloc = 0
// memory = 0
// alloc = 0
// memory = 0
// alloc = 0
// 1 2 3 4

#include <cstddef>

#include <cassert>

namespace short_alloc {

template <std::size_t N, std::size_t alignment = alignof(std::max_align_t)>
class arena {
  alignas(alignment) char buf_[N];
  char *ptr_;

public:
  ~arena() { ptr_ = nullptr; }
  arena() noexcept : ptr_(buf_) {}
  arena(const arena &) = delete;
  arena &operator=(const arena &) = delete;

  template <std::size_t ReqAlign> char *allocate(std::size_t n);
  void deallocate(char *p, std::size_t n) noexcept;

  static constexpr std::size_t size() noexcept { return N; }
  std::size_t used() const noexcept {
    return static_cast<std::size_t>(ptr_ - buf_);
  }
  void reset() noexcept { ptr_ = buf_; }

private:
  static std::size_t align_up(std::size_t n) noexcept {
    return (n + (alignment - 1)) & ~(alignment - 1);
  }

  bool pointer_in_buffer(char *p) noexcept {
    return buf_ <= p && p <= buf_ + N;
  }
};

template <std::size_t N, std::size_t alignment>
template <std::size_t ReqAlign>
char *arena<N, alignment>::allocate(std::size_t n) {
  static_assert(ReqAlign <= alignment, "alignment is too small for this arena");
  assert(pointer_in_buffer(ptr_) && "short_alloc has outlived arena");
  auto const aligned_n = align_up(n);
  if (static_cast<decltype(aligned_n)>(buf_ + N - ptr_) >= aligned_n) {
    char *r = ptr_;
    ptr_ += aligned_n;
    return r;
  }

  static_assert(alignment <= alignof(std::max_align_t),
                "you've chosen an "
                "alignment that is larger than alignof(std::max_align_t), and "
                "cannot be guaranteed by normal operator new");
  return static_cast<char *>(::operator new(n));
}

template <std::size_t N, std::size_t alignment>
void arena<N, alignment>::deallocate(char *p, std::size_t n) noexcept {
  assert(pointer_in_buffer(ptr_) && "short_alloc has outlived arena");
  if (pointer_in_buffer(p)) {
    n = align_up(n);
    if (p + n == ptr_)
      ptr_ = p;
  } else {
    ::operator delete(p);
  }
}

template <class T, std::size_t N, std::size_t Align = alignof(std::max_align_t)>
class short_alloc {
public:
  using value_type = T;
  static auto constexpr alignment = Align;
  static auto constexpr size = N;
  using arena_type = arena<size, alignment>;

private:
  arena_type &a_;

public:
  short_alloc(const short_alloc &) = default;
  short_alloc &operator=(const short_alloc &) = delete;

  // Purposefully not explicit
  short_alloc(arena_type &a) noexcept : a_(a) {
    static_assert(size % alignment == 0,
                  "size N needs to be a multiple of alignment Align");
  }
  template <class U>
  explicit short_alloc(const short_alloc<U, N, alignment> &a) noexcept
      : a_(a.a_) {}

  template <class _Up> struct rebind {
    using other = short_alloc<_Up, N, alignment>;
  };

  T *allocate(std::size_t n) {
    return reinterpret_cast<T *>(
        a_.template allocate<alignof(T)>(n * sizeof(T)));
  }
  void deallocate(T *p, std::size_t n) noexcept {
    a_.deallocate(reinterpret_cast<char *>(p), n * sizeof(T));
  }

  template <class T1, std::size_t N1, std::size_t A1, class U, std::size_t M,
            std::size_t A2>
  friend bool operator==(const short_alloc<T1, N1, A1> &x,
                         const short_alloc<U, M, A2> &y) noexcept;

  template <class U, std::size_t M, std::size_t A> friend class short_alloc;
};

template <class T, std::size_t N, std::size_t A1, class U, std::size_t M,
          std::size_t A2>
inline bool operator==(const short_alloc<T, N, A1> &x,
                       const short_alloc<U, M, A2> &y) noexcept {
  return N == M && A1 == A2 && &x.a_ == &y.a_;
}

template <class T, std::size_t N, std::size_t A1, class U, std::size_t M,
          std::size_t A2>
inline bool operator!=(const short_alloc<T, N, A1> &x,
                       const short_alloc<U, M, A2> &y) noexcept {
  return !(x == y);
}

} // namespace short_alloc

#endif // SRC_DATA_STRUCTURES_SHORT_ALLOC_H
