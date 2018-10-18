/* Copyright (C) 2018 New Joy - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the GPLv3
 *
 *
 * You should have received a copy of the GPLv3 license with
 * this file. If not, please visit https://www.gnu.org/licenses/gpl-3.0.en.html
 *
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Author: yangqp5@outlook.com (New Joy)
 *
 */

#ifndef _MATHJOY_SIEVE_HPP_
#define _MATHJOY_SIEVE_HPP_

#include <algorithm>

namespace mathjoy {

  template <typename RandomAccessIterator, typename Integer>
  void mark_sieve(RandomAccessIterator first, RandomAccessIterator last,
                  Integer factor) {
    *first = false;
    while (last - first > factor) {
      first = first + factor;
      *first = false;
    }
  }

  // template<typename Integer> inline
  // Integer index_square(Integer n) {
  //   return 2*n*n + 6*n + 3;
  // }

  template <typename RandomAccessIterator, typename Integer>
  void sift0(RandomAccessIterator first, Integer n) {
    std::fill(first, first + n, true);
    Integer i(0);
    Integer index_square(3);

    while (index_square < n) {
      // invariant: index_square 2i^2 + 6i + 3
      if (first[i]) {
        mark_sieve(first + index_square,
                   first + n,
                   i + i + 3); // factor
      }
      ++i;
      //index_square = 2 * i * i + 6 * i + 3;
      index_square = 2 * i * (i + 3) + 3;
    }

  }

  template <typename RandomAccessIterator, typename Integer>
  void sift1(RandomAccessIterator first, Integer n) {
    Integer i(0);
    Integer index_square(3);
    RandomAccessIterator last(first + n);
    Integer factor(3);
    std::fill(first, last, true);

    while (index_square < n) {
      if (first[i]) {
        mark_sieve(first, last, factor);
      }
      ++i;
      index_square = 2 * i * (i + 3) + 3;
      factor = i + i + 3;
    }
  }

  template <typename RandomAccessIterator, typename Integer>
  void sift(RandomAccessIterator first, Integer n) {
    Integer i(0);
    Integer index_square(3);
    RandomAccessIterator last(first + n);
    Integer factor(3);
    std::fill(first, last, true);

    while (index_square < n) {
      if (first[i]) {
        mark_sieve(first, last, factor);
      }
      ++i;
      index_square += factor;
      factor += Integer(2); // using constructor here is actually a great detail
      index_square += factor;
    }
  }
}

#endif /* _MATHJOY_SIEVE_HPP_ */
