#pragma once

#ifndef __MATH_JOY_SORT__
#define __MATH_JOY_SORT__

/*
 * https://stackoverflow.com/questions/18022192/insertion-sort-with-binary-search
 * Insert sort complexity is O(n^2) in most implementation
 * Even complexity for bianry search is O(nlogn),
 * but the insertion would actually be O(n^2), with std::vector
 * If with a std::list the binary search is not possible to achieve O(logn)
 */
/*
 * The algorithm is pretty striaght forward
 * It finds the upper bound of the to be inserted element n
 * in the range of [first, n-1] <=> [first, n), and insert element n
 * by rotating [to_be_inserted, n-1] to the back
 */
template <typename ForwardIt>
ForwardIt insertionSort(ForwardIt first, ForwardIt last) {
  ForwardIt iter = first;

  while ( iter != last ) {
    ForwardIt to_be_inserted = std::upper_bound(first, iter, *iter);
    std::rotate(to_be_inserted, iter, iter+1);
    ++iter;
  }

  // returns the first iterator
  return first;
}

#endif
