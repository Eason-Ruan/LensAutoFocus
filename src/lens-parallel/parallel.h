#ifndef _PRL2_PARALLEL_H
#define _PRL2_PARALLEL_H

#include <functional>

#include "../ext/ThreadPool/ThreadPool.h"

namespace Prl2 {

class Parallel {
 public:
  Parallel();

  // Execute a for loop in parallel
  // job: The functi。on to be parallelized
  // nChunks: Number of loop divisions
  // n: Number of iterations
  void parallelFor1D(const std::function<void(unsigned int)>& job,
                     unsigned int nChunks, unsigned int n);

  // Execute a nested for loop in parallel
  // job: The function to be parallelized
  // nChunks_x: Number of loop divisions in X direction
  // nChunks_y: Number of loop divisions in Y direction
  // nx: Number of iterations in X direction
  // ny: Number of iterations in Y direction
  void parallelFor2D(const std::function<void(unsigned int, unsigned int)>& job,
                     unsigned int nChunks_x, unsigned int nChunks_y,
                     unsigned int nx, unsigned int ny);

 private:
  ThreadPool pool;
};

}  // namespace Prl2

#endif