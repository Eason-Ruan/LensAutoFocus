#ifndef SAMPLER_H
#define SAMPLER_H

#include <cstdint>
#include <memory>

namespace Prl2 {

// Class for generating random numbers
class Sampler {
 public:
  Sampler() = default;
  virtual ~Sampler() = default;

  // Set the seed value
  virtual void setSeed(uint64_t seed) = 0;

  // Get a random number for the next dimension
  virtual Real getNext() = 0;

  // Get two random numbers for the next dimension
  virtual CGL::Vector2D getNext2D() = 0;

  // Clone the sampler
  virtual std::unique_ptr<Sampler> clone(uint64_t seed) = 0;
};

}  // namespace Prl2

#endif
