#ifndef RANDOM_H
#define RANDOM_H
#include "rng.h"
#include "sampler.h"
#include "type.h"
#include "vector2D.h"

namespace Prl2 {

// Sampler that simply returns random numbers
class RandomSampler final : public Sampler {
 public:
  RandomSampler(){};
  explicit RandomSampler(const uint64_t seed) : rng(RNG(seed)){};

  void setSeed(const uint64_t seed) override { rng.setSeed(seed); };

  CGL::Real getNext() override { return rng.uniformReal(); }
  CGL::Vector2D getNext2D() override {
    return CGL::Vector2D(rng.uniformReal(), rng.uniformReal());
  };

  std::unique_ptr<Sampler> clone(const uint64_t seed) override {
    return std::unique_ptr<Sampler>(new RandomSampler(seed));
  };

 private:
  RNG rng;
};

}  // namespace Prl2

#endif
