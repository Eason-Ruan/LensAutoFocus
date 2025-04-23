#include "rng.h"


namespace Prl2 {

RNG::RNG() {
  state.state = PCG32_DEFAULT_STATE;
  state.inc = PCG32_DEFAULT_STREAM;
  uniformUInt32();
}

RNG::RNG(const uint64_t seed) {
  state.state = seed;
  state.inc = PCG32_DEFAULT_STREAM;
  uniformUInt32();
}

void RNG::setSeed(const uint64_t seed) {
  state.state = seed;
  state.inc = PCG32_DEFAULT_STREAM;
  uniformUInt32();
  uniformUInt32();
}

uint32_t RNG::uniformUInt32() { return pcg32_random_r(&state); }

CGL::Real RNG::uniformReal() {
  return std::min(uniformUInt32() * 0x1p-32, 1 - EPS_D);
}

}  // namespace Prl2