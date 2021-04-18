#include <RandomNumber.h>

RandomNumberDouble::RandomNumberDouble(const double _a,
                                       const double _b,
                                       const bool   _same_on_all_ranks)
  : a(_a)
  , b(_b)
  , same_on_all_ranks(_same_on_all_ranks)
  , uniform_distribution(a, b)
  , timeSeed(
      (same_on_all_ranks ?
         0.0 :
         std::chrono::high_resolution_clock::now().time_since_epoch().count()))
  , seed_sequence{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)}
{
  rng.seed(seed_sequence);
}


void
RandomNumberDouble::reinit()
{
  // re-initialize the random number generator with time-dependent seed
  timeSeed =
    (same_on_all_ranks ?
       0.0 :
       std::chrono::high_resolution_clock::now().time_since_epoch().count());
  std::seed_seq seed_sequence{uint32_t(timeSeed & 0xffffffff),
                              uint32_t(timeSeed >> 32)};
  rng.seed(seed_sequence);
}

double
RandomNumberDouble::generate()
{
  return uniform_distribution(rng);
}