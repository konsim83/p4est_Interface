#ifndef __RANDOMNUMBER_H__
#define __RANDOMNUMBER_H__

// STL
#include <chrono>
#include <random>

/*
 * Generate a random double number between [a,b)
 */
class RandomNumberDouble
{
public:
  RandomNumberDouble() = delete;
  RandomNumberDouble(const double _a,
                     const double _b,
                     const bool   _same_on_all_ranks);

  void
  reinit();

  double
  generate();

private:
  double a, b;

  bool same_on_all_ranks;

  std::uniform_real_distribution<double> uniform_distribution;

  uint64_t timeSeed;

  std::seed_seq seed_sequence;

  std::mt19937_64 rng;
};


#endif // __RANDOMNUMBER_H__