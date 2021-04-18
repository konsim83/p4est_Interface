#ifndef __UTILITIES_H__
#define __UTILITIES_H__

// Deal.ii
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

// STL
#include <string>

// My headers
#include <RandomNumber.h>

///////////////////////////////////////////

template <int dim>
void
fill_points_randomly(std::vector<dealii::Point<dim>> &points,
                     double                           a,
                     double                           b,
                     bool                             same_on_all_ranks);

// Extern template instantiations
extern template void
  fill_points_randomly<2>(std::vector<dealii::Point<2>> &points,
                          double                         a,
                          double                         b,
                          bool                           same_on_all_ranks);

extern template void
  fill_points_randomly<3>(std::vector<dealii::Point<3>> &points,
                          double                         a,
                          double                         b,
                          bool                           same_on_all_ranks);

///////////////////////////////////////////

template <int dim>
void
write_mesh(
  const dealii::parallel::distributed::Triangulation<dim> &triangulation,
  const std::string &                                      filename);

// Extern template instantiations
extern template void
write_mesh<2>(
  const dealii::parallel::distributed::Triangulation<2> &triangulation,
  const std::string &                                    filename);

extern template void
write_mesh<3>(
  const dealii::parallel::distributed::Triangulation<3> &triangulation,
  const std::string &                                    filename);

///////////////////////////////////////////

#endif // __UTILITIES_H__