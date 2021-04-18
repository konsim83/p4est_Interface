#ifndef __DISTRIBUTEDFEFIELDFUNCTION_H__
#define __DISTRIBUTEDFEFIELDFUNCTION_H__

// Deal.ii
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/p4est_wrappers.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

// Boost
#include <boost/optional.hpp>

// C
#include <stdio.h>

// STL
#include <fstream>
#include <iostream>
#include <map>
#include <typeinfo>
#include <vector>

// My headers
#include <RandomNumber.h>
#include <p4estInterface.h>

/*
 * Class finds the MPI rank of the partition of the distributed triangulation
 * object that owns points.
 */
template <int dim>
class DistributedFEFieldFunction
{
public:
  DistributedFEFieldFunction();

  DistributedFEFieldFunction(const DistributedFEFieldFunction<dim> &other) =
    delete;

  DistributedFEFieldFunction<dim> &
  operator=(const DistributedFEFieldFunction<dim> &other) = delete;

  void
  generate_triangualtion(const unsigned int n_refine);

  void
  write_mesh(const std::string &filename);

  double
  value(const dealii::Point<dim> &p);

  void
  value_list(const std::vector<dealii::Point<dim>> &points,
             std::vector<double>                    values);

private:
  boost::optional<dealii::Point<dim>>
  get_reference_coordinates(
    const typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
    const dealii::Point<dim> &                                    point) const;

  MPI_Comm mpi_communicator;

  dealii::parallel::distributed::Triangulation<dim> triangulation;

  dealii::FE_Q<dim> fe;

  dealii::DoFHandler<dim> dof_handler;

  dealii::ConditionalOStream pcout;

  dealii::TimerOutput computing_timer;

  const dealii::Mapping<dim> &mapping;

  typename dealii::DoFHandler<dim>::active_cell_iterator cell_hint;

  const bool is_periodic;

  bool is_initialized;
};

// Extern template instantiations
extern template class DistributedFEFieldFunction<2>;
extern template class DistributedFEFieldFunction<3>;

#endif // __DISTRIBUTEDFEFIELDFUNCTION_H__