#ifndef __PARTITIONSEARCH_H__
#define __PARTITIONSEARCH_H__

// Deal.ii
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/p4est_wrappers.h>
#include <deal.II/distributed/tria.h>

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
#include <Utilities.h>
#include <p4estInterface.h>

/*
 * Class finds the MPI rank of the partition of the distributed triangulation
 * object that owns points.
 */
template <int dim>
class PartitionSearch
{
public:
  /*
   * Convenience typedef
   */
  using ForrestType = typename dealii::internal::p4est::types<dim>::forest;

  PartitionSearch();

  PartitionSearch(const PartitionSearch<dim> &other) = delete;

  PartitionSearch<dim> &
  operator=(const PartitionSearch<dim> &other) = delete;

  void
  generate_triangualtion(const unsigned int n_refine);

  int
  find_owner_rank_p4est(const dealii::Point<dim> &p);

private:
  /*******************************/
  /********** Callbacks **********/
  /*******************************/
  static int
  my_local_quadrant_fn(typename internal::types<dim>::forest *  p4est,
                       typename internal::types<dim>::topidx    which_tree,
                       typename internal::types<dim>::quadrant *quadrant,
                       int                                      pfirst,
                       int                                      plast,
                       void *                                   point);

  static int
  my_local_point_fn(typename internal::types<dim>::forest *  p4est,
                    typename internal::types<dim>::topidx    which_tree,
                    typename internal::types<dim>::quadrant *quadrant,
                    int                                      pfirst,
                    int                                      plast,
                    void *                                   point);
  /*******************************/
  /*******************************/

  MPI_Comm mpi_communicator;

  dealii::parallel::distributed::Triangulation<dim> triangulation;

  dealii::ConditionalOStream pcout;

  dealii::TimerOutput computing_timer;

  static const bool verbose_in_point_fn    = true;
  static const bool verbose_in_quadrant_fn = false;
};

// Extern template instantiations
extern template class PartitionSearch<2>;
extern template class PartitionSearch<3>;

#endif // __PARTITIONSEARCH_H__