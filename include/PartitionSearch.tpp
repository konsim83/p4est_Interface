#ifndef __PARTITIONSEARCH_TPP__
#define __PARTITIONSEARCH_TPP__

#include <PartitionSearch.h>


template <int dim>
PartitionSearch<dim>::PartitionSearch()
  : mpi_communicator(MPI_COMM_WORLD)
  , triangulation(mpi_communicator,
                  typename dealii::Triangulation<dim>::MeshSmoothing(
                    dealii::Triangulation<dim>::smoothing_on_refinement |
                    dealii::Triangulation<dim>::smoothing_on_coarsening))
  , pcout(std::cout,
          (dealii::Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
  , computing_timer(mpi_communicator,
                    pcout,
                    dealii::TimerOutput::summary,
                    dealii::TimerOutput::wall_times)
{}



template <int dim>
void
PartitionSearch<dim>::generate_triangualtion(const unsigned int n_refine)
{
  dealii::TimerOutput::Scope t(computing_timer, "mesh generation");

  dealii::GridGenerator::hyper_L(triangulation, 0, 1);

  // dealii::GridGenerator::hyper_cube(triangulation,
  //                                   0,
  //                                   1,
  //                                   /* colorize */ true);

  triangulation.refine_global(n_refine);

  // mesh outputfor dignostic purposes
  write_mesh(triangulation, "my_triangulation");
}



/*******************************/
/********** Callbacks **********/
/*******************************/
template <int dim>
int
PartitionSearch<dim>::my_local_quadrant_fn(
  typename internal::types<dim>::forest *  p4est,
  typename internal::types<dim>::topidx    which_tree,
  typename internal::types<dim>::quadrant *quadrant,
  int                                      pfirst,
  int                                      plast,
  void *                                   point)
{
  // User pointer is a parallel triangualtion
  dealii::parallel::distributed::Triangulation<dim> *this_triangualtion_ptr =
    static_cast<dealii::parallel::distributed::Triangulation<dim> *>(
      p4est->user_pointer);

  // Check some things
  Assert(this_triangualtion_ptr != nullptr, dealii::ExcInternalError());
  Assert(this_triangualtion_ptr == this_triangualtion_ptr,
         dealii::ExcInternalError());
  //    Assert (0 <= pfirst && pfirst <= plast && plast < g->mpisize,
  //    dealii::ExcInternalError());
  Assert(point == nullptr,
         dealii::ExcInternalError()); // point must be nullptr here

  MPI_Comm this_mpi_communicator = this_triangualtion_ptr->get_communicator();

  // if (verbose_in_quadrant_fn)

  /* we are not trying to find local spheres */
  //    if (pfirst == plast && pfirst == g->mpirank)
  //      {
  //        return 0;
  //      }

  //    std::cout << "[Rank "
  //              <<
  //              dealii::Utilities::MPI::this_mpi_process(this_mpi_communicator)
  //              << "]   "
  //              << "Quadrant function called in tree: " << which_tree << "
  //              ("
  //              << pfirst << ", " << plast << ")   "
  //              << "Quadrant length   " <<
  //              P4EST_QUADRANT_LEN(quadrant->level)
  //              << "   on level   " << quadrant->level
  //              << "   coord:   " << quadrant->x << " " << quadrant->y << "
  //              ";
  //    if (dim == 3)
  //      {
  //        std::cout << quadrant->z << " ";
  //      }
  //    std::cout << std::endl;

  return /* true */ 1;
}



template <int dim>
int
PartitionSearch<dim>::my_local_point_fn(
  typename internal::types<dim>::forest *  p4est,
  typename internal::types<dim>::topidx    which_tree,
  typename internal::types<dim>::quadrant *quadrant,
  int                                      pfirst,
  int                                      plast,
  void *                                   point)
{
  // User pointer is a parallel triangualtion
  dealii::parallel::distributed::Triangulation<dim> *this_triangualtion_ptr =
    static_cast<dealii::parallel::distributed::Triangulation<dim> *>(
      p4est->user_pointer);

  // Check some things
  Assert(this_triangualtion_ptr != nullptr, dealii::ExcInternalError());
  Assert(this_triangualtion_ptr == this_triangualtion_ptr,
         dealii::ExcInternalError());
  //    Assert (0 <= pfirst && pfirst <= plast && plast < g->mpisize,
  //    dealii::ExcInternalError());
  Assert(point != nullptr,
         dealii::ExcInternalError()); // point must NOT be nullptr here

  MPI_Comm this_mpi_communicator = this_triangualtion_ptr->get_communicator();

  /*
   * Do some point checks for the point to return true (nonzero) or false
   * (zero).
   */

  const auto quad_length_on_level =
    1 << (static_cast<int>(P4EST_MAXLEVEL) - static_cast<int>(quadrant->level));
  // P4EST_QUADRANT_LEN(quadrant->level);

  double *point_double = (double *)point;

  if (verbose_in_point_fn)
    {
      std::cout << "[Rank "
                << dealii::Utilities::MPI::this_mpi_process(
                     this_mpi_communicator)
                << "]" << std::endl;

      std::cout << "   "
                << "physical point:   ( ";
      for (unsigned int d = 0; d < dim; ++d)
        {
          if (d < dim - 1)
            std::cout << point_double[d] << " | ";
          else
            std::cout << point_double[d] << " )" << std::endl;
        }

      std::cout << "   "
                << "Point function called in tree   : " << which_tree
                << std::endl
                << "   "
                << "pfirst = " << pfirst << std::endl
                << "   "
                << "plast = " << plast << std::endl
                << "   "
                << "level = " << static_cast<int>(quadrant->level) << std::endl
                << "   "
                << "element coords:   (" << quadrant->x << " | " << quadrant->y;

      // std::cout << " | " << quadrant->z << ")" << std::endl;
      std::cout << ")" << std::endl;
    }

  /* we may be up in the tree branches */
  if (pfirst < plast)
    {
      if (verbose_in_point_fn)
        std::cout << "   "
                  << "---> continuing recursion (approx match found)"
                  << std::endl
                  << std::endl;

      /* an optimistic match is good enough when we're walking the tree */
      return /* true */ 1;
    }

  if (verbose_in_point_fn)
    std::cout << "   "
              << "---> stopping recursion (pfirst==plast)" << std::endl
              << "   "
              << "---> pfirst    : " << pfirst << std::endl
              << "   "
              << "---> owner rank: ?" << std::endl
              << std::endl;

  return /* false */ 0;
}
/*******************************/
/*******************************/



template <int dim>
int
PartitionSearch<dim>::find_owner_rank_p4est(const dealii::Point<dim> &p)
{
  dealii::Timer timer;
  timer.restart();

  /*
   * Get access to some p4est internals
   */
  ForrestType *forrest = const_cast<ForrestType *>(triangulation.get_p4est());

  int mpi_rank = -1;

  // assign pointer to an array of points.
  sc_array_t *single_point;

  // allocate memory for a single dim-dimensional point
  single_point = sc_array_new_count(sizeof(double[dim]), 1);

  // provide a pointer to access the point's coordinates
  double *my_point = (double *)sc_array_index_int(single_point, 0);

  // now assigne the actual value
  for (unsigned int d = 0; d < dim; ++d)
    {
      my_point[d] = p(d);
    }

  // Call the search partition function
  internal::functions<dim>::search_partition(
    forrest, 0, my_local_quadrant_fn, my_local_point_fn, single_point);

  timer.stop();

  ////////////////////////////////

  std::cout << "[MPI Rank "
            << dealii::Utilities::MPI::this_mpi_process(mpi_communicator)
            << "] searched for point   (" << p << ")" << std::endl
            << "              -> "
            << "done in   " << timer.cpu_time() << " sec" << std::endl
            << "              -> "
            << "owner rank   " << mpi_rank << std::endl;

  my_point = nullptr;
  // Free the allocated memory
  sc_array_destroy_null(&single_point);

  return mpi_rank;
}

#endif // __PARTITIONSEARCH_TPP__