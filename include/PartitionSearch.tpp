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
PartitionSearch<dim>::generate_triangulation(const unsigned int n_refine)
{
  dealii::TimerOutput::Scope t(computing_timer, "mesh generation");

  dealii::GridGenerator::hyper_L(triangulation, 0, 1);

  // dealii::GridGenerator::hyper_cube(triangulation,
  //                                   0,
  //                                   1,
  //                                   /* colorize */ true);

  triangulation.refine_global(n_refine);

  auto cell = triangulation.begin_active();
  cell->set_refine_flag();
  ++cell;
  ++cell;
  cell->set_refine_flag();
  // triangulation.prepare_coarsening_and_refinement();
  triangulation.execute_coarsening_and_refinement();

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
  // User pointer is a parallel triangulation
  dealii::parallel::distributed::Triangulation<dim> *this_triangulation_ptr =
    static_cast<dealii::parallel::distributed::Triangulation<dim> *>(
      p4est->user_pointer);

  // Check some things
  Assert(this_triangulation_ptr != nullptr, dealii::ExcInternalError());
  Assert(this_triangulation_ptr == this_triangulation_ptr,
         dealii::ExcInternalError());
  //    Assert (0 <= pfirst && pfirst <= plast && plast < g->mpisize,
  //    dealii::ExcInternalError());
  Assert(point == nullptr,
         dealii::ExcInternalError()); // point must be nullptr here

  MPI_Comm this_mpi_communicator = this_triangulation_ptr->get_communicator();

  // if (verbose_in_quadrant_fn)



  part_global_t *g = (part_global_t *)p4est->user_pointer;
  /* compute coordinate range of this quadrant */
  loopquad(g, which_tree, quadrant, g->lxyz, g->hxyz, g->dxyz);
  /* always return 1 to search particles individually */
  return 1;



  static void loopquad(part_global_t * g,
                       p4est_topidx_t tt,
                       p4est_quadrant_t * quad,
                       double lxyz[3],
                       double hxyz[3],
                       double dxyz[3])
  {
    int            i;
    p4est_qcoord_t qh;

    qh = P4EST_QUADRANT_LEN(quad->level);
    p4est_qcoord_to_vertex(g->conn, tt, quad->x, quad->y, quad->z, lxyz);
    p4est_qcoord_to_vertex(
      g->conn, tt, quad->x + qh, quad->y + qh, quad->z + qh, hxyz);

    for (i = 0; i < 3; ++i)
      {
        lxyz[i] /= g->bricklength;
        hxyz[i] /= g->bricklength;
        dxyz[i] = hxyz[i] - lxyz[i];
      }
  }

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
  // User pointer is a parallel triangulation
  dealii::parallel::distributed::Triangulation<dim> *this_triangulation_ptr =
    static_cast<dealii::parallel::distributed::Triangulation<dim> *>(
      p4est->user_pointer);

  // Check some things
  Assert(this_triangulation_ptr != nullptr, dealii::ExcInternalError());
  Assert(this_triangulation_ptr == this_triangulation_ptr,
         dealii::ExcInternalError());
  //    Assert (0 <= pfirst && pfirst <= plast && plast < g->mpisize,
  //    dealii::ExcInternalError());
  Assert(point != nullptr,
         dealii::ExcInternalError()); // point must NOT be nullptr here

  MPI_Comm this_mpi_communicator = this_triangulation_ptr->get_communicator();

  /*
   * Do some point checks for the point to return true (nonzero) or false
   * (zero).
   */

  const int p4est_maxlevel = P4EST_MAXLEVEL;

  const auto quad_length_on_level =
    1 << (static_cast<int>(p4est_maxlevel) - static_cast<int>(quadrant->level));
  // P4EST_QUADRANT_LEN(quadrant->level);

  double *point_double = (double *)point;

  if (verbose_in_point_fn)
    {
      std::cout << "[Rank "
                << dealii::Utilities::MPI::this_mpi_process(
                     this_mpi_communicator)
                << "] p = ( ";
      for (unsigned int d = 0; d < dim; ++d)
        {
          if (d < dim - 1)
            std::cout << point_double[d] << " | ";
          else
            std::cout << point_double[d] << " )";
        }

      std::cout << " | Point_fnc called | tree = " << which_tree
                << " | pfirst = " << pfirst << " | plast = " << plast
                << " | level = " << static_cast<int32_t>(quadrant->level)
                << " | element coords = (" << static_cast<int32_t>(quadrant->x)
                << " | " << static_cast<int32_t>(quadrant->y);

      // std::cout << " | " << quadrant->z << ")" << std::endl;
      std::cout << ")";
    }

  /* we may be up in the tree branches */
  if (pfirst < plast)
    {
      if (verbose_in_point_fn)
        std::cout << " ---> continuing ..." << std::endl << std::endl;

      /* an optimistic match is good enough when we're walking the tree */
      return /* true */ 1;
    }

  if (verbose_in_point_fn)
    std::cout << " ---> stopping!" << std::endl << std::endl;

  return /* false */ 1;
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
            << "] searched for point   (" << p << ")"
            << " | done in   " << timer.cpu_time() << " sec"
            << " | owner rank   " << mpi_rank << std::endl;

  my_point = nullptr;
  // Free the allocated memory
  sc_array_destroy_null(&single_point);

  return mpi_rank;
}

#endif // __PARTITIONSEARCH_TPP__