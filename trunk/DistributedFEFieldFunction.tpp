#ifndef __DISTRIBUTEDFEFIELDFUNCTION_TPP__
#define __DISTRIBUTEDFEFIELDFUNCTION_TPP__

#include <DistributedFEFieldFunction.h>

template <int dim>
DistributedFEFieldFunction<dim>::DistributedFEFieldFunction()
  : mpi_communicator(MPI_COMM_WORLD)
  , triangulation(mpi_communicator,
                  typename dealii::Triangulation<dim>::MeshSmoothing(
                    dealii::Triangulation<dim>::smoothing_on_refinement |
                    dealii::Triangulation<dim>::smoothing_on_coarsening))
  , fe(1)
  , dof_handler(triangulation)
  , pcout(std::cout,
          (dealii::Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
  , computing_timer(mpi_communicator,
                    pcout,
                    dealii::TimerOutput::summary,
                    dealii::TimerOutput::wall_times)
  , mapping(dealii::StaticMappingQ1<dim>::mapping)
  , cell_hint()
  , is_periodic(false)
  , is_initialized(false)
{}



template <int dim>
void
DistributedFEFieldFunction<dim>::generate_triangualtion(
  const unsigned int n_refine)
{
  dealii::TimerOutput::Scope t(computing_timer, "mesh generation");

  dealii::GridGenerator::hyper_cube(triangulation,
                                    0,
                                    1,
                                    /* colorize */ true);
  if (is_periodic)
    {
      std::vector<dealii::GridTools::PeriodicFacePair<
        typename dealii::parallel::distributed::Triangulation<
          dim>::cell_iterator>>
        periodicity_vector;

      for (unsigned int d = 0; d < dim; ++d)
        {
          dealii::GridTools::collect_periodic_faces(triangulation,
                                                    /*b_id1*/ 2 * (d + 1) - 2,
                                                    /*b_id2*/ 2 * (d + 1) - 1,
                                                    /*direction*/ d,
                                                    periodicity_vector);
        }

      triangulation.add_periodicity(periodicity_vector);
    } // if


  triangulation.refine_global(n_refine);

  dof_handler.distribute_dofs(fe);

  //  cell_hint = triangulation.begin_active();
  cell_hint = dof_handler.begin_active();

  is_initialized = true;
}



template <int dim>
void
DistributedFEFieldFunction<dim>::write_mesh(const std::string &filename)
{
  Assert(is_initialized, dealii::ExcNotInitialized());

  dealii::TimerOutput::Scope t(computing_timer, "write mesh");

  pcout << std::endl
        << "*** Writing mesh ***" << std::endl
        << "*** MPI ranks used       "
        << dealii::Utilities::MPI::n_mpi_processes(mpi_communicator)
        << std::endl
        << "*** Dimension:           " << dim << std::endl
        << "*** No. of cells:        " << triangulation.n_active_cells()
        << std::endl;

  /*
   * Print some general mesh info
   */
  {
    std::map<dealii::types::boundary_id, unsigned int> boundary_count;
    for (auto &cell : triangulation.active_cell_iterators())
      {
        for (unsigned int face = 0;
             face < dealii::GeometryInfo<dim>::faces_per_cell;
             ++face)
          {
            if (cell->face(face)->at_boundary())
              boundary_count[cell->face(face)->boundary_id()]++;
          }
      }

    pcout << "*** Boundary indicators: ";
    for (const std::pair<const dealii::types::boundary_id, unsigned int> &pair :
         boundary_count)
      {
        pcout << pair.first << "(" << pair.second << " times) ";
      }
    pcout << std::endl;
  }

  dealii::GridOut grid_out;

  grid_out.write_mesh_per_processor_as_vtu(triangulation,
                                           filename,
                                           /* view_levels */ false,
                                           /* include_artificials */ false);

  pcout << "*** Written to:          " << filename << ".pvtu" << std::endl
        << std::endl;
}



template <int dim>
double
DistributedFEFieldFunction<dim>::value(const dealii::Point<dim> &p)
{
  Assert(is_initialized, dealii::ExcNotInitialized());

  dealii::Timer timer;
  timer.restart();

  /*
   * Get access to some p4est internals
   */
  ForrestType *forrest = const_cast<ForrestType *>(triangulation.get_p4est());

  double value    = -1;
  int    mpi_rank = -1;

  typename dealii::DoFHandler<dim>::active_cell_iterator cell = cell_hint;
  if (cell == dof_handler.end())
    {
      cell = dof_handler.begin_active();
    }

  boost::optional<dealii::Point<dim>> qp =
    get_reference_coordinates(cell_hint, p);
  if (qp)
    {
      /*
       * If point is found to be in the (locally owned) cell = cell_hint
       * return this_mpi_rank
       */
      const dealii::CellId &cell_id = cell->id();

      mpi_rank = dealii::Utilities::MPI::this_mpi_process(mpi_communicator);
    }
  else
    {
      /*
       * If point is not found to in the (locally owned) cell = cell_hint
       * search for it on this processor.
       */
      const std::pair<
        typename dealii::internal::
          ActiveCellIterator<dim, dim, dealii::DoFHandler<dim>>::type,
        dealii::Point<dim>>
        my_pair = dealii::GridTools::find_active_cell_around_point(mapping,
                                                                   dof_handler,
                                                                   p);

      cell = my_pair.first;
      qp   = my_pair.second;
    }

  if (cell->is_locally_owned())
    {
      /*
       * If the cell found on this processor is neither ghost not artificial
       * then return this_mpi_rank.
       */
      mpi_rank = dealii::Utilities::MPI::this_mpi_process(mpi_communicator);
    }
  else
    {
      /*
       * Point is not owned by a local cell.
       */
      mpi_rank = -1;
    }

  timer.stop();
  std::cout << "---> MPI rank   "
            << dealii::Utilities::MPI::this_mpi_process(mpi_communicator)
            << "   search for point   " << p << " ....."
            << " done in   " << timer.cpu_time() << "   seconds.   "
            << " Found cell_id   " << cell->id() << "   and owner rank   "
            << mpi_rank << std::endl;

  return value;
}


template <int dim>
void
DistributedFEFieldFunction<dim>::value_list(
  const std::vector<dealii::Point<dim>> &points,
  std::vector<double>                    values)
{
  Assert(is_initialized, dealii::ExcNotInitialized());

  dealii::TimerOutput::Scope t(
    computing_timer,
    "finding owner's MPI rank (" +
      dealii::Utilities::int_to_string(points.size()) + " point)");

  //
  // Fill with life
  //
}


template <int dim>
boost::optional<dealii::Point<dim>>
DistributedFEFieldFunction<dim>::get_reference_coordinates(
  const typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
  const dealii::Point<dim> &                                    point) const
{
  try
    {
      dealii::Point<dim> qp = mapping.transform_real_to_unit_cell(cell, point);
      if (dealii::GeometryInfo<dim>::is_inside_unit_cell(qp))
        return qp;
      else
        return boost::optional<dealii::Point<dim>>();
    }
  catch (const typename dealii::Mapping<dim>::ExcTransformationFailed &)
    {
      // transformation failed, so
      // assume the point is
      // outside
      return boost::optional<dealii::Point<dim>>();
    }
}

#endif // __DISTRIBUTEDFEFIELDFUNCTION_TPP__