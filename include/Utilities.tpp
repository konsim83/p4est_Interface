#ifndef __UTILITIES_TPP__
#define __UTILITIES_TPP__

#include <Utilities.h>

template <int dim>
void
fill_points_randomly(std::vector<dealii::Point<dim>> &points,
                     double                           a,
                     double                           b,
                     bool                             same_on_all_ranks)
{
  RandomNumberDouble random_double_generator(a, b, same_on_all_ranks);

  for (auto &p : points)
    {
      for (unsigned int d = 0; d < dim; ++d)
        {
          p(d) = random_double_generator.generate();
        }
    }
}



template <int dim>
void
write_mesh(
  const dealii::parallel::distributed::Triangulation<dim> &triangulation,
  const std::string &                                      filename)
{
  MPI_Comm mpi_communicator = triangulation.get_communicator();
  dealii::ConditionalOStream pcout(std::cout,
                                   (dealii::Utilities::MPI::this_mpi_process(
                                      mpi_communicator) == 0));

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
        pcout << pair.first << " (" << pair.second << " times) ";
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

#endif // __UTILITIES_TPP__