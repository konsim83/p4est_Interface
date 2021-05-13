// Deal.ii
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_q.h>

// Boost
#include <boost/optional.hpp>

// C
#include <stdio.h>

// STL
#include <fstream>
#include <iostream>

// My headers
#include <PartitionSearch.h>
#include <Utilities.h>


int
main(int argc, char *argv[])
{
  try
    {
      // dimension and global refinement
      const int          dim      = 2;
      const unsigned int n_refine = 2;

      // If we want to limit the search on one rank only
      const bool         limit_search_to_single_rank = true;
      const unsigned int do_search_on_rank           = 0;
      bool               same_points_on_all_ranks    = true;

      // MPI_Init
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, dealii::numbers::invalid_unsigned_int);

      // Setup partition search class
      PartitionSearch<dim> partition_search;
      partition_search.generate_triangulation(n_refine);

      // points whose owner ranks we want to find
      // generate a number of points
      const unsigned int              n_points = 1;
      std::vector<dealii::Point<dim>> test_points(n_points,
                                                  dealii::Point<dim>());

      if (limit_search_to_single_rank)
        {
          if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
              do_search_on_rank)
            {
              std::cout << std::endl
                        << "**************************************" << std::endl
                        << "Searching partition only on rank "
                        << do_search_on_rank << ":" << std::endl
                        << "**************************************" << std::endl
                        << std::endl;

              fill_points_randomly(test_points,
                                   /* min */ 0,
                                   /* max */ 1,
                                   same_points_on_all_ranks);

              for (auto &p : test_points)
                partition_search.find_owner_rank_p4est(p);
            }
        }
      else
        {
          // here all points
          fill_points_randomly(test_points,
                               /* min */ 0,
                               /* max */ 1,
                               same_points_on_all_ranks);

          for (auto &p : test_points)
            partition_search.find_owner_rank_p4est(p);
        }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
