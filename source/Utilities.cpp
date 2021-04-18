#include <Utilities.h>
#include <Utilities.tpp>

///////////////////////////////////////////

// Template instantiations
template void
fill_points_randomly<2>(std::vector<dealii::Point<2>> &points,
                        double                         a,
                        double                         b,
                        bool                           same_on_all_ranks);
template void
fill_points_randomly<3>(std::vector<dealii::Point<3>> &points,
                        double                         a,
                        double                         b,
                        bool                           same_on_all_ranks);

///////////////////////////////////////////

// Template instantiations
template void
write_mesh<2>(
  const dealii::parallel::distributed::Triangulation<2> &triangulation,
  const std::string &                                    filename);

template void
write_mesh<3>(
  const dealii::parallel::distributed::Triangulation<3> &triangulation,
  const std::string &                                    filename);

///////////////////////////////////////////