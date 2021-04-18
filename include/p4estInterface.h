#ifndef __P4ESTINTERFACE_H__
#define __P4ESTINTERFACE_H__

// p4est
#include <p4est.h>
#include <p4est_search.h>
#include <p8est.h>
#include <p8est_search.h>

/*
 * We deliberately stay as close as possible to the existing interface in
 * Deal.II
 */

namespace internal
{
  /*
   * A unified interface for types in p4est
   */
  template <int dim>
  struct types;

  template <>
  struct types<2>
  {
    using connectivity     = p4est_connectivity_t;
    using forest           = p4est_t;
    using tree             = p4est_tree_t;
    using quadrant         = p4est_quadrant_t;
    using topidx           = p4est_topidx_t;
    using locidx           = p4est_locidx_t;
    using gloidx           = p4est_gloidx_t;
    using balance_type     = p4est_connect_type_t;
    using ghost            = p4est_ghost_t;
    using transfer_context = p4est_transfer_context_t;

    // Callback types (function pointers)
    using search_partition_callback = p4est_search_partition_t;
  };

  template <>
  struct types<3>
  {
    using connectivity     = p8est_connectivity_t;
    using forest           = p8est_t;
    using tree             = p8est_tree_t;
    using quadrant         = p8est_quadrant_t;
    using topidx           = p4est_topidx_t;
    using locidx           = p4est_locidx_t;
    using gloidx           = p4est_gloidx_t;
    using balance_type     = p8est_connect_type_t;
    using ghost            = p8est_ghost_t;
    using transfer_context = p8est_transfer_context_t;

    // Callback types (function pointers)
    using search_partition_callback = p8est_search_partition_t;
  };

  /////////////////////////////////////
  /////////////////////////////////////

  /*
   * Declare a unified interface for p4est functions
   */
  template <int dim>
  struct functions;

  template <>
  struct functions<2>
  {
    static void (&search_partition)(
      types<2>::forest *                  p4est,
      int                                 call_post,
      types<2>::search_partition_callback quadrant_fn,
      types<2>::search_partition_callback point_fn,
      sc_array_t *                        points);
  };

  template <>
  struct functions<3>
  {
    static void (&search_partition)(
      types<3>::forest *                  p4est,
      int                                 call_post,
      types<3>::search_partition_callback quadrant_fn,
      types<3>::search_partition_callback point_fn,
      sc_array_t *                        points);
  };

  /////////////////////////////////////
  /////////////////////////////////////

  /*
   * Now assign the declared interface the actual function in p4est
   */
  void (&functions<2>::search_partition)(
    types<2>::forest *                  p4est,
    int                                 call_post,
    types<2>::search_partition_callback quadrant_fn,
    types<2>::search_partition_callback point_fn,
    sc_array_t *                        points) = p4est_search_partition;


  void (&functions<3>::search_partition)(
    types<3>::forest *                  p4est,
    int                                 call_post,
    types<3>::search_partition_callback quadrant_fn,
    types<3>::search_partition_callback point_fn,
    sc_array_t *                        points) = p8est_search_partition;
  /////////////////////////////////////
  /////////////////////////////////////

} // namespace internal

#endif // __P4ESTINTERFACE_H__