#pragma once

#include "generic.h"
#include "meshes/triangle_mesh.template.h"

#include "modified_volume_constraint_elements_with_integrals.h"

//==start_of_problem_class============================================
/// Problem class to simulate inviscid bubble propagating along 2D channel
//====================================================================
template <class ELEMENT> class BubbleInChannelProblem : public Problem {
  /// Make a copy for using in bifurcation tracking

public:
  double Q_jack;

  /// Constructor
  BubbleInChannelProblem();

  Problem *make_copy();

  /// Destructor
  ~BubbleInChannelProblem();

  /// Actions before adapt: Wipe the mesh of free surface elements
  void actions_before_adapt();

  /// Actions after adapt: Rebuild the mesh of free surface elements
  void actions_after_adapt();

  double global_temporal_error_norm();

  /// Update the after solve (empty)
  void actions_after_newton_solve() {}

  /// Update the problem specs before solve
  void actions_before_newton_solve() {}

  /// \short Set boundary conditions and complete the build of all elements
  void complete_problem_setup();

  // Custom function to perform continuation in 1 parameter
  void jack_solve(double Q, bool jack_arc, unsigned n_iter, double sign);

  // Custom function to perform continuation in 2 parameters
  void jack_solve(double Q, double V, bool jack_arc, unsigned n_iter, int sign);

  void dump_it(ofstream &dump_file);

  /// Doc the solution
  void doc_solution(const std::string &comment = "");

  /// Compute the error estimates and assign to elements for plotting
  void compute_error_estimate(double &max_err, double &min_err);

  void reset_lagrangian_coordinates();

  double get_bubble_pressure();

  void create_inflow_elements();

  void delete_inflow_elements();

  void delete_outflow_elements();

  void pin_all_nodal_positions();

  void unpin_all_nodal_positions();

  void pin_bubble_pressure();

  void pin_all_pressure();

  void unpin_all_pressure();

  void pin_curvature_equations();

  void pin_tangential_lagrange();

  void unpin_curvature_equations();

  void reset_interface_equations();

  void redo_equation_numbering();

  //// Operations on parameters

  /// Volume
  void set_V(double new_V);
  double get_V();
  void increment_V(double dV);

  /// Frame speed (if pinned, frame does not move with bubble)
  void set_U(double new_U);
  double get_U();
  void pin_U();
  void unpin_U();

  /// Bubble pressure  (volume constraint is not enforced if P is pinned)
  void set_P(double new_P);
  double get_P();
  void pin_P();
  void unpin_P();

  /// Imposed flow rate
  void set_Q(double new_Q);
  double get_Q();
  void increment_Q(double d_Q);
  double get_Q_inv();

  /// We measure CoM_Y through an integral constraint.
  double get_CoM_Y();
  void set_CoM_Y(double new_CoM_Y);

  /// Channel geometry parameters
  void set_asymmetry(double new_asymmetry);
  void set_h(double new_h);
  void set_w(double new_w);
  void set_alpha(double new_alpha);
  double get_asymmetry();
  double get_h();
  double get_w();
  double get_maj_rad();
  double get_min_rad();
  double get_xcenter();
  double get_ycenter();
  unsigned get_circpts();
  double get_alpha();
  void increment_h(double dh);

  void output_jacobian();

  void set_steady();

private:
  /// \short Create free surface elements
  void create_free_surface_elements();

  /// \short Delete free surface elements
  void delete_free_surface_elements(); // end of delete_free_surface_elements

  void create_volume_constraint_elements();
  void delete_volume_constraint_elements();

  void create_CoM_X_constraint_elements();
  void delete_CoM_X_constraint_elements();

  void create_CoM_Y_constraint_elements();
  void delete_CoM_Y_constraint_elements();

  /// Pointers to mesh of free surface elements
  Mesh *Free_surface_mesh_pt;

  // /// Pointer to mesh containing elements that impose volume constraint
  Mesh *Volume_constraint_mesh_pt;
  Mesh *CoM_X_constraint_mesh_pt;
  Mesh *CoM_Y_constraint_mesh_pt;
  Mesh *Integral_measures_mesh_pt;

  /// Pointer to Fluid_mesh
  RefineableSolidTriangleMesh<ELEMENT> *Fluid_mesh_pt;

  /// Vector storing pointer to the bubble polygons
  Vector<TriangleMeshPolygon *> Bubble_polygon_pt;

  /// Triangle mesh polygon for outer boundary
  TriangleMeshPolygon *Outer_boundary_polyline_pt;

  Vector<TriangleMeshPolygon *> Obstacle_polygon_pt;

  /// Pointer to a global bubble pressure datum
  Data *Bubble_pressure_data_pt;

  Data *St_data_pt;
  Data *Alpha_data_pt;
  Data *Obstacle_height_data_pt;
  Data *Obstacle_width_data_pt;
  Data *G_data_pt;
  ///---------------------
public:
  Data *U_data_pt;
  Data *Q_inv_data_pt;
  ///---------------------
  Data *Asymmetry_data_pt;
  Data *CoM_Y_data_pt;

  Data *Integral_measures_data_pt;

  Mesh *Inflow_mesh_pt;
  Mesh *Outflow_mesh_pt;
  //// /// Pointer to element that imposes volume constraint for bubble
  VolumeConstraintElement *Vol_constraint_el_pt;
  VolumeConstraintElement *CoM_X_constraint_el_pt;
  //    VolumeConstraintElement* CoM_Y_constraint_el_pt;
  SelfReferentialVolumeConstraintElement *CoM_Y_constraint_el_pt;
  SelfReferentialVolumeConstraintElement *Integral_measures_el_pt;

  /// Enumeration of mesh boundaries
  enum {
    Inflow_boundary_id = 0,
    Upper_wall_boundary_id = 1,
    Outflow_boundary_id = 2,
    Bottom_wall_boundary_id = 3,
    First_bubble_boundary_id = 4,
    Second_bubble_boundary_id = 5,
  };

  ///-------------------------------------------------------------------------------

  void set_steady_constraints();

  void set_steady_constraints_free_Q();

  void pin_G();
  void pin_Q();

  void unpin_G();
  void unpin_Q();

  bool get_jacobian_sign_change();
  ///-------------------------------------------------------------------------------
}; // end_of_problem_class

#include "hele_shaw_bubble_problem_with_integrals.tpp"
