//==start_of_namespace==============================
/// Namespace for Problem Parameter
//==================================================
namespace Problem_Parameter
{
  /// Doc info object
  DocInfo Doc_info;

  /// Pseudo-solid Poisson ratio
  double Nu = 0.3;

  // double Radius=0.7;


  bool Minisolve_state = false;
  bool ignore_height_effects_in_dynamic_bc = false;
  bool ignore_height_effects_in_viscous_terms = false;

  bool Set_steady_conditions_after_adapt = false;

  double* global_G_pt = 0;
  double* global_G_upper_outflow_pt = 0;
  double* global_Q_inv_pt = 0;
  double* global_Frame_speed_pt = 0;
  double* global_Obstacle_height_pt = 0;
  double* global_Obstacle_width_pt = 0;
  double* global_alpha_pt = 0;
  double* global_Asymmetry_pt = 0;
  double* global_Bubble_pressure_pt = 0;
  double* global_St_pt = 0;
  double* global_delta_pt = 0;
  double* global_CoM_Y_pt = 0;
  double* global_measure_pt = 0;
  int global_Q_zero = 0;

  double Major_Radius;
  double Minor_Radius;
  double xcenter;
  double ycenter;
  unsigned circpts;

  // Jack
  double hessian_tolerance;
  double tressian_tolerance;
  int hessian_order;
  int tressian_order;
  int eps_number = 0;
  int theta_number = 0;
  ofstream Trace_weakly_nonlinear;

  unsigned n_integral_measures = 12;

  // Node* Tip_node_pt;
  bool ignore_bulk_effects = false;

  bool fill_in_solid_jacobian_interface_by_fd = false;

  bool hopf1;

  /// \short Volume of the bubble (negative because it's outside the
  /// fluid!)
  double Volume = -MathematicalConstants::Pi * Major_Radius * Minor_Radius;
  double Centre_of_mass = 0;

  /// \short Scaling factor for inflow velocity (allows it to be switched off
  /// to do hydrostatics)
  double Inflow_veloc_magnitude = 0.0;

  /// \short Length of the channel
  double Length = 4.0;

  /// Constitutive law used to determine the mesh deformation
  ConstitutiveLaw* Constitutive_law_pt = 0;

  /// Trace file
  ofstream Trace_file;

  /// \short File to document the norm of the solution (for validation
  /// purposes -- triangle doesn't give fully reproducible results so
  /// mesh generation/adaptation may generate slightly different numbers
  /// of elements on different machines!)
  ofstream Norm_file;

  ofstream OccluHeight_file;

  void channel_height_function(const Vector<double>& x, double& b)
  {
    /// This function should have obstacle width and height, asymmetry and
    /// possibly sharpness.

    double height = 0.1;
    double width = 0.33;
    double asymmetry = 0;
    double sharpness = 40; /// Quite smooth for now.
    sharpness = 40; /// ALICE!
    if (global_Obstacle_height_pt != 0)
    {
      height = *global_Obstacle_height_pt;
    }
    if (global_Obstacle_width_pt != 0)
    {
      width = *global_Obstacle_width_pt;
    }
    if (global_Asymmetry_pt != 0)
    {
      asymmetry = *global_Asymmetry_pt;
    }
    double y = x[1];
    double tanh_plus = std::tanh(sharpness * (y + width));
    double tanh_minus = std::tanh(sharpness * (y - width));
    ///--- Introduce oscillatory roughness to the occlusion 17/09/2015 ---
    //      double w1 = 50.0;   double w2 = 100.0;   double w3 = 200.0;
    //      double A1 = 0.05;   double A2 = 0.1;     double A3 = 0.2;
    //      double rough = A1*sin(w1*y)+A2*sin(w2*y)+A3*sin(w3*y);
    ///
    //      double A = 0.01; double B = 0.01;
    //      double w = 50;
    //      double a1 = 50;   double a2 = 1000;
    //      double roughCell = ( A*(sin(w*y)+1) + B*exp(-a1*y*y) )*exp(-a2*y*y);
    ///------ wave cross section ------
    //      double n = 50;
    //      double C = 0.0005;  /// worked with 0.005
    //      double waveCell = C*sin(n*M_PI*y);
    ///--- Introduce random roughness to the occlusion 18/09/2015 ---
    int sdr = abs((int)(y * 100000));
    srand(sdr);
    double rd = 0.01 * ((double)rand() / (double)(RAND_MAX));
    ///-------------------------------------------------------

    b = 1 - height * 0.5 * (tanh_plus - tanh_minus); /// current profile
    //     b = 1 - ( height * 0.5 *(tanh_plus - tanh_minus) * (1+rough) +
    //     waveCell ); b = 1 - ( height * 0.5 *(tanh_plus - tanh_minus) +
    //     waveCell );
    double c = 1 - height * 0.5 * (tanh_plus - tanh_minus) * (1 + rd);
    b += -y * asymmetry;

    double Q = 1 / (*global_Q_inv_pt);
    if (Q > 0.0100409 && Q < 0.0101489)
    {
      OccluHeight_file << y << "   " << b << "   " << sdr << "   " << rd
                       << "   " << c << std::endl;
    }
  }

  void channel_height_derivative(const Vector<double>& x, Vector<double>& dbdx)
  {
    // std::cout<<"I am here"<<std::endl;
    double y = x[1];
    double height = 0.1;
    double width = 0.33;
    double asymmetry = 0;
    double sharpness = 40; /// Quite smooth for now.
    sharpness = 40; /// ALICE!
    if (global_Obstacle_height_pt != 0)
    {
      height = *global_Obstacle_height_pt;
    }
    if (global_Obstacle_width_pt != 0)
    {
      width = *global_Obstacle_width_pt;
    }
    if (global_Asymmetry_pt != 0)
    {
      asymmetry = *global_Asymmetry_pt;
    }

    double cosh_plus = pow(std::cosh(sharpness * (y + width)), 2);
    double cosh_minus = pow(std::cosh(sharpness * (y - width)), 2);

    dbdx[1] = -sharpness * height * 0.5 * (1.0 / cosh_plus - 1.0 / cosh_minus);
    dbdx[1] += -asymmetry;
    dbdx[0] = 0.0;
  }


  ofstream UpperWall_file;

  void upper_wall_flux_fct(const Vector<double>& x,
                           double& b,
                           double& dbdt,
                           Vector<double>& dbdx,
                           Vector<double>& d_dbdt_dx)
  {
    channel_height_function(x, b);
    if (ignore_bulk_effects)
    {
      b = 1;
      UpperWall_file << b << std::endl;
    }
    dbdt = 0.0;

    channel_height_derivative(x, dbdx);

    d_dbdt_dx[0] = 0.0;
    d_dbdt_dx[1] = 0.0;
  }

  void upper_wall_fct(const Vector<double>& x, double& b, double& dbdt)
  {
    channel_height_function(x, b);
    if (ignore_bulk_effects)
    {
      b = 1;
      UpperWall_file << b << std::endl;
    }
    dbdt = 0;
  }

  void frame_speed_fct(const Vector<double>& x, Vector<double>& U_frame)
  {
    if (global_Frame_speed_pt != 0)
    {
      U_frame[0] = *global_Frame_speed_pt;
      U_frame[1] = 0;
    }
    else
    {
      U_frame[0] = 0;
      U_frame[1] = 0;
    }
  }

  void bubble_pressure_function(const Vector<double>& x, double& pressure)
  {
    /// For slight simplicity, we can put the transverse curvature term as a
    /// spatially varying bubble pressure.

    double height = 1.0;
    if (ignore_height_effects_in_dynamic_bc == false)
    {
      channel_height_function(x, height);
    }
    double alpha = *global_alpha_pt;
    double Q = 1 / (*global_Q_inv_pt);
    pressure = -1.0 / (3.0 * Q * alpha * height) + (*global_Bubble_pressure_pt);
  }

  void bubble_pressure_gradient_function(const Vector<double>& x,
                                         Vector<double>& grad_pb)
  {
    /// ALICE_FLAG
    std::cout << "Evaluating bubble pressure gradient" << std::endl;
    std::cout << "Missing height dependence!" << std::endl;
  }

  void normal_flux_ahead_of_bubble(const Vector<double>& x, double& flux)
  {
    /// Ahead of the bubble we have the boundary condition dp/dx = -G.
    /// We need to supply the function b^3 n.grad p = -G b^3.
    double b;
    channel_height_function(x, b);
    if (ignore_bulk_effects)
    {
      b = 1;
    }
    double G;
    if (global_G_pt != 0)
    {
      G = *global_G_pt;
    }
    else
    {
      G = 1;
    }
    double p_x = -1 * G;
    flux = p_x * (b * b * b);
  }

  void normal_flux_behind_bubble(const Vector<double>& x, double& flux)
  {
    /// Ahead of the bubble we have the boundary condition dp/dx = -G.
    /// We need to supply the function b^3 n.grad p = -G b^3.
    double b;
    channel_height_function(x, b);
    if (ignore_bulk_effects)
    {
      b = 1;
    }
    double G;
    if (global_G_pt != 0)
    {
      G = *global_G_pt;
    }
    else
    {
      G = 1;
    }
    double p_x = G;
    flux = p_x * (b * b * b);
  }

  Vector<double> vector_of_eigenvalues_rp(10, 2.0);
  Vector<double> vector_of_eigenvalues_ip(10, 2.0);

  ofstream M_file;
} // namespace Problem_Parameter


namespace oomph
{
  class MyNewElement
    : public virtual ProjectableHeleShawElement<
        PseudoSolidNodeUpdateElement<THeleShawElement<3>, TPVDElement<2, 3>>>
  {
  private:
    /// Storage for elemental error estimate -- used for post-processing
    double Error;

  public:
    /// Constructor initialise error
    MyNewElement()
    {
      Error = 0.0;
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      /// The bulk element has no time derivative terms.
      fill_in_contribution_to_jacobian(residuals, jacobian);
    }

    /// Set error value for post-processing
    void set_error(const double& error)
    {
      Error = error;
    }

    void get_error(double& error)
    {
      error = Error;
    }

    /// Overload output function
    void output(std::ostream& outfile, const unsigned& nplot)
    {
      // Vector of local coordinates and velocity
      Vector<double> s(2);
      Vector<double> x(2);
      Vector<double> velocity(2);
      Vector<double> pressure_gradient(2);

      double h, dhdt;
      Vector<double> dhdx(2), d_dhdt_dx(2);
      // Tecplot header info
      outfile << tecplot_zone_string(nplot);

      // Loop over plot points
      unsigned num_plot_points = nplot_points(nplot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates and velocity at plot point
        get_s_plot(iplot, nplot, s);
        get_velocity(s, velocity);
        get_pressure_gradient(s, pressure_gradient);

        for (unsigned i = 0; i < 2; i++)
        {
          x[i] = interpolated_x(s, i);
          outfile << interpolated_x(s, i) << " ";
        }

        get_upper_wall_flux_data(0, x, h, dhdt, dhdx, d_dhdt_dx);

        outfile << velocity[0] << " " << velocity[1] << " "
                << interpolated_p_hele_shaw(s) << " " << pressure_gradient[0]
                << " " << pressure_gradient[1] << " " << Error << " " << size()
                << " " << h << " "
                << "\n";
      }

      // Write tecplot footer (e.g. FE connectivity lists)
      write_tecplot_zone_footer(outfile, nplot);
    }


    /// Overload output function
    void output(std::ostream& outfile)
    {
      unsigned nplot = 3;
      output(outfile, nplot);
    }
  };


  template<>
  class FaceGeometry<MyNewElement> : public virtual SolidTElement<1, 3>
  {
  public:
    FaceGeometry() : SolidTElement<1, 3>() {}
  };

  template<>
  class FaceGeometry<FaceGeometry<MyNewElement>>
    : public virtual SolidPointElement
  {
  public:
    FaceGeometry() : SolidPointElement() {}
  };


} // namespace oomph
