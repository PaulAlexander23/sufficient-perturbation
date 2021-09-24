#include <string>

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
