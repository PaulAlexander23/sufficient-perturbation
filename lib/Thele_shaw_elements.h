// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//           Version 0.90. August 3, 2009.
// LIC//
// LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================

#pragma once

// Header file for THeleShaw elements
#ifndef OOMPH_THELE_SHAW_ELEMENTS_HEADER
#define OOMPH_THELE_SHAW_ELEMENTS_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// OOMPH-LIB headers
// hierher uncomment these
/* #include "../generic/nodes.h" */
/* #include "../generic/oomph_utilities.h" */
/* #include "../generic/Telements.h" */

#include "hele_shaw_elements.h"

namespace oomph
{
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  // THeleShawElement
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //======================================================================
  /// THeleShawElement<DIM,NNODE_1D> elements are isoparametric triangular
  /// DIM-dimensional HeleShaw elements with  NNODE_1D nodal points along each
  /// element edge. Inherits from TElement and HeleShawEquations
  //======================================================================
  template<unsigned NNODE_1D>
  class THeleShawElement : public virtual TElement<2, NNODE_1D>,
                           public virtual HeleShawEquations,
                           public virtual ElementWithZ2ErrorEstimator
  {
  public:
    ///\short  Constructor: Call constructors for TElement and
    /// HeleShaw equations
    THeleShawElement() : TElement<2, NNODE_1D>(), HeleShawEquations() {}


    /// Broken copy constructor
    THeleShawElement(const THeleShawElement<NNODE_1D>& dummy)
    {
      BrokenCopy::broken_copy("THeleShawElement");
    }

    /// Broken assignment operator
    void operator=(const THeleShawElement<NNODE_1D>&)
    {
      BrokenCopy::broken_assign("THeleShawElement");
    }

    /// \short  Access function for Nvalue: # of `values' (pinned or dofs)
    /// at node n (always returns the same value at every node, 1)
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return Initial_Nvalue;
    }

    /// \short Output function:
    ///  x,y,u   or    x,y,z,u
    void output(std::ostream& outfile)
    {
      HeleShawEquations::output(outfile);
    }

    ///  \short Output function:
    ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      HeleShawEquations::output(outfile, n_plot);
    }


    /// \short C-style output function:
    ///  x,y,u   or    x,y,z,u
    void output(FILE* file_pt)
    {
      HeleShawEquations::output(file_pt);
    }


    ///  \short C-style output function:
    ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      HeleShawEquations::output(file_pt, n_plot);
    }


    /// \short Output function for an exact solution:
    ///  x,y,u_exact
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
    {
      HeleShawEquations::output_fct(outfile, n_plot, exact_soln_pt);
    }


    /// \short Output function for a time-dependent exact solution.
    ///  x,y,u_exact (calls the steady version)
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      HeleShawEquations::output_fct(outfile, n_plot, time, exact_soln_pt);
    }

    //    void fill_in_contribution_to_jacobian_and_mass_matrix(
    //        Vector<double> &residuals,
    //        DenseMatrix<double> &jacobian,
    //        DenseMatrix<double> &mass_matrix
    //        )
    //        {
    ////            std::cout << "Calling THeleshaw" << std::endl;
    //            fill_in_contribution_to_jacobian(residuals,jacobian);
    //        }


  protected:
    /// Shape, test functions & derivs. w.r.t. to global coords. Return
    /// Jacobian.
    inline double dshape_and_dtest_eulerian_hele_shaw(const Vector<double>& s,
                                                      Shape& psi,
                                                      DShape& dpsidx,
                                                      Shape& test,
                                                      DShape& dtestdx) const;


    /// Shape, test functions & derivs. w.r.t. to global coords. Return
    /// Jacobian.
    inline double dshape_and_dtest_eulerian_at_knot_hele_shaw(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const;

    /// \short Shape/test functions and derivs w.r.t. to global coords at
    /// integration point ipt; return Jacobian of mapping (J). Also compute
    /// derivatives of dpsidx, dtestdx and J w.r.t. nodal coordinates.
    inline double dshape_and_dtest_eulerian_at_knot_hele_shaw(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      RankFourTensor<double>& d_dpsidx_dX,
      Shape& test,
      DShape& dtestdx,
      RankFourTensor<double>& d_dtestdx_dX,
      DenseMatrix<double>& djacobian_dX) const;

    /// \short Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return (NNODE_1D - 1);
    }

    /// Number of 'flux' terms for Z2 error estimation
    unsigned num_Z2_flux_terms()
    {
      return 2;
    }

    //// /// Get 'flux' for Z2 error recovery: Pressure gradient
    //// void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    ////  {
    ////////      std::cout << "Get Z2 flux" << std::endl;
    ////    Vector<double> pressure_gradient(2,0.0);
    //////    , velocity(2,0.0);
    ////      this->get_pressure_gradient(s,pressure_gradient);
    //////     this->get_velocity(s,velocity);
    ////     flux[0] = pressure_gradient[0];
    //////      - 2*velocity[0];
    //////      - 2*velocity[0];
    ////     flux[1] = pressure_gradient[1];
    //////     - 2*velocity[1];
    //////        std::cout << flux[0] << " " << flux[1] << std::endl;
    ////}

    /// Get 'flux' for Z2 error recovery: Pressure gradient
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
      Vector<double> pressure_gradient(2, 0.0), velocity(2, 0.0);
      this->get_pressure_gradient(s, pressure_gradient);
      this->get_velocity(s, velocity);
      flux[0] = 0 * pressure_gradient[0] - 2 * velocity[0];
      flux[1] = 0 * pressure_gradient[1] - 2 * velocity[1];
    }


    /// \short Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return TElement<2, NNODE_1D>::nvertex_node();
    }

    /// \short Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return TElement<2, NNODE_1D>::vertex_node_pt(j);
    }


  private:
    /// Static unsigned that holds the (same) number of variables at every node
    static const unsigned Initial_Nvalue;
  };


  // Inline functions:


  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned NNODE_1D>
  double THeleShawElement<NNODE_1D>::dshape_and_dtest_eulerian_hele_shaw(
    const Vector<double>& s,
    Shape& psi,
    DShape& dpsidx,
    Shape& test,
    DShape& dtestdx) const
  {
    unsigned n_node = this->nnode();

    // Call the geometrical shape functions and derivatives
    double J = this->dshape_eulerian(s, psi, dpsidx);

    // Loop over the test functions and derivatives and set them equal to the
    // shape functions
    for (unsigned i = 0; i < n_node; i++)
    {
      test[i] = psi[i];
      dtestdx(i, 0) = dpsidx(i, 0);
      dtestdx(i, 1) = dpsidx(i, 1);
    }

    // Return the jacobian
    return J;
  }


  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned NNODE_1D>
  double THeleShawElement<NNODE_1D>::
    dshape_and_dtest_eulerian_at_knot_hele_shaw(const unsigned& ipt,
                                                Shape& psi,
                                                DShape& dpsidx,
                                                Shape& test,
                                                DShape& dtestdx) const
  {
    // Call the geometrical shape functions and derivatives
    double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidx);

    // Set the pointers of the test functions
    test = psi;
    dtestdx = dpsidx;

    // Return the jacobian
    return J;
  }

  /// Define the shape functions (psi) and test functions (test) and
  /// their derivatives w.r.t. global coordinates (dpsidx and dtestdx)
  /// and return Jacobian of mapping (J). Additionally compute the
  /// derivatives of dpsidx, dtestdx and J w.r.t. nodal coordinates.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned NNODE_1D>
  double THeleShawElement<NNODE_1D>::
    dshape_and_dtest_eulerian_at_knot_hele_shaw(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      RankFourTensor<double>& d_dpsidx_dX,
      Shape& test,
      DShape& dtestdx,
      RankFourTensor<double>& d_dtestdx_dX,
      DenseMatrix<double>& djacobian_dX) const
  {
    // Call the geometrical shape functions and derivatives
    const double J = this->dshape_eulerian_at_knot(
      ipt, psi, dpsidx, djacobian_dX, d_dpsidx_dX);

    // Set the pointers of the test functions
    test = psi;
    dtestdx = dpsidx;
    d_dtestdx_dX = d_dpsidx_dX;

    // Return the jacobian
    return J;
  }


  //=======================================================================
  /// Face geometry for the THeleShawElement elements: The spatial
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<THeleShawElement<NNODE_1D>>
    : public virtual TElement<1, NNODE_1D>
  {
  public:
    /// \short Constructor: Call the constructor for the
    /// appropriate lower-dimensional TElement
    FaceGeometry() : TElement<1, NNODE_1D>() {}
  };


} // namespace oomph

#endif


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


// hierher break

// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//           Version 0.90. August 3, 2009.
// LIC//
// LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
// Non-inline functions for HeleShaw elements

// hierher uncomment
//#include "Thele_shaw_elements.h"


namespace oomph
{
  /////////////////////////////////////////////////////////////////////////
  // THeleShawElement
  /////////////////////////////////////////////////////////////////////////


  //======================================================================
  // Set the data for the number of Variables at each node, always 1
  //======================================================================
  template<unsigned NNODE_1D>
  const unsigned THeleShawElement<NNODE_1D>::Initial_Nvalue = 1;

  //====================================================================
  // Force build of templates
  //====================================================================
  template class THeleShawElement<2>;
  template class THeleShawElement<3>;
  template class THeleShawElement<4>;

} // namespace oomph