/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/scalar/CovariantLaplacianFermion.h

Copyright (C) 2018

Author: Guido Cossu <guido.cossu@ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */

#ifndef COVARIANT_LAPLACIAN_FERMION_H
#define COVARIANT_LAPLACIAN_FERMION_H

namespace Grid
{
namespace QCD
{

////////////////////////////////////////////////////////////
// Laplacian operator L on fields in the fundamental rep
//
// phi: field in the fundamental rep
//
// L phi(x) = Sum_mu [ U_mu(x) phi(x+mu) + U_mu(x-mu) phi(x-mu) - 2ND phi(x)]
//
// Operator designed to be encapsulated by
// an HermitianLinearOperator<.. , ..>
////////////////////////////////////////////////////////////

// has to inherit from an implementation
template <class Impl>
class LaplacianFundamental
{
public:
  INHERIT_IMPL_TYPES(Impl);

  // add a bool to smear only in the spatial directions
  LaplacianFundamental(GridBase *grid, bool spatial = false)
    : U(Nd, grid), spatial_laplacian(spatial){};
  
  void ImportGauge(const GaugeField &_U)
  {
    for (int mu = 0; mu < Nd; mu++)
      U[mu] = PeekIndex<LorentzIndex>(_U, mu);
  }
  
  void M(const FermionField &in, FermionField &out)
  {
    int dims = spatial_laplacian ? (Nd - 1) : Nd;
    
    out = -2.0 * dims * in;
    // eventually speed up with the stencil operator, if necessary
    for (int mu = 0; mu < dims; mu++)
      out += Impl::CovShiftForward(U[mu], mu, in) + Impl::CovShiftBackward(U[mu], mu, in);
  }
  
private:
  bool spatial_laplacian;
  std::vector<GaugeLinkField> U;
}; // Laplacian

} // QCD
} // Grid
#endif
