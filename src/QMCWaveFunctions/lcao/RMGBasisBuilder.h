//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_RMGBASIS_BUILDER_H
#define QMCPLUSPLUS_RMGBASIS_BUILDER_H


#include "Message/MPIObjectBase.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "io/hdf_archive.h"

namespace qmcplusplus
{
/** atomic basisset builder
   * @tparam COT, CenteredOrbitalType = SoaAtomicBasisSet<RF,SH>
   *
   * Reimplement AtomiSPOSetBuilder.h
   */
template<typename VALT>
class AOBasisBuilder<RMGBasisSet<VALT>> : public MPIObjectBase
{

//Cij need to be the same format as LCAO
//orbitals can be defined by itself. better to be hdf5 format inwhich can be folder files .. 
//eigenvector and eigenvalues need to have the same format. 
private:
  using COT = RMGBasisSet<VALT>;
public:
  AOBasisBuilder(const std::string& eName, Communicate* comm)
      : MPIObjectBase(comm)
  { }

  bool putH5(hdf_archive& hin)
  { }

  COT* createAOSetH5(hdf_archive& hin)
  { 

      int num_orb;
      hin.read(num_orb, "NumOrbThiscenter");
      std::vector<double> hgrid(3);
      std::vector<int> grid_dim(3), grid_start(3);
      hin.read(hgrid, "hgrid");

      hin.read(grid_dim, "grid_dim");
      hin.read(grid_start, "grid_start");
      double r0[3];
      r0[0] = grid_start[0] * hgrid[0];
      r0[1] = grid_start[1] * hgrid[1];
      r0[2] = grid_start[2] * hgrid[2];

      COT* aoBasis = new COT(num_orb, r0, hgrid.data(), grid_dim.data());
      for(int st = 0; st < num_orb; st++)
      {
          std::string  orb= "orbital_" + std::to_string(st);
          std::string  orb_x= "orbital_x_" + std::to_string(st);
          std::string  orb_y= "orbital_y_" + std::to_string(st);
          std::string  orb_z= "orbital_z_" + std::to_string(st);
          std::string  orb_L= "orbital_L_" + std::to_string(st);
          hin.read(aoBasis->Vgrids[st]  ,   orb);
          hin.read(aoBasis->Vgrids_x[st], orb_x);
          hin.read(aoBasis->Vgrids_y[st], orb_y);
          hin.read(aoBasis->Vgrids_z[st], orb_z);
          hin.read(aoBasis->Vgrids_L[st], orb_L);
      }
      return aoBasis;
  }
};

} // namespace qmcplusplus
#endif
