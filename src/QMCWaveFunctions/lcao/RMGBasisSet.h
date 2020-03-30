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

/** @file RMGBasisSet.h
 */
#ifndef QMCPLUSPLUS_RMG_BASISSET_H
#define QMCPLUSPLUS_RMG_BASISSET_H
//need a unit test
namespace qmcplusplus
{
/* A basis set for RMG localized orbitals
   *
   * @tparam VALT : value type
   */
template<typename VALT>
class RMGBasisSet
{
//all orbitals on the same center.

private:
  ///size of the basis set
// numnumber of borbitals. 
  int BasisSetSize;
  VALT r0[3]; //left-bottom corner of the orbitals
  VALT r1[3]; //right up corner of the orbitals
  VALT hgrid[3]; // gridspacing in bohr
  int num_grid[3]; // number of grid in x, y, z director
public:
  std::vector<std::vector<VALT>> Vgrids;
  std::vector<std::vector<VALT>> Vgrids_x;
  std::vector<std::vector<VALT>> Vgrids_y;
  std::vector<std::vector<VALT>> Vgrids_z;
  std::vector<std::vector<VALT>> Vgrids_L;

  RMGBasisSet<VALT>* makeClone() const
  {
    return new RMGBasisSet(*this);
  }

  RMGBasisSet<VALT>(int num_orbitals_thiscenter, double r0_in[3], double hgrid_in[3], int num_grid_in[3])
  {
    for(int i = 0; i < 3; i++)
    {
      r0[i] = r0_in[i];
      r1[i] = r0_in[i] + num_grid_in[i] * hgrid_in[i];
      hgrid[i] = hgrid_in[i];
      num_grid[i] = num_grid_in[i]; 
    }


    BasisSetSize = num_orbitals_thiscenter;
    Vgrids.resize(BasisSetSize);
    Vgrids_x.resize(BasisSetSize);
    Vgrids_y.resize(BasisSetSize);
    Vgrids_z.resize(BasisSetSize);
    Vgrids_L.resize(BasisSetSize);

    for(int ib = 0; ib < BasisSetSize; ib++)
    {
        Vgrids[ib].resize  (num_grid[0] * num_grid[1] * num_grid[2]);
        Vgrids_x[ib].resize(num_grid[0] * num_grid[1] * num_grid[2]);
        Vgrids_y[ib].resize(num_grid[0] * num_grid[1] * num_grid[2]);
        Vgrids_z[ib].resize(num_grid[0] * num_grid[1] * num_grid[2]);
        Vgrids_L[ib].resize(num_grid[0] * num_grid[1] * num_grid[2]);
    }

  }
  inline int getBasisSetSize() const
  {
      return BasisSetSize;
  }

  void checkInVariables(opt_variables_type& active)
  {
  }

  void checkOutVariables(const opt_variables_type& active)
  {
  }

  void resetParameters(const opt_variables_type& active)
  {
  }

  /** evaluate VGL

   */

  template<typename LAT, typename T, typename PosType, typename VGL>
      inline void evaluateVGL(const LAT& lattice, const T r, const PosType& dr, const size_t offset, VGL& vgl,PosType Tv)
      {

          //   values, gradients, Laplacian

          double Val[3], disp[3]{0.0, 0.0, 0.0};
          int grid_pos[3]{-1,-1,-1};
         //check where the real space point dr map to the orbital's real space grid and displacment.
         /*   int i(0), j(0), k(0);
                  {
                      Val[0] = dr[0] + i * lattice.R(0, 0) + j * lattice.R(1, 0) + k * lattice.R(2, 0) ;
                      Val[1] = dr[1] + i * lattice.R(0, 1) + j * lattice.R(1, 1) + k * lattice.R(2, 1) ;
                      Val[2] = dr[2] + i * lattice.R(0, 2) + j * lattice.R(1, 2) + k * lattice.R(2, 2) ;
                      if(Val[0] > r0[0] && Val[0] < r1[0] && Val[1] > r0[1] && Val[1] < r1[1] && Val[2] > r0[2] && Val[2] < r1[2])
                      {
                          grid_pos[0] = (Val[0]-r0[0]) / hgrid[0];
                          disp[0] = Val[0]-r0[0] - grid_pos[0] * hgrid[0];
                          grid_pos[1] = (Val[1]-r0[1]) / hgrid[1];
                          disp[1] = Val[1]-r0[1] - grid_pos[1] * hgrid[1];
                          grid_pos[2] = (Val[2]-r0[2]) / hgrid[2];
                          disp[2] = Val[2]-r0[2] - grid_pos[2] * hgrid[2];
                       //   break;
                      }
                  }
o*/
          grid_pos[0] = (dr[0]-r0[0]) / hgrid[0];
          disp[0] = dr[0]-r0[0] - grid_pos[0] * hgrid[0];
          grid_pos[1] = (dr[1]-r0[1]) / hgrid[1];
          disp[1] = dr[1]-r0[1] - grid_pos[1] * hgrid[1];
          grid_pos[2] = (dr[2]-r0[2]) / hgrid[2];
          disp[2] = dr[2]-r0[2] - grid_pos[2] * hgrid[2];
          // set oiuters to vgl datatype

          auto* restrict phi      = vgl.data(0) + offset;
          auto* restrict dphi_x   = vgl.data(1) + offset;
          auto* restrict dphi_y   = vgl.data(2) + offset;
          auto* restrict dphi_z   = vgl.data(3) + offset;
          auto* restrict d2phi    = vgl.data(4) + offset;

          for(int ib = 0; ib < BasisSetSize; ib++)
          {
              phi[ib] = 0.0;
              dphi_x[ib] = 0.0;
              dphi_y[ib] = 0.0;
              dphi_z[ib] = 0.0;
              d2phi[ib] = 0.0;
          }

          int incx = num_grid[1] * num_grid[2];

          if(grid_pos[0] > 1 && grid_pos[1] > 1 && grid_pos[2] > 1 && 
                  grid_pos[0] < num_grid[0] - 3 && grid_pos[1] < num_grid[1] - 3 &&
                  grid_pos[2] < num_grid[2] - 3)
          {
              int idx000 = grid_pos[0] * incx + grid_pos[1] * num_grid[2] + grid_pos[2];
              int idx001 = grid_pos[0] * incx + grid_pos[1] * num_grid[2] + grid_pos[2] +1;
              int idx010 = grid_pos[0] * incx + (grid_pos[1]+1) * num_grid[2] + grid_pos[2];
              int idx100 = (grid_pos[0]+1) * incx + grid_pos[1] * num_grid[2] + grid_pos[2];
              int idx101 = (grid_pos[0]+1) * incx + grid_pos[1] * num_grid[2] + grid_pos[2] +1;
              int idx011 = grid_pos[0] * incx + (grid_pos[1]+1) * num_grid[2] + grid_pos[2] +1;
              int idx110 = (grid_pos[0]+1) * incx + (grid_pos[1]+1) * num_grid[2] + grid_pos[2];
              int idx111 = (grid_pos[0]+1) * incx + (grid_pos[1]+1) * num_grid[2] + grid_pos[2]+1;
              double coef000 = (1.0 - disp[0]) * (1.0 - disp[1]) * (1.0 - disp[2]);
              double coef001 = (1.0 - disp[0]) * (1.0 - disp[1]) * disp[2];
              double coef010 = (1.0 - disp[0]) * disp[1] * (1.0 - disp[2]);
              double coef100 = disp[0] * (1.0 - disp[1]) * (1.0 - disp[2]);
              double coef101 = disp[0] * (1.0 - disp[1]) * disp[2];
              double coef011 = (1.0 - disp[0]) * disp[1] * disp[2];
              double coef110 = disp[0] * disp[1] * (1.0 - disp[2]);
              double coef111 = disp[0] * disp[1] * disp[2];
              for(int ib = 0; ib < BasisSetSize; ib++)
              {
                  phi[ib] = Vgrids[ib][idx000] * coef000 
                      + Vgrids[ib][idx001] * coef001 
                      + Vgrids[ib][idx010] * coef010 
                      + Vgrids[ib][idx100] * coef100 
                      + Vgrids[ib][idx101] * coef101 
                      + Vgrids[ib][idx011] * coef011 
                      + Vgrids[ib][idx110] * coef110 
                      + Vgrids[ib][idx111] * coef111;
                  dphi_x[ib] = Vgrids_x[ib][idx000] * coef000 
                      + Vgrids_x[ib][idx001] * coef001 
                      + Vgrids_x[ib][idx010] * coef010 
                      + Vgrids_x[ib][idx100] * coef100 
                      + Vgrids_x[ib][idx101] * coef101 
                      + Vgrids_x[ib][idx011] * coef011 
                      + Vgrids_x[ib][idx110] * coef110 
                      + Vgrids_x[ib][idx111] * coef111;
                  dphi_y[ib] = Vgrids_y[ib][idx000] * coef000 
                      + Vgrids_y[ib][idx001] * coef001 
                      + Vgrids_y[ib][idx010] * coef010 
                      + Vgrids_y[ib][idx100] * coef100 
                      + Vgrids_y[ib][idx101] * coef101 
                      + Vgrids_y[ib][idx011] * coef011 
                      + Vgrids_y[ib][idx110] * coef110 
                      + Vgrids_y[ib][idx111] * coef111;
                  dphi_z[ib] = Vgrids_z[ib][idx000] * coef000 
                      + Vgrids_z[ib][idx001] * coef001 
                      + Vgrids_z[ib][idx010] * coef010 
                      + Vgrids_z[ib][idx100] * coef100 
                      + Vgrids_z[ib][idx101] * coef101 
                      + Vgrids_z[ib][idx011] * coef011 
                      + Vgrids_z[ib][idx110] * coef110 
                      + Vgrids_z[ib][idx111] * coef111;
                  d2phi[ib] = Vgrids_L[ib][idx000] * coef000 
                      + Vgrids_L[ib][idx001] * coef001 
                      + Vgrids_L[ib][idx010] * coef010 
                      + Vgrids_L[ib][idx100] * coef100 
                      + Vgrids_L[ib][idx101] * coef101 
                      + Vgrids_L[ib][idx011] * coef011 
                      + Vgrids_L[ib][idx110] * coef110 
                      + Vgrids_L[ib][idx111] * coef111;
                  //      interpolate_value(grid_pos, disp, Vgrids[ib], phi[ib]); 
                  //      interpolate_value(grid_pos, disp, Vgrids_x[ib], dphi_x[ib]); 
                  //     interpolate_value(grid_pos, disp, Vgrids_y[ib], dphi_y[ib]); 
                  //    interpolate_value(grid_pos, disp, Vgrids_z[ib], dphi_z[ib]); 
                  //   interpolate_value(grid_pos, disp, Vgrids_L[ib], d2phi[ib]); 
              }


          }


          // must be implemented. interface arguments needs change
      }

  template<typename LAT, typename T, typename PosType, typename VGH>
      inline void evaluateVGH(const LAT& lattice, const T r, const PosType& dr, const size_t offset, VGH& vgh)
      {
      }

  template<typename LAT, typename T, typename PosType, typename VGHGH>
      inline void evaluateVGHGH(const LAT& lattice, const T r, const PosType& dr, const size_t offset, VGHGH& vghgh)
      {
      }

  /** evaluate V
   */
  template<typename LAT, typename T, typename PosType, typename VT>
      inline void evaluateV(const LAT& lattice, const T r, const PosType& dr, VT* restrict phi,PosType Tv)
      {
          // must be implemented. interface arguments needs change
          double Val[3], disp[3]{0.0, 0.0, 0.0};
          int grid_pos[3]{-1,-1,-1};
          //check where the real space point dr map to the orbital's real space grid and displacment.
          //for(int i = -1; i <= 1; i++)
          //    for(int j = -1; j <= 1; j++)
          //        for(int k = -1; k <= 1; k++)
          /*
             int i(0), j(0), k(0);
             {
          //     Val[0] = dr[0] + i * lattice.R(0, 0) + j * lattice.R(1, 0) + k * lattice.R(2, 0) ;
          //     Val[1] = dr[1] + i * lattice.R(0, 1) + j * lattice.R(1, 1) + k * lattice.R(2, 1) ;
          //     Val[2] = dr[2] + i * lattice.R(0, 2) + j * lattice.R(1, 2) + k * lattice.R(2, 2) ;
          // if(Val[0] > r0[0] && Val[0] < r1[0] && Val[1] > r0[1] && Val[1] < r1[1] && Val[2] > r0[2] && Val[2] < r1[2])
          {
          grid_pos[0] = (dr[0]-r0[0]) / hgrid[0];
          //     disp[0] = Val[0]-r0[0] - grid_pos[0] * hgrid[0];
          grid_pos[1] = (dr[1]-r0[1]) / hgrid[1];
          //     disp[1] = Val[1]-r0[1] - grid_pos[1] * hgrid[1];
          grid_pos[2] = (dr[2]-r0[2]) / hgrid[2];
          //     disp[2] = Val[2]-r0[2] - grid_pos[2] * hgrid[2];
          //   break;
          }
          }
           */
          grid_pos[0] = (dr[0]-r0[0]) / hgrid[0];
          disp[0] = dr[0]-r0[0] - grid_pos[0] * hgrid[0];
          grid_pos[1] = (dr[1]-r0[1]) / hgrid[1];
          disp[1] = dr[1]-r0[1] - grid_pos[1] * hgrid[1];
          grid_pos[2] = (dr[2]-r0[2]) / hgrid[2];
          disp[2] = dr[2]-r0[2] - grid_pos[2] * hgrid[2];

          //std::cout<<  grid_pos[0] << " " << grid_pos[1] << " "<< grid_pos[2]<< " " <<dr[0] << " " <<dr[1] <<" " <<dr[2]<< " "<<disp[0]<<std::endl;
          // set oiuters to vgl datatype

          for(int ib = 0; ib < BasisSetSize; ib++)
          {
              phi[ib] = 0.0;
          }


          int incx = num_grid[1] * num_grid[2];
          if(grid_pos[0] > 1 && grid_pos[1] > 1 && grid_pos[2] > 1 && 
                  grid_pos[0] < num_grid[0] - 3 && grid_pos[1] < num_grid[1] - 3 &&
                  grid_pos[2] < num_grid[2] - 3)
          {
              int idx000 = grid_pos[0] * incx + grid_pos[1] * num_grid[2] + grid_pos[2];
              int idx001 = grid_pos[0] * incx + grid_pos[1] * num_grid[2] + grid_pos[2] +1;
              int idx010 = grid_pos[0] * incx + (grid_pos[1]+1) * num_grid[2] + grid_pos[2];
              int idx100 = (grid_pos[0]+1) * incx + grid_pos[1] * num_grid[2] + grid_pos[2];
              int idx101 = (grid_pos[0]+1) * incx + grid_pos[1] * num_grid[2] + grid_pos[2] +1;
              int idx011 = grid_pos[0] * incx + (grid_pos[1]+1) * num_grid[2] + grid_pos[2] +1;
              int idx110 = (grid_pos[0]+1) * incx + (grid_pos[1]+1) * num_grid[2] + grid_pos[2];
              int idx111 = (grid_pos[0]+1) * incx + (grid_pos[1]+1) * num_grid[2] + grid_pos[2]+1;
              double coef000 = (1.0 - disp[0]) * (1.0 - disp[1]) * (1.0 - disp[2]);
              double coef001 = (1.0 - disp[0]) * (1.0 - disp[1]) * disp[2];
              double coef010 = (1.0 - disp[0]) * disp[1] * (1.0 - disp[2]);
              double coef100 = disp[0] * (1.0 - disp[1]) * (1.0 - disp[2]);
              double coef101 = disp[0] * (1.0 - disp[1]) * disp[2];
              double coef011 = (1.0 - disp[0]) * disp[1] * disp[2];
              double coef110 = disp[0] * disp[1] * (1.0 - disp[2]);
              double coef111 = disp[0] * disp[1] * disp[2];
              for(int ib = 0; ib < BasisSetSize; ib++)
              {
                  phi[ib] = Vgrids[ib][idx000] * coef000 
                      + Vgrids[ib][idx001] * coef001 
                      + Vgrids[ib][idx010] * coef010 
                      + Vgrids[ib][idx100] * coef100 
                      + Vgrids[ib][idx101] * coef101 
                      + Vgrids[ib][idx011] * coef011 
                      + Vgrids[ib][idx110] * coef110 
                      + Vgrids[ib][idx111] * coef111;
              }


          }


          // must be implemented. interface arguments needs change
      }

  inline void interpolate_value(int grid_pos[3], double disp[3], std::vector<VALT> vgrids, VALT &val)
  {
      int idx = grid_pos[0] * num_grid[1] * num_grid[2] + grid_pos[1] * num_grid[2] + grid_pos[2];
      val = vgrids[idx];
      //   return;
      //   // at the boundary all values are zero
      //   if(grid_pos[0] < 2 || grid_pos[0] > num_grid[0] -3 ) val =  0.0;
      //   if(grid_pos[1] < 2 || grid_pos[1] > num_grid[1] -3 ) val =  0.0;
      //   if(grid_pos[2] < 2 || grid_pos[2] > num_grid[2] -3 ) val =  0.0;

      //   // cubic interpolation
      //   double cc[4][3];
      //   for(int i = 0; i < 3; i++)
      //   {
      //       double frac = disp[i]/hgrid[i];
      //       cc[0][i] = -frac * (1.0 - frac) * (2.0 - frac) / 6.0;
      //       cc[1][i] = (1.0 + frac) * (1.0 - frac) * (2.0 - frac) / 2.0;
      //       cc[2][i] = (1.0 + frac) * frac * (2.0 - frac) / 2.0;
      //       cc[3][i] = -(1.0 + frac) * frac * (1.0 - frac) / 6.0;
      //   }


      //   VALT result = 0.0;
      //   for(int i = 0; i < 4; i++)
      //       for(int j = 0; j < 4; j++)
      //           for(int k = 0; k < 4; k++)
      //           {
      //               int ix = grid_pos[0] +i -1;
      //               int iy = grid_pos[1] +j -1;
      //               int iz = grid_pos[2] +k -1;
      //               int idx = ix * num_grid[1] * num_grid[2] + iy * num_grid[2] + iz;
      //               double coeff = cc[i][0] * cc[j][1] * cc[k][2];
      //               result += coeff * vgrids[idx];
      //           }

      //   val = result;

  }
};

} // namespace qmcplusplus
#endif
