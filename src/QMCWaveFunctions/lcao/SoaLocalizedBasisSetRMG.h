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


/** @file SoaLocalizedBasisSetRMG.h
 * @brief specialization for RMGBasisSet
 *
 */
#ifndef QMCPLUSPLUS_SOA_LOCALIZEDBASISSET_RMG_H
#define QMCPLUSPLUS_SOA_LOCALIZEDBASISSET_RMG_H

//need a unit test. 
namespace qmcplusplus
{
/** specialization for SoaLocalizedBasisSet<RMGBasisSet<T>, ORBT>
 */
template<typename T, typename ORBT>
struct SoaLocalizedBasisSet<RMGBasisSet<T>, ORBT> : public SoaBasisSetBase<ORBT>
{
  using COT = RMGBasisSet<ORBT>;
  using BaseType = SoaBasisSetBase<ORBT>;
  typedef typename BaseType::vgl_type vgl_type;
  typedef typename BaseType::vgh_type vgh_type;
  typedef typename BaseType::vghgh_type vghgh_type;
  typedef typename ParticleSet::PosType PosType;

  using BaseType::BasisSetSize;

  ///number of centers, e.g., ions
  size_t NumCenters;
  ///number of quantum particles
  size_t NumTargets;
  ///ion particle set
  const ParticleSet& ions_;
  ///number of quantum particles
  const int myTableIndex;
  ///Global Coordinate of Supertwist read from HDF5
  PosType SuperTwist;

  /** container to store the offsets of the basis functions
   *
   * the number of basis states for center J is BasisOffset[J+1]-Basis[J]
   */
  aligned_vector<size_t> BasisOffset;

  /** container of the unique pointers to the Atomic Orbitals
   *
   * size of LOBasisSet = number  of unique centers
   */
  aligned_vector<COT*> LOBasisSet;

  /** constructor
   * @param ions ionic system
   * @param els electronic system
   */
  SoaLocalizedBasisSet(ParticleSet& ions, ParticleSet& els)
      : ions_(ions), myTableIndex(els.addTable(ions, DT_SOA)), SuperTwist(0.0)
  {
    NumCenters = ions.getTotalNum();
    NumTargets = els.getTotalNum();
    LOBasisSet.resize(NumCenters);
    BasisOffset.resize(NumCenters + 1);
    BasisSetSize = 0;
  }

  /** copy constructor */
  SoaLocalizedBasisSet(const SoaLocalizedBasisSet& a) = default;

  /** makeClone */
  //SoaLocalizedBasisSet<COT>* makeClone() const
  BaseType* makeClone() const
  {
    SoaLocalizedBasisSet<COT, ORBT>* myclone = new SoaLocalizedBasisSet<COT, ORBT>(*this);
    for (int i = 0; i < LOBasisSet.size(); ++i)
      myclone->LOBasisSet[i] = LOBasisSet[i]->makeClone();
    return myclone;
  }

  void setBasisSetSize(int nbs)
  {
    const auto& IonID(ions_.GroupID);
    if (BasisSetSize > 0 && nbs == BasisSetSize)
      return;

    //evaluate the total basis dimension and offset for each center
    BasisOffset[0] = 0;
    for (int c = 0; c < NumCenters; c++)
    {
      std::cout << "aaa " << c <<"  " <<BasisOffset[c] << std::endl;
      BasisOffset[c + 1] = BasisOffset[c] + LOBasisSet[c]->getBasisSetSize();
    }
    BasisSetSize = BasisOffset[NumCenters];
  }

  void setPBCParams(const TinyVector<int, 3>& PBCImages,
                    const TinyVector<double, 3> Sup_Twist,
                    const std::vector<QMCTraits::ValueType>& phase_factor)
  {
  }
  /** compute VGL 
   * @param P quantum particleset
   * @param iat active particle
   * @param vgl Matrix(5,BasisSetSize)
   * @param trialMove if true, use getTempDists()/getTempDispls()
   */
  inline void evaluateVGL(const ParticleSet& P, int iat, vgl_type& vgl)
  {
    const auto& IonID(ions_.GroupID);
    const auto& coordR  = P.activeR(iat);
    const auto& d_table = P.getDistTable(myTableIndex);
    const auto& dist    = (P.activePtcl == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
    const auto& displ   = (P.activePtcl == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);

    PosType Tv;
    for (int c = 0; c < NumCenters; c++)
    {
        Tv[0] = (ions_.R[c][0] - coordR[0]) - displ[c][0];
        Tv[1] = (ions_.R[c][1] - coordR[1]) - displ[c][1];
        Tv[2] = (ions_.R[c][2] - coordR[2]) - displ[c][2];
        LOBasisSet[c]->evaluateVGL(P.Lattice, dist[c], displ[c], BasisOffset[c], vgl, Tv);
    }


    // Ye: needs implemenation
  }


  /** compute VGH 
   * @param P quantum particleset
   * @param iat active particle
   * @param vgl Matrix(10,BasisSetSize)
   * @param trialMove if true, use getTempDists()/getTempDispls()
   */
  inline void evaluateVGH(const ParticleSet& P, int iat, vgh_type& vgh)
  {
  }

  /** compute VGHGH 
   * @param P quantum particleset
   * @param iat active particle
   * @param vghgh Matrix(20,BasisSetSize)
   * @param trialMove if true, use getTempDists()/getTempDispls()
   */
  inline void evaluateVGHGH(const ParticleSet& P, int iat, vghgh_type& vghgh)
  {
  }

  /** compute values for the iat-paricle move
   *
   * Always uses getTempDists() and getTempDispls()
   * Tv is a translation vector; In PBC, in order to reduce the number
   * of images that need to be summed over when generating the AO the 
   * nearest image displacement, dr, is used. Tv corresponds to the 
   * translation that takes the 'general displacement' (displacement
   * between ion position and electron position) to the nearest image 
   * displacement. We need to keep track of Tv because it must be add
   * as a phase factor, i.e., exp(i*k*Tv).
   */
  inline void evaluateV(const ParticleSet& P, int iat, ORBT* restrict vals)
  {
    const auto& IonID(ions_.GroupID);
    const auto& coordR  = P.activeR(iat);
    const auto& d_table = P.getDistTable(myTableIndex);
    const auto& dist    = (P.activePtcl == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
    const auto& displ   = (P.activePtcl == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);

    PosType Tv;
    for (int c = 0; c < NumCenters; c++)
    {
        Tv[0] = (ions_.R[c][0] - coordR[0]) - displ[c][0];
        Tv[1] = (ions_.R[c][1] - coordR[1]) - displ[c][1];
        Tv[2] = (ions_.R[c][2] - coordR[2]) - displ[c][2];
        LOBasisSet[c]->evaluateV(P.Lattice, dist[c], displ[c], vals + BasisOffset[c], Tv);
    }
      // Ye: needs implemenation
  }

  inline void evaluateGradSourceV(const ParticleSet& P, int iat, const ParticleSet& ions, int jion, vgl_type& vgl)
  {
  }

  inline void evaluateGradSourceVGL(const ParticleSet& P, int iat, const ParticleSet& ions, int jion, vghgh_type& vghgh)
  {
  }

  /** add a new set of Centered Atomic Orbitals
   * @param icenter the index of the center
   * @param aos a set of Centered Atomic Orbitals
   */
  void add(int icenter, COT* aos) { LOBasisSet[icenter] = aos; }
};
} // namespace qmcplusplus
#endif
