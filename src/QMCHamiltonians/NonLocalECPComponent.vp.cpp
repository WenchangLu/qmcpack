//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim and Simone Chiesa
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/NonLocalECPComponent.h"

namespace qmcplusplus
{

/** evaluate the non-local potential of the iat-th ionic center
 * @param W electron configuration
 * @param iat ionic index
 * @param psi trial wavefunction
 * @param return the non-local component
 *
 * Currently, we assume that the ratio-only evaluation does not change the state
 * of the trial wavefunction and do not call psi.rejectMove(ieL).
 */
NonLocalECPComponent::RealType
NonLocalECPComponent::evaluateVP(const ParticleSet& W, int iat, TrialWaveFunction& psi)
{
  RealType esum=0.0;
  RealType pairpot;
  ParticleSet::ParticlePos_t deltarV(nknot);
  const DistanceTableData* myTable = W.DistTables[myTableIndex];

  for(int nn=myTable->M[iat],iel=0; nn<myTable->M[iat+1]; nn++,iel++)
  {
    register RealType r(myTable->r(nn));
    if(r>Rmax) continue;
    register RealType rinv(myTable->rinv(nn));
    register PosType  dr(myTable->dr(nn));

    // Compute ratio of wave functions
    for (int j=0; j < nknot ;j++) deltarV[j]=r*rrotsgrid_m[j]-dr;
    VP->makeMoves(iel,deltarV);
    psi.evaluateRatios(*VP,psiratio);
    for(int j=0; j<nknot; ++j) psiratio[j]*=sgridweight_m[j];

    for(int ip=0; ip< nchannel; ip++) vrad[ip]=nlpp_m[ip]->splint(r)*wgt_angpp_m[ip];

    // Compute spherical harmonics on grid
    for (int j=0, jl=0; j<nknot ; j++)
    {
      RealType zz=dot(dr,rrotsgrid_m[j])*rinv;
      // Forming the Legendre polynomials
      lpol[0]=1.0;
      RealType lpolprev=0.0;
      for (int l=0 ; l< lmax ; l++)
      {
        //Not a big difference
        //lpol[l+1]=(2*l+1)*zz*lpol[l]-l*lpolprev;
        //lpol[l+1]/=(l+1);
        lpol[l+1]=Lfactor1[l]*zz*lpol[l]-l*lpolprev;
        lpol[l+1]*=Lfactor2[l];
        lpolprev=lpol[l];
      }
      for(int l=0; l <nchannel; l++,jl++)
        Amat[jl]=lpol[ angpp_m[l] ];
    }
    if(nchannel==1)
    {
      pairpot = vrad[0]*BLAS::dot(nknot, &Amat[0],&psiratio[0]);
    }
    else
    {
      BLAS::gemv(nknot, nchannel, &Amat[0], &psiratio[0], &wvec[0]);
      pairpot = BLAS::dot(nchannel, &vrad[0], &wvec[0]);
    }
    if(tracing_particle_quantities)
    {
      (*Vi_sample)(iat) += .5*pairpot;
      (*Ve_sample)(iel) += .5*pairpot;
    }
    esum += pairpot;
  }   /* end loop over electron */
  return esum;
}


NonLocalECPComponent::RealType
NonLocalECPComponent::evaluateVP(const ParticleSet& W, int iat, TrialWaveFunction& psi,vector<NonLocalData>& Txy)
{
  const DistanceTableData* myTable = W.DistTables[myTableIndex];
  RealType esum=0.0;
  ParticleSet::ParticlePos_t deltarV(nknot);
  for(int nn=myTable->M[iat],iel=0; nn<myTable->M[iat+1]; nn++,iel++)
  {
    register RealType r(myTable->r(nn));
    if(r>Rmax) continue;
    register RealType rinv(myTable->rinv(nn));
    register PosType  dr(myTable->dr(nn));

    for (int j=0; j < nknot ;j++) deltarV[j]=r*rrotsgrid_m[j]-dr;
    VP->makeMoves(iel,deltarV);
    psi.evaluateRatios(*VP,psiratio);
    for(int j=0; j<nknot; ++j) psiratio[j]*=sgridweight_m[j];

    // Compute radial potential
    for(int ip=0; ip< nchannel; ip++)
      vrad[ip]=nlpp_m[ip]->splint(r)*wgt_angpp_m[ip];

    RealType pairpot=0; 
    // Compute spherical harmonics on grid
    for (int j=0, jl=0; j<nknot ; j++)
    {
      RealType zz=dot(dr,rrotsgrid_m[j])*rinv;
      // Forming the Legendre polynomials
      lpol[0]=1.0;
      RealType lpolprev=0.0;
      for (int l=0 ; l< lmax ; l++)
      {
        //Not a big difference
        //lpol[l+1]=(2*l+1)*zz*lpol[l]-l*lpolprev;
        //lpol[l+1]/=(l+1);
        lpol[l+1]=Lfactor1[l]*zz*lpol[l]-l*lpolprev;
        lpol[l+1]*=Lfactor2[l];
        lpolprev=lpol[l];
      }

      RealType lsum=0;
      for(int l=0; l <nchannel; l++)
        lsum += vrad[l]*lpol[ angpp_m[l] ];
      lsum *= psiratio[j];
      Txy.push_back(NonLocalData(iel,lsum,deltarV[j]));
      pairpot+=lsum;
    }

    if(tracing_particle_quantities)
    {
      (*Vi_sample)(iat) += .5*pairpot;
      (*Ve_sample)(iel) += .5*pairpot;
    }
    esum += pairpot;
  }   /* end loop over electron */
  return esum;
}


}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 6181 $   $Date: 2014-01-22 11:58:36 -0500 (Wed, 22 Jan 2014) $
 * $Id: NonLocalECPComponent.cpp 6181 2014-01-22 16:58:36Z jnkim $
 ***************************************************************************/
