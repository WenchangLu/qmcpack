#ifndef QMCPLUSPLUS_AFQMC_ENERGYESTIMATOR_H
#define QMCPLUSPLUS_AFQMC_ENERGYESTIMATOR_H

#include<Message/MPIObjectBase.h>
#include"AFQMC/config.h"
#include<vector>
#include<queue>
#include<string>
#include<iostream>
#include<fstream>

#include "io/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"
#include "Utilities/NewTimer.h"

#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/Walkers/WalkerConfig.hpp"

namespace qmcplusplus
{
namespace afqmc
{

class EnergyEstimator: public EstimatorBase
{

  public:

  EnergyEstimator(afqmc::TaskGroup_& tg_, AFQMCInfo info, xmlNodePtr cur,
        Wavefunction& wfn, bool impsamp_=true, bool timer=true):
            EstimatorBase(info),TG(tg_),wfn0(wfn),importanceSampling(impsamp_)
  {

    data.resize(2);
  }

  ~EnergyEstimator() {}

  void accumulate_step(WalkerSet& wlks, std::vector<ComplexType>& curData) {}

  void accumulate_block(WalkerSet& wset)
  {
    AFQMCTimers[energy_timer]->start();
    size_t nwalk = wset.size();
    if(eloc.size(0) != nwalk || eloc.size(1) != 3)
      eloc.reextent({nwalk,3});
    if(ovlp.size(0) != nwalk)
      ovlp.reextent(iextensions<1u>{nwalk});
    if(wprop.size(0) != 4 || wprop.size(1) != nwalk)
      wprop.reextent({4,nwalk});

    ComplexType dum, et;
    wfn0.Energy(wset,eloc,ovlp);
    // in case GPU 
    ComplexMatrix<std::allocator<ComplexType>> eloc_(eloc);
    ComplexVector<std::allocator<ComplexType>> ovlp_(ovlp);
    if(TG.TG_local().root()) {
      wset.getProperty(WEIGHT,wprop[0]);  
      wset.getProperty(OVLP,wprop[1]);  
      wset.getProperty(PHASE,wprop[2]);  
      data[0] = data[1] = std::complex<double>(0,0);
      for(int i=0; i<nwalk; i++) {
        if(std::isnan(real(wprop[0][i]))) continue;
        if(importanceSampling) {
          dum = (wprop[0][i])*ovlp_[i]/(wprop[1][i]);
        } else {
          dum = (wprop[0][i])*ovlp_[i]*(wprop[2][i]);
        }
        et = eloc_[i][0]+eloc_[i][1]+eloc_[i][2];
        if( (!std::isfinite(real(dum))) || (!std::isfinite(real(et*dum))) ) continue;
        data[1] += dum;
        data[0] += et*dum;
      }
      TG.TG_heads().all_reduce_in_place_n(data.begin(),data.size(),std::plus<>());
    }
    AFQMCTimers[energy_timer]->stop();

  }

  void tags(std::ofstream& out)
  {
    if(TG.Global().root()) {
      out<<"EnergyEstim_" <<name <<"_nume_real  EnergyEstim_" <<name <<"_nume_imag "
         <<"EnergyEstim_" <<name <<"_deno_real  EnergyEstim_" <<name <<"_deno_imag "
         <<"EnergyEstim_" <<name <<"_timer ";
    }
  }

  void print(std::ofstream& out,WalkerSet& wset)
  {
    if(TG.Global().root()) {
     int n = wset.get_global_target_population();
      out<< data[0].real()/n << " " << data[0].imag()/n << " "
         << data[1].real()/n << " " << data[1].imag()/n << " "
         <<AFQMCTimers[energy_timer]->get_total() <<" ";
      AFQMCTimers[energy_timer]->reset();
    }
  }

  private:

  std::string name;

  TaskGroup_& TG;

  Wavefunction& wfn0;

  ComplexMatrix<device_allocator<ComplexType>> eloc;
  ComplexVector<device_allocator<ComplexType>> ovlp;
  ComplexMatrix<std::allocator<ComplexType>> wprop;

  std::vector<std::complex<double> > data;

  bool importanceSampling;

};
}
}

#endif
