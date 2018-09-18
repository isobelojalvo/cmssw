

#include "DataFormats/CaloClusters/interface/L1CaloCluster.h"

using std::ostream;
using std::endl;
using std::hex;
using std::dec;

typedef std::vector<L1CaloCluster> L1CaloClusterCollection;

// default constructor
L1CaloCluster::L1CaloCluster() : m_data(0) { }



// destructor
L1CaloCluster::~L1CaloCluster() { }


// print to stream
ostream& operator << (ostream& os, const L1CaloCluster& clus) {
  os << "L1CaloCluster:";
  os << " Reco -> ET = " <<clus.p4().Pt();
          os <<" Eta = " <<clus.p4().Eta();
	  os <<" Phi = " <<clus.p4().Phi() <<std::endl;

  os << " Et = "      << clus.et();
  os << " towerEta = "<< clus.towerEta();
  os << " towerPhi = "<< clus.towerPhi()<<std::endl;

  os << " ecalEnergy = "<< clus.ecalEnergy();
  os << " hcalEnergy = "<< clus.hcalEnergy();
  os << " caloEnergy = "<< clus.caloEnergy()<<std::endl;

  os << " EoH = "<<clus.EoH();
  os << " HoE = "<<clus.HoE()<<std::endl;

  os << " raw = "<<clus.raw()<<std::endl;

  if(clus.isPhoton())
    os << " Photon = True "<<std::endl;
  else
    os << " Photon = False "<<std::endl;

  if(clus.isPi0())
    os << " Pi0 = True "<<std::endl;
  else
    os << " Pi0 = False "<<std::endl;

  if(clus.isNeutralHadron())
    os << " NeutralHadron = True "<<std::endl;
  else
    os << " NeutralHadron = False "<<std::endl;

  if(!clus.isPhoton() && !clus.isNeutralHadron())
    os << " Not Photon or Neutral Hadron"<<std::endl;
  return os;

}

