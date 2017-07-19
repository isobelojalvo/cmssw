

#include "DataFormats/L1CaloTrigger/interface/L1PFObject.h"

using std::ostream;
using std::endl;
using std::hex;
using std::dec;

// default constructor
L1PFObject::L1PFObject() : m_data(0) { }



// destructor
L1PFObject::~L1PFObject() { }


// print to stream
ostream& operator << (ostream& os, const L1PFObject& pf) {
  os << "L1PFObject:";
  os << " Reco -> ET = " <<pf.p4().Pt();
          os <<" Eta = " <<pf.p4().Eta();
	  os <<" Phi = " <<pf.p4().Phi() <<std::endl;

  os << " Et = "      << pf.et();
  os << " towerEta = "<< pf.towerEta();
  os << " towerPhi = "<< pf.towerPhi()<<std::endl;

  os << " ecalEnergy = "<< pf.ecalEnergy();
  os << " hcalEnergy = "<< pf.hcalEnergy();
  os << " caloEnergy = "<< pf.caloEnergy()<<std::endl;

  os << " EoH = "<<pf.EoH();
  os  <<" HoE = "<<pf.HoE()<<std::endl;

  os << " raw = "<<pf.raw()<<std::endl;

  if(pf.isElectron())
    os << " Electron = True "<<std::endl;

  if(pf.isChargedHadron())
    os << " ChargedHadron = True "<<std::endl;

  if(!pf.isElectron() && !pf.isChargedHadron())
    os << " Not Photon or Neutral Hadron"<<std::endl;
  return os;
}

