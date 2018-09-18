#ifndef L1CALOCLUSTER_H
#define L1CALOCLUSTER_H

#include <ostream>
#include "DataFormats/Candidate/interface/LeafCandidate.h"
//#include <TLorentzVector.h>

/*
 * 28 bit number
 * |iphi
 * |         |ieta sign
 * |         | 
 * |         ||ieta
 * |         ||         |et
 * 0000 0000 0000 0000 0000 0000 0000
 */



class L1CaloCluster :  public reco::LeafCandidate
{
public:

  /// default constructor
  L1CaloCluster();

  /// destructor
  ~L1CaloCluster();

  // get/set methods for the data

  /// reset the data content (not position id!)
  void reset() { m_data = 0; }

  /// get raw data
  uint32_t raw() const { return m_data; }

  /// get Et
  //unsigned et() const { return (m_et); }
  int towerEtaSide() const{ return (m_towerEtaSide); }
  unsigned towerEta() const{ return (m_towerEta); }
  unsigned towerPhi() const{ return (m_towerPhi); }
  unsigned crystalEta() const{ return (m_crystalEta); }
  unsigned crystalPhi() const{ return (m_crystalPhi); }

  //only initialized for neutral clusters
  bool isPhoton()  const{ return (m_isPhoton); }
  bool isPi0()  const{ return (m_isPi0); }
  bool isNeutralHadron()  const{ return (m_isNeutralHadron); }

  void setEt(unsigned inputEt) { 
    m_et = inputEt; }

  void setEt(float inputEt) { 
    if(inputEt>200)
      inputEt = 200;
    uint32_t uInputEt = (uint32_t) inputEt;
    uInputEt = uInputEt*10;
    m_et = (unsigned) uInputEt;
    m_data = m_data&0xFFFFF800;//first clear bottom rows
    m_data = m_data|(uInputEt&0x7FF);
  }

  void setCrystalEta(unsigned inputEta ) { 
    m_crystalEta = inputEta;
  }

  void setCrystalPhi(unsigned inputPhi ) { 
    m_crystalPhi = inputPhi;
  }

  void setTowerEta(unsigned inputEta ) { 
    (m_towerEta = inputEta); 
    m_data = m_data&0xFFF807FF;//first clear eta values
    m_data = m_data|((inputEta&0x7F)<<11);
  }
  void setTowerEtaSide(int inputEtaSide ) { 
    (m_towerEtaSide = inputEtaSide); 
    m_data = m_data&0xFFF7FFFF;//first clear eta sign value
    m_data = m_data|((inputEtaSide&0x1)<<19);    
}
  void setTowerPhi(unsigned inputPhi ) { 
    (m_towerPhi = inputPhi); 
    m_data = m_data&0x00FFFFFF;//first clear phi values
    m_data = m_data|((inputPhi&0xFF)<<20);

  }
  void setEoH(unsigned inputEoH) {(m_EoH = inputEoH);}
  void setHoE(unsigned inputHoE) {(m_HoE = inputHoE);}

  /// set data
  void setRawData(uint32_t data) { m_data = data; }


  // reco level quantities to be set manually, temporary aid for algo development
  //LorentzVector p4() const {return m_p4;};
  float ecalEnergy() const{return m_ecalEnergy;}
  float hcalEnergy() const{return m_hcalEnergy;}
  float caloEnergy() const{return (m_ecalEnergy + m_hcalEnergy);}
  int EoH() const{return (m_EoH);}
  int HoE() const{return (m_HoE);}

  //void setp4(LorentzVector input) {m_p4 = input;};
  //void setPtEtaPhiE(float pt, float eta, float phi, float et){m_p4.SetPtEtaPhiE(pt,eta,phi,et);};
  void setEcalEnergy(float input){ m_ecalEnergy = input;};
  void setHcalEnergy(float input){ m_hcalEnergy = input;};
  void setIsPhoton(bool input){ m_isPhoton = input;};
  void setIsPi0(bool input){ m_isPi0 = input;};
  void setIsNeutralHadron(bool input){ m_isNeutralHadron = input;};

  /// is there any information in the candidate
  bool empty() const { return (m_data == 0); }

  /// print to stream
  friend std::ostream& operator << (std::ostream& os, const L1CaloCluster& reg);

 private:

  uint32_t m_data;
  //TLorentzVector m_p4;
  //for temporary use
  float m_ecalEnergy;
  float m_hcalEnergy;

  int m_crystalEta;
  int m_crystalPhi;
  unsigned m_towerEta;
  int m_towerEtaSide;
  unsigned m_towerPhi;
  unsigned m_maxCrystalEta;
  unsigned m_maxCrystalPhi;
  unsigned m_et;
  unsigned m_EoH;
  unsigned m_HoE;
  bool m_isPhoton;
  bool m_isPi0;
  bool m_isNeutralHadron;

};

typedef std::vector<L1CaloCluster> L1CaloClusterCollection;

#endif /*L1CALOCLUSTER_H*/
