#ifndef L1PFTAU_H
#define L1PFTAU_H

#include <ostream>
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include <TLorentzVector.h>



class L1PFTau
{



public:

  /// default constructor
  L1PFTau();

  //L1PFTau(const L1PFTau &input){  
  //setEt( input.et() );
  //setPtEtaPhiE( input.p4().Pt(), input.p4().Eta(), input.p4().Phi(), input.p4().Pt());
  //};

  /// destructor
  ~L1PFTau();


  // get/set methods for the data

  /// reset the data content (not position id!)
  void reset() { m_data = 0; }

  /// get raw data
  uint16_t raw() const { return m_data; }

  /// get Et
  unsigned et() const { return (m_et); }
  unsigned towerEta() const{ return (m_towerEta); }
  unsigned towerPhi() const{ return (m_towerPhi); }
  TTTrack< Ref_Phase2TrackerDigi_ > trackRef() const{return (m_trackRef);}

  void setEt(unsigned inputEt) { m_et = inputEt;}
  void setTowerEta(unsigned inputEta ) { (m_towerEta = inputEta); }
  void setTowerEtaSide(unsigned inputEtaSide ) { (m_towerEtaSide = inputEtaSide); }
  void setTowerPhi(unsigned inputPhi ) { (m_towerPhi = inputPhi); }
  void setEoH(unsigned inputEoH) {(m_EoH = inputEoH);}
  void setHoE(float inputHoE){m_HoE = inputHoE;};  

  /// set data
  void setRawData(uint32_t data) { m_data = data; }

  // reco level quantities to be set manually, temporary aid for algo development
  TLorentzVector p4() const {return m_p4;};
  float ecalEnergy() const{return m_ecalEnergy;}
  float hcalEnergy() const{return m_hcalEnergy;}
  float caloEnergy() const{return m_caloEnergy;}
  float HoE() const{return m_HoE;}
  float EoH() const{return m_EoH;}
  float relIso() const{return m_relativeIsolation;}
  float rawIso() const{return m_rawIsolation;}
  int tauType()const{return m_tauType;}

  void setp4(TLorentzVector input) {m_p4 = input;};
  void setPtEtaPhiE(float pt, float eta, float phi, float et){m_p4.SetPtEtaPhiE(pt,eta,phi,et);};
  void setEcalEnergy(float input){ m_ecalEnergy = input;};
  void setHcalEnergy(float input){ m_hcalEnergy = input;};
  void setCaloEnergy(float input){ m_caloEnergy = input;};
  void setTrackRef(TTTrack< Ref_Phase2TrackerDigi_ > trackRef){m_trackRef = trackRef;};
  void setRelIso(float inputIso){m_relativeIsolation = inputIso;};

  void setTauType(float input){ m_tauType = input;};
  void setRawIso(float inputIso){m_rawIsolation = inputIso;};

  /// is there any information in the candidate
  bool empty() const { return (m_data == 0); }
  /*
   * Strip tracks i.e. e/g cands
   */



  /// print to stream
  friend std::ostream& operator << (std::ostream& os, const L1PFTau& reg);

 private:

  uint16_t m_data;
  TLorentzVector m_p4;
  //for temporary use
  float m_ecalEnergy;
  float m_hcalEnergy;
  float m_caloEnergy;
  float m_relativeIsolation;

  unsigned m_towerEta;
  unsigned m_towerEtaSide;
  unsigned m_towerPhi;
  unsigned m_maxCrystalEta;
  unsigned m_maxCrystalPhi;
  unsigned m_et;
  unsigned m_EoH;
  unsigned m_HoE;

  unsigned m_tauType;
  float m_rawIsolation;


  TTTrack< Ref_Phase2TrackerDigi_ > m_trackRef;
};

typedef std::vector<L1PFTau> L1PFTauCollection;

#endif /*L1PFTAU_H*/
