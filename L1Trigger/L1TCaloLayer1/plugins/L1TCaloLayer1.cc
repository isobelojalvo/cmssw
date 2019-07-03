// -*- C++ -*-
//
// Package:    L1Trigger/L1TCaloLayer1
// Class:      L1TCaloLayer1
// 
/**\class L1TCaloLayer1 L1TCaloLayer1.cc L1Trigger/L1TCaloLayer1/plugins/L1TCaloLayer1.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Sridhara Rao Dasu
//         Created:  Thu, 08 Oct 2015 09:20:16 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "L1Trigger/L1TCaloLayer1/src/UCTParameters.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTGeometry.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTLogging.hh"

#include "L1Trigger/L1TCaloLayer1/src/UCTLayer1.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTCrate.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTCard.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTRegion.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTTower.hh"

#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"

#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"

#include "L1Trigger/L1TCaloLayer1/src/L1TCaloLayer1FetchLUTs.hh"

using namespace l1t;
using namespace l1tcalo;

//
// class declaration
//

class L1TCaloLayer1 : public edm::EDProducer {
public:
  explicit L1TCaloLayer1(const edm::ParameterSet&);
  ~L1TCaloLayer1();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
      
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;

  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------

  edm::EDGetTokenT<EcalTrigPrimDigiCollection> ecalTPSource;
  std::string ecalTPSourceLabel;
  edm::EDGetTokenT<HcalTrigPrimDigiCollection> hcalTPSource;
  std::string hcalTPSourceLabel;
  
  std::vector< std::vector< std::vector < uint32_t > > > ecalLUT;
  std::vector< std::vector< std::vector < uint32_t > > > hcalLUT;
  std::vector< std::vector< uint32_t > > hfLUT;

  std::vector< UCTTower* > twrList;

  bool useLSB;
  bool useCalib;
  bool useECALLUT;
  bool useHCALLUT;
  bool useHFLUT;
  bool verbose;
  bool unpackHcalMask;
  bool unpackEcalMask;

  UCTParameters uctParameters;
  UCTLayer1 *layer1;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
L1TCaloLayer1::L1TCaloLayer1(const edm::ParameterSet& iConfig) :
  ecalTPSource(consumes<EcalTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("ecalToken"))),
  ecalTPSourceLabel(iConfig.getParameter<edm::InputTag>("ecalToken").label()),
  hcalTPSource(consumes<HcalTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("hcalToken"))),
  hcalTPSourceLabel(iConfig.getParameter<edm::InputTag>("hcalToken").label()),
  ecalLUT(28, std::vector< std::vector<uint32_t> >(2, std::vector<uint32_t>(256))),
  hcalLUT(28, std::vector< std::vector<uint32_t> >(2, std::vector<uint32_t>(256))),
  hfLUT(12, std::vector < uint32_t >(256)),
  useLSB(iConfig.getParameter<bool>("useLSB")),
  useCalib(iConfig.getParameter<bool>("useCalib")),
  useECALLUT(iConfig.getParameter<bool>("useECALLUT")),
  useHCALLUT(iConfig.getParameter<bool>("useHCALLUT")),
  useHFLUT(iConfig.getParameter<bool>("useHFLUT")),
  verbose(iConfig.getParameter<bool>("verbose")), 
  unpackHcalMask(iConfig.getParameter<bool>("unpackHcalMask")),
  unpackEcalMask(iConfig.getParameter<bool>("unpackEcalMask")),
  uctParameters(iConfig.getParameter<double>("activityFraction"), 
		iConfig.getParameter<double>("ecalActivityFraction"), 
		iConfig.getParameter<double>("miscActivityFraction"))
{
  produces<CaloTowerBxCollection>();
  produces<L1CaloRegionCollection>();
  produces<std::vector<UCTRegion> >();
  layer1 = new UCTLayer1(&uctParameters);
  vector<UCTCrate*> crates = layer1->getCrates();
  for(uint32_t crt = 0; crt < crates.size(); crt++) {
    vector<UCTCard*> cards = crates[crt]->getCards();
    for(uint32_t crd = 0; crd < cards.size(); crd++) {
      vector<UCTRegion*> regions = cards[crd]->getRegions();
      for(uint32_t rgn = 0; rgn < regions.size(); rgn++) {
	vector<UCTTower*> towers = regions[rgn]->getTowers();
	for(uint32_t twr = 0; twr < towers.size(); twr++) {
	  twrList.push_back(towers[twr]);
	}
      }
    }
  }

  // This sort corresponds to the sort condition on
  // the output CaloTowerBxCollection
  std::sort(twrList.begin(), twrList.end(), [](UCTTower* a, UCTTower* b) {
      return CaloTools::caloTowerHash(a->caloEta(), a->caloPhi()) < CaloTools::caloTowerHash(b->caloEta(), b->caloPhi());
      });
}

L1TCaloLayer1::~L1TCaloLayer1() {
  if(layer1 != 0) delete layer1;
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
L1TCaloLayer1::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<EcalTrigPrimDigiCollection> ecalTPs;
  iEvent.getByToken(ecalTPSource, ecalTPs);
  edm::Handle<HcalTrigPrimDigiCollection> hcalTPs;
  iEvent.getByToken(hcalTPSource, hcalTPs);

  std::unique_ptr<CaloTowerBxCollection> towersColl (new CaloTowerBxCollection);
  std::unique_ptr<L1CaloRegionCollection> rgnCollection (new L1CaloRegionCollection);
  //std::unique_ptr<std::vector<UCTRegion> > uctRgnCollection (new std::vector<UCTRegion>);

  uint32_t expectedTotalET = 0;
  if(!layer1->clearEvent()) {
    LOG_ERROR << "UCT: Failed to clear event" << std::endl;
    return;
  }

  for ( const auto& ecalTp : *ecalTPs ) {
    if ( unpackEcalMask && ((ecalTp.sample(0).raw()>>13) & 0x1) ) continue;
    int caloEta = ecalTp.id().ieta();
    int caloPhi = ecalTp.id().iphi();
    int et = ecalTp.compressedEt();
    bool fgVeto = ecalTp.fineGrain();
    UCTTowerIndex t = UCTTowerIndex(caloEta, caloPhi);
    if(!layer1->setECALData(t,fgVeto,et)) {
      LOG_ERROR << "UCT: Failed loading an ECAL tower" << std::endl;
      return;
    }
    expectedTotalET += et;
  }

  for ( const auto& hcalTp : *hcalTPs ) {
    if ( unpackHcalMask && ((hcalTp.sample(0).raw()>>13) & 0x1) ) continue;
    int caloEta = hcalTp.id().ieta();
    uint32_t absCaloEta = abs(caloEta);
    // Tower 29 is not used by Layer-1
    if(absCaloEta == 29) {
      continue;
    }
    // Prevent usage of HF TPs with Layer-1 emulator if HCAL TPs are old style
    else if(hcalTp.id().version() == 0 && absCaloEta > 29) {
      continue;
    }
    else if(absCaloEta <= 41) {
      int caloPhi = hcalTp.id().iphi();
      int et = hcalTp.SOI_compressedEt();
      bool fg  = hcalTp.t0().fineGrain(0);
      bool fg2 = hcalTp.t0().fineGrain(1);
      if(caloPhi <= 72) {
        UCTTowerIndex t = UCTTowerIndex(caloEta, caloPhi);
        uint32_t featureBits = 0;
        if(fg)  featureBits |= 0b01;
        // fg2 should only be set for HF
        if(absCaloEta > 29 && fg2) featureBits |= 0b10;
        if(!layer1->setHCALData(t, featureBits, et)) {
          LOG_ERROR << "caloEta = " << caloEta << "; caloPhi =" << caloPhi << std::endl;
          LOG_ERROR << "UCT: Failed loading an HCAL tower" << std::endl;
          return;
        }
        expectedTotalET += et;
      }
      else {
	LOG_ERROR << "Illegal Tower: caloEta = " << caloEta << "; caloPhi =" << caloPhi << "; et = " << et << std::endl;	
      }
    }
    else {
      LOG_ERROR << "Illegal Tower: caloEta = " << caloEta << std::endl;
    }
  }
  
   //Process
  if(!layer1->process()) {
    LOG_ERROR << "UCT: Failed to process layer 1" << std::endl;
  }

  int theBX = 0; // Currently we only read and process the "hit" BX only

  for(uint32_t twr = 0; twr < twrList.size(); twr++) {
    CaloTower caloTower;
    caloTower.setHwPt(twrList[twr]->et());               // Bits 0-8 of the 16-bit word per the interface protocol document
    caloTower.setHwEtRatio(twrList[twr]->er());          // Bits 9-11 of the 16-bit word per the interface protocol document
    caloTower.setHwQual(twrList[twr]->miscBits());       // Bits 12-15 of the 16-bit word per the interface protocol document
    caloTower.setHwEta(twrList[twr]->caloEta());         // caloEta = 1-28 and 30-41
    caloTower.setHwPhi(twrList[twr]->caloPhi());         // caloPhi = 1-72
    caloTower.setHwEtEm(twrList[twr]->getEcalET());      // This is provided as a courtesy - not available to hardware
    caloTower.setHwEtHad(twrList[twr]->getHcalET());     // This is provided as a courtesy - not available to hardware
    towersColl->push_back(theBX, caloTower);
  }

  iEvent.put(std::move(towersColl), "UCTTowers");

  UCTGeometry g;
  vector<UCTCrate*> crates = layer1->getCrates();
  for(uint32_t crt = 0; crt < crates.size(); crt++) {
    vector<UCTCard*> cards = crates[crt]->getCards();
    for(uint32_t crd = 0; crd < cards.size(); crd++) {
      vector<UCTRegion*> regions = cards[crd]->getRegions();
      for(uint32_t rgn = 0; rgn < regions.size(); rgn++) {
	uint32_t rawData = regions[rgn]->rawData();
	uint32_t regionData = rawData & 0x0000FFFF;
	uint32_t crate = regions[rgn]->getCrate();
	uint32_t card = regions[rgn]->getCard();
	uint32_t region = regions[rgn]->getRegion();
	bool negativeEta = regions[rgn]->isNegativeEta();
	uint32_t rPhi = g.getUCTRegionPhiIndex(crate, card);
	// We want to reuse L1CaloRegion and L1CaloRegionDetID
	// We do not want to change those classes too much
	// We want comparison to legacy for Barrel and Endcap to work transparently
	// Noting that rEta is packed in 5 bits of L1CaloRegionDetID, we have a scheme!
	// We store the Barrel and Endcap regions in the same location as done for RCT
	// HF has changed in the upgrade, 6x2 HF regions instead of 4x2 in case of RCT
	// Note that for upgrade region numbers range 0-6 for Barrel/Endcap and 7-12 for HF
	// So, the scheme used for rEta for upgrade is:
	// rEta= 0- 3 for -HF regions 7-10
	// rEta= 4-10 for -B/E regions 0-6
	// rEta=11-17 for +B/E regions 0-6
	// rEta=18-23 for +HF regions 7-12
	// rEta=30 for -HF region 11
	// rEta=31 for -HF region 12
	uint32_t rEta = 10 - region;
	if(negativeEta && region == 11) rEta = 30;
	if(negativeEta && region == 12) rEta = 31;
	if(!negativeEta) rEta = 11 + region; // Positive eta portion is offset by 11
	rgnCollection->push_back(L1CaloRegion((uint16_t) regionData, (unsigned) rEta, (unsigned) rPhi, (int16_t) 0));
	//uctRgnCollection->push_back(*regions[rgn]);
      }
    }
  }  
  iEvent.put(std::move(rgnCollection),"L1Regions");
  //iEvent.put(std::move(uctRgnCollection),"UCTRegions");

}



// ------------ method called once each job just before starting event loop  ------------
void 
L1TCaloLayer1::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1TCaloLayer1::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void
L1TCaloLayer1::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  if(!L1TCaloLayer1FetchLUTs(iSetup, ecalLUT, hcalLUT, hfLUT, useLSB, useCalib, useECALLUT, useHCALLUT, useHFLUT)) {
    LOG_ERROR << "L1TCaloLayer1::beginRun: failed to fetch LUTS - using unity" << std::endl;
  }
  for(uint32_t twr = 0; twr < twrList.size(); twr++) {
    twrList[twr]->setECALLUT(&ecalLUT);
    twrList[twr]->setHCALLUT(&hcalLUT);
    twrList[twr]->setHFLUT(&hfLUT);
  }
}


// ------------ method called when ending the processing of a run  ------------
/*
  void
  L1TCaloLayer1::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
  void
  L1TCaloLayer1::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void
  L1TCaloLayer1::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TCaloLayer1::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TCaloLayer1);
/* vim: set ts=8 sw=2 tw=0 et :*/
