// -*- C++ -*-
//
// Package:    L1Trigger/mtdTPGenerator
// Class:      mtdTPGenerator
// 
/**\class mtdTPGenerator mtdTPGenerator.cc L1Trigger/mtdTPGenerator/plugins/mtdTPGenerator.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Isobel Ojalvo
//         Created:  Mon, 18 Feb 2019 16:25:13 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/ForwardDetId/interface/MTDDetId.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/FTLDigi/interface/FTLDigiCollections.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

#include "DataFormats/FTLDigi/interface/FTLDigiCollections.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHit.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"

#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"

using namespace edm;
using namespace std;

struct simple_tpg{
  float eta;
  float phi;
};

struct mtd_tpg{
  float eta;
  float phi; 
  float time;
  float x;
  float y;
  float z;
  float energy;
};

struct MTDinfo {

  float sim_energy;
  float sim_time;
  float sim_x;
  float sim_y;
  float sim_z;

  uint32_t digi_row[2];
  uint32_t digi_col[2];
  uint32_t digi_charge[2];
  uint32_t digi_time1[2];
  uint32_t digi_time2[2];

  float ureco_charge[2];
  float ureco_time[2];

  float reco_energy;
  float reco_time;

};

//
// class declaration
//

class mtdTPGenerator : public edm::stream::EDProducer<> {
   public:
      explicit mtdTPGenerator(const edm::ParameterSet&);
      ~mtdTPGenerator();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      void initializeMTDMenuBased(vector<l1t::PFCluster> emCands, vector<l1t::PFCluster> hcCands, vector<l1t::Muon> muCands, std::vector<mtd_tpg> &mtd_tpgs_out);
      void initializeMTDROIBased(vector<l1t::PFCluster> emCands, vector<l1t::PFCluster> hcCands, vector<l1t::Muon> muCands, std::vector<mtd_tpg> &mtd_tpgs_out);
      void initializeMTDHits(edm::Event& iEvent, const edm::EventSetup& iSetup);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

  const MTDGeometry* geom_; 

  EDGetTokenT< vector<l1t::PFCandidate> > L1PFToken_;

  EDGetTokenT< BTLDigiCollection >    btlDigisToken_;
  EDGetTokenT< ETLDigiCollection >    etlDigisToken_;
  
  EDGetTokenT< FTLRecHitCollection >  btlRecHitToken_;
  EDGetTokenT< FTLRecHitCollection >  etlRecHitToken_;

  EDGetTokenT< FTLClusterCollection > btlClusterToken_;
  EDGetTokenT< FTLClusterCollection > etlClusterToken_;

  std::vector<mtd_tpg> mtd_tpgs;
  std::unordered_map<uint32_t, MTDinfo> btl_hits;
  std::unordered_map<uint32_t, MTDinfo> etl_hits[2];

  edm::EDGetTokenT<l1t::MuonBxCollection> muCands_;
  
  std::vector<edm::EDGetTokenT<l1t::PFClusterCollection>> emCands_;
  std::vector<edm::EDGetTokenT<l1t::PFClusterCollection>> hadCands_;

  float emPtCut_, hadPtCut_;
  bool useMenu_;
  float menuMuPtCut_,menuEmPtCut_,menuHcPtCut_;
  bool useROI_;
  float roiMuPtCut_, roiEmPtCut_, roiHcPtCut_;
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
mtdTPGenerator::mtdTPGenerator(const edm::ParameterSet& cfg):
  L1PFToken_(          consumes< vector<l1t::PFCandidate> >(cfg.getParameter<edm::InputTag>("L1PFObjects"))),
  btlClusterToken_(    consumes< FTLClusterCollection > ( cfg.getParameter<InputTag>("mtdClusterBarrel"))),
  etlClusterToken_(    consumes< FTLClusterCollection > ( cfg.getParameter<InputTag>("mtdClusterEndcap"))),
  btlDigisToken_(      consumes< BTLDigiCollection >   ( cfg.getParameter<InputTag>("FTLBarrel"))),
  etlDigisToken_(      consumes< ETLDigiCollection >   ( cfg.getParameter<InputTag>("FTLEndcap"))),
  btlRecHitToken_(     consumes< FTLRecHitCollection > ( cfg.getParameter<InputTag>("recHitBarrel"))), // finish me
  etlRecHitToken_(     consumes< FTLRecHitCollection > ( cfg.getParameter<InputTag>("recHitEndcap"))), // finish me
  muCands_(            consumes< l1t::MuonBxCollection> ( cfg.getParameter<edm::InputTag>("muons"))),
  emPtCut_(            cfg.getParameter<double>("emPtCut")),
  hadPtCut_(           cfg.getParameter<double>("hadPtCut")),
  //menu minimum mu pt, ecal calo pt, minimum hcal calo pt
  useMenu_(            cfg.getParameter<bool>("useMenu")),
  menuMuPtCut_(        cfg.getParameter<double>("menuMuPtCut")),
  menuEmPtCut_(        cfg.getParameter<double>("menuEmPtCut")),
  menuHcPtCut_(        cfg.getParameter<double>("menuHcPtCut")),
  //ROI 
  useROI_(             cfg.getParameter<bool>("useROI")),
  roiMuPtCut_(         cfg.getParameter<double>("roiMuPtCut")),
  roiEmPtCut_(         cfg.getParameter<double>("roiEmPtCut")),
  roiHcPtCut_(         cfg.getParameter<double>("roiHcPtCut"))
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/   //now do what ever other initialization is needed

  produces<l1t::PFCandidateCollection>("Time");

  for (auto & tag : cfg.getParameter<std::vector<edm::InputTag>>("emClusters")) {
    emCands_.push_back(consumes<l1t::PFClusterCollection>(tag));
  }
  for (auto & tag : cfg.getParameter<std::vector<edm::InputTag>>("hadClusters")) {
    hadCands_.push_back(consumes<l1t::PFClusterCollection>(tag));
  }
  
 
}


mtdTPGenerator::~mtdTPGenerator()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
mtdTPGenerator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   //Make the product

   std::vector<mtd_tpg> mtd_tpgs_selected;

   vector<l1t::PFCluster> emCands;
   vector<l1t::PFCluster> hadCands;
   vector<l1t::Muon>      muCands;

       // ------ READ CALOS -----
    edm::Handle<l1t::PFClusterCollection> caloHandle;
    for (const auto & tag : emCands_) {
        iEvent.getByToken(tag, caloHandle);
        const auto & calos = *caloHandle;
        for (unsigned int ic = 0, nc = calos.size(); ic < nc; ++ic) {
            const auto & calo = calos[ic];
            //if (debugR_ > 0 && deltaR(calo.eta(),calo.phi(),debugEta_,debugPhi_) > debugR_) continue;
            if (calo.pt() > emPtCut_) emCands.push_back(calo);
        }
    }

    for (const auto & tag : hadCands_) {
        iEvent.getByToken(tag, caloHandle);
        const auto & calos = *caloHandle;
        for (unsigned int ic = 0, nc = calos.size(); ic < nc; ++ic) {
            const auto & calo = calos[ic];
            //if (debugR_ > 0 && deltaR(calo.eta(),calo.phi(),debugEta_,debugPhi_) > debugR_) continue;
            if (calo.pt() > hadPtCut_) hadCands.push_back(calo);

        }
    }
  
    /// ------ READ MUONS ----
    edm::Handle<l1t::MuonBxCollection> muons;
    iEvent.getByToken(muCands_, muons);
    for (auto it = muons->begin(0), ed = muons->end(0); it != ed; ++it) {
        const l1t::Muon & mu = *it;
        //if (debugR_ > 0 && deltaR(mu.eta(),mu.phi(),debugEta_,debugPhi_) > debugR_) continue;
        muCands.push_back(mu);
    }
    
    /// ------ READ L1 PFCANDS ------

    edm::Handle< l1t::PFCandidateCollection > l1PFCandidates;
    iEvent.getByToken( L1PFToken_, l1PFCandidates);

   //Get the MTD hits, put in terms of pt, eta, phi
   initializeMTDHits(iEvent, iSetup);

   std::vector<mtd_tpg> mtdCands;

   //vector MTD TPGs
   if(useMenu_==true){
     initializeMTDMenuBased(emCands, hadCands, muCands, mtdCands);
   }
   else if(useROI_==true){
     initializeMTDROIBased( emCands, hadCands, muCands, mtdCands);
   }

   //std::cout<<"Number of MTDCands "<<mtdCands.size()<<std::endl;

   //std::cout<<"creating the new L1PFCands"<<std::endl;
   std::unique_ptr<l1t::PFCandidateCollection> newL1PFCands(new l1t::PFCandidateCollection);

   //assign time to pfcandidates
   for(auto l1PFCand : *l1PFCandidates){

     l1t::PFCandidate newPFCandidate = l1PFCand;
     //std::cout<<"new PFCand Pt: "<<l1PFCand.pt()<<std::endl;
     newPFCandidate.setTime(0);
     for(auto mtdCand : mtdCands){
       if(fabs(mtdCand.eta-l1PFCand.eta())+fabs(mtdCand.phi-l1PFCand.phi()) < 0.05 ){

	 //if(mtdCand.time < 0 || mtdCand.time > 19)
	   //std::cout<<"match to MTD found but time is strange, pt, eta, phi, energy, time: "<< mtdCand.eta<<", "<< mtdCand.phi<< ", "<< mtdCand.energy<< ", "<<mtdCand.time<<std::endl;

	 newPFCandidate.setTime(mtdCand.time);

	 break;
       }
     }
     newL1PFCands->push_back(newPFCandidate);
   }

   // transfer the product
   iEvent.put(std::move(newL1PFCands),"Time");
   
}

void 
mtdTPGenerator::initializeMTDMenuBased(vector<l1t::PFCluster> emCands, vector<l1t::PFCluster> hcCands, vector<l1t::Muon> muCands, std::vector<mtd_tpg> &mtd_tpgs_out){
  bool passMenu = false;
  //std::cout<<"initialize the menu "<<std::endl;
  //check to see if there is one calo object above pt threshold
  for(auto emCand:emCands){
    if(emCand.pt() > menuEmPtCut_)
      passMenu = true;
  }
  //std::cout<<" Finished EM Cands"<<std::endl;
  for(auto hcCand:hcCands){
    if(hcCand.pt() > menuHcPtCut_)
      passMenu = true;
  }
  //std::cout<<" Finished HC Cands"<<std::endl;
  //check to see if there is one muon above pt thershold
  for(auto muCand:muCands){
    if(muCand.pt() > menuMuPtCut_)
      passMenu = true;
  }
  //std::cout<<" Finished MU Cands"<<std::endl;
  //transfer all MTD digis from btl to L1 vector
  if(passMenu){
    std::cout<<"Passed the menu requirements"<<std::endl;
    for(auto mtd : mtd_tpgs)
      mtd_tpgs_out.push_back(mtd);
  }


}

void 
mtdTPGenerator::initializeMTDROIBased(vector<l1t::PFCluster> emCands, vector<l1t::PFCluster> hcCands, vector<l1t::Muon> muCands, std::vector<mtd_tpg> &mtd_tpgs_out){
  vector<simple_tpg> tpgs_above_threshold;

  //Find Calo Objects above threshold, match in eta/phi to MTD hits, fill vector
  for(auto emCand:emCands){
    if(emCand.pt() > roiEmPtCut_){
      simple_tpg temp;
      tpgs_above_threshold.push_back(temp);
    }
  }
  for(auto hcCand:hcCands){
    if(hcCand.pt() > roiHcPtCut_){
      simple_tpg temp;
      tpgs_above_threshold.push_back(temp);
    }
  }

  //Find Muon Objects above threshold, match in eta/phi to MTD hits, fill vector
  for(auto muCand:muCands){
    if(muCand.pt() > roiMuPtCut_){
      simple_tpg temp;
      tpgs_above_threshold.push_back(temp);
    }
  }

  //transfer all MTD digis from btl to L1 vector
  std::vector<mtd_tpg>::iterator iter_mtd;
  for(iter_mtd = mtd_tpgs.begin(); iter_mtd != mtd_tpgs.end(); ){
    for(auto tpg : tpgs_above_threshold){
      if(fabs(iter_mtd->eta - tpg.eta)<0.4 && fabs(iter_mtd->phi-tpg.phi)<0.4){
	mtd_tpgs_out.push_back(*iter_mtd);
	//remove the mtd from the list to avoid duplicating this process if there is ecal and hcal and muon
	iter_mtd = mtd_tpgs.erase(iter_mtd);
      }
      else
	iter_mtd++;
    }
  }
}


void mtdTPGenerator::initializeMTDHits(edm::Event& iEvent, const edm::EventSetup& iSetup){
  //std::cout<<"initialize the MTD Hits "<<std::endl;
  edm::ESHandle<MTDGeometry> geom;
  iSetup.get<MTDDigiGeometryRecord>().get(geom);
  geom_ = geom.product();

  mtd_tpgs.clear();

  edm::Handle<BTLDigiCollection> h_BTL_digi;
  iEvent.getByToken(btlDigisToken_,h_BTL_digi);

  edm::Handle<ETLDigiCollection> h_ETL_digi;
  iEvent.getByToken(etlDigisToken_,h_ETL_digi);

  edm::Handle< FTLRecHitCollection > btlRecHits;
  iEvent.getByToken(btlRecHitToken_,btlRecHits);

  edm::Handle< FTLRecHitCollection > etlRecHits;
  iEvent.getByToken(etlRecHitToken_,etlRecHits);

  unsigned int n_digi_btl[2] = {0,0};
  btl_hits.clear();
  
  if (h_BTL_digi->size() > 0 ) {

    //std::cout << " ----------------------------------------" << std::endl;
    //std::cout << " BTL DIGI collection:" << std::endl;
    
    for (const auto& dataFrame: *h_BTL_digi) {
      // in case print outs are needed for debugging this can be uncommented
      /*
      // --- detector element ID:
      std::cout << "   det ID:  det = " << dataFrame.id().det() 
		<< "  subdet = "  << dataFrame.id().mtdSubDetector() 
		<< "  side = "    << dataFrame.id().mtdSide() 
		<< "  rod = "     << dataFrame.id().mtdRR() 
		<< "  mod = "     << dataFrame.id().module() 
		<< "  type = "    << dataFrame.id().modType() 
		<< "  crystal = " << dataFrame.id().crystal() 
		<< std::endl;


      // --- loop over the dataFrame samples
      for (int isample = 0; isample<dataFrame.size(); ++isample){

	const auto& sample = dataFrame.sample(isample);

	std::cout << "       sample " << isample << ":"; 
	if ( sample.data()==0 && sample.toa()==0 ) {
	  std::cout << std::endl;
	  continue;
	}
	std::cout << "  amplitude = " << sample.data() 
		  << "  time1 = " <<  sample.toa() 
		  << "  time2 = " <<  sample.toa2() << std::endl;
      */


      DetId id =  dataFrame.id();
      
      const auto& sample_L = dataFrame.sample(0);
      const auto& sample_R = dataFrame.sample(1);

      btl_hits[id.rawId()].digi_row[0] = sample_L.row();
      btl_hits[id.rawId()].digi_row[1] = sample_R.row();
      btl_hits[id.rawId()].digi_col[0] = sample_L.column();
      btl_hits[id.rawId()].digi_col[1] = sample_R.column();

      btl_hits[id.rawId()].digi_charge[0] = sample_L.data();
      btl_hits[id.rawId()].digi_charge[1] = sample_R.data();
      btl_hits[id.rawId()].digi_time1[0]  = sample_L.toa();
      btl_hits[id.rawId()].digi_time1[1]  = sample_R.toa();
      btl_hits[id.rawId()].digi_time2[0]  = sample_L.toa2();
      btl_hits[id.rawId()].digi_time2[1]  = sample_R.toa2();

      if ( sample_L.data() > 0 )
	n_digi_btl[0]++;

      if ( sample_R.data() > 0 )
	n_digi_btl[1]++;
      
      
    } // digi loop

  } // if (h_BTL_digi->size() > 0 )

  // --- ETL

  unsigned int n_digi_etl[2] = {0,0};

  if ( h_ETL_digi->size() > 0 ) {

    for (const auto& dataFrame: *h_ETL_digi) {

      ETLDetId id =  dataFrame.id();
      int idet = (id.zside()+1)/2;

      // --- loop over the dataFrame samples
      for (int isample = 0; isample<dataFrame.size(); ++isample){

	const auto& sample = dataFrame.sample(isample);

	if ( sample.data()!=0 && sample.toa()!=0 ) {

	  // on-time sample
	  if ( isample == 2 ) {

	    etl_hits[idet][id.rawId()].digi_row[0] = sample.row();
	    etl_hits[idet][id.rawId()].digi_col[0] = sample.column();

	    etl_hits[idet][id.rawId()].digi_charge[0] = sample.data();
	    etl_hits[idet][id.rawId()].digi_time1[0]  = sample.toa();

	    n_digi_etl[idet]++;

	  }

	}

      } // isaple loop

    } // dataFrame loop

  } // if ( h_ETL_digi->size() > 0 )
  else
    std::cout<<"WARNING ETL DIGI SIZE IS 0!!"<<std::endl;

  unsigned int n_reco_btl = 0;

  if (btlRecHits->size() > 0 ) {

    for (const auto& recHit: *btlRecHits) {

      DetId id = recHit.id();

      btl_hits[id.rawId()].reco_energy = recHit.energy();
      btl_hits[id.rawId()].reco_time   = recHit.time();

      if ( recHit.energy() > 0. )
	n_reco_btl++;


    } // recHit loop

  } // if ( h_BTL_reco->size() > 0 )

    // --- ETL

  unsigned int n_reco_etl[2] = {0,0};

  if ( etlRecHits->size() > 0 ) {

    for (const auto& recHit: *etlRecHits) {

      ETLDetId id = recHit.id();
      int idet = (id.zside()+1)/2;

      etl_hits[idet][id.rawId()].reco_energy = recHit.energy();
      etl_hits[idet][id.rawId()].reco_time   = recHit.time();

      if ( recHit.energy() > 0. )
	n_reco_etl[idet]++;


    } // recHit loop

  } // if ( h_ETL_reco->size() > 0 )




  //
   for (auto const& hit: btl_hits) {

      BTLDetId detId(hit.first); 
      DetId geoId = BTLDetId(detId.mtdSide(),detId.mtdRR(),detId.module()+14*(detId.modType()-1),0,1);
      // this seg faults from time to time for specific events! 
      // it appears to originate in the external MTD geometry function...
      // this should be solved with Lindsey and other developers.
      // this portion is needed to convert the detID to eta/phi

      const MTDGeomDet* thedet           = geom_->idToDet(geoId);
      const ProxyMTDTopology& topoproxy  = static_cast<const ProxyMTDTopology&>(thedet->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());

      if ( (hit.second).reco_energy < 0.1 ) continue;

      Local3DPoint simscaled( 0.1*(hit.second).sim_x, 
			      0.1*(hit.second).sim_y, 
			      0.1*(hit.second).sim_z );

      simscaled              = topo.pixelToModuleLocalPoint( simscaled, detId.row( topo.nrows()), detId.column(topo.nrows()) );

      // Finally get an object in terms of eta, phi
      const auto& global_pos = thedet->toGlobal(simscaled);

      mtd_tpg temptpg;
      temptpg.eta    = (double)global_pos.eta();
      temptpg.phi    = (double)global_pos.phi();
      temptpg.time   = (double)(hit.second).reco_time;
      temptpg.x      = (double)global_pos.x();
      temptpg.y      = (double)global_pos.y();
      temptpg.z      = (double)global_pos.z();
      temptpg.energy = (double)(hit.second).reco_energy;

      // get rid of the nonsense
      if(temptpg.time < 19 && temptpg.time > 0)
	mtd_tpgs.push_back(temptpg);

      //if(temptpg.time < 0 )
      //std::cout<<"error!! time is less than 0, eta, phi, energy, time: "<< temptpg.eta<<","<< temptpg.phi <<","<< temptpg.energy<<","<<temptpg.time<<std::endl;
   }
   
   for (int idet=0; idet<2; ++idet){
     //std::cout<<"found etl hits size: "<<etl_hits[idet].size()<<std::endl;
     for (auto const& hit: etl_hits[idet]) {

       ETLDetId detId(hit.first);
       DetId geoId = ETLDetId(detId.mtdSide(),detId.mtdRR(),detId.module(),0);
       const MTDGeomDet* thedet = geom_->idToDet(geoId);
       const PixelTopology& topo = static_cast<const PixelTopology&>(thedet->topology());

       // --- SIM

       if ( (hit.second).reco_time != 0. ) {
	 // Get the SIM hit global position

	 Local3DPoint simscaled(0.1*(hit.second).sim_x,0.1*(hit.second).sim_y,0.1*(hit.second).sim_z);
	 const auto& global_pos = thedet->toGlobal(simscaled);

	 mtd_tpg temptpg;
	 temptpg.eta    = (double)global_pos.eta();
	 temptpg.phi    = (double)global_pos.phi();
	 temptpg.time   = (double)(hit.second).reco_time;
	 temptpg.x      = (double)global_pos.x();
	 temptpg.y      = (double)global_pos.y();
	 temptpg.z      = (double)global_pos.z();
	 temptpg.energy = (double)(hit.second).reco_energy;
	 //std::cout<<"mtd time "<<temptpg.time<<std::endl;
	 // get rid of the nonsense
	 if(temptpg.time < 19 && temptpg.time > 0)
	   mtd_tpgs.push_back(temptpg);

	 //if(temptpg.time < 0 )
	   //std::cout<<"error!! time is less than 0, eta, phi, energy, time: "<< temptpg.eta<<","<< temptpg.phi <<","<< temptpg.energy<<","<<temptpg.time<<std::endl;
       }
     }
   }
   std::cout<<"Finish the MTD Hits "<<std::endl;

}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
mtdTPGenerator::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
mtdTPGenerator::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
mtdTPGenerator::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
mtdTPGenerator::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
mtdTPGenerator::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
mtdTPGenerator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
mtdTPGenerator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(mtdTPGenerator);
