// -*- C++ -*-
//
// Package:    L1TkPartAnalyzer
// Class:      TkTriggerParticleAnalzer
// 
/**\class TkTriggerParticleAnalzer TkTriggerParticleAnalzer.cc SLHCUpgradeSimulations/TkTriggerParticleAnalzer/src/TkTriggerParticleAnalzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Emmanuelle Perez,40 1-A28,+41227671915,
//         Created:  Thu Nov 14 11:22:13 CET 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"


// Gen-level stuff:
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkPrimaryVertex.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEtMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEtMissParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkJetParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkJetParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkHTMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkHTMissParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace l1extra;


//
// class declaration
//

class TkTriggerParticleAnalzer : public edm::EDAnalyzer {
public:

  typedef TTTrack< Ref_PixelDigi_ >  L1TkTrackType;
  typedef std::vector< L1TkTrackType >  L1TkTrackCollectionType;

  explicit TkTriggerParticleAnalzer(const edm::ParameterSet&);
  ~TkTriggerParticleAnalzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  // ----------member data ---------------------------
  // float tp_pt, tp_eta, tp_phi, tp_z;
  std::vector<float>tp_pt; std::vector<float>tp_eta;  std::vector<float>tp_phi; std::vector<float>tp_z;
  std::vector<float>tk_pt; std::vector<float>tk_eta;  std::vector<float>tk_phi; std::vector<float>tk_z;
  std::vector<float>tk_chi2; std::vector<int>tk_truth;
  TTree*AmMuons;
  edm::InputTag GenPartInputTag;
  edm::InputTag TrackPartInputTag;
  edm::InputTag TTTracksInputTag;
  edm::InputTag TTTracksAssocInputTag;
  int PDG;
};



void TkTriggerParticleAnalzer::beginJob()
{
  edm::Service<TFileService> fs;
  AmMuons = fs->make<TTree>("AmMuons", "");
  // AmMuons->Branch("tp_eta", &tp_eta, "tp_eta/F");
  // AmMuons->Branch("tp_phi", &tp_phi, "tp_phi/F");
  // AmMuons->Branch("tp_pt", &tp_pt, "tp_pt/F");
  // AmMuons->Branch("tp_z", &tp_z, "tp_z/F");
  AmMuons->Branch("tp_eta", &tp_eta);
  AmMuons->Branch("tp_phi", &tp_phi);
  AmMuons->Branch("tp_pt", &tp_pt);
  AmMuons->Branch("tp_z", &tp_z);

  AmMuons->Branch("tk_eta", &tk_eta);
  AmMuons->Branch("tk_phi", &tk_phi);
  AmMuons->Branch("tk_pt", &tk_pt);
  AmMuons->Branch("tk_z", &tk_z);
  AmMuons->Branch("tk_chi2", &tk_chi2);
  AmMuons->Branch("tk_truth", &tk_truth);
}




TkTriggerParticleAnalzer::TkTriggerParticleAnalzer(const edm::ParameterSet& iConfig)

{
  GenPartInputTag   = iConfig.getParameter<edm::InputTag>("GenPartInputTag");
  TrackPartInputTag =iConfig.getParameter<edm::InputTag>("TrackPartTag");
  TTTracksInputTag  =iConfig.getParameter<edm::InputTag>("TTTracksInputTag");
  TTTracksAssocInputTag=iConfig.getParameter<edm::InputTag>("inputTagMC");
  PDG=iConfig.getParameter< int >("ParticleType");
}


TkTriggerParticleAnalzer::~TkTriggerParticleAnalzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
TkTriggerParticleAnalzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle;
  iEvent.getByLabel(TrackPartInputTag, TrackingParticleHandle);

  edm::Handle< std::vector< TTTrack< Ref_PixelDigi_ > > > TTTrackHandle;
  iEvent.getByLabel(TTTracksInputTag,TTTrackHandle);

  edm::Handle< TTTrackAssociationMap< Ref_PixelDigi_ > > MCTruthTTTrackHandle;
  iEvent.getByLabel(TTTracksAssocInputTag,MCTruthTTTrackHandle);

  // tp_eta=999;
  // tp_phi=999;
  // tp_pt=-1;
  // tp_z=999;
  tp_eta.resize(0);
  tp_phi.resize(0);
  tp_pt.resize(0);
  tp_z.resize(0);

 
  tk_eta.resize(0);
  tk_phi.resize(0);
  tk_pt.resize(0);
  tk_z.resize(0);
  tk_chi2.resize(0);
  tk_truth.resize(0);

  if(TrackingParticleHandle.isValid()) {
    std::vector< TrackingParticle >::const_iterator iterTP;
    for(iterTP=TrackingParticleHandle->begin(); iterTP!=TrackingParticleHandle->end(); ++iterTP){
      // if(abs(iterTP->pdgId())!=PDG) continue;
      //pt, eta, phi, z
      if (iterTP->pt() < 3) continue;
      // tp_eta=iterTP->eta();
      // tp_phi=iterTP->phi();
      // tp_pt=iterTP->pt();
      // tp_z=iterTP->vz();
      tp_eta.push_back(iterTP->eta());
      tp_phi.push_back(iterTP->phi());
      tp_pt.push_back(iterTP->pt());
      tp_z.push_back(iterTP->vz());
      // break;
    }
  }

  if(TTTrackHandle.isValid()) {
    std::vector< TTTrack< Ref_PixelDigi_ > >::const_iterator iterL1Track;

    int this_l1track = 0;
    for ( iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++ ) {
      edm::Ptr< TTTrack< Ref_PixelDigi_ > > l1track_ptr(TTTrackHandle, this_l1track);
      ++this_l1track;
      tk_eta.push_back(iterL1Track->getMomentum(5).eta());
      tk_phi.push_back(iterL1Track->getMomentum(5).phi());
      tk_pt.push_back(iterL1Track->getMomentum(5).perp());
      tk_z.push_back(iterL1Track->getPOCA(5).z());
      tk_chi2.push_back(iterL1Track->getChi2(5));
      int tmp_trk_genuine = 0;

      if (MCTruthTTTrackHandle->isGenuine(l1track_ptr)) tmp_trk_genuine = 1;
      else tmp_trk_genuine=0;
      tk_truth.push_back(tmp_trk_genuine);
    }
  }
  AmMuons->Fill();
}


// ------------ method called once each job just after ending the event loop  ------------
void 
TkTriggerParticleAnalzer::endJob() 
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TkTriggerParticleAnalzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TkTriggerParticleAnalzer);
