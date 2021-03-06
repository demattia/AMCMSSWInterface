#ifndef NTupleTools_AMTrackProducer_h_


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/PatternMatcher.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackFitter.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "SLHCUpgradeSimulations/L1TrackTrigger/interface/StubPtConsistency.h"
#include "SLHCL1TrackTriggerSimulations/NTupleTools/interface/MapTTStubs.h" 
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TriggerTowerMap.h"
//#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
//#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include<TFile.h>
#include<TTree.h>
#include<TH2F.h>
#include<map>
class AMTrackProducer : public edm::EDProducer {
  public:
    explicit AMTrackProducer(const edm::ParameterSet&);

  private:
    virtual void beginJob();
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void initOptions(ProgramOption &option, int PBIndex);
double genTrackDistanceLongitudinal(const double &z0, const double &cotTheta, const double &pt, const double &d0,
const int charge, const double &B, const double &R, const double &z);
double genTrackDistanceTransverse(const double &pt, const double &phi0, const double &d0,
                                                 const int charge, const double &B, const double &phi, const double &R);
int PatternBankMap(float eta, float phi);

 virtual void endJob();

ProgramOption option_;
PatternMatcher* matcher_;

std::vector<float>tp_pt, tp_eta, tp_phi;
std::vector<float>tp_vz, tp_chgOverPt;
std::vector<float>matchtp_pt, matchtp_eta, matchtp_phi;
std::vector<float>matchtp_vz, matchtp_chgOverPt;

std::vector<float>trk_pt, trk_eta, trk_phi;
std::vector<float>trk_vz, trk_chgOverPt;
std::vector<float>trk_chi2, trk_cotheta;
std::vector<float>trk_rinv;
std::vector<int>trk_class;
std::vector<int>trk_ghost;
int NRoads;int NStubs;
int NStubs0, NStubs1, NStubs2, NStubs3, NStubs4, NStubs5;
int NgoodStubs0, NgoodStubs1, NgoodStubs2, NgoodStubs3, NgoodStubs4, NgoodStubs5;
std::vector<int> StubsinRoad;
std::vector<int> goodStubsinRoad;
TTree*Amtracks;
edm::InputTag StubsTag_;
TH2F*TrigTowerMap;
TriggerTowerMap * ttmap_;
TFile*fin;
};
void AMTrackProducer::beginJob(){
   ttmap_ = new TriggerTowerMap();
   ttmap_->read("/fdata/hepx/store/user/rish/AMSIMULATION/CMSSW_6_2_0_SLHC25_patch3/src/SLHCL1TrackTriggerSimulations/AMSimulation/data/"); 
   fin=new TFile("TrigTowerMap.root", "READ");
   TrigTowerMap=(TH2F*)fin->Get("h2TrigMap");
   edm::Service<TFileService> fs;


    Amtracks = fs->make<TTree>("Amtracks", "");
    Amtracks->Branch("tp_pt", &tp_pt);
    Amtracks->Branch("tp_eta", &tp_eta);
    Amtracks->Branch("tp_phi", &tp_phi);
    Amtracks->Branch("tp_vz", &tp_vz);
    Amtracks->Branch("tp_chgOverPt", &tp_chgOverPt);
    Amtracks->Branch("matchtp_pt", &matchtp_pt);
    Amtracks->Branch("matchtp_eta", &matchtp_eta);
    Amtracks->Branch("matchtp_phi", &matchtp_phi);
    Amtracks->Branch("matchtp_vz", &matchtp_vz);
    Amtracks->Branch("matchtp_chgOverPt", &matchtp_chgOverPt);
    Amtracks->Branch("trk_pt", &trk_pt);
    Amtracks->Branch("trk_eta", &trk_eta);
    Amtracks->Branch("trk_phi", &trk_phi);
    Amtracks->Branch("trk_vz", &trk_vz);
    Amtracks->Branch("trk_chi2", &trk_chi2);
    Amtracks->Branch("trk_rinv", &trk_rinv);
    Amtracks->Branch("trk_cotheta", &trk_cotheta);
    Amtracks->Branch("trk_class", &trk_class);
    Amtracks->Branch("trk_ghost", &trk_ghost);
    Amtracks->Branch("NRoads", &NRoads, "NRoads/I");
    Amtracks->Branch("NStubs", &NStubs, "NStubs/I");
    Amtracks->Branch("NStubs0", &NStubs0, "NStubs0/I");
    Amtracks->Branch("NStubs1", &NStubs1, "NStubs1/I");
    Amtracks->Branch("NStubs2", &NStubs2, "NStubs2/I");
    Amtracks->Branch("NStubs3", &NStubs3, "NStubs3/I");
    Amtracks->Branch("NStubs4", &NStubs4, "NStubs4/I");
    Amtracks->Branch("NStubs5", &NStubs5, "NStubs5/I");
    Amtracks->Branch("NgoodStubs0", &NgoodStubs0, "NgoodStubs0/I");
    Amtracks->Branch("NgoodStubs1", &NgoodStubs1, "NgoodStubs1/I");
    Amtracks->Branch("NgoodStubs2", &NgoodStubs2, "NgoodStubs2/I");
    Amtracks->Branch("NgoodStubs3", &NgoodStubs3, "NgoodStubs3/I");
    Amtracks->Branch("NgoodStubs4", &NgoodStubs4, "NgoodStubs4/I");
    Amtracks->Branch("NgoodStubs5", &NgoodStubs5, "NgoodStubs5/I");
    Amtracks->Branch("StubsinRoad", &StubsinRoad);
    Amtracks->Branch("goodStubsinRoad", &goodStubsinRoad);


}
void AMTrackProducer::endJob(){
delete ttmap_;
delete fin;
}
AMTrackProducer::AMTrackProducer(const edm::ParameterSet& iConfig) {
     StubsTag_ =(iConfig.getParameter<edm::InputTag>("inputTagStub"));
     produces< std::vector< TTTrack< Ref_PixelDigi_ > > >( "Level1TTTracks" ).setBranchAlias("Level1TTTracks");

}

void AMTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    //std::auto_ptr<std::vector<float> > dummies(new std::vector<float>());
//  typedef TTStubAssociationMap<Ref_PixelDigi_> assocmap_stub;

  std::auto_ptr< std::vector< TTTrack< Ref_PixelDigi_ > > > L1TkTracksForOutput( new std::vector< TTTrack< Ref_PixelDigi_ > > );
	 edm::Handle<edmNew::DetSetVector<TTStub<Ref_PixelDigi_> > > pixelDigiTTStubs;
        iEvent.getByLabel(StubsTag_, pixelDigiTTStubs);

    std::auto_ptr<std::vector< TTTrack<Ref_PixelDigi_> > > TTTrackVector(new std::vector<TTTrack<Ref_PixelDigi_> >());
    ProgramOption option;
/*
  edm::ESHandle<StackedTrackerGeometry> stackedGeometryHandle;
    iSetup.get<StackedTrackerGeometryRecord>().get(stackedGeometryHandle);
 const StackedTrackerGeometry *theStackedGeometry = stackedGeometryHandle.product();
edm::ESHandle<MagneticField> magneticFieldHandle;
iSetup.get<IdealMagneticFieldRecord>().get(magneticFieldHandle);
const MagneticField* theMagneticField = magneticFieldHandle.product();
double mMagneticFieldStrength = theMagneticField->inTesla(GlobalPoint(0,0,0)).z();
*/
    TTStubPlusTPReader reader(0);
    //std::cout<<"Just declare constructor "<<std::endl;
    if(!pixelDigiTTStubs.isValid())return;
    reader.init(iEvent);//just to get tracking particles
    //std::cout<<"After Init of TTStubPlusTPReader "<<std::endl; 
    NStubs= reader.vb_r->size();
    NStubs0=0;
    NStubs1=0;
    NStubs2=0;
    NStubs3=0;
    NStubs4=0;
    NStubs5=0;
    NgoodStubs0=0;
    NgoodStubs1=0;
    NgoodStubs2=0;
    NgoodStubs3=0;
    NgoodStubs4=0;
    NgoodStubs5=0;
    for (unsigned int istub=0; istub<reader.vb_r->size(); ++istub) {
	float r=reader.vb_r->at(istub);
          // unsigned moduleId = reader.vb_modId   ->at(istub);
//	    int tpID=reader.vb_tpId->at(istub);
	
/*
	    for( int t=0; t<48; ++t){
	    std::map<unsigned, bool> ttrmap = ttmap_->getTriggerTowerReverseMap(t);
            if (ttrmap.find(moduleId) != ttrmap.end()){
			if(tpID>-1)std::cout<<"Module ID "<<tpID<<" Tower "<<t<<std::endl;	
			break;
		}
*/
	
	if(r>20 && r<30)++NStubs0;
        if(r>20 && r<30 && reader.vb_tpId->at(istub)>=0)++NgoodStubs0;
        if(r>30 && r<40)++NStubs1;
        if(r>30 && r<40 && reader.vb_tpId->at(istub)>=0)++NgoodStubs1;
        if(r>40 && r<60)++NStubs2;
        if(r>40 && r<60 && reader.vb_tpId->at(istub)>=0)++NgoodStubs2;
        if(r>60 && r<70)++NStubs3;
        if(r>60 && r<70 && reader.vb_tpId->at(istub)>=0)++NgoodStubs3;
        if(r>70 && r<95)++NStubs4;
        if(r>70 && r<95 && reader.vb_tpId->at(istub)>=0)++NgoodStubs4;
        if(r>100)++NStubs5;
	if(r>100  && reader.vb_tpId->at(istub)>=0)++NgoodStubs5;
 
    }


    option.datadir = "/fdata/hepx/store/user/rish/AMSIMULATION/Forked/CMSSW_6_2_0_SLHC25_patch3/"; //std::getenv("CMSSW_BASE");

    option.datadir += "/src/SLHCL1TrackTriggerSimulations/AMSimulation/data/";

    tp_pt.resize(0); tp_eta.resize(0); tp_phi.resize(0);
    tp_vz.resize(0); tp_chgOverPt.resize(0);
    int PBIndex=-1;

for(unsigned int tp=0; tp<reader.vp2_pt->size(); ++tp){
	if(abs(reader.vp2_pdgId->at(tp))!=13)continue;
	PBIndex=PatternBankMap(reader.vp2_eta->at(tp),(reader.vp2_phi->at(tp)));
	//std::cout<<"eta, phi "<<reader.vp2_eta->at(tp)<<", "<<reader.vp2_phi->at(tp)<<"MY pattern Bank "<<PBIndex<<std::endl;
	
	tp_pt.push_back(reader.vp2_pt->at(tp));	
        tp_eta.push_back(reader.vp2_eta->at(tp));
        tp_phi.push_back(reader.vp2_phi->at(tp));
        tp_vz.push_back(reader.vp2_vz->at(tp));
	tp_chgOverPt.push_back(reader.vp2_charge->at(tp)/reader.vp2_pt->at(tp));
}
    if( fabs(tp_eta[0])>1.42)return;
    float rmax=0;
  unsigned int inTower=0;
   TH1F*OtherTowers=new TH1F("OtherTowers", "",48, 1, 49);
   // std::vector<int>OtherTower;
    for (unsigned int istub=0; istub<reader.vb_r->size(); ++istub) {
	float r=reader.vb_r->at(istub);
            unsigned moduleId = reader.vb_modId   ->at(istub);
	    int tpID=reader.vb_tpId->at(istub);
	    if(tpID<0)continue;
	  
	  //  for( int t=0; t<48; ++t){
	    std::map<unsigned, bool> ttrmap = ttmap_->getTriggerTowerReverseMap(PBIndex);
            if (ttrmap.find(moduleId) != ttrmap.end()){
			//std::cout<<"Module ID "<<r<<" Tower "<<PBIndex<<std::endl;	
			if(r>rmax){rmax=r; ++inTower;}
			//float dz=genTrackDistanceLongitudinal( tp_vz[0], 
		}
	    else{
		for( int t=0; t<48; ++t){
			if(t==PBIndex)continue;
			ttrmap = ttmap_->getTriggerTowerReverseMap(t);	
			if(ttrmap.find(moduleId) != ttrmap.end()){
				
				if(r>rmax){rmax=r; OtherTowers->Fill(t);}		
				std::cout<<"Module ID "<<r<<"Other Tower "<<t<<std::endl;	
			}
		}
	}
}

    if(OtherTowers->Integral()>inTower){
	   int binmax = OtherTowers->GetMaximumBin();
   double x = OtherTowers->GetXaxis()->GetBinLowEdge(binmax);
	std::cout<<"Other Tower with Hits "<<x<<std::endl;
	PBIndex=x;
//	for(unsigned i=0; i<OtherTower.size(); ++i)
		//if(i>0 && OtherTower[i]===OtherTower[i-1]){
	//		++frequency;
	//`	}
//		std::cout<<" Other Towers "<<OtherTower[i]<<std::endl;
    }	
      
    std::cout<<"In Tower "<<inTower<<" Other Tower "<<OtherTowers->Integral()<<std::endl;
    std::cout<<"PBIndex "<<PBIndex<<std::endl;

    if(PBIndex<1)return;
    initOptions(option,PBIndex);
    
    matcher_=new PatternMatcher(option);
 
    matcher_->loadPatterns(option.bankfile);
  
   // std::vector<TTRoad> roads=matcher_->makeRoads(iEvent);
    delete matcher_;
    delete OtherTowers;
    //NRoads=roads.size();
   // std::cout<<"Roads "<<roads.size()<<std::endl; 
 /*   
 TrackFitter trackFit(option);
 //std::vector<TTTrack2>ttracks;
   std::vector<TTTrack2>ttracks=trackFit.makeTracks(iEvent, roads);
    StubsinRoad.resize(0);
    goodStubsinRoad.resize(0);
    for (unsigned iroad=0; iroad<roads.size(); ++iroad) {
        StubsinRoad.push_back(roads[iroad].stubRefs.size());
   	int good=0; 
        for(unsigned istub=0; istub<roads[iroad].stubRefs.size(); ++istub){
	   if(reader.vb_tpId->at(istub)>=0)++good;
	}
	goodStubsinRoad.push_back(good);
}
    
//now just take L1 ttracks2 and build the standard collection:
    for(unsigned int t=0; t<ttracks.size();++t){
	TTTrack<Ref_PixelDigi_> aTrack;
	 std::vector<unsigned>StubsRefIndex=ttracks[t].stubRefs();
	unsigned sindex = 0;
	edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >::const_iterator inputIter;
	edmNew::DetSet< TTStub< Ref_PixelDigi_ > >::const_iterator stubIter;
	for ( inputIter = pixelDigiTTStubs->begin(); inputIter != pixelDigiTTStubs->end(); ++inputIter ){
    		for ( stubIter = inputIter->begin(); stubIter != inputIter->end(); ++stubIter ){
			std::vector<unsigned>::iterator it=std::find(StubsRefIndex.begin(), StubsRefIndex.end(),sindex);
			if(it!=StubsRefIndex.end()){
			edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > tempStubRef = edmNew::makeRefTo( pixelDigiTTStubs, stubIter);
			  aTrack.addStubRef(tempStubRef);
		}
		++sindex;
   	   }
	}
	GlobalVector p3(GlobalVector::Cylindrical(ttracks[t].pt(), ttracks[t].phi0(),ttracks[t].pt()*sinh(ttracks[t].eta()) ));
	aTrack.setMomentum(p3,4);
	aTrack.setRInv(ttracks[t].rinv(),4);
        aTrack.setChi2(ttracks[t].chi2()/ttracks[t].ndof(),4);
	GlobalPoint POCA(0, 0,ttracks[t].vz());
	float consistency4par = StubPtConsistency::getConsistency(aTrack, theStackedGeometry, mMagneticFieldStrength, 4);
	aTrack.setPOCA(POCA,4);
	aTrack.setStubPtConsistency(consistency4par,4);
   	//if(ttracks[t].chi2Red()<5 && aTrack.getStubRefs().size()>=4 && !ttracks[t].isGhost())
	L1TkTracksForOutput->push_back( aTrack);
   }
   std::cout<<"Tracks "<<L1TkTracksForOutput->size()<<std::endl;
//iEvent.put( L1TkTracksForOutput, "Level1TTTracks");
   //track particle loop
matchtp_pt.resize(0); matchtp_eta.resize(0); matchtp_phi.resize(0);
matchtp_vz.resize(0); matchtp_chgOverPt.resize(0);
trk_pt.resize(0); trk_eta.resize(0); trk_phi.resize(0);
trk_vz.resize(0); trk_chgOverPt.resize(0);trk_chi2.resize(0);
trk_cotheta.resize(0);trk_rinv.resize(0);
trk_class.resize(0);trk_ghost.resize(0);

    for(unsigned int t=0; t<ttracks.size(); ++t){//track loop
	trk_pt.push_back(ttracks[t].pt());	
	trk_eta.push_back(ttracks[t].eta());
	trk_phi.push_back(ttracks[t].phi());
	trk_vz.push_back(ttracks[t].vz());
	trk_chi2.push_back(ttracks[t].chi2Red());
	trk_rinv.push_back(ttracks[t].rinv());
	trk_cotheta.push_back(ttracks[t].cottheta());
	trk_class.push_back(ttracks[t].synTpId());

	if(ttracks[t].isGhost())trk_ghost.push_back(1);
	else trk_ghost.push_back(0);

//trk_chgOverPt.push_back(ttracks[t].z0());

	if(ttracks[t].tpId()>=0){
		matchtp_pt.push_back(reader.vp2_pt->at(ttracks[t].tpId()));
		matchtp_eta.push_back(reader.vp2_eta->at(ttracks[t].tpId()));
                matchtp_phi.push_back(reader.vp2_phi->at(ttracks[t].tpId()));
                matchtp_vz.push_back(reader.vp2_vz->at(ttracks[t].tpId()));
		matchtp_chgOverPt.push_back(reader.vp2_charge->at(ttracks[t].tpId())/reader.vp2_pt->at(ttracks[t].tpId()));
	}
    }
*/
	Amtracks->Fill();
}

double AMTrackProducer::genTrackDistanceTransverse(const double &pt, const double &phi0, const double &d0,
                                                 const int charge, const double &B, const double &phi, const double &R) 
{
  //std::cout<<"Transverse input: "<<pt<<", "<<phi0<<", "<<d0<<", "<<charge<<", "<<B<<", "<<phi<<", "<<R<<std::endl;
    double rho = charge*pt/(B*0.003);
    double phiGen = phi0 - asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0)));
    double deltaPhi = (phi - phiGen);
    if (deltaPhi > M_PI) deltaPhi -= M_PI;
    else if (deltaPhi < -M_PI) deltaPhi += M_PI;
  //            //std::cout<<"Transverse output: "<<deltaPhi<<std::endl;
    return deltaPhi;
  }
double AMTrackProducer::genTrackDistanceLongitudinal(const double &z0, const double &cotTheta, const double &pt, const double &d0,
const int charge, const double &B, const double &R, const double &z) {
  //                                                                   //std::cout<<"Longitudinal input: "<<z0<<", "<<cotTheta<<", "<<pt<<", "<<d0<<", "<<charge<<", "<<B<<", "<<R<<", "<<z<<std::endl;
   double rho = charge*pt/(B*0.003);
 double zGen = z0 + charge*rho*cotTheta*acos((pow(rho,2) + pow(d0+rho,2) - pow(R,2))/(2*rho*(rho+d0)));
  //                                                                         //std::cout<<"Longitudinal output: "<<z-zGen<<std::endl;
  return (z - zGen);
  }

int  AMTrackProducer::PatternBankMap(float eta, float phi){
int pb=0;
int ib=TrigTowerMap->FindBin(eta,phi);
if(ib>0)pb=TrigTowerMap->GetBinContent(ib);
return pb;
}

void AMTrackProducer::initOptions(ProgramOption &option, int PBIndex){
//    option.bankfile=TString::Format("/fdata/hepx/store/user/rish/AMSIMULATION/CMSSW_6_2_0_SLHC25_patch3/src/SLHCL1TrackTriggerSimulations/AMSimulation/data/Patterns/patternBank_tt%d_sf1_nz1_pt3_100M.root", PBIndex);
    option.bankfile=TString::Format("/fdata/hepx/store/user/rish/AMSIMULATION/Forked/CMSSW_6_2_0_SLHC25_patch3/src/patternBank_tt%d_sf1_nz4_pt3_5M.root", PBIndex);
    option.verbose=2;
    option.tower=PBIndex;
    option.superstrip="sf1_nz4";
 //   option.maxPatterns=150327;
    option.nLayers=6;
    option.minFrequency=1;
     option.maxStubs=99999;
    option.maxMisses=0;
    option.maxRoads=99999;

    option.algo="LTF";
   option.maxChi2=99999;
    option.minNdof=1;
    option.maxCombs=99999;
    option.maxTracks=99999;
    option.hitBits=0;
    option.view="XYZ";
    //option.matrix="
}
DEFINE_FWK_MODULE(AMTrackProducer);

#endif
