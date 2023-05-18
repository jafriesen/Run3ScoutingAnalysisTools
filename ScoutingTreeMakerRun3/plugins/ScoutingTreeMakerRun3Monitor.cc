// -*- C++ -*-
//
// Package:    Run3ScoutingAnalysisTools/ScoutingTreeMakerRun3Monitor
// Class:      ScoutingTreeMakerRun3Monitor
//
/**\class ScoutingTreeMakerRun3Monitor ScoutingTreeMakerRun3Monitor.cc Run3ScoutingAnalysisTools/ScoutingTreeMakerRun3/plugins/ScoutingTreeMakerRun3Monitor.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  David Sperka
//         Created:  Sat, 11 Feb 2023 14:15:08 GMT
//
//

// system include files
#include <memory>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

#include "DataFormats/PatCandidates/interface/Muon.h"

//
// class declaration
//

class ScoutingTreeMakerRun3Monitor : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit ScoutingTreeMakerRun3Monitor(const edm::ParameterSet&);
  ~ScoutingTreeMakerRun3Monitor() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    const edm::InputTag triggerResultsTag;
    const edm::EDGetTokenT<edm::TriggerResults>                 triggerResultsToken;
    const edm::EDGetTokenT<std::vector<Run3ScoutingMuon> >      muonsToken;
    const edm::EDGetTokenT<std::vector<pat::Muon> >               offlineMuonsToken;
    const edm::EDGetTokenT<std::vector<Run3ScoutingElectron> >  electronsToken;
    //const edm::EDGetTokenT<std::vector<Run3ScoutingVertex> >    verticesToken;
    const edm::EDGetTokenT<double>                              rhoToken;
    const edm::EDGetTokenT<std::vector<Run3ScoutingPhoton> >    photonsToken;
    const edm::EDGetTokenT<std::vector<Run3ScoutingParticle> >  pfcandsToken;
    const edm::EDGetTokenT<std::vector<Run3ScoutingPFJet> >     pfjetsToken;
    const edm::EDGetTokenT<std::vector<Run3ScoutingTrack> >     tracksToken;

    std::vector<std::string> triggerPathsVector;
    std::map<std::string, int> triggerPathsMap;

    bool doL1;
    triggerExpression::Data triggerCache_;

    edm::InputTag                algInputTag_;
    edm::InputTag                extInputTag_;
    edm::EDGetToken              algToken_;
    std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
    std::vector<std::string>     l1Seeds_;
    std::vector<std::string>     l1MonitorSeeds_;
    std::vector<bool>            l1Result_;

    TTree* tree;
    float mass;
    float pt;
    float dr;
    float pt1;
    float pt2;
    float eta1;
    float eta2;
    int id1;
    int id2;
    float rho;
    int nScoutingMuons;

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
ScoutingTreeMakerRun3Monitor::ScoutingTreeMakerRun3Monitor(const edm::ParameterSet& iConfig):
    triggerResultsTag        (iConfig.getParameter<edm::InputTag>("triggerresults")),
    triggerResultsToken      (consumes<edm::TriggerResults>                    (triggerResultsTag)),
    muonsToken               (consumes<std::vector<Run3ScoutingMuon> >             (iConfig.getParameter<edm::InputTag>("muons"))),
    offlineMuonsToken        (consumes<std::vector<pat::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>("offlineMuons"))),
    electronsToken           (consumes<std::vector<Run3ScoutingElectron> >         (iConfig.getParameter<edm::InputTag>("electrons"))),
    //verticesToken            (consumes<std::vector<Run3ScoutingVertex> >           (iConfig.getParameter<edm::InputTag>("vertices"))),
    rhoToken                 (consumes<double>                                 (iConfig.getParameter<edm::InputTag>("rho"))), 
    photonsToken           (consumes<std::vector<Run3ScoutingPhoton> >         (iConfig.getParameter<edm::InputTag>("photons"))),
    pfcandsToken             (consumes<std::vector<Run3ScoutingParticle> >         (iConfig.getParameter<edm::InputTag>("pfcands"))),
    pfjetsToken              (consumes<std::vector<Run3ScoutingPFJet> >            (iConfig.getParameter<edm::InputTag>("pfjets"))),
    tracksToken              (consumes<std::vector<Run3ScoutingTrack> >            (iConfig.getParameter<edm::InputTag>("tracks"))),
    doL1                     (iConfig.existsAs<bool>("doL1")               ?    iConfig.getParameter<bool>  ("doL1")            : false)
{
    usesResource("TFileService");
    if (doL1) {
        algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
        extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
        algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
        l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
        l1MonitorSeeds_ = iConfig.getParameter<std::vector<std::string> >("l1MonitorSeeds");
        l1GtUtils_ = std::make_unique<l1t::L1TGlobalUtil>(iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);
    }
    else {
        l1Seeds_ = std::vector<std::string>();
        l1MonitorSeeds_ = std::vector<std::string>();
        l1GtUtils_ = 0;
    }
}

ScoutingTreeMakerRun3Monitor::~ScoutingTreeMakerRun3Monitor() {
}

//
// member functions
//

// ------------ method called for each event  ------------
void ScoutingTreeMakerRun3Monitor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;

  /*
  Handle<vector<ScoutingVertex> > verticesH;
  iEvent.getByToken(verticesToken, verticesH);
  */

  Handle<vector<pat::Muon> > offlineMuonsH;
  iEvent.getByToken(offlineMuonsToken, offlineMuonsH);
  
  if (offlineMuonsH->size()<2) return;

  int nOfflineMuons=0;
  vector<int> idx;
  int j=0;
  for (auto muons_iter = offlineMuonsH->begin(); muons_iter != offlineMuonsH->end(); ++muons_iter) {
      //cout<<"offline muon pt: "<<muons_iter->pt()<<" id: "<<muons_iter->isMediumMuon()<<" eta: "<<muons_iter->eta()<<endl;
      if (muons_iter->pt()>4 && muons_iter->isMediumMuon() and abs(muons_iter->eta())<1.9) {
          nOfflineMuons+=1;
          idx.push_back(j);
      }
      j+=1;
  }
  
  if (idx.size()<2) {/*cout<<"failed offline muons"<<endl;*/ return;}

  edm::Handle<edm::TriggerResults> triggerResultsH;
  iEvent.getByToken(triggerResultsToken, triggerResultsH);

  bool passDST=false;
  for (size_t i = 0; i < triggerPathsVector.size(); i++) {
      if (triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) passDST=true;
  }

  if (!passDST) {/*cout<<"failed DST"<<endl;*/ return;}

  l1Result_.clear();
  nScoutingMuons=0;

  bool passMonitor=false;
  if (doL1) {
      l1GtUtils_->retrieveL1(iEvent,iSetup,algToken_);
      /*for(unsigned int r = 0; r<100; r++){
        string name ("empty");
                bool algoName_ = false;
                        algoName_ = l1GtUtils_->getAlgNameFromBit(i,name);
                        cout << "getAlgNameFromBit = " << algoName_  << endl;
                        cout << "L1 bit number = " << i << " ; L1 bit name = " << name << endl;
                        }*/
      for( unsigned int iseed = 0; iseed < l1Seeds_.size(); iseed++ ) {
          bool l1htbit = 0;
          l1GtUtils_->getFinalDecisionByName(string(l1Seeds_[iseed]), l1htbit);
          l1Result_.push_back( l1htbit );
      }
      for( unsigned int iseed = 0; iseed < l1MonitorSeeds_.size(); iseed++ ) {
          bool l1htbit = 0;
          l1GtUtils_->getFinalDecisionByName(string(l1MonitorSeeds_[iseed]), l1htbit);
          if (l1htbit) passMonitor=true;
      }
  }

  if (!passMonitor) {/*cout<<"failed L1 seed"<<endl;*/ return;}
    
  Handle<double> rhoH;
  iEvent.getByToken(rhoToken, rhoH);

  cout<<"passed all selections"<<endl;

  pt1=offlineMuonsH->at(idx[0]).pt();
  pt2=offlineMuonsH->at(idx[1]).pt();
  eta1=offlineMuonsH->at(idx[0]).eta();
  eta2=offlineMuonsH->at(idx[1]).eta();
  float phi1=offlineMuonsH->at(idx[0]).phi();
  float phi2=offlineMuonsH->at(idx[1]).phi();
  id1=offlineMuonsH->at(idx[0]).pdgId();
  id2=offlineMuonsH->at(idx[1]).pdgId();

  TLorentzVector mu1;
  mu1.SetPtEtaPhiM(pt1,eta1,phi1,0.105658);

  TLorentzVector mu2;
  mu2.SetPtEtaPhiM(pt2,eta2,phi2,0.105658);

  TLorentzVector dimu = mu1+mu2;
  mass=dimu.M();
  pt=dimu.Pt();
  dr=mu1.DeltaR(mu2);

  rho=*rhoH;

  Handle<vector<Run3ScoutingMuon> > muonsH;
  iEvent.getByToken(muonsToken, muonsH);

  nScoutingMuons=0;
  for (auto muons_iter = muonsH->begin(); muons_iter != muonsH->end(); ++muons_iter) {
      if (muons_iter->pt()>4 and abs(muons_iter->eta())<1.9) {
          nScoutingMuons+=1;
      }
  }

  tree->Fill();

}

// ------------ method called once each job just before starting event loop  ------------
void ScoutingTreeMakerRun3Monitor::beginJob() {
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("tree"       , "tree");
    tree->Branch("mass"                 , &mass                          , "mass/F");
    tree->Branch("pt"                 , &pt                          , "pt/F");
    tree->Branch("dr"                 , &dr                          , "dr/F");
    tree->Branch("pt1"                 , &pt1                          , "pt1/F");
    tree->Branch("pt2"                 , &pt2                          , "pt2/F");
    tree->Branch("eta1"                 , &eta1                          , "eta1/F");
    tree->Branch("eta2"                 , &eta2                          , "eta2/F");
    tree->Branch("id1"                 , &id1                          , "id1/I");
    tree->Branch("id2"                 , &id2                          , "id2/I");
    tree->Branch("rho"                 , &rho                          , "rho/F");
    tree->Branch("l1Result", "std::vector<bool>"             ,&l1Result_, 32000, 0);
    tree->Branch("nScoutingMuons", &nScoutingMuons             , "nScoutingMuons/I");

}

// ------------ method called once each job just after ending the event loop  ------------
void ScoutingTreeMakerRun3Monitor::endJob() {
  // please remove this method if not needed
}

void ScoutingTreeMakerRun3Monitor::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
    // HLT paths
    triggerPathsVector.push_back("DST_Run3_PFScoutingPixelTracking_v*");

    HLTConfigProvider hltConfig;
    bool changedConfig = false;
    hltConfig.init(iRun, iSetup, triggerResultsTag.process(), changedConfig);

    for (size_t i = 0; i < triggerPathsVector.size(); i++) {
        triggerPathsMap[triggerPathsVector[i]] = -1;
    }

    for(size_t i = 0; i < triggerPathsVector.size(); i++){
        TPRegexp pattern(triggerPathsVector[i]);
        for(size_t j = 0; j < hltConfig.triggerNames().size(); j++){
            std::string pathName = hltConfig.triggerNames()[j];
            if(TString(pathName).Contains(pattern)){
                triggerPathsMap[triggerPathsVector[i]] = j;
            }
        }
    }
}


void ScoutingTreeMakerRun3Monitor::endRun(edm::Run const&, edm::EventSetup const&) {
}

void ScoutingTreeMakerRun3Monitor::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void ScoutingTreeMakerRun3Monitor::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ScoutingTreeMakerRun3Monitor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ScoutingTreeMakerRun3Monitor);
