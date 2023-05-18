// system include files
#include <memory>
#include <TTree.h>
#include <TLorentzVector.h>

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

#include "DataFormats/Scouting/interface/ScoutingElectron.h"
#include "DataFormats/Scouting/interface/ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/ScoutingVertex.h"
#include "DataFormats/Scouting/interface/ScoutingTrack.h"
#include "DataFormats/Scouting/interface/ScoutingMuon.h"
#include "DataFormats/Scouting/interface/ScoutingParticle.h"

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

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class TestScoutingTreeMakerRun2 : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TestScoutingTreeMakerRun2(const edm::ParameterSet&);
  ~TestScoutingTreeMakerRun2() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
    //void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //void endRun(edm::Run const&, edm::EventSetup const&) override;
  //void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;   

  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<edm::TriggerResults>             triggerResultsToken;
  const edm::EDGetTokenT<std::vector<ScoutingMuon> >      muonsToken;
  const edm::EDGetTokenT<std::vector<ScoutingElectron> >  electronsToken;
  const edm::EDGetTokenT<std::vector<ScoutingVertex> >    primaryVerticesToken;
  const edm::EDGetTokenT<std::vector<ScoutingVertex> >    verticesToken;
  const edm::EDGetTokenT<std::vector<ScoutingPhoton> >  photonsToken;
  const edm::EDGetTokenT<std::vector<ScoutingTrack> >  tracksToken;

  std::vector<std::string> triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;

  bool doL1;
  triggerExpression::Data triggerCache_;

  edm::InputTag                algInputTag_;
  edm::InputTag                extInputTag_;
  edm::EDGetToken              algToken_;
  std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
  std::vector<std::string>     l1Seeds_;
  std::vector<bool>            l1Result_;

  TTree* tree;
  float mass;
  float pt;
  float dr;
  float pt1;
  float pt2;
  float eta1;
  float eta2;

  int nMuonsID;
  bool signCheck;

  float avgPrimaryX;
  float avgPrimaryY;

  float lxy;
  float vtxErr;
  float IP;

  float lxy2;
  float vtxErr2;
  float IP2;

  int nvtx;
  int ndvtx;
  bool vtxMatch;

  std::vector<float> vtxX;
  std::vector<float> vtxY;
  std::vector<float> vtxXError;
  std::vector<float> vtxYError;

  const float BSx = 0.0965;
  const float BSy = -0.062;
  const float phi0 = -0.571;
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
TestScoutingTreeMakerRun2::TestScoutingTreeMakerRun2(const edm::ParameterSet& iConfig):
    triggerResultsTag        (iConfig.getParameter<edm::InputTag>("triggerresults")),
    triggerResultsToken      (consumes<edm::TriggerResults>                    (triggerResultsTag)),
    muonsToken               (consumes<std::vector<ScoutingMuon> >             (iConfig.getParameter<edm::InputTag>("muons"))),
    electronsToken           (consumes<std::vector<ScoutingElectron> >         (iConfig.getParameter<edm::InputTag>("electrons"))),
    primaryVerticesToken     (consumes<std::vector<ScoutingVertex> >           (iConfig.getParameter<edm::InputTag>("primaryVertices"))),
    verticesToken            (consumes<std::vector<ScoutingVertex> >           (iConfig.getParameter<edm::InputTag>("displacedVertices"))),
    photonsToken             (consumes<std::vector<ScoutingPhoton> >         (iConfig.getParameter<edm::InputTag>("photons"))),
    tracksToken              (consumes<std::vector<ScoutingTrack> >            (iConfig.getParameter<edm::InputTag>("tracks"))),
    doL1                     (iConfig.existsAs<bool>("doL1")               ?    iConfig.getParameter<bool>  ("doL1")            : false)
{
    usesResource("TFileService");
    if (doL1) {
        algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
        extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
        algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
        l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
        l1GtUtils_ = std::make_unique<l1t::L1TGlobalUtil>(iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);
    }
    else {
        l1Seeds_ = std::vector<std::string>();
        l1GtUtils_ = 0;
    }
}

TestScoutingTreeMakerRun2::~TestScoutingTreeMakerRun2() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void TestScoutingTreeMakerRun2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;

  Handle<vector<ScoutingMuon> > muonsH;
  iEvent.getByToken(muonsToken, muonsH);

  if (muonsH->size()<2) return;

  int nMuons=0;
  nMuonsID=0;
  vector<int> idx;

  int j=0;
  for (auto muons_iter = muonsH->begin(); muons_iter != muonsH->end(); ++muons_iter) {
      //std::cout<<"pt: "<<muons_iter->pt()<<std::endl;                                             
      //std::cout<<"trkiso: "<<muons_iter->trackIso()<<" pix hits: "<< muons_iter->nValidPixelHits()<<" layers: "<<muons_iter->nTrackerLayersWithMeasurement()<<" trk chi2: "<<muons_iter->trk_chi2()<<std::endl; 

      if (muons_iter->pt()>4) {
          nMuons+=1;
          if ((muons_iter->trackIso()<0.15) &&
              (muons_iter->nValidPixelHits()>0) &&
              (muons_iter->nTrackerLayersWithMeasurement()>5)&&
              (muons_iter->chi2()<10)) {
              nMuonsID+=1;
          }
      }

      idx.push_back(j);
      j+=1;
  }

  std::cout<<std::endl<<idx.size()<<std::endl;

  if (idx.size()>1) {
      if ((muonsH->at(idx[0]).charge())*(muonsH->at(idx[1]).charge()) > 0) {
        std::cout<<"Fail sign check"<<std::endl;
        return;
      }

      pt1=muonsH->at(idx[0]).pt();
      pt2=muonsH->at(idx[1]).pt();
      signCheck=((muonsH->at(idx[0]).charge())*(muonsH->at(idx[1]).charge()) < 0);

      eta1=muonsH->at(idx[0]).eta();
      eta2=muonsH->at(idx[1]).eta();      
      float phi1=muonsH->at(idx[0]).phi();
      float phi2=muonsH->at(idx[1]).phi();
      
      TLorentzVector mu1;
      mu1.SetPtEtaPhiM(pt1,eta1,phi1,0.105658);

      TLorentzVector mu2;
      mu2.SetPtEtaPhiM(pt2,eta2,phi2,0.105658);

      TLorentzVector dimu = mu1+mu2;
      mass=dimu.M();
      pt=dimu.Pt();
      dr=mu1.DeltaR(mu2);

      std::cout<<"pt1: "<<pt1<<" pt2: "<<pt2<<" nMuonsID: "<<nMuonsID<<std::endl;
      std::cout<<"mass: "<<mass<<" pt: "<<pt<<" charge0: "<<(muonsH->at(idx[0]).charge())<<" charge1: "<<(muonsH->at(idx[1]).charge())<<" charge product: "<<(muonsH->at(idx[0]).charge())*(muonsH->at(idx[1]).charge())<<std::endl;

      Handle<vector<ScoutingVertex> > primaryVerticesH;
      iEvent.getByToken(primaryVerticesToken, primaryVerticesH);

      vtxX.clear();
      vtxY.clear();
      vtxXError.clear();
      vtxYError.clear();

      nvtx = 0;
      for (auto vtx_iter = primaryVerticesH->begin(); vtx_iter != primaryVerticesH->end(); ++vtx_iter) {
        //std::cout<<"primary x: "<<vtx_iter->x() <<" y: "<<  vtx_iter->y()<<" ex: "<<vtx_iter->xError()<<" ey: "<<vtx_iter->yError()<<std::endl;
        nvtx++;
        vtxX.push_back(vtx_iter->x());   
        vtxY.push_back(vtx_iter->y());        
        vtxXError.push_back(vtx_iter->xError());        
        vtxYError.push_back(vtx_iter->yError());           
      }

      avgPrimaryX = 0;
      if (!vtxX.empty()) {
        avgPrimaryX = std::reduce(vtxX.begin(), vtxX.end(), 0.0) / vtxX.size();
      }

      avgPrimaryY = 0;
      if (!vtxY.empty()) {
        avgPrimaryY = std::reduce(vtxY.begin(), vtxY.end(), 0.0) / vtxY.size();
      }

      std::cout<<"avgPrimaryX: "<<avgPrimaryX<<" avgPrimaryY: "<<avgPrimaryY<<std::endl;

      Handle<vector<ScoutingVertex> > verticesH;
      iEvent.getByToken(verticesToken, verticesH);


      std::vector<int> vtxIndx1 = (muonsH->at(idx[0])).vtxIndx();
      std::vector<int> vtxIndx2 = (muonsH->at(idx[1])).vtxIndx();

      std::cout<<"vtxIndx1 size: "<<vtxIndx1.size()<<" vtxIndx2 size: "<<vtxIndx2.size()<<" num vtx: "<<verticesH->size()<<std::endl;

      lxy = 0;
      vtxErr = 0;
      IP = 0;

      lxy2 = 0;
      vtxErr2 = 0;
      IP2 = 0;

      ndvtx = verticesH->size();
      vtxMatch = vtxIndx1.size() > 0 && vtxIndx2.size() > 0 && vtxIndx1[0] == 0 && vtxIndx2[0] == 0 && ndvtx > 0;

      if (vtxMatch) {

        std::cout<<"vtxIndx12: "<<(muonsH->at(idx[0])).vtxIndx()[0]<<" "<<(muonsH->at(idx[1])).vtxIndx()[0]<<endl;

        auto vtx = verticesH->begin();
        double dx = (vtx->x()) - avgPrimaryX;
        double dy = (vtx->y()) - avgPrimaryY;

        lxy = sqrt(dx*dx + dy*dy);
        vtxErr = sqrt(dx*dx*(vtx->xError())*(vtx->xError()) + dy*dy*(vtx->yError())*(vtx->yError())) / lxy;
        IP = lxy/vtxErr;

        double dx2 = (vtx->x()) - BSx;
        double dy2 = (vtx->y()) - BSy;

        lxy2 = sqrt(dx2*dx2 + dy2*dy2);
        vtxErr2 = sqrt(dx2*dx2*(vtx->xError())*(vtx->xError()) + dy2*dy2*(vtx->yError())*(vtx->yError())) / lxy2;
        IP2 = lxy2/vtxErr2;
      }
      
      std::cout << "lxy: " << lxy << " lxy2: " << lxy2 << std::endl;

      std::cout<<"tree filling with mass: "<<mass<<", pt: "<<pt<<std::endl;
      tree->Fill();
  }
}

// ------------ method called once each job just before starting event loop  ------------
void TestScoutingTreeMakerRun2::beginJob() {
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("tree"      , "tree");
    tree->Branch("mass"                , &mass                        , "mass/F"    );
    tree->Branch("pt"                  , &pt                          , "pt/F"      );
    tree->Branch("dr"                  , &dr                          , "dr/F"      );
    tree->Branch("pt1"                 , &pt1                         , "pt1/F"     );
    tree->Branch("pt2"                 , &pt2                         , "pt2/F"     );
    tree->Branch("eta1"                , &eta1                        , "eta1/F"    );
    tree->Branch("eta2"                , &eta2                        , "eta2/F"    );
    tree->Branch("nvtx"                , &nvtx                        , "nvtx/I"    );
    tree->Branch("ndvtx"               , &ndvtx                       , "ndvtx/I"   );
    tree->Branch("vtxMatch"            , &vtxMatch                    , "vtxMatch/B");
    tree->Branch("lxy"                 , &lxy                         , "lxy/F"     );
    tree->Branch("vtxErr"              , &vtxErr                      , "vtxErr/F"  );
    tree->Branch("IP"                  , &IP                          , "IP/F"      );
    tree->Branch("lxy2"                 , &lxy2                        , "lxy2/F"     );
    tree->Branch("vtxErr2"              , &vtxErr2                      , "vtxErr2/F"  );
    tree->Branch("IP2"                  , &IP2                          , "IP2/F"      );
    tree->Branch("nMuonsID"            , &nMuonsID                    , "nMuonsID/I");

}

// ------------ method called once each job just after ending the event loop  ------------
void TestScoutingTreeMakerRun2::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TestScoutingTreeMakerRun2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(TestScoutingTreeMakerRun2);
