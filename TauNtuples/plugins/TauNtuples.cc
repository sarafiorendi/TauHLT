// -*- C++ -*-
//
// Package:    tauHLT/TauNtuples
// Class:      TauNtuples
//
/**\class TauNtuples TauNtuples.cc tauHLT/TauNtuples/plugins/TauNtuples.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Sara Fiorendi
//         Created:  Thu, 18 Nov 2021 21:51:35 GMT
//
//

// system include files
#include <memory>
#include <vector>       

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Common/interface/TriggerNames.h"
// #include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/TauReco/interface/PFTauTransverseImpactParameterAssociation.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TTree.h"

#include "TauHLT/TauNtuples/src/tauTree.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class TauNtuples : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
 
    using TauIPCollection = edm::AssociationVector<reco::PFTauRefProd,
                                                   std::vector<reco::PFTauTransverseImpactParameterRef>>;

    explicit TauNtuples(const edm::ParameterSet&);
    ~TauNtuples();

    virtual void beginEvent() ;
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    void fillL1Muons(const edm::Handle<l1t::MuonBxCollection> &,
                     const edm::Event   &
                    );

    void fillL1Taus(const edm::Handle<l1t::TauBxCollection> &,  
                    const edm::Event   &
                   );

    void fillHltMuons(const edm::Handle<reco::RecoChargedCandidateCollection> &,
                      const edm::Event      &
                     );
    void fillHltTaus(const edm::Handle<reco::PFTauCollection> &,
                     const edm::Event      &,
                     const edm::EDGetTokenT<TauIPCollection> &,
                     const edm::EDGetTokenT<reco::PFTauDiscriminator> &,
                     const edm::EDGetTokenT<reco::PFTauDiscriminator> &
                    );
    void fillHlt(const edm::Handle<edm::TriggerResults> &, 
                 const edm::Handle<trigger::TriggerEvent> &,
                 const edm::TriggerNames &,
                 const edm::Event &,
                 bool 
                );
  private:
    void beginJob() override;
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void endJob() override;
      

  // ----------member data ---------------------------
    edm::Service<TFileService> outfile_;

    edm::EDGetTokenT<l1t::MuonBxCollection> l1mucandToken_; 
    edm::EDGetTokenT<l1t::TauBxCollection> l1taucandToken_; 

    edm::EDGetTokenT<reco::RecoChargedCandidateCollection> hltmucandToken_; 
    edm::EDGetTokenT<reco::PFTauCollection> hlttaucandToken_; 
    edm::EDGetTokenT<TauIPCollection> hlttauIPToken_; 
    edm::EDGetTokenT<reco::PFTauDiscriminator> hlttauIsoToken_; 
    edm::EDGetTokenT<reco::PFTauDiscriminator> hlttauIsoValueToken_; 

    edm::EDGetTokenT<edm::TriggerResults>   triggerResultToken_;
    edm::EDGetTokenT<trigger::TriggerEvent> triggerSummToken_;

    edm::EDGetTokenT<LumiScalersCollection> lumiScalerToken_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken_;
    
    TauEvent event_;
    std::map<std::string,TTree*> tree_;

};



TauNtuples::TauNtuples(const edm::ParameterSet& iConfig)
    : 
      l1mucandToken_ (consumes<l1t::MuonBxCollection>(iConfig.getUntrackedParameter<edm::InputTag>("L1MuonCandidates"))),
      l1taucandToken_(consumes<l1t::TauBxCollection>(iConfig.getUntrackedParameter<edm::InputTag>("L1TauCandidates"))),
      hltmucandToken_(consumes<reco::RecoChargedCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("hltMuCandidates"))),
      hlttaucandToken_(consumes<reco::PFTauCollection>(iConfig.getUntrackedParameter<edm::InputTag>("hltTauCandidates"))),
      hlttauIPToken_(consumes<TauIPCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tauIP"))),
      hlttauIsoToken_(consumes<reco::PFTauDiscriminator>(iConfig.getUntrackedParameter<edm::InputTag>("tauIso"))),
      hlttauIsoValueToken_(consumes<reco::PFTauDiscriminator>(iConfig.getUntrackedParameter<edm::InputTag>("tauIsoValue"))),
      triggerResultToken_(consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("triggerResult"))),
      triggerSummToken_(consumes<trigger::TriggerEvent>(iConfig.getUntrackedParameter<edm::InputTag>("triggerSummary"))),
      lumiScalerToken_(consumes<LumiScalersCollection>(iConfig.getUntrackedParameter<edm::InputTag>("lumiScalerTag"))),
      puToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getUntrackedParameter<edm::InputTag>("puInfoTag")))
     {}

TauNtuples::~TauNtuples(){}



// ------------ method called for each event  ------------
void TauNtuples::analyze(const edm::Event& event, const edm::EventSetup& iSetup) {
  
  using namespace edm;
  
  beginEvent();

  // Fill general info
  event_.runNumber             = event.id().run();
  event_.luminosityBlockNumber = event.id().luminosityBlock();
  event_.eventNumber           = event.id().event();

  // Handle the L1 collections and fill info
  edm::Handle<l1t::MuonBxCollection> mul1cands;
  if (event.getByToken(l1mucandToken_, mul1cands))
    fillL1Muons(mul1cands, event);
  else
    edm::LogWarning("") << "Online L1 muon collection not found !!!";

  edm::Handle<l1t::TauBxCollection> taul1cands;
  if (event.getByToken(l1taucandToken_, taul1cands))
    fillL1Taus(taul1cands, event);
  else
    edm::LogWarning("") << "Online L1 tau collection not found !!!";
  
  // Handle the hlt muon collection and fill online muons
  edm::Handle<reco::RecoChargedCandidateCollection> l3mucands;
  if (event.getByToken(hltmucandToken_, l3mucands))
    fillHltMuons(l3mucands, event); 
  else
    edm::LogWarning("") << "Online L3 muons collection not found !!!";

  // Handle the hlt taus collection and fill online taus
  edm::Handle<reco::PFTauCollection> hlttaucands;
  if (event.getByToken(hlttaucandToken_, hlttaucands))
    fillHltTaus(hlttaucands, event, hlttauIPToken_, hlttauIsoToken_, hlttauIsoValueToken_); 
  else
    edm::LogWarning("") << "Online PF taus collection not found !!!";

  // Fill trigger information for the new trigger
  edm::Handle<edm::TriggerResults>   triggerResults;
  edm::Handle<trigger::TriggerEvent> triggerEvent;

  if (event.getByToken(triggerResultToken_, triggerResults) &&
      event.getByToken(triggerSummToken_  , triggerEvent)) {
      
    edm::TriggerNames triggerNames_ = event.triggerNames(*triggerResults);
    fillHlt(triggerResults, triggerEvent, triggerNames_, event, false);
  }
  else 
    edm::LogError("") << "New trigger collection not found !!!";

  // Fill bx and inst lumi info
  if (event.isRealData()) {
    event_.bxId  = event.bunchCrossing();
    edm::Handle<LumiScalersCollection> lumiScaler;
    if (event.getByToken(lumiScalerToken_, lumiScaler) && lumiScaler->begin() != lumiScaler->end())
        event_.instLumi = lumiScaler->begin()->instantLumi(); 
  }
  
  // Fill PU info
  if (!event.isRealData()) {
    edm::Handle<std::vector< PileupSummaryInfo > > puInfo;
    if (event.getByToken(puToken_, puInfo)){
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      for(PVI = puInfo->begin(); PVI != puInfo->end(); ++PVI) {
        if(PVI->getBunchCrossing()==0) {
          event_.trueNI   = PVI->getTrueNumInteractions();
          continue;
        }
      }
    } 
    else  
      edm::LogError("") << "PU collection not found !!!";
  }


  // endEvent();
  tree_["tauTree"] -> Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void TauNtuples::beginJob() {
  tree_["tauTree"] = outfile_-> make<TTree>("tauTree","tauTree");
  tree_["tauTree"] -> Branch("event" ,&event_, 64000,2);
}

// ------------ method called once each job just after ending the event loop  ------------
void TauNtuples::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TauNtuples::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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

// ---------------------------------------------------------------------
void TauNtuples::fillHltTaus( const edm::Handle<reco::PFTauCollection> & taucands ,
                              const edm::Event                         & event     ,
                              const edm::EDGetTokenT<TauIPCollection>  & tauIPToken_,
                              const edm::EDGetTokenT<reco::PFTauDiscriminator>  & tauIsoToken_,
                              const edm::EDGetTokenT<reco::PFTauDiscriminator>  & tauIsoValueToken_
                               )
{

  edm::Handle<TauIPCollection> tauIPPars;
  edm::Handle<reco::PFTauDiscriminator> tauIso;
  edm::Handle<reco::PFTauDiscriminator> tauIsoValue;
  for (size_t idx = 0; idx < taucands->size(); ++idx) {
    reco::PFTauRef PFTauRef(taucands, idx);
    
    HLTTauCand theTau;

    theTau.pt      = PFTauRef->pt();
    theTau.eta     = PFTauRef->eta();
    theTau.phi     = PFTauRef->phi();
    theTau.charge  = PFTauRef->charge();
    theTau.decaymode = PFTauRef->decayMode();  
// //     theL3Mu.dz = candref -> vz();

    // tau components
    theTau.nChargedHad = PFTauRef->signalChargedHadrCands().size() ;
    theTau.nGamma      = PFTauRef->signalGammaCands().size() ;
    theTau.leadChargedCandPdgId = PFTauRef->leadChargedHadrCand()->pdgId(); 
    theTau.leadChargedCandPt = PFTauRef->leadChargedHadrCand()->pt(); 
    theTau.leadCandPdgId = PFTauRef->leadCand()->pdgId(); 
    theTau.leadCandPt    = PFTauRef->leadCand()->pt(); 

    float sum_pt_charged = 0;
    for (auto ich : PFTauRef->signalChargedHadrCands())
        sum_pt_charged += ich->pt();
    theTau.sum_pt_charged = sum_pt_charged;

    float sum_pt_neutral = 0;
    for (auto ineu : PFTauRef->signalGammaCands())
        sum_pt_neutral += ineu->pt();
    theTau.sum_pt_neutral = sum_pt_neutral;

//   Float_t maxHCALPFClusterEt;
    
    // tau IP related info
    if (event.getByToken(tauIPToken_, tauIPPars)){
    
      const reco::PFTauTransverseImpactParameter& tauIPInfo = *(*tauIPPars)[PFTauRef];
      theTau.dxy    = tauIPInfo.dxy() ;
      theTau.dxyerr = tauIPInfo.dxy_error() ;

      theTau.ip3d       = tauIPInfo.ip3d() ;
      theTau.sigmaip3d  = tauIPInfo.ip3d_error() ;
      theTau.hasSV      = tauIPInfo.hasSecondaryVertex() ;
      theTau.sv_x       = tauIPInfo.secondaryVertexPos().x() ;
      theTau.sv_y       = tauIPInfo.secondaryVertexPos().y() ;
      theTau.sv_z       = tauIPInfo.secondaryVertexPos().z() ;
      theTau.l_x        = tauIPInfo.flightLength().x() ;
      theTau.l_y        = tauIPInfo.flightLength().y() ;
      theTau.l_z        = tauIPInfo.flightLength().z() ;
      theTau.l_sigma    = tauIPInfo.flightLengthSig() ;
    }
    else{
      theTau.dxy       = -9999;
      theTau.dxyerr    = -9999;
      theTau.ip3d      = -9999;
      theTau.sigmaip3d = -9999;
      theTau.hasSV     = -9999;
      theTau.sv_x      = -9999;
      theTau.sv_y      = -9999;
      theTau.sv_z      = -9999;
      theTau.l_x       = -9999;
      theTau.l_y       = -9999;
      theTau.l_z       = -9999;
      theTau.l_sigma   = -9999;
    }

    // tau iso info
    theTau.passChargedIso = -1;
    if (event.getByToken(tauIsoToken_, tauIso)){
        if ((*tauIso)[PFTauRef] < 0.5) 
            theTau.passChargedIso = 1;
        else    
            theTau.passChargedIso = 0;
    }
    theTau.chargedIso = -9999.;
    if (event.getByToken(tauIsoToken_, tauIsoValue)){            
        theTau.chargedIso = (*tauIso)[PFTauRef];
    }
            
    event_.hlttaus   .push_back(theTau);
  }
}

// ---------------------------------------------------------------------
void TauNtuples::fillHltMuons(const edm::Handle<reco::RecoChargedCandidateCollection> & l3cands ,
                               const edm::Event                                        & event
                               )
{

//   edm::Handle<reco::IsoDepositMap> trkDepMap;
//   edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMap;
//   edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMap;

  for( unsigned int il3 = 0; il3 < l3cands->size(); ++il3) {
    HLTMuonCand theL3Mu;

    reco::RecoChargedCandidateRef candref(l3cands, il3);
    theL3Mu.pt      = candref -> pt();
    theL3Mu.eta     = candref -> eta();
    theL3Mu.phi     = candref -> phi();
    theL3Mu.charge  = candref -> charge();
//     theL3Mu.dxy     = candref -> dxy();  // to be checked
    
    theL3Mu.dz = candref -> vz();

    event_.hltmuons   .push_back(theL3Mu);
  }
}

// ---------------------------------------------------------------------
void TauNtuples::fillL1Muons(const edm::Handle<l1t::MuonBxCollection> & l1cands ,
                              const edm::Event                         & event    
                              )
{

  for (int ibx = l1cands->getFirstBX(); ibx <= l1cands->getLastBX(); ++ibx) {
    if (ibx != 0) continue;
    for (auto it = l1cands->begin(ibx); it != l1cands->end(ibx); it++) {

      l1t::MuonRef muon(l1cands, distance(l1cands->begin(l1cands->getFirstBX()),it) );

      L1MuonCand theL1Mu;
      theL1Mu.pt       = muon -> pt();
      theL1Mu.eta      = muon -> eta();
      theL1Mu.phi      = muon -> phi();
      theL1Mu.charge   = muon -> charge();
      theL1Mu.quality  = muon -> hwQual();

      event_.L1muons.push_back(theL1Mu);
    }
  }
}

//---------------------------------------------------

void TauNtuples::fillL1Taus(const edm::Handle<l1t::TauBxCollection> & l1cands ,
                            const edm::Event                         & event    
                           )
{

  for (int ibx = l1cands->getFirstBX(); ibx <= l1cands->getLastBX(); ++ibx) {
    if (ibx != 0) continue;
    for (auto it = l1cands->begin(ibx); it != l1cands->end(ibx); it++) {

      l1t::TauRef tau(l1cands, distance(l1cands->begin(l1cands->getFirstBX()),it) );

      L1TauCand theL1Tau;
      theL1Tau.pt       = tau -> pt();
      theL1Tau.eta      = tau -> eta();
      theL1Tau.phi      = tau -> phi();
      theL1Tau.charge   = tau -> charge();
      theL1Tau.iso      = tau -> hwIso();

      event_.L1taus.push_back(theL1Tau);
    }
  }
}

//---------------------------------------------------


void TauNtuples::fillHlt(const edm::Handle<edm::TriggerResults>    & triggerResults, 
                          const edm::Handle<trigger::TriggerEvent> & triggerEvent  ,
                          const edm::TriggerNames                  & triggerNames  ,
                          const edm::Event                         & event         ,
                          bool                                       isTag         )
{    
   
  for (unsigned int itrig=0; itrig < triggerNames.size(); ++itrig) {
    LogDebug ("triggers") << triggerNames.triggerName(itrig) ;
    if (triggerResults->accept(itrig)){
      std::string pathName = triggerNames.triggerName(itrig);
      if ( pathName.find ("HLT_IsoMu"  ) !=std::string::npos ||
           pathName.find ("HLT_Mu45"   ) !=std::string::npos ||
           pathName.find ("HLT_Mu"     ) !=std::string::npos ||
           pathName.find ("HLT_TkMu"   ) !=std::string::npos ||
           pathName.find ("HLT_Mu17"   ) !=std::string::npos ||
           pathName.find ("HLT_L2Mu"   ) !=std::string::npos ||
           pathName.find ("HLT_DoubleDisplaced" ) !=std::string::npos ||
           pathName.find ("HLT_DoubleTight") !=std::string::npos ||
           pathName.find ("HLT_DoubleMedium"  ) !=std::string::npos
      ){
        if (isTag) event_.hltTag.triggers.push_back(pathName);
        else       event_.hlt   .triggers.push_back(pathName);
      }
    }
  }
     
     
  const trigger::size_type nFilters(triggerEvent->sizeFilters());
  for (trigger::size_type iFilter=0; iFilter!=nFilters; ++iFilter) {
    std::string filterTag = triggerEvent->filterTag(iFilter).encode();

    if ( ( filterTag.find ("sMu"     ) !=std::string::npos ||
           filterTag.find ("SingleMu") !=std::string::npos ||
           filterTag.find ("SingleMu") !=std::string::npos ||
           filterTag.find ("DiMuon"  ) !=std::string::npos 
           ) &&
           filterTag.find ("Tau"       ) ==std::string::npos   &&
           filterTag.find ("EG"        ) ==std::string::npos   &&
           filterTag.find ("MultiFit"  ) ==std::string::npos
       )
    {
      std::string filterTag = triggerEvent->filterTag(iFilter).encode();
  
      trigger::Keys objectKeys = triggerEvent->filterKeys(iFilter);
      const trigger::TriggerObjectCollection& triggerObjects(triggerEvent->getObjects());
      
      for (trigger::size_type iKey=0; iKey<objectKeys.size(); ++iKey) {  
        trigger::size_type objKey = objectKeys.at(iKey);
        const trigger::TriggerObject& triggerObj(triggerObjects[objKey]);
        
        HLTObjCand hltObj;
        
        hltObj.filterTag = filterTag;
  
        hltObj.pt  = triggerObj.pt();
        hltObj.eta = triggerObj.eta();
        hltObj.phi = triggerObj.phi();
        
        if (isTag)       event_.hltTag.objects.push_back(hltObj);
        else             event_.hlt   .objects.push_back(hltObj);
        
      }
    }         
  }
  
  // fill hlt rho information
//   if (!isTag){
//     edm::Handle <double>  hltRhoCollection;
//     if (event.getByToken(rhoCorrectionToken_, hltRhoCollection) && hltRhoCollection.isValid()){
//        event_.hlt   .rho = *(hltRhoCollection.product());
//     } 
//   }
}

//---------------------------------------------------

void TauNtuples::beginEvent()
{

  event_.hlt.triggers.clear();
  event_.hlt.objects.clear();
  event_.hlt.rho        = -1;

  event_.hltTag.triggers.clear();
  event_.hltTag.objects.clear();
  event_.hltTag.rho = -1;
// 
//   event_.genParticles.clear();
//   event_.muons.clear();
  event_.hlttaus.clear();
  event_.hltmuons.clear();
  event_.L1taus.clear();
  event_.L1muons.clear();

//   event_.L2muons.clear();

//   for (unsigned int ix=0; ix<3; ++ix) {
//     event_.primaryVertex[ix] = 0.;
//     for (unsigned int iy=0; iy<3; ++iy) {
//       event_.cov_primaryVertex[ix][iy] = 0.;
//     }
//   }

  event_.nVtx       = -1;
  event_.trueNI     = -1;
  event_.rho        = -1;
  event_.bxId       = -1;
  event_.instLumi   = -1;
  
//   nGoodVtx = 0; 
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauNtuples);
