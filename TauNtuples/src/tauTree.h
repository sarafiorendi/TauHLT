#ifndef  tauTree_h
#define  tauTree_h

#include "TROOT.h"
#include "TMath.h"
#include <vector>
#include <string>



// class GenParticleCand {
// public:
//   Int_t   pdgId; 
//   Int_t   status; 
//   Float_t energy; 
//   Float_t pt; 
//   Float_t eta; 
//   Float_t phi; 
//   std::vector<Int_t>  pdgMother; 
//   std::vector<Int_t>  pdgRealMother; 
// 
//   GenParticleCand(){};
//   virtual ~GenParticleCand(){};
//   
//   ClassDef(GenParticleCand,1)
// };
// 


class HLTTauCand {

public:

  Float_t pt;           
  Float_t eta;          
  Float_t phi;          
  Float_t dxy;
  Float_t dxyerr;
  Int_t   charge;         
  Int_t   decaymode;         

  Int_t leadChargedCandPdgId;
  Int_t leadCandPdgId;
  Int_t nChargedHad;
  Int_t nGamma;
  
  Float_t leadChargedCandPt;
  Float_t leadCandPt;
  Float_t sum_pt_charged;
  Float_t sum_pt_neutral;
  Float_t maxHCALPFClusterEt;

  Float_t ip3d      ;
  Float_t sigmaip3d ;
  Float_t sv_x      ;
  Float_t sv_y      ;
  Float_t sv_z      ;
  Float_t l_x       ;
  Float_t l_y       ;
  Float_t l_z       ;
  Float_t l_sigma   ;

  Int_t hasSV       ;
  
  Int_t passChargedIso;
  Float_t chargedIso;


  HLTTauCand(){};
  virtual ~HLTTauCand(){};

  ClassDef(HLTTauCand,1);   
};

class HLTMuonCand {
public:

  Float_t pt;           
  Float_t eta;          
  Float_t phi;          
  Float_t dxy;          
  Float_t dz;
//   Float_t trkpt;         // pt of the track of the hlt muon [GeV]
  Int_t   charge;         // pt of the track of the hlt muon [GeV]
  
//   Float_t ecalDep;
//   Float_t hcalDep;
//   Float_t trkDep;


  HLTMuonCand(){};
  virtual ~HLTMuonCand(){};

  ClassDef(HLTMuonCand,1);
};
// 
// 
// class L2MuonCand {
// public:
// 
//   Float_t pt;           
//   Float_t eta;          
//   Float_t phi;          
//   Int_t   charge;      
// 
//   Int_t   minStations;      
//   Int_t   minHits    ;      
//   
//   L2MuonCand(){};
//   virtual ~L2MuonCand(){};
// 
//   ClassDef(L2MuonCand,1)
// 
// };
// 

class L1MuonCand {
public:

  Float_t pt;           
  Float_t eta;          
  Float_t phi;          
  Int_t   charge;      
  Int_t   quality;      
  
  L1MuonCand(){};
  virtual ~L1MuonCand(){};

  ClassDef(L1MuonCand,1)

};

class L1TauCand {
public:

  Float_t pt;           
  Float_t eta;          
  Float_t phi;          
  Int_t   charge;      
  Int_t   iso;      
  
  L1TauCand(){};
  virtual ~L1TauCand(){};

  ClassDef(L1TauCand,1)

};



class HLTObjCand {
public:

  std::string filterTag; // name of filter passed by the object
  Float_t pt;            // pt of the object passing the filter [GeV]
  Float_t eta;           // eta of the object passing the filter
  Float_t phi;           // phi of the object passing the filter
  
  HLTObjCand(){};
  virtual ~HLTObjCand(){};

  ClassDef(HLTObjCand,1)

};





class HLTInfo {
public:
  std::vector<std::string>  triggers;  
  std::vector<HLTObjCand>   objects;   
  double                    rho;     // ecal+hcal(m3)


  HLTInfo(){};
  virtual ~HLTInfo(){};
  bool match( const std::string & path ) {
	if (  std::find (  triggers.begin(), triggers.end(), path ) != triggers.end() )  return true;
//     if (! iname.compare("HLT_Mu20_v1") == 0) continue;
	return false;
  }

  bool find( const std::string & path ) {
	for ( std::vector<std::string>::const_iterator it = triggers.begin(); it != triggers.end(); ++it ) {
//       std::cout << *it << std::endl;
      if ( it-> compare(path) == 0) return true;
//       if ( it->find ( path ) != std::string::npos ) return true;
	}
	return false;
  }

  ClassDef(HLTInfo,1)

};


class TauEvent {
public:

  Int_t   runNumber;             
  Int_t   luminosityBlockNumber; 
  Int_t   eventNumber;           

  Int_t   nVtx;                    
  Float_t primaryVertex[3];        
//   Float_t cov_primaryVertex[3][3]; 

  Float_t rho; 

  Float_t bxId;
  Float_t instLumi; 
  Float_t trueNI;   

//   std::vector <GenParticleCand> genParticles; 
//   std::vector <L2MuonCand>      L2muons;      

  std::vector <L1TauCand>       L1taus;      
  std::vector <L1MuonCand>      L1muons;      

  std::vector <HLTTauCand>      hlttaus;      
  std::vector <HLTMuonCand>     hltmuons;      

  HLTInfo                       hlt;           
  HLTInfo                       hltTag;            

  TauEvent(){};
  virtual ~TauEvent(){};

  ClassDef(TauEvent,1)
};


#endif
