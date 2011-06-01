#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//@@ dwjang // #include "FWCore/Utilities/interface/Exception.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/Vector3D.h"
#include "Math/LorentzVector.h"
#include <vector>
#include <string>

/////////////////////////////////////////////////////////////////////////
JetCorrectionUncertainty::JetCorrectionUncertainty () 
{
  mJetEta = -9999;
  mJetPt  = -9999;
  mJetPhi = -9999;
  mJetE   = -9999;
  mJetEMF = -9999;
  mLepPx  = -9999;
  mLepPy  = -9999;
  mLepPz  = -9999;
  mIsJetEset   = false;
  mIsJetPtset  = false;
  mIsJetPhiset = false;
  mIsJetEtaset = false;
  mIsJetEMFset = false;
  mIsLepPxset  = false;
  mIsLepPyset  = false;
  mIsLepPzset  = false;
  mAddLepToJet = false;
  mUncertainty = new SimpleJetCorrectionUncertainty();
}
/////////////////////////////////////////////////////////////////////////
JetCorrectionUncertainty::JetCorrectionUncertainty(const std::string& fDataFile)  
{
  mJetEta = -9999;
  mJetPt  = -9999;
  mJetPhi = -9999;
  mJetE   = -9999;
  mJetEMF = -9999;
  mLepPx  = -9999;
  mLepPy  = -9999;
  mLepPz  = -9999;
  mIsJetEset   = false;
  mIsJetPtset  = false;
  mIsJetPhiset = false;
  mIsJetEtaset = false;
  mIsJetEMFset = false;
  mIsLepPxset  = false;
  mIsLepPyset  = false;
  mIsLepPzset  = false;
  mAddLepToJet = false;
  mUncertainty = new SimpleJetCorrectionUncertainty(fDataFile);
}
/////////////////////////////////////////////////////////////////////////
JetCorrectionUncertainty::JetCorrectionUncertainty(const JetCorrectorParameters& fParameters)  
{
  mJetEta = -9999;
  mJetPt  = -9999;
  mJetPhi = -9999;
  mJetE   = -9999;
  mJetEMF = -9999;
  mLepPx  = -9999;
  mLepPy  = -9999;
  mLepPz  = -9999;
  mIsJetEset   = false;
  mIsJetPtset  = false;
  mIsJetPhiset = false;
  mIsJetEtaset = false;
  mIsJetEMFset = false;
  mIsLepPxset  = false;
  mIsLepPyset  = false;
  mIsLepPzset  = false;
  mAddLepToJet = false;
  mUncertainty = new SimpleJetCorrectionUncertainty(fParameters);
}
/////////////////////////////////////////////////////////////////////////
JetCorrectionUncertainty::~JetCorrectionUncertainty () 
{
  delete mUncertainty;
}
/////////////////////////////////////////////////////////////////////////
void JetCorrectionUncertainty::setParameters(const std::string& fDataFile) 
{
  //---- delete the mParameters pointer before setting the new address ---
  delete mUncertainty; 
  mUncertainty = new SimpleJetCorrectionUncertainty(fDataFile);
}
/////////////////////////////////////////////////////////////////////////
float JetCorrectionUncertainty::getUncertainty(bool fDirection) 
{
  float result;
  std::vector<float> vx,vy;
  vx = fillVector(mUncertainty->parameters().definitions().binVar());
  vy = fillVector(mUncertainty->parameters().definitions().parVar());
  result = mUncertainty->uncertainty(vx,vy[0],fDirection);
  mIsJetEset   = false;
  mIsJetPtset  = false;
  mIsJetPhiset = false;
  mIsJetEtaset = false;
  mIsJetEMFset = false;
  mIsLepPxset  = false;
  mIsLepPyset  = false;
  mIsLepPzset  = false;
  return result;
}
//------------------------------------------------------------------------ 
//--- Reads the parameter names and fills a vector of floats -------------
//------------------------------------------------------------------------
std::vector<float> JetCorrectionUncertainty::fillVector(const std::vector<std::string>& fNames)
{
  std::vector<float> result;
  for(unsigned i=0;i<fNames.size();i++)
    {
      if (fNames[i] == "JetEta")
        {
	  //@@ dwjang // //          if (!mIsJetEtaset)
	    //@@ dwjang // //            throw cms::Exception("JetCorrectionUncertainty::")<<" jet eta is not set";
          result.push_back(mJetEta);
        }
      else if (fNames[i] == "JetPt")
        {
	  //@@ dwjang // //          if (!mIsJetPtset)
	    //@@ dwjang // //            throw cms::Exception("JetCorrectionUncertainty::")<<" jet pt is not set";  
          result.push_back(mJetPt);
        }
      else if (fNames[i] == "JetPhi")
        {
	  //@@ dwjang // //          if (!mIsJetPhiset)
	    //@@ dwjang // //            throw cms::Exception("JetCorrectionUncertainty::")<<" jet phi is not set";  
          result.push_back(mJetPt);
        }
      else if (fNames[i] == "JetE")
        {
	  //@@ dwjang // //          if (!mIsJetEset)
	    //@@ dwjang // //            throw cms::Exception("JetCorrectionUncertainty::")<<" jet energy is not set";
          result.push_back(mJetE);
        }
      else if (fNames[i] == "JetEMF")
        {
	  //@@ dwjang // //          if (!mIsJetEMFset)
	    //@@ dwjang // //            throw cms::Exception("JetCorrectionUncertainty::")<<" jet emf is not set";
          result.push_back(mJetEMF);
        } 
      else if (fNames[i] == "LepPx")
        {
	  //@@ dwjang // //          if (!mIsLepPxset)
	    //@@ dwjang // //            throw cms::Exception("JetCorrectionUncertainty::")<<" lepton px is not set";  
          result.push_back(mLepPx);
        }
      else if (fNames[i] == "LepPy")
        {
	  //@@ dwjang // //          if (!mIsLepPyset)
	    //@@ dwjang // //            throw cms::Exception("JetCorrectionUncertainty::")<<" lepton py is not set";  
          result.push_back(mLepPy);
        }
      else if (fNames[i] == "LepPz")
        {
	  //@@ dwjang // //          if (!mIsLepPzset)
	    //@@ dwjang // //            throw cms::Exception("JetCorrectionUncertainty::")<<" lepton pz is not set";  
          result.push_back(mLepPz);
        }
      else {
	//@@ dwjang // //        throw cms::Exception("JetCorrectionUncertainty::")<<" unknown parameter "<<fNames[i];
	;
      }
    }     
  return result;      
}
//------------------------------------------------------------------------ 
//--- Calculate the PtRel (needed for the SLB) ---------------------------
//------------------------------------------------------------------------
float JetCorrectionUncertainty::getPtRel()
{
  typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > PtEtaPhiELorentzVector;
  typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<float> > XYZVector;
  PtEtaPhiELorentzVector jet;
  XYZVector lep;
  jet.SetPt(mJetPt);
  jet.SetEta(mJetEta);
  jet.SetPhi(mJetPhi);
  jet.SetE(mJetE);
  lep.SetXYZ(mLepPx,mLepPy,mLepPz);
  float lj_x = (mAddLepToJet) ? lep.X()+jet.Px() : jet.Px();
  float lj_y = (mAddLepToJet) ? lep.Y()+jet.Py() : jet.Py();
  float lj_z = (mAddLepToJet) ? lep.Z()+jet.Pz() : jet.Pz();
  // absolute values squared
  float lj2  = lj_x*lj_x+lj_y*lj_y+lj_z*lj_z;
//@@ dwjang // //   if (!(lj2 > 0))
//@@ dwjang // //     throw cms::Exception("JetCorrectionUncertainty")<<" not positive lepton-jet momentum: "<<lj2;
  float lep2 = lep.X()*lep.X()+lep.Y()*lep.Y()+lep.Z()*lep.Z();
  // projection vec(mu) to lepjet axis
  float lepXlj = lep.X()*lj_x+lep.Y()*lj_y+lep.Z()*lj_z;
  // absolute value squared and normalized
  float pLrel2 = lepXlj*lepXlj/lj2;
  // lep2 = pTrel2 + pLrel2
  float pTrel2 = lep2-pLrel2;
  return (pTrel2 > 0) ? std::sqrt(pTrel2) : 0.0;
}
//------------------------------------------------------------------------ 
//--- Setters ------------------------------------------------------------
//------------------------------------------------------------------------
void JetCorrectionUncertainty::setJetEta(float fEta)
{
  mJetEta = fEta;
  mIsJetEtaset = true;
}
//------------------------------------------------------------------------
void JetCorrectionUncertainty::setJetPt(float fPt)
{
  mJetPt = fPt;
  mIsJetPtset  = true;
}
//------------------------------------------------------------------------
void JetCorrectionUncertainty::setJetPhi(float fPhi)
{
  mJetPhi = fPhi;
  mIsJetPhiset  = true;
}
//------------------------------------------------------------------------
void JetCorrectionUncertainty::setJetE(float fE)
{
  mJetE = fE;
  mIsJetEset   = true;
}
//------------------------------------------------------------------------
void JetCorrectionUncertainty::setJetEMF(float fEMF)
{
  mJetEMF = fEMF;
  mIsJetEMFset = true;
}
//------------------------------------------------------------------------
void JetCorrectionUncertainty::setLepPx(float fPx)
{
  mLepPx = fPx;
  mIsLepPxset  = true;
}
//------------------------------------------------------------------------
void JetCorrectionUncertainty::setLepPy(float fPy)
{
  mLepPy = fPy;
  mIsLepPyset  = true;
}
//------------------------------------------------------------------------
void JetCorrectionUncertainty::setLepPz(float fPz)
{
  mLepPz = fPz;
  mIsLepPzset  = true;
}
//------------------------------------------------------------------------
