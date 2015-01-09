// Recipe taken from the RazorBoost gurus N. Strobbe, S. Sekmen
//
//  https://github.com/nstrobbe/RazorBoost/blob/master/analyzer/utils.h
//

#include <vector>
#include "TLorentzVector.h"

#include "../interface/Data.h"

using namespace std;

// Hemispheres:
vector<TLorentzVector> CombineJets(vector<TLorentzVector> myjets)
{
  vector<TLorentzVector> mynewjets;
  TLorentzVector j1, j2;
  //bool foundGood = false;
  int N_comb = 1;
  for(unsigned int i = 0; i < myjets.size(); i++){
    N_comb *= 2;
  }
  double M_min = 9999999999.0;
  int j_count;
  for(int i = 1; i < N_comb-1; i++){
    TLorentzVector j_temp1, j_temp2;
    int itemp = i;
    j_count = N_comb/2;
    int count = 0;
    while(j_count > 0){
      if(itemp/j_count == 1){
	j_temp1 += myjets[count];
      } else {
	j_temp2 += myjets[count];
      }
      itemp -= j_count*(itemp/j_count);
      j_count /= 2;
      count++;
    }
    double M_temp = j_temp1.M2()+j_temp2.M2();
    // smallest mass
    if(M_temp < M_min){
      M_min = M_temp;
      j1 = j_temp1;
      j2 = j_temp2;
    }
  }
  if(j2.Pt() > j1.Pt()){
    TLorentzVector temp = j1;
    j1 = j2;
    j2 = temp;
  }
  mynewjets.push_back(j1);
  mynewjets.push_back(j2);
  return mynewjets;
}

// MR
double CalcMR(TLorentzVector ja, TLorentzVector jb){
  double A = ja.P();
  double B = jb.P();
  double az = ja.Pz();
  double bz = jb.Pz();
  TVector3 jaT, jbT;
  jaT.SetXYZ(ja.Px(),ja.Py(),0.0);
  jbT.SetXYZ(jb.Px(),jb.Py(),0.0);
  double ATBT = (jaT+jbT).Mag2();
  double temp = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
		     (jbT.Dot(jbT)-jaT.Dot(jaT))*(jbT.Dot(jbT)-jaT.Dot(jaT))/(jaT+jbT).Mag2());
  double mybeta = (jbT.Dot(jbT)-jaT.Dot(jaT))/sqrt(ATBT*((A+B)*(A+B)-(az+bz)*(az+bz)));
  double mygamma = 1./sqrt(1.-mybeta*mybeta);
  //gamma times MRstar
  temp *= mygamma;
  return temp;
}

// MTR
double CalcMTR(TLorentzVector ja, TLorentzVector jb, TVector3 met){
  double temp = met.Mag()*(ja.Pt()+jb.Pt()) - met.Dot(ja.Vect()+jb.Vect());
  temp /= 2.;
  temp = sqrt(temp);
  return temp;
}

// MT
double CalcMT(TLorentzVector lepton, TLorentzVector pfmet){
  return sqrt( 2 * lepton.Pt() * pfmet.Pt() * ( 1 - cos( pfmet.Phi() - lepton.Phi() ) ) );
}

// Calculate Razor variables here
void calcRazorAK4(Data& d) {
  
  // Select the best pair of jets (AK4, pt>40, |eta| < 3.0)
  std::vector<TLorentzVector> sjetl;
  while(d.jetAK4.Loop()) {
    if (!(d.jetAK4.Pt > 40) ) continue;
    if (!(fabs(d.jetAK4.Eta) < 3) ) continue;
    TLorentzVector jl;
    jl.SetPtEtaPhiE(d.jetAK4.Pt, d.jetAK4.Eta,
		    d.jetAK4.Phi, d.jetAK4.E);
    sjetl.push_back(jl);
  }
  std::vector<TLorentzVector> hemis = CombineJets(sjetl);
  
  // ---------------------
  // -- Razor variables --
  // ---------------------
  
  TVector3 metl;
  metl.SetPtEtaPhi(d.met.metPt[0], 0, d.met.metPhi[0]);
  
  if (hemis.size() == 2) {
    d.jetAK4.MR = CalcMR(hemis[0], hemis[1]);
    d.jetAK4.MTR = CalcMTR(hemis[0], hemis[1], metl);
    d.jetAK4.R = d.jetAK4.MTR / d.jetAK4.MR;
    d.jetAK4.R2 = pow(d.jetAK4.R, 2);
  }
}

void calcRazorAK8(Data& d) {
  
  // Select the best pair of jets (AK8, pt>40, |eta| < 3.0)
  std::vector<TLorentzVector> sjetl;
  while(d.jetAK8.Loop()) {
    if (!(d.jetAK8.Pt > 40) ) continue;
    if (!(fabs(d.jetAK8.Eta) < 3) ) continue;
    TLorentzVector jl;
    jl.SetPtEtaPhiE(d.jetAK8.Pt, d.jetAK8.Eta,
		    d.jetAK8.Phi, d.jetAK8.E);
    sjetl.push_back(jl);
  }
  std::vector<TLorentzVector> hemis = CombineJets(sjetl);
  
  // ---------------------
  // -- Razor variables --
  // ---------------------
  
  TVector3 metl;
  metl.SetPtEtaPhi(d.met.metPt[0], 0, d.met.metPhi[0]);
  
  if (hemis.size() == 2) {
    d.jetAK8.MR = CalcMR(hemis[0], hemis[1]);
    d.jetAK8.MTR = CalcMTR(hemis[0], hemis[1], metl);
    d.jetAK8.R = d.jetAK8.MTR / d.jetAK8.MR;
    d.jetAK8.R2 = pow(d.jetAK8.R, 2);
  }
}
