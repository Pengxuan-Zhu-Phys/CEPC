// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Tools/Cutflow.hh"
#include "HepMC3/GenParticle.h"
#include <tuple>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <cmath>

class Vector3D {
public:
    double x, y, z;

    // Make vector 3 
    Vector3D(double x = 0.0, double y = 0.0, double z = 0.0) : x(x), y(y), z(z) {}
    
    double getX() const { return x; }
    double getY() const { return y; }
    double getZ() const { return z; }
    // add  
    Vector3D operator+(const Vector3D& v) const{
      return Vector3D(x + v.x, y + v.y, z + v.z); 
    }

    // minus
    Vector3D operator-(const Vector3D& v) const{
      return Vector3D(x - v.x, y - v.y, z - v.z);
    }

    // multiple 
    Vector3D operator*(const double scalar) const{
      return Vector3D(x * scalar, y * scalar, z * scalar);
    }
    // Friend function for scalar-vector multiplication
    friend Vector3D operator*(double scalar, const Vector3D& vec) {
        return vec * scalar;  // Reuse the vector-scalar multiplication
    }
    // calculate the length 
    double length() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    // dot 
    double dot(const Vector3D& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    // cross 
    Vector3D cross(const Vector3D& v) const {
      return {
        y * v.z - z * v.y, 
        z * v.x - x * v.z, 
        x * v.y - y * v.x
      };
    }

    // return unit vector 
    Vector3D normalize() const {
        double len = length();
        if (len != 0) {
            return *this * (1.0 / len);
        }
        return *this; // 返回自身，如果长度为0，保持不变
    }

};


namespace Rivet
{
  class DataFrame
  {
  private:
    std::vector<std::unordered_map<std::string, double>> rows;

  public:
    void addRow(const std::unordered_map<std::string, double> &row)
    {
      rows.push_back(row);
    }

    void toCSV(const std::string &filename)
    {
      std::ofstream file(filename);

      if (!file.is_open())
      {
        std::cerr << "Failed to open the file." << std::endl;
        return;
      }

      // Writing header
      auto header_iter = rows[0].begin();
      while (header_iter != rows[0].end())
      {
        file << header_iter->first;
        if (++header_iter != rows[0].end())
          file << ",";
      }
      file << "\n";

      // Writing data
      for (const auto &row : rows)
      {
        auto row_iter = row.begin();
        while (row_iter != row.end())
        {
          file << row_iter->second;
          if (++row_iter != row.end())
            file << ",";
        }
        file << "\n";
      }

      file.close();
      std::cout << "Data written to " << filename << std::endl;
    }
  };

  /// @brief Add a short analysis description here
  class CEPC_mRC : public Analysis
  {
  public:
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CEPC_mRC);

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init()
    {
      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      dfname = getOption("DFNAME");
      dfname.append(".csv");

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      declare(Beam(), "Beams");

      const FinalState fs(Cuts::abseta < 3.0);
      declare(VisibleFinalState(Cuts::abseta < 3.0), "vfs");

      Cut lepton_cuts = ((Cuts::abseta < 3.0) && (Cuts::E > 0.5 * GeV));

      FinalState bare_elec(Cuts::abspid == PID::ELECTRON && lepton_cuts);
      FinalState bare_muon(Cuts::abspid == PID::MUON && lepton_cuts);
      // PromptFinalState bare_elec(Cuts::abspid == PID::ELECTRON && lepton_cuts);
      // PromptFinalState bare_muon(Cuts::abspid == PID::MUON && lepton_cuts);

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");

      // FinalState of prompt photons and bare muons and electrons in the event
      // FinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);

      // Dress the prompt bare leptons with prompt photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons
      // DressedLeptons elec_drs(photons, bare_elec, 0.1, lepton_cuts);
      // DressedLeptons muon_drs(photons, bare_muon, 0.1, lepton_cuts);
      DressedLeptons elec_drs(photons, bare_elec, 0.1, lepton_cuts);
      DressedLeptons muon_drs(photons, bare_muon, 0.1, lepton_cuts);
      declare(muon_drs, "muons");
      declare(elec_drs, "elecs");
      // declare(bare_elec, "elecs");
      // declare(bare_muon, "muons");

      // Missing momentum
      declare(MissingMomentum(fs), "MET");

      // Book histograms
      // specify custom binning
      book(_hist_Wmin, "hist_Wmass_mRCmin", 50, 75, 85);
      book(_hist_Wmax, "hist_Wmass_mRCmax", 50, 75, 85);
      book(_Norm_hist_Wmin, "Nhist_Wmass_mRCmin", 50, 75, 85);
      book(_Norm_hist_Wmax, "Nhist_Wmass_mRCmax", 50, 75, 85);
      book(_hist_mRCmin, "hist_mRCmin", 240, 0., 125.);
      book(_hist_mRCmax, "hist_mRCmax", 240, 0., 125.);
      book(_Norm_hist_mRCmin, "Nhist_mRCmin", 250, 0., 125.);
      book(_Norm_hist_mRCmax, "Nhist_mRCmax", 250, 0., 125.);
    }

    /// Perform the per-event analysis
    void analyze(const Event &event)
    {

      FourMomentum pMiss;
      for (const Particle &p : apply<VisibleFinalState>(event, "vfs").particles())
      {
        pMiss -= p.momentum();
      }
      const Particles &muonFS = apply<FinalState>(event, "muons").particlesByPt();
      const int nmuon = muonFS.size();

      const Particles &elecFS = apply<FinalState>(event, "elecs").particlesByPt();
      const int nelec = elecFS.size();

      // const FastJets &jetsAntiKt4 = apply<FastJets>(event, "jets");
      // const Jets &jets = jetsAntiKt4.jetsByPt(10.0 * GeV);
      // const int njet = jets.size();

      const ParticlePair &beams = apply<Beam>(event, "Beams").beams();
      const double eC1 = beams.first.E();
      const double eC2 = beams.second.E();
      const Vector3 axis = (beams.first.charge() > 0) ? beams.first.momentum().p3().unit() : beams.second.momentum().p3().unit();

      FourMomentum P_Sum;
      P_Sum.setE(eC1 + eC2);
      P_Sum.setPx(0.);
      P_Sum.setPy(0.);
      P_Sum.setPz(0.);

      pMiss.setE(eC1 + eC2 + pMiss.E());
      // const double Emiss = pMiss.E();

      bool wmassPre = false;
      if (nmuon == 1 && nelec == 1 && pMiss.E() >= 1.0)
      {
        if (muonFS[0].charge() * elecFS[0].charge() < 0)
        {
          wmassPre = true;
        }
      }
      // Pass Preselection
      if (wmassPre)
      {
        FourMomentum P_ISR = P_Sum - pMiss - muonFS[0].momentum() - elecFS[0].momentum();
        FourMomentum P_recoil = P_Sum - muonFS[0].momentum() - elecFS[0].momentum();
        const double mRecoil = P_recoil.mass();
        const double Emiss = pMiss.E();

        const double mVV = (muonFS[0].momentum() + elecFS[0].momentum()).mass();
        vector<FourMomentum> mrc = solve_semi_pIa(pMiss, muonFS[0].mom(), elecFS[0].mom(), P_ISR, 0.);
        // const double mRCmin = mrc[0];
        // const double mRCmax = mrc[1];
        // const double mLSPmax = mrc[2];
        // const double mRCLSP = mrc[3];

        FourMomentum P_Wa;
        FourMomentum P_Wb;

        const HepMC3::GenEvent *genEvent = event.genEvent();
        MSG_INFO("New Events found particles ");
        for (HepMC3::ConstGenParticlePtr particle : genEvent->particles())
        {
          // MSG_INFO("Found W particle tree  -> " << particle->pid() << "\t" << particle->momentum());
          if (particle->pid() == 24 || particle->pid() == -24)
          {
            if (hasChild(particle, 13))
            {
              P_Wa = particle->momentum();
              MSG_INFO("P_Wa is " << P_Wa);
            }
            if (hasChild(particle, 11))
            {
              P_Wb = particle->momentum();
              MSG_INFO("P_Wb is " << P_Wb);
            }
          }
        }
        MSG_INFO("New Events end loop particles ");

        const double mWa = P_Wa.mass();
        // MSG_INFO("mWa is -> " << mWa);
        const double mWb = P_Wb.mass();
        _hist_Wmin->fill(mRCmin);
        _hist_Wmax->fill(mRCmax);
        _Norm_hist_Wmin->fill(mRCmin);
        _Norm_hist_Wmax->fill(mRCmax);
        _hist_mRCmin->fill(mRCmin);
        _hist_mRCmax->fill(mRCmax);
        _Norm_hist_mRCmin->fill(mRCmin);
        _Norm_hist_mRCmax->fill(mRCmax);

        df.addRow({{"mRCmin", mRCmin},
                   {"mRCmax", mRCmax},
                   {"mRCLSP", mRCLSP},
                   {"mLSPmax", mLSPmax},
                   {"mRecoil", mRecoil},
                   {"pVa", muonFS[0].p3().mod()},
                   {"pVb", elecFS[0].p3().mod()},
                   {"mVV", mVV},
                   {"EMiss", Emiss},
                   {"pMiss", pMiss.p3().mod()},
                   {"mWa", mWa},
                   {"mWb", mWb}});
      }
      else
        return;
    }

    /// Normalise histograms etc., after the run
    void finalize()
    {
      const double sf = 409.2825;
      const double Lint = 5.05e3;
      double norm = sf * Lint;
      MSG_INFO("data file is " << dfname);
      df.toCSV(dfname);

      MSG_INFO("Norm is " << norm);
      MSG_INFO("Total Cross section is " << crossSection() / femtobarn << " femtobarn!");
      scale(_hist_mRCmin, sf * Lint / sumW()); // norm to generated cross-section in pb (after cuts)
      scale(_hist_mRCmax, sf * Lint / sumW()); // norm to generated cross-section in pb (after cuts)
      scale(_hist_Wmin, sf * Lint / sumW());   // norm to generated cross-section in pb (after cuts)
      scale(_hist_Wmax, sf * Lint / sumW());   // norm to generated cross-section in pb (after cuts)

      normalize(_Norm_hist_Wmin);
      normalize(_Norm_hist_Wmax);
      normalize(_Norm_hist_mRCmin);
      normalize(_Norm_hist_mRCmax);
    }

    vector<FourMomentum> solve_semi_pIa(const FourMomentum &met, const FourMomentum &Va, const FourMomentum &Vb, const FourMomentum &ISR, const double &mI = 0.0) const
    {
      // MSG_INFO("Tag 1");
      FourMomentum CM = met + Va + Vb;

      FourMomentum bmet;
      FourMomentum bVa;
      FourMomentum bVb;
      FourMomentum pISR;

      LorentzTransform LT = LorentzTransform::mkFrameTransformFromBeta(CM.betaVec());
      LorentzTransform ivsLT = LT.inverse();


      pISR = LT.transform(ISR);
      bmet = LT.transform(met);
      bVa = LT.transform(Va);
      bVb = LT.transform(Vb);

      const double ss = (bmet.E() + bVa.E() + bVb.E()) / 2.;
      const double EVa = bVa.E();
      const double EVb = bVb.E();
      const double EIa = ss - EVa;
      const double EIb = ss - EVb;
      
      const Vector3D p3I(bmet.px(), bmet.py(), bmet.pz());
      const double lI = p3I.length();
      const Vector3D nI = p3I.normalize();

      const Vector3D p3Va(bVa.px(), bVa.py(), bVa.pz());
      const double lVa = p3Va.length();
      const Vector3D nVa = p3Va.normalize();

      const Vector3D p3Vb(bVb.px(), bVb.py(), bVb.pz());
      const double lVb = p3Vb.length();
      const Vector3D nVb = p3Vb.normalize(); 

      // const double pI1max = sqrt(EVa * EVa - mI * mI);
      // const double pI2max = sqrt(EVb * EVb - mI * mI);

      const double costhetaIVa = nI.dot(nVa);
      const double costhetaIVb = nI.dot(nVb);
      const double costhetaVab = nVa.dot(nVb);

      // Define the basis vector 
      Vector3D h3V = nI * lVa * costhetaIVa - p3Va; 
      Vector3D nhV = h3V.normalize(); 

      Vector3D r3V = nI.cross(nhV); 
      Vector3D rV  = r3V.normalize(); 
      // basis vector : {nI, nhV, rV}

      // Solve the visible triangle 
      const vector<double> pos_V = solveXY(lVa, lVb, lI); 
      const vector<double> pos_I = solveXY(EIa, EIb, lI); 

      const Vector3D p3O = pos_I[0] * nI; 
      const Vector3D p3A = pos_I[0] * nI - pos_I[1] * nhV; 
      const Vector3D p3B = pos_I[0] * nI + pos_I[1] * nhV; 
      const Vector3D p3M(0., 0., 0.); 
      const Vector3D p3N = lI * nI; 
      const Vector3D p3P = pos_V[0] * nI + pos_V[1] * nhV; 

      const Vector3D p3Pa_A = p3P - p3A; 
      const Vector3D p3Pa_B = p3P - p3B; 
      const Vector3D p3Pa_O = p3P - p3O; 

      const FourMomentum pPa_O(ss, p3Pa_O.getX(), p3Pa_O.getY(), p3Pa_O.getZ()); 
      const FourMomentum pPa_A(ss, p3Pa_A.getX(), p3Pa_A.getY(), p3Pa_A.getZ());
      const FourMomentum pPa_B(ss, p3Pa_B.getX(), p3Pa_B.getY(), p3Pa_B.getZ()); 

      const FourMomentum pIa_O(EIa, p3O.getX(), p3O.getY(), p3O.getZ());
      const FourMomentum pIa_A(EIa, p3A.getX(), p3A.getY(), p3A.getZ());
      const FourMomentum pIa_B(EIa, p3B.getX(), p3B.getY(), p3B.getZ());
      
      const FourMomentum pIaO = ivsLT.transform(pIa_O);
      const FourMomentum pIaA = ivsLT.transform(pIa_A);
      const FourMomentum pIaB = ivsLT.transform(pIa_A);
      const double mImax = sqrt(EIa * EIa - p3O.dot(p3O)); 
      const double mPmax = sqrt(EVa * EVa - p3Pa_A.dot(p3Pa_A)); 
      const double mPmin = sqrt(EVa * EVa - p3Pa_B.dot(p3Pa_B)); 
      const double mPLSP = sqrt(EVa * EVa - p3Pa_O.dot(p3Pa_O)); 

      const vector<FourMomentum> pIa = {pIaO, pIaA, pIaB};
      return pIa;

      


      // if (pI1max > 0. && pI2max > 0. && pI1max + pI2max > pMiss && pL1 + pL2 > pMiss && abs(pI1max - pI2max) < pMiss && abs(pL1 - pL2) < pMiss)
      // {
        

      //   const vector<double> pos_C = solveXY(pI1max, pI2max, pMiss);
      //   const vector<double> pos_B = solveXY(pL1, pL2, pMiss);

      //   const double pMax2 = pow(pos_B[0] - pos_C[0], 2) + pow(pos_B[1] + pos_C[1], 2);
      //   const double pMinX2 = pow(pos_B[0] - pos_C[0], 2);
      //   const double pMinY2 = pos_B[1] > pos_C[1] ? pow(pos_B[1] - pos_C[1], 2) : 0.0;
      //   const double pLSP2 = pow(pos_B[0] - pos_C[0], 2) + pow(pos_B[1], 2);

      //   const double mYmax = sqrt(ss * ss - pMinX2 - pMinY2);
      //   const double mYmin = sqrt(ss * ss - pMax2);
      //   const double mImax = sqrt(EI1 * EI1 - pos_C[0] * pos_C[0]);
      //   const double mYLSP = sqrt(ss * ss - pLSP2);

      //   // return mYmin;
      //   const vector<double> mrc = {mYmin, mYmax, mImax, mYLSP};
      //   return mrc;
      // }
      // else
      // {
      //   MSG_INFO("Erroring in reconstructing mass");
      //   const vector<double> mYmax = {-1.0, -1.0};
      //   return mYmax;
      // }
    }

    vector<double> solveXY(const double &p1, const double &p2, const double &pMiss) const
    {
      const double x = 0.5 / pMiss * (pow(p1, 2) - pow(p2, 2) + pow(pMiss, 2));
      const double numinator = 2.0 * (pow(p1, 2) * pow(p2, 2) + pow(p1, 2) * pow(pMiss, 2) + pow(p2, 2) * pow(pMiss, 2)) - pow(p1, 4) - pow(p2, 4) - pow(pMiss, 4);
      const double y = 0.5 / pMiss * sqrt(numinator);
      const vector<double> pos = {x, y};
      return pos;
    }

    bool hasChild(HepMC3::ConstGenParticlePtr particle, const int &cid)
    {
      if (particle->end_vertex())
      {
        for (auto child : particle->end_vertex()->particles_out())
        {
          if (child->pdg_id() == cid || child->pdg_id() == -cid)
          {
            return true; // Found an electron or positron among the children
          }
        }
      }
      return false; // No electron or positron found among the children
    }

    /// @}

    /// @name Histograms
    /// @{
    std::string dfname;
    DataFrame df;

    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    /// @}
    Histo1DPtr _hist_Wmin;
    Histo1DPtr _hist_Wmax;
    Histo1DPtr _Norm_hist_Wmin;
    Histo1DPtr _Norm_hist_Wmax;
    Histo1DPtr _hist_mRCmin;
    Histo1DPtr _hist_mRCmax;
    Histo1DPtr _Norm_hist_mRCmin;
    Histo1DPtr _Norm_hist_mRCmax;
    /// @}
  };

  RIVET_DECLARE_PLUGIN(CEPC_mRC);

}
