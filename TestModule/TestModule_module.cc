////////////////////////////////////////////////////////////////////////
// Class:       TestModule
// Plugin Type: analyzer (art v2_05_01)
// File:        TestModule_module.cc
//
// Generated at Mon Mar 12 16:28:34 2018 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/GTruth.h"

#include "TTree.h"

class TestModule;


class TestModule : public art::EDAnalyzer {
  public:
    explicit TestModule(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    TestModule(TestModule const &) = delete;
    TestModule(TestModule &&) = delete;
    TestModule & operator = (TestModule const &) = delete;
    TestModule & operator = (TestModule &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    // selected optional functions
    void beginJob() override;

    static const int kMaxMCParticles = 100;

    // Metadata
    int run;
    int subrun;
    int event;

    // MCFlux
    int MCFlux_evtno;
    double MCFlux_NuPosX, MCFlux_NuPosY, MCFlux_NuPosZ;
    double MCFlux_NuMomX, MCFlux_NuMomY, MCFlux_NuMomZ, MCFlux_NuMomE;
    double MCFlux_genx, MCFlux_geny, MCFlux_genz;
    int MCFlux_ntype;
    int MCFlux_ptype;
    double MCFlux_nimpwt;
    double MCFlux_dk2gen;
    double MCFlux_nenergyn;
    double MCFlux_tpx, MCFlux_tpy, MCFlux_tpz;
    double MCFlux_vx, MCFlux_vy, MCFlux_vz;
    int MCFlux_tptype;

    // MCTruth
    int MCTruth_NParticles;
    int MCTruth_particles_TrackId[kMaxMCParticles];
    int MCTruth_particles_PdgCode[kMaxMCParticles];
    int MCTruth_particles_Mother[kMaxMCParticles];
    int MCTruth_particles_StatusCode[kMaxMCParticles];
    int MCTruth_particles_NumberDaughters[kMaxMCParticles];
    int MCTruth_particles_Daughters[100][100]; // only on multi-dimensional array allowed... thanks root.
    double MCTruth_particles_Gvx[kMaxMCParticles];
    double MCTruth_particles_Gvy[kMaxMCParticles];
    double MCTruth_particles_Gvz[kMaxMCParticles];
    double MCTruth_particles_Gvt[kMaxMCParticles];
    double MCTruth_particles_px0[kMaxMCParticles];
    double MCTruth_particles_py0[kMaxMCParticles];
    double MCTruth_particles_pz0[kMaxMCParticles];
    double MCTruth_particles_e0[kMaxMCParticles];
    int MCTruth_particles_Rescatter[kMaxMCParticles];
    double MCTruth_particles_polx[kMaxMCParticles];
    double MCTruth_particles_poly[kMaxMCParticles];
    double MCTruth_particles_polz[kMaxMCParticles];
    int MCTruth_neutrino_CCNC;
    int MCTruth_neutrino_mode;
    int MCTruth_neutrino_interactionType;
    int MCTruth_neutrino_target;
    int MCTruth_neutrino_nucleon;
    int MCTruth_neutrino_quark;
    double MCTruth_neutrino_QSqr;
    double MCTruth_neutrino_W;
    double MCTruth_neutrino_X;
    double MCTruth_neutrino_Y;

    // GTruth
    bool GTruth_IsSeaQuark;
    int GTruth_tgtPDG;
    double GTruth_weight;
    double GTruth_probability;
    double GTruth_Xsec;
    double GTruth_DiffXsec;
    double GTruth_vertexX;
    double GTruth_vertexY;
    double GTruth_vertexZ;
    double GTruth_vertexT;
    int GTruth_Gscatter;
    int GTruth_Gint;
    int GTruth_ResNum;
    int GTruth_NumPiPlus;
    int GTruth_NumPi0;
    int GTruth_NumPiMinus;
    int GTruth_NumProton;
    int GTruth_NumNeutron;
    bool GTruth_IsCharm;
    double GTruth_gX;
    double GTruth_gY;
    //double GTruth_gZ;
    double GTruth_gT;
    double GTruth_gW;
    double GTruth_gQ2;
    double GTruth_gq2;
    int GTruth_ProbePDG;
    double GTruth_ProbeP4x;
    double GTruth_ProbeP4y;
    double GTruth_ProbeP4z;
    double GTruth_ProbeP4E;
    double GTruth_HitNucP4x;
    double GTruth_HitNucP4y;
    double GTruth_HitNucP4z;
    double GTruth_HitNucP4E;
    double GTruth_FShadSystP4x;
    double GTruth_FShadSystP4y;
    double GTruth_FShadSystP4z;
    double GTruth_FShadSystP4E;

  private:
    TTree* tree;
    art::ServiceHandle< art::TFileService > tfs;

    bool isDebug;


};


TestModule::TestModule(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
    // More initializers here.
{

  isDebug = p.get<bool> ("IsDebug");

}

void TestModule::beginJob()
{

  tree = tfs->make<TTree>("tree", "tree");
  // Metadata
  tree->Branch("run"     , &run);
  tree->Branch("subrun"  , &subrun);
  tree->Branch("event"   , &event);

  // MCFlux
  tree->Branch("MCFlux_evtno"   , &MCFlux_evtno);
  tree->Branch("MCFlux_NuPosX"  , &MCFlux_NuPosX);
  tree->Branch("MCFlux_NuPosY"  , &MCFlux_NuPosY);
  tree->Branch("MCFlux_NuPosZ"  , &MCFlux_NuPosZ);
  tree->Branch("MCFlux_NuMomX"  , &MCFlux_NuMomX);
  tree->Branch("MCFlux_NuMomY"  , &MCFlux_NuMomY);
  tree->Branch("MCFlux_NuMomZ"  , &MCFlux_NuMomZ);
  tree->Branch("MCFlux_NuMomE"  , &MCFlux_NuMomE);
  tree->Branch("MCFlux_genx"    , &MCFlux_genx);
  tree->Branch("MCFlux_geny"    , &MCFlux_geny);
  tree->Branch("MCFlux_genz"    , &MCFlux_genz);
  tree->Branch("MCFlux_ntype"   , &MCFlux_ntype);
  tree->Branch("MCFlux_ptype"   , &MCFlux_ptype);
  tree->Branch("MCFlux_nimpwt"  , &MCFlux_nimpwt);
  tree->Branch("MCFlux_dk2gen"  , &MCFlux_dk2gen);
  tree->Branch("MCFlux_nenergyn", &MCFlux_nenergyn);
  tree->Branch("MCFlux_tpx"     , &MCFlux_tpx);
  tree->Branch("MCFlux_tpy"     , &MCFlux_tpy);
  tree->Branch("MCFlux_tpz"     , &MCFlux_tpz);
  tree->Branch("MCFlux_tptype"  , &MCFlux_tptype);
  tree->Branch("MCFlux_vx"      , &MCFlux_vx);
  tree->Branch("MCFlux_vy"      , &MCFlux_vy);
  tree->Branch("MCFlux_vz"      , &MCFlux_vz);

  // MCTruth
  tree->Branch("MCTruth_NParticles", &MCTruth_NParticles);
  tree->Branch("MCTruth_particles_TrackId", MCTruth_particles_TrackId, "MCTruth_particles_TrackId[MCTruth_NParticles]/I");
  tree->Branch("MCTruth_particles_PdgCode", &MCTruth_particles_PdgCode, "MCTruth_particles_PdgCode[MCTruth_NParticles]/I");
  tree->Branch("MCTruth_particles_Mother", &MCTruth_particles_Mother, "MCTruth_particles_Mother[MCTruth_NParticles]/I"); 
  tree->Branch("MCTruth_particles_StatusCode", MCTruth_particles_StatusCode, "MCTruth_particles_StatusCode[MCTruth_NParticles]/I"); 
  tree->Branch("MCTruth_particles_NumberDaughters", MCTruth_particles_NumberDaughters, "MCTruth_particles_NumberDaughters[MCTruth_NParticles]/I");
  tree->Branch("MCTruth_particles_Daughters", MCTruth_particles_Daughters, "MCTruth_particles_Daughters[MCTruth_NParticles][100]/I");
  tree->Branch("MCTruth_particles_Gvx", MCTruth_particles_Gvx, "MCTruth_particles_Gvx[MCTruth_NParticles]/D");
  tree->Branch("MCTruth_particles_Gvy", MCTruth_particles_Gvy, "MCTruth_particles_Gvy[MCTruth_NParticles]/D");
  tree->Branch("MCTruth_particles_Gvz", MCTruth_particles_Gvz, "MCTruth_particles_Gvz[MCTruth_NParticles]/D");
  tree->Branch("MCTruth_particles_Gvt", MCTruth_particles_Gvt, "MCTruth_particles_Gvt[MCTruth_NParticles]/D");
  tree->Branch("MCTruth_particles_px0", MCTruth_particles_px0, "MCTruth_particles_px0[MCTruth_NParticles]/D");
  tree->Branch("MCTruth_particles_py0", MCTruth_particles_py0, "MCTruth_particles_py0[MCTruth_NParticles]/D");
  tree->Branch("MCTruth_particles_pz0", MCTruth_particles_pz0, "MCTruth_particles_pz0[MCTruth_NParticles]/D");
  tree->Branch("MCTruth_particles_e0", MCTruth_particles_e0, "MCTruth_particles_e0[MCTruth_NParticles]/D");
  tree->Branch("MCTruth_particles_Rescatter", MCTruth_particles_Rescatter, "MCTruth_particles_Rescatter[MCTruth_NParticles]/I");
  tree->Branch("MCTruth_particles_polx", MCTruth_particles_polx, "MCTruth_particles_polx[MCTruth_NParticles]/D");
  tree->Branch("MCTruth_particles_poly", MCTruth_particles_poly, "MCTruth_particles_poly[MCTruth_NParticles]/D");
  tree->Branch("MCTruth_particles_polz", MCTruth_particles_polz, "MCTruth_particles_polz[MCTruth_NParticles]/D");
  tree->Branch("MCTruth_neutrino_CCNC", &MCTruth_neutrino_CCNC);
  tree->Branch("MCTruth_neutrino_mode", &MCTruth_neutrino_mode);
  tree->Branch("MCTruth_neutrino_interactionType", &MCTruth_neutrino_interactionType);
  tree->Branch("MCTruth_neutrino_target", &MCTruth_neutrino_target);
  tree->Branch("MCTruth_neutrino_nucleon", &MCTruth_neutrino_nucleon);
  tree->Branch("MCTruth_neutrino_quark", &MCTruth_neutrino_quark);
  tree->Branch("MCTruth_neutrino_W", &MCTruth_neutrino_W);
  tree->Branch("MCTruth_neutrino_X", &MCTruth_neutrino_X);
  tree->Branch("MCTruth_neutrino_Y", &MCTruth_neutrino_Y);
  tree->Branch("MCTruth_neutrino_QSqr",&MCTruth_neutrino_QSqr);

  //GTruth
  tree->Branch("GTruth_IsSeaQuark"   , &GTruth_IsSeaQuark);
  tree->Branch("GTruth_tgtPDG"       , &GTruth_tgtPDG);
  tree->Branch("GTruth_weight"       , &GTruth_weight);
  tree->Branch("GTruth_probability"  , &GTruth_probability);
  tree->Branch("GTruth_Xsec"         , &GTruth_Xsec);
  tree->Branch("GTruth_DiffXsec"     , &GTruth_DiffXsec);
  tree->Branch("GTruth_vertexX"      , &GTruth_vertexX);
  tree->Branch("GTruth_vertexY"      , &GTruth_vertexY);
  tree->Branch("GTruth_vertexZ"      , &GTruth_vertexZ);
  tree->Branch("GTruth_vertexT"      , &GTruth_vertexT);
  tree->Branch("GTruth_Gscatter"     , &GTruth_Gscatter);
  tree->Branch("GTruth_Gint"         , &GTruth_Gint);
  tree->Branch("GTruth_ResNum"       , &GTruth_ResNum);
  tree->Branch("GTruth_NumPiPlus"    , &GTruth_NumPiPlus);
  tree->Branch("GTruth_NumPi0"       , &GTruth_NumPiMinus);
  tree->Branch("GTruth_NumPiMinus"   , &GTruth_NumPiMinus);
  tree->Branch("GTruth_NumProton"    , &GTruth_NumProton);
  tree->Branch("GTruth_NumNeutron"   , &GTruth_NumNeutron);
  tree->Branch("GTruth_IsCharm"      , &GTruth_IsCharm);
  tree->Branch("GTruth_gX"           , &GTruth_gX);
  tree->Branch("GTruth_gY"           , &GTruth_gY);
  tree->Branch("GTruth_gT"           , &GTruth_gT);
  tree->Branch("GTruth_gW"           , &GTruth_gW);
  tree->Branch("GTruth_gQ2"          , &GTruth_gQ2);
  tree->Branch("GTruth_gq2"          , &GTruth_gq2);
  tree->Branch("GTruth_ProbePDG"     , &GTruth_ProbePDG);
  tree->Branch("GTruth_ProbeP4x"     , &GTruth_ProbeP4x);
  tree->Branch("GTruth_ProbeP4y"     , &GTruth_ProbeP4y);
  tree->Branch("GTruth_ProbeP4z"     , &GTruth_ProbeP4z);
  tree->Branch("GTruth_ProbeP4E"     , &GTruth_ProbeP4E);
  tree->Branch("GTruth_HitNucP4x"    , &GTruth_HitNucP4x);
  tree->Branch("GTruth_HitNucP4y"    , &GTruth_HitNucP4y);
  tree->Branch("GTruth_HitNucP4z"    , &GTruth_HitNucP4z);
  tree->Branch("GTruth_HitNucP4E"    , &GTruth_HitNucP4E);
  tree->Branch("GTruth_FShadSystP4x" , &GTruth_FShadSystP4x);
  tree->Branch("GTruth_FShadSystP4y" , &GTruth_FShadSystP4y);
  tree->Branch("GTruth_FShadSystP4z" , &GTruth_FShadSystP4z);
  tree->Branch("GTruth_FShadSystP4E" , &GTruth_FShadSystP4E);



}

void TestModule::analyze(art::Event const & e)
{

  art::Handle< std::vector<simb::MCFlux> > mcFluxHandle;
  e.getByLabel("generator", mcFluxHandle);
  if (!mcFluxHandle.isValid()) return;
  std::vector< art::Ptr<simb::MCFlux> > mcFluxVec;
  art::fill_ptr_vector(mcFluxVec, mcFluxHandle);
  if (mcFluxVec.size() == 0){
    std::cout << ">> No MCFlux information" << std::endl;
    return;
  }

  art::Handle< std::vector<simb::MCTruth> > mcTruthHandle;
  e.getByLabel("generator", mcTruthHandle);
  if (!mcTruthHandle.isValid()) return;
  std::vector< art::Ptr<simb::MCTruth> > mcTruthVec;
  art::fill_ptr_vector(mcTruthVec, mcTruthHandle);
  if (mcTruthVec.size() == 0){
    std::cout << ">> No MCTruth information" << std::endl;
    return;
  }

  art::Handle< std::vector< simb::GTruth > > gTruthHandle;
  e.getByLabel("generator", gTruthHandle);
  if (!gTruthHandle.isValid()) return;
  std::vector< art::Ptr<simb::GTruth> > gTruthVec;
  art::fill_ptr_vector(gTruthVec, gTruthHandle);
  if (gTruthVec.size() == 0){
    std::cout << ">> No GTruth information" << std::endl;
    return;
  }

  const art::Ptr<simb::MCFlux> mcFlux = mcFluxVec.at(0);
  const art::Ptr<simb::MCTruth> mcTruth = mcTruthVec.at(0);
  const simb::MCParticle& nu = mcTruth->GetNeutrino().Nu();
  const art::Ptr<simb::GTruth> gTruth = gTruthVec.at(0);

  run = e.run();
  subrun = e.subRun();
  event = e.event();

  // possibly the wrong variables, but let's see for now...
  MCFlux_evtno     = mcFlux->fevtno;
  MCFlux_NuPosX    = nu.Vx();
  MCFlux_NuPosY    = nu.Vy();
  MCFlux_NuPosZ    = nu.Vz();
  MCFlux_NuMomX    = nu.Px(); 
  MCFlux_NuMomY    = nu.Py(); 
  MCFlux_NuMomZ    = nu.Pz(); 
  MCFlux_NuMomE    = nu.E();
  MCFlux_genx      = mcFlux->fgenx;
  MCFlux_geny      = mcFlux->fgeny;
  MCFlux_genz      = mcFlux->fgenz;
  MCFlux_ntype     = mcFlux->fntype;
  MCFlux_ptype     = mcFlux->fptype;
  MCFlux_nimpwt    = mcFlux->fnimpwt;
  MCFlux_dk2gen    = mcFlux->fdk2gen;
  MCFlux_nenergyn  = mcFlux->fnenergyn;
  MCFlux_tpx       = mcFlux->ftpx;
  MCFlux_tpy       = mcFlux->ftpy;
  MCFlux_tpz       = mcFlux->ftpz;
  MCFlux_tptype    = mcFlux->ftptype;
  MCFlux_vx        = mcFlux->fvx;
  MCFlux_vy        = mcFlux->fvy;
  MCFlux_vz        = mcFlux->fvz;

  // loop MCParticle info for MCTruth object

  MCTruth_NParticles = mcTruth->NParticles();

  for (int i = 0; i < MCTruth_NParticles; i++){

    const simb::MCParticle& mcParticle = mcTruth->GetParticle(i);

    MCTruth_particles_TrackId[i] = mcParticle.TrackId();
    MCTruth_particles_PdgCode[i] = mcParticle.PdgCode();
    MCTruth_particles_Mother[i]  = mcParticle.Mother();
    MCTruth_particles_StatusCode[i] = mcParticle.StatusCode();
    MCTruth_particles_NumberDaughters[i] = mcParticle.NumberDaughters();

    for (int j = 0; j < MCTruth_particles_NumberDaughters[i]; j++){

      const simb::MCParticle& daughterMcParticle = mcTruth->GetParticle(j);
      MCTruth_particles_Daughters[i][j] = daughterMcParticle.TrackId();

    }

    MCTruth_particles_Gvx[i] = mcParticle.Gvx();
    MCTruth_particles_Gvy[i] = mcParticle.Gvy();
    MCTruth_particles_Gvz[i] = mcParticle.Gvz();
    MCTruth_particles_Gvt[i] = mcParticle.Gvt();
    MCTruth_particles_px0[i] = mcParticle.Px(0);
    MCTruth_particles_py0[i] = mcParticle.Py(0);
    MCTruth_particles_pz0[i] = mcParticle.Pz(0);
    MCTruth_particles_e0[i] = mcParticle.E(0);
    MCTruth_particles_Rescatter[i] = mcParticle.Rescatter();
    MCTruth_particles_polx[i] = mcParticle.Polarization().X();
    MCTruth_particles_poly[i] = mcParticle.Polarization().Y();
    MCTruth_particles_polz[i] = mcParticle.Polarization().Z();
  }

  const simb::MCNeutrino& mcNeutrino = mcTruth->GetNeutrino();

  MCTruth_neutrino_CCNC = mcNeutrino.CCNC();
  MCTruth_neutrino_mode = mcNeutrino.Mode();
  MCTruth_neutrino_interactionType = mcNeutrino.InteractionType();
  MCTruth_neutrino_target = mcNeutrino.Target();
  MCTruth_neutrino_nucleon = mcNeutrino.HitNuc();
  MCTruth_neutrino_quark = mcNeutrino.HitQuark();
  MCTruth_neutrino_W = mcNeutrino.W();
  MCTruth_neutrino_X = mcNeutrino.X();
  MCTruth_neutrino_Y = mcNeutrino.Y();
  MCTruth_neutrino_QSqr = mcNeutrino.QSqr();

  GTruth_IsSeaQuark = gTruth->fIsSeaQuark;
  GTruth_tgtPDG = gTruth->ftgtPDG;
  GTruth_weight = gTruth->fweight;
  GTruth_probability = gTruth->fprobability;
  GTruth_Xsec = gTruth->fXsec;
  GTruth_DiffXsec = gTruth->fDiffXsec;
  GTruth_vertexX = gTruth->fVertex.X();
  GTruth_vertexY = gTruth->fVertex.Y();
  GTruth_vertexZ = gTruth->fVertex.Z();
  GTruth_vertexT = gTruth->fVertex.T();
  GTruth_Gscatter = gTruth->fGscatter;
  GTruth_Gint = gTruth->fGint;
  GTruth_ResNum = gTruth->fResNum;
  GTruth_NumPiPlus = gTruth->fNumPiPlus;
  GTruth_NumPi0 = gTruth->fNumPi0;
  GTruth_NumPiMinus = gTruth->fNumPiMinus;
  GTruth_NumProton = gTruth->fNumProton;
  GTruth_NumNeutron = gTruth->fNumNeutron;
  GTruth_IsCharm = gTruth->fIsCharm;
  GTruth_gX = gTruth->fgX;
  GTruth_gY = gTruth->fgY;
  GTruth_gT = gTruth->fgT;
  GTruth_gW = gTruth->fgW;
  GTruth_gQ2 = gTruth->fgQ2;
  GTruth_gq2 = gTruth->fgq2;
  GTruth_ProbePDG = gTruth->fProbePDG;
  GTruth_ProbeP4x = gTruth->fProbeP4.X();
  GTruth_ProbeP4y = gTruth->fProbeP4.Y();
  GTruth_ProbeP4z = gTruth->fProbeP4.Z();
  GTruth_ProbeP4E = gTruth->fProbeP4.E();
  GTruth_HitNucP4x = gTruth->fHitNucP4.X();
  GTruth_HitNucP4y = gTruth->fHitNucP4.Y();
  GTruth_HitNucP4z = gTruth->fHitNucP4.Z();
  GTruth_HitNucP4E = gTruth->fHitNucP4.E();
  GTruth_FShadSystP4x = gTruth->fFShadSystP4.X();
  GTruth_FShadSystP4y = gTruth->fFShadSystP4.Y();
  GTruth_FShadSystP4z = gTruth->fFShadSystP4.Z();
  GTruth_FShadSystP4E = gTruth->fFShadSystP4.E();
  
  tree->Fill();

}

DEFINE_ART_MODULE(TestModule)
