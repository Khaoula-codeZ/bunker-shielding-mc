// BunkerMC - Step 1 of the showcase: validate Geant4 against NCRP 151 TVLs
// for a 6 MV photon beam in ordinary concrete (2.35 g/cm3).
// Geometry: 200 cm concrete slab split into 20 x 10 cm layers.
// Variance reduction: geometric importance splitting (importance 2^i per layer).
// v4 (02/10/2026): air volume behind the slab carries the last-layer importance and the
// slab spans the full world width, so no importance jump > 2 remains (fixes GeomBias1001).
// Set BUNKER_ANALOG=1 in the environment to run without biasing (bias check).
// Scoring: photon fluence spectra (per cm2, per source photon) on each 10 cm plane,
// written as CSV histograms; H*(10) and TVLs are computed in analyze.py.
// v5 (03/10/2026): per-history H*(10) estimator on every plane (sum and sum of squares of the
// history score, written to hstat.csv). Uses the same ICRP 74 coefficients evaluated at the
// histogram bin centre, so its mean reproduces the histogram result exactly; its variance
// accounts for the correlation between split copies of one history (DV01 open action 5).
#include "G4RunManagerFactory.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "G4GeometrySampler.hh"
#include "G4ImportanceBiasing.hh"
#include "G4IStore.hh"
#include "G4UImanager.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4Run.hh"
#include "G4Gamma.hh"
#include "G4Event.hh"
#include <algorithm>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <iomanip>

namespace cfg {
  constexpr G4int    nLayers   = 20;          // 20 x 10 cm = 200 cm concrete
  constexpr G4double layer     = 10 * cm;
  constexpr G4double halfXY    = 150 * cm;    // slab lateral half-size
  constexpr G4double scoreHalf = 25 * cm;     // central 50x50 cm scoring window
  constexpr G4int    nBins     = 70;
  constexpr G4double eMax      = 7 * MeV;
}

namespace hstar {  // ICRP 74 photon H*(10)/fluence, pSv cm2 (log-log interpolation, clamped)
  const G4double E[25] = {0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.08,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.8,1,1.5,2,3,4,5,6,8,10};
  const G4double H[25] = {0.061,0.83,1.05,0.81,0.64,0.55,0.51,0.53,0.61,0.89,1.20,1.80,2.38,2.93,3.44,4.38,5.20,6.90,8.60,11.1,13.4,15.5,17.6,21.6,25.6};
  inline G4double h(G4double e) {
    if (e <= E[0]) return H[0];
    if (e >= E[24]) return H[24];
    int i = 0; while (E[i + 1] < e) ++i;
    const G4double t = std::log(e / E[i]) / std::log(E[i + 1] / E[i]);
    return std::exp(std::log(H[i]) + t * std::log(H[i + 1] / H[i]));
  }
}

namespace tally {  // per-history H*(10) score per plane (pSv per source photon, unnormalised by N)
  std::vector<G4double> ev(cfg::nLayers + 1, 0.), s(cfg::nLayers + 1, 0.), s2(cfg::nLayers + 1, 0.);
}

class Detector : public G4VUserDetectorConstruction {
 public:
  Detector() { Build(); }                      // built eagerly: sampler needs the world
  G4VPhysicalVolume* Construct() override { return fWorld; }
  G4VPhysicalVolume* World() const { return fWorld; }

  void CreateImportanceStore() {
    G4IStore* is = G4IStore::GetInstance();
    is->AddImportanceGeometryCell(1., *fWorld);
    for (G4int i = 0; i < cfg::nLayers; ++i)
      is->AddImportanceGeometryCell(std::pow(2., i), *fLayers[i], i);
    is->AddImportanceGeometryCell(std::pow(2., cfg::nLayers - 1), *fBack, 0);
  }

 private:
  void Build() {
    auto nist = G4NistManager::Instance();
    auto air  = nist->FindOrBuildMaterial("G4_AIR");
    auto conc = nist->BuildMaterialWithNewDensity("Concrete_235", "G4_CONCRETE", 2.35 * g / cm3);

    // World lateral size = slab size: no low-importance air beside the slab
    auto worldS = new G4Box("World", cfg::halfXY, cfg::halfXY, 300 * cm);
    auto worldL = new G4LogicalVolume(worldS, air, "World");
    fWorld = new G4PVPlacement(nullptr, {}, worldL, "World", nullptr, false, 0);

    auto layS = new G4Box("Layer", cfg::halfXY, cfg::halfXY, cfg::layer / 2);
    for (G4int i = 0; i < cfg::nLayers; ++i) {
      auto layL = new G4LogicalVolume(layS, conc, "Layer" + std::to_string(i));
      fLayers.push_back(new G4PVPlacement(nullptr, G4ThreeVector(0, 0, (i + 0.5) * cfg::layer),
                                          layL, "Layer" + std::to_string(i), worldL, false, i));
    }
    // Air region behind the slab (z = 200 to 300 cm) with the last-layer importance
    const G4double zEnd = cfg::nLayers * cfg::layer;
    auto backS = new G4Box("Back", cfg::halfXY, cfg::halfXY, (300 * cm - zEnd) / 2);
    auto backL = new G4LogicalVolume(backS, air, "Back");
    fBack = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, (300 * cm + zEnd) / 2),
                              backL, "Back", worldL, false, 0);
  }
  G4VPhysicalVolume* fWorld = nullptr;
  std::vector<G4VPhysicalVolume*> fLayers;
  G4VPhysicalVolume* fBack = nullptr;
};

class Primary : public G4VUserPrimaryGeneratorAction {
 public:
  void GeneratePrimaries(G4Event* e) override { fGPS.GeneratePrimaryVertex(e); }
 private:
  G4GeneralParticleSource fGPS;   // configured in run6MV.mac
};

class Stepping : public G4UserSteppingAction {
 public:
  void UserSteppingAction(const G4Step* s) override {
    auto post = s->GetPostStepPoint();
    if (post->GetStepStatus() != fGeomBoundary) return;
    if (s->GetTrack()->GetDefinition() != G4Gamma::Definition()) return;
    const G4ThreeVector p = post->GetPosition();
    if (std::abs(p.x()) > cfg::scoreHalf || std::abs(p.y()) > cfg::scoreHalf) return;
    const G4double k = p.z() / cfg::layer;
    const G4int plane = static_cast<G4int>(std::lround(k));
    if (std::abs(k - plane) > 1e-6 || plane < 0 || plane > cfg::nLayers) return;
    G4double cosz = std::abs(post->GetMomentumDirection().z());
    if (cosz < 0.05) cosz = 0.05;                                   // grazing-angle cap
    const G4double area = std::pow(2 * cfg::scoreHalf / cm, 2);     // cm2
    // Pre-step weight = weight while crossing, unaffected by split/roulette at the boundary
    const G4double w = s->GetPreStepPoint()->GetWeight() / cosz / area;
    const G4double e = post->GetKineticEnergy() / MeV;
    G4AnalysisManager::Instance()->FillH1(plane, e, w);
    // Per-history tally: same bin-centre H*(10) coefficient as analyze.py; under/overflow skipped
    const G4double binW = cfg::eMax / MeV / cfg::nBins;
    if (e >= 0. && e < cfg::eMax / MeV) {
      const G4int b = static_cast<G4int>(e / binW);
      tally::ev[plane] += w * hstar::h((b + 0.5) * binW);
    }
  }
};

class EventAction : public G4UserEventAction {
 public:
  void BeginOfEventAction(const G4Event*) override { std::fill(tally::ev.begin(), tally::ev.end(), 0.); }
  void EndOfEventAction(const G4Event*) override {
    for (std::size_t i = 0; i < tally::ev.size(); ++i) {
      tally::s[i] += tally::ev[i];
      tally::s2[i] += tally::ev[i] * tally::ev[i];
    }
  }
};

class RunAction : public G4UserRunAction {
 public:
  RunAction() {
    auto am = G4AnalysisManager::Instance();
    am->SetDefaultFileType("csv");
    am->SetVerboseLevel(0);
    for (G4int i = 0; i <= cfg::nLayers; ++i)
      am->CreateH1("plane" + std::to_string(i), "photon fluence spectrum (cm-2 per bin)",
                   cfg::nBins, 0., cfg::eMax / MeV);
  }
  void BeginOfRunAction(const G4Run*) override { G4AnalysisManager::Instance()->OpenFile("bunker"); }
  void EndOfRunAction(const G4Run* r) override {
    auto am = G4AnalysisManager::Instance();
    am->Write();
    am->CloseFile();
    const G4int N = r->GetNumberOfEvent();
    G4cout << "NPRIMARIES " << N << G4endl;
    std::ofstream f("hstat.csv");
    f << "plane,depth_cm,N,sum,sum2\n" << std::setprecision(17);
    for (std::size_t i = 0; i < tally::s.size(); ++i) {
      f << i << "," << i * cfg::layer / cm << "," << N << "," << tally::s[i] << "," << tally::s2[i] << "\n";
      const G4double mean = tally::s[i] / N;
      const G4double var = (N > 1) ? (tally::s2[i] / N - mean * mean) / (N - 1) : 0.;
      G4cout << "HSTAT " << i * cfg::layer / cm << " " << mean << " "
             << (mean > 0 ? std::sqrt(std::max(var, 0.)) / mean : -1) << G4endl;
    }
  }
};

class Actions : public G4VUserActionInitialization {
 public:
  void Build() const override {
    SetUserAction(new Primary);
    SetUserAction(new RunAction);
    SetUserAction(new Stepping);
    SetUserAction(new EventAction);
  }
};

int main(int argc, char** argv) {
  auto rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  auto det = new Detector;
  rm->SetUserInitialization(det);

  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = factory.GetReferencePhysList("QBBC");
  phys->SetDefaultCutValue(1 * cm);                     // electrons irrelevant for H*(10) here
  const bool analog = std::getenv("BUNKER_ANALOG") != nullptr;
  G4GeometrySampler sampler(det->World(), "gamma");
  if (!analog) phys->RegisterPhysics(new G4ImportanceBiasing(&sampler));
  G4cout << (analog ? "MODE ANALOG (no biasing)" : "MODE IMPORTANCE SPLITTING") << G4endl;
  rm->SetUserInitialization(phys);
  rm->SetUserInitialization(new Actions);
  rm->Initialize();
  if (!analog) det->CreateImportanceStore();

  if (argc > 1) G4UImanager::GetUIpointer()->ApplyCommand("/control/execute " + G4String(argv[1]));
  delete rm;
  return 0;
}
