// BunkerMC - Step 1 of the showcase: validate Geant4 against NCRP 151 TVLs
// for a 6 MV photon beam in ordinary concrete (2.35 g/cm3).
// Geometry: 200 cm concrete slab split into 20 x 10 cm layers.
// Variance reduction: geometric importance splitting (importance 2^i per layer).
// v4 (02/10/2026): air volume behind the slab carries the last-layer importance and the
// slab spans the full world width, so no importance jump > 2 remains (fixes GeomBias1001).
// Set BUNKER_ANALOG=1 in the environment to run without biasing (bias check).
// Scoring: photon fluence spectra (per cm2, per source photon) on each 10 cm plane,
// written as CSV histograms; H*(10) and TVLs are computed in analyze.py.
#include "G4RunManagerFactory.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserRunAction.hh"
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
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>

namespace cfg {
  constexpr G4int    nLayers   = 20;          // 20 x 10 cm = 200 cm concrete
  constexpr G4double layer     = 10 * cm;
  constexpr G4double halfXY    = 150 * cm;    // slab lateral half-size
  constexpr G4double scoreHalf = 25 * cm;     // central 50x50 cm scoring window
  constexpr G4int    nBins     = 70;
  constexpr G4double eMax      = 7 * MeV;
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
    G4AnalysisManager::Instance()->FillH1(plane, post->GetKineticEnergy() / MeV, w);
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
    G4cout << "NPRIMARIES " << r->GetNumberOfEvent() << G4endl;
  }
};

class Actions : public G4VUserActionInitialization {
 public:
  void Build() const override {
    SetUserAction(new Primary);
    SetUserAction(new RunAction);
    SetUserAction(new Stepping);
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
