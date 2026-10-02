// MazeMC v2 - 6 MV bunker maze door, audit-driven version.
//  * 40 cm water phantom at isocentre; absorbed dose scored in a 4x4x4 cm water voxel at isocentre
//    -> door dose normalised per Gy at isocentre in phantom (NCRP workload convention).
//  * Realistic door opening 1.2 m x 2.1 m in the maze end wall; optional lead door (argv[2], mm).
//  * H*(10) scored by track length in a 10 cm deep volume 30 cm beyond the door:
//    HFULL = averaged over the door area, HCENT = central 50x50 cm at 1.05 m above floor.
//  * Importance splitting along the maze; per-history statistics.
// Geometry (m, x east, y north, z up, isocentre at origin): room x[-3.5,3.5] y[-4,2.5] z[-1.5,1.5];
// inner maze wall x[-3.5,1.5] y[2.5,2.5+tW]; maze corridor x[-3.5,3.5] y[2.5+tW,4.0+tW]; door at x = -3.5.
// v3: inner maze wall thickness tW (argv[3], cm, default 50). Corridor width (1.5 m) and outer
//     maze wall (1.2 m) are kept; the block is extended north by (tW - 0.5 m).
#include "G4RunManagerFactory.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserEventAction.hh"
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
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include <cmath>
#include <vector>
#include <utility>
#include <cstdlib>
#include <string>
#include <algorithm>

namespace cfg {
  G4double tPb = 0.;                                   // lead door thickness
  G4double tW = 0.5 * m;                               // inner maze wall thickness
  G4double doorY = 3.75 * m;                           // door centre y = 3.25 m + tW (set in main)
  const G4double doorZ = -0.45 * m;                    // opening z[-1.5,0.6]
  const G4double hdY = 0.6 * m, hdZ = 1.05 * m;
  const G4double scoreHalfX = 5 * cm;                  // 10 cm deep scoring slab
  G4double vFull = 0., vCent = 0.;                     // cm3
  const G4double voxHalf = 2 * cm;                     // 4x4x4 cm dose voxel
}

namespace hstar {  // ICRP 74 photon H*(10)/fluence, pSv cm2
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

namespace tally {  // per-history: 0 HFULL (pSv), 1 HCENT (pSv), 2 DISO, 3 DMAX, 4 D10 (Gy)
  const int n = 5;
  G4double ev[n] = {0}, s[n] = {0}, s2[n] = {0};
}

class Detector : public G4VUserDetectorConstruction {
 public:
  Detector() { Build(); }
  G4VPhysicalVolume* Construct() override { return fWorld; }
  G4VPhysicalVolume* World() const { return fWorld; }
  void CreateImportanceStore() {
    G4IStore* is = G4IStore::GetInstance();
    for (auto& c : fCells) is->AddImportanceGeometryCell(c.second, *c.first, 0);
  }

 private:
  G4VPhysicalVolume* Box(const G4String& n, G4ThreeVector cW, G4ThreeVector half, G4Material* mat,
                         G4LogicalVolume* mother, G4ThreeVector motherW, G4double imp) {
    auto l = new G4LogicalVolume(new G4Box(n, half.x(), half.y(), half.z()), mat, n);
    auto p = new G4PVPlacement(nullptr, cW - motherW, l, n, mother, false, 0, true);
    fCells.push_back({p, imp});
    return p;
  }
  void Build() {
    auto nist = G4NistManager::Instance();
    auto air = nist->FindOrBuildMaterial("G4_AIR");
    auto water = nist->FindOrBuildMaterial("G4_WATER");
    auto conc = nist->BuildMaterialWithNewDensity("Concrete_235", "G4_CONCRETE", 2.35 * g / cm3);
    auto pb = nist->FindOrBuildMaterial("G4_Pb");

    auto wL = new G4LogicalVolume(new G4Box("World", 7 * m, 7 * m, 4 * m), air, "World");
    fWorld = new G4PVPlacement(nullptr, {}, wL, "World", nullptr, false, 0);
    fCells.push_back({fWorld, 1.});

    const G4ThreeVector O, cB(0, cfg::tW / 2, 0), cR(0, -0.75 * m, 0);
    auto bL = Box("Concrete", cB, {5.5 * m, 5.2 * m + cfg::tW / 2, 3.0 * m}, conc, wL, O, 1.)->GetLogicalVolume();
    auto rL = Box("Room", cR, {3.5 * m, 3.25 * m, 1.5 * m}, air, bL, cB, 1.)->GetLogicalVolume();
    auto phL = Box("Phantom", O, {20 * cm, 20 * cm, 20 * cm}, water, rL, cR, 1.)->GetLogicalVolume();
    Box("DoseVoxel", O, {cfg::voxHalf, cfg::voxHalf, cfg::voxHalf}, water, phL, O, 1.);
    // Beam-axis voxels for gantry 0 (beam enters phantom top at z = +20 cm, SSD 80 cm):
    Box("DoseDmax", {0, 0, 18.5 * cm}, {2 * cm, 2 * cm, 0.5 * cm}, water, phL, O, 1.);   // 1.5 cm depth
    Box("Dose10",   {0, 0, 10.0 * cm}, {2 * cm, 2 * cm, 1.0 * cm}, water, phL, O, 1.);   // 10 cm depth
    auto reg = new G4Region("PhantomRegion");                // fine cuts only where dose is scored
    reg->AddRootLogicalVolume(phL);
    auto pc = new G4ProductionCuts; pc->SetProductionCut(1 * mm); reg->SetProductionCuts(pc);

    const G4double yC = 3.25 * m + cfg::tW;                            // corridor centre
    Box("MazeGap", {2.5 * m, 2.5 * m + cfg::tW / 2, 0}, {1.0 * m, cfg::tW / 2, 1.5 * m},  air, bL, cB, 2.);
    Box("MazeC1",  {2.5 * m, yC, 0},   {1.0 * m, 0.75 * m, 1.5 * m},  air, bL, cB, 4.);
    Box("MazeC2",  {0.25 * m, yC, 0},  {1.25 * m, 0.75 * m, 1.5 * m}, air, bL, cB, 16.);
    Box("MazeC3",  {-2.25 * m, yC, 0}, {1.25 * m, 0.75 * m, 1.5 * m}, air, bL, cB, 64.);

    const G4ThreeVector cP(-4.5 * m, cfg::doorY, cfg::doorZ);           // door opening through W wall
    auto pL = Box("Passage", cP, {1.0 * m, cfg::hdY, cfg::hdZ}, air, bL, cB, 64.)->GetLogicalVolume();
    if (cfg::tPb > 0)
      Box("LeadDoor", {-3.5 * m - cfg::tPb / 2, cfg::doorY, cfg::doorZ}, {cfg::tPb / 2, cfg::hdY - 1 * um, cfg::hdZ - 1 * um},
          pb, pL, cP, 64.);
    const G4ThreeVector cS(-3.5 * m - cfg::tPb - 30 * cm - cfg::scoreHalfX, cfg::doorY, cfg::doorZ);
    const G4ThreeVector hS(cfg::scoreHalfX, cfg::hdY - 1 * mm, cfg::hdZ - 1 * mm);
    auto sL = Box("ScoreFull", cS, hS, air, pL, cP, 64.)->GetLogicalVolume();
    const G4ThreeVector hC(cfg::scoreHalfX - 0.1 * mm, 25 * cm, 25 * cm);
    Box("ScoreCentral", cS, hC, air, sL, cS, 64.);
    cfg::vFull = 8 * hS.x() * hS.y() * hS.z() / cm3;
    cfg::vCent = 8 * hC.x() * hC.y() * hC.z() / cm3;
  }
  G4VPhysicalVolume* fWorld = nullptr;
  std::vector<std::pair<G4VPhysicalVolume*, G4double>> fCells;
};

class Primary : public G4VUserPrimaryGeneratorAction {
 public:
  void GeneratePrimaries(G4Event* e) override { fGPS.GeneratePrimaryVertex(e); }
 private:
  G4GeneralParticleSource fGPS;
};

class Stepping : public G4UserSteppingAction {
 public:
  void UserSteppingAction(const G4Step* s) override {
    const auto pre = s->GetPreStepPoint();
    const G4String& vol = pre->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
    const G4double w = pre->GetWeight();
    if (vol == "DoseVoxel" || vol == "DoseDmax" || vol == "Dose10") {
      const int i = (vol == "DoseVoxel") ? 2 : (vol == "DoseDmax") ? 3 : 4;
      const G4double mass = (i == 2 ? 64. : i == 3 ? 16. : 32.) * 1e-3;                 // kg (water)
      tally::ev[i] += w * s->GetTotalEnergyDeposit() / joule / mass;                    // Gy
      return;
    }
    const bool full = (vol == "ScoreFull"), cent = (vol == "ScoreCentral");
    if (!full && !cent) return;
    if (s->GetTrack()->GetDefinition() != G4Gamma::Definition()) return;
    const G4double L = s->GetStepLength() / cm, e = pre->GetKineticEnergy() / MeV, hh = hstar::h(e);
    tally::ev[0] += w * L / cfg::vFull * hh;                                            // full-area average
    G4AnalysisManager::Instance()->FillH1(0, e, w * L / cfg::vFull);
    if (cent) tally::ev[1] += w * L / cfg::vCent * hh;
  }
};

class EventAction : public G4UserEventAction {
 public:
  void BeginOfEventAction(const G4Event*) override { for (int i = 0; i < tally::n; ++i) tally::ev[i] = 0; }
  void EndOfEventAction(const G4Event*) override {
    for (int i = 0; i < tally::n; ++i) { tally::s[i] += tally::ev[i]; tally::s2[i] += tally::ev[i] * tally::ev[i]; }
  }
};

class RunAction : public G4UserRunAction {
 public:
  RunAction() {
    auto am = G4AnalysisManager::Instance();
    am->SetDefaultFileType("csv"); am->SetVerboseLevel(0); am->SetFileName("door");
    am->CreateH1("door", "photon fluence 30 cm beyond door (cm-2 per bin)", 70, 0., 7.);
  }
  void BeginOfRunAction(const G4Run*) override {
    auto am = G4AnalysisManager::Instance();
    am->SetFileName(am->GetFileName() + "_Pb" + std::to_string((int)std::lround(cfg::tPb / mm))
                   + "_W" + std::to_string((int)std::lround(cfg::tW / cm)));
    am->OpenFile();
  }
  void EndOfRunAction(const G4Run* r) override {
    auto am = G4AnalysisManager::Instance(); am->Write(); am->CloseFile();
    const G4double N = r->GetNumberOfEvent();
    const char* tag[5] = {"HFULL", "HCENT", "DISO", "DMAX", "D10"};
    G4cout << "NPRIMARIES " << r->GetNumberOfEvent() << G4endl;
    for (int i = 0; i < tally::n; ++i) {
      const G4double mean = tally::s[i] / N, var = (tally::s2[i] / N - mean * mean) / (N - 1);
      G4cout << tag[i] << " " << mean << " " << (mean > 0 ? std::sqrt(std::max(var, 0.)) / mean : -1) << G4endl;
    }
  }
};

class Actions : public G4VUserActionInitialization {
 public:
  void Build() const override {
    SetUserAction(new Primary); SetUserAction(new RunAction);
    SetUserAction(new Stepping); SetUserAction(new EventAction);
  }
};

int main(int argc, char** argv) {
  if (argc > 2) cfg::tPb = std::atof(argv[2]) * mm;
  if (argc > 3) cfg::tW = std::atof(argv[3]) * cm;
  cfg::doorY = 3.25 * m + cfg::tW;
  G4cout << "INNERWALL_CM " << cfg::tW / cm << " PB_MM " << cfg::tPb / mm << G4endl;
  auto rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  auto det = new Detector;
  rm->SetUserInitialization(det);
  G4PhysListFactory factory;
  auto phys = factory.GetReferencePhysList("QBBC");
  phys->SetDefaultCutValue(1 * cm);
  G4GeometrySampler sampler(det->World(), "gamma");
  phys->RegisterPhysics(new G4ImportanceBiasing(&sampler));
  rm->SetUserInitialization(phys);
  rm->SetUserInitialization(new Actions);
  rm->Initialize();
  det->CreateImportanceStore();
  if (argc > 1) G4UImanager::GetUIpointer()->ApplyCommand("/control/execute " + G4String(argv[1]));
  delete rm;
  return 0;
}
