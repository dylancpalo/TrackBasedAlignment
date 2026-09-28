// Single-panel KalSeed analyzer for tracker selection and hit/POCA diagnostics.

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/RecoDataProducts/inc/KalIntersection.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/DataProducts/inc/StrawId.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TVector3.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

namespace mu2e {

class AlignmentTrackAnalyzer : public art::EDAnalyzer {

    public:
    // Runtime configuration read from the FHiCL file.
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> kalSeedsTag{Name("kalSeedsTag"), Comment("Input tag for a KalSeedCollection.")};
      fhicl::Atom<int>              maxPrint{Name("maxPrint"),    Comment("Maximum number of events to print.")};
      fhicl::Atom<double>               tmin{Name("tmin"),        Comment("Fiducial time cut.")};

      // Track/hit selection parameters.  Defaults are supplied in FHiCL.
      fhicl::Atom<double> maxChi2Dof{Name("maxChi2Dof"), Comment("Maximum track chi2/ndof.")};
      fhicl::Atom<double> minAbsFitDOCA{Name("minAbsFitDOCA"), Comment("Minimum absolute fitted DOCA for a good hit [mm].")};
      fhicl::Atom<int> minEvenStrawsPerPanel{Name("minEvenStrawsPerPanel"), Comment("Minimum number of even straws for a quality panel intersection.")};
      fhicl::Atom<int> minOddStrawsPerPanel{Name("minOddStrawsPerPanel"), Comment("Minimum number of odd straws for a quality panel intersection.")};
      fhicl::Atom<int> minGoodHitsPerQualityPanel{Name("minGoodHitsPerQualityPanel"), Comment("Minimum number of good hits in a quality panel intersection.")};
      fhicl::Atom<int> minGoodHits{Name("minGoodHits"), Comment("Minimum number of good hits on the track.")};
      fhicl::Atom<int> minQualityPanelIntersections{Name("minQualityPanelIntersections"), Comment("Minimum number of quality panel intersections on the track.")};
      fhicl::Atom<int> maxQualityPanelIntersections{Name("maxQualityPanelIntersections"), Comment("Maximum number of quality panel intersections on the track; use -1 for no upper bound.")};
      fhicl::Atom<double> minXYSpan{Name("minXYSpan"), Comment("Minimum maximum pairwise XY span of quality hits in the global frame [mm].")};
      fhicl::Atom<double> minZSpan{Name("minZSpan"), Comment("Minimum Z span of quality hits in the global frame [mm].")};
    };
    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit AlignmentTrackAnalyzer(const Parameters& conf);

    TVector3 ComputePOCA2D(
    double DOCA,
    const TVector3& trackDir,
    int leftRight,
    const TVector3& wireStart,
    const TVector3& wireEnd
    );

    TVector3 LocalToGlobal(const TVector3 &pLocal,
    const TVector3 &wireEnd1,
    const TVector3 &wireEnd2);

    double computeQualityScore(const std::vector<double>& values, double maxPossibleSpan, double k);

    bool PointInsidePanelTrapezoid(double x, double y, int panelIndex) const;
    void InitializePanelTrapezoids(const Tracker& trk);

    void beginJob() override;

    void analyze( const art::Event& event) override;

    private:
    Config _conf;

    // Configurable selection values.
    art::ProductToken<KalSeedCollection>   _kalSeedsToken;
    double                                 _tmin;
    int                                    _maxPrint;
    double                                 _maxChi2Dof;
    double                                 _minAbsFitDOCA;
    int                                    _minEvenStrawsPerPanel;
    int                                    _minOddStrawsPerPanel;
    int                                    _minGoodHitsPerQualityPanel;
    int                                    _minGoodHits;
    int                                    _minQualityPanelIntersections;
    int                                    _maxQualityPanelIntersections;
    double                                 _minXYSpan;
    double                                 _minZSpan;

    int _nEvents = 0;

    // Four-point trapezoidal approximation to each of the 216 tracker panels.
    // A/B are the two endpoints of straw 0; C/D are the two endpoints of straw 94.
    bool _panelTrapezoidsInitialized = false;
    TVector3 _panelA[216];
    TVector3 _panelB[216];
    TVector3 _panelC[216];
    TVector3 _panelD[216];
    double _panelZ[216] = {0.0};

    // Histograms and profiles.
    TH1F* _hNTracks    = nullptr;
    TH1F* _hHasCalo    = nullptr;
    TH1F* _hnDOF       = nullptr;
    TH1F* _ht0         = nullptr;
    TH1F* _hp          = nullptr;
    TH1F* _hpErr       = nullptr;
    TH1F* _hnSkip      = nullptr;
    TProfile* _hTD     = nullptr;
    TProfile *_hDocaErrorRDrift = nullptr;
    TProfile *_hDocaErrorCDrift = nullptr;
    TProfile *_hDocaErrorUDOCA = nullptr;

    
    TH1F* _hNHitsQuality         = nullptr;
    TH1F* _hNPanelsQuality          = nullptr;
    TH1F* _hTrackPhi       = nullptr;
    TH1F* _hTrackTheta      = nullptr;

    

    TH1F* _hDOCAError      = nullptr;
    TH1F* _hXPOCAError      = nullptr;
    TH1F* _hYPOCAError      = nullptr;
    TH1F* _hZPOCAError      = nullptr;
    TH1F* _hXPOCAErrorRatio = nullptr;
    TH1F* _hYPOCAErrorRatio = nullptr;
    TH1F* _hZPOCAErrorRatio = nullptr;
    TH1F* _hFullSpan = nullptr;
    TH1F* _hXYSpan = nullptr;
    TH1F* _hZSpan = nullptr;

    TH1F* _hDOCAErrorCleaned      = nullptr;

    TH1F* _hdriftTime      = nullptr;
    TH1F* _hRDrift      = nullptr;
    TH1F* _hCDrift      = nullptr;
    TH1F* _hUDOCA      = nullptr;

    TH1F* _hdriftTimePanel[12*18];
    TH1F* _hEDepPanel[12*18];
    TH1F* _hDOCAErrorPanel[12*18];
    TH1F* _hXPOCAErrorPanel[12*18];
    TH1F* _hYPOCAErrorPanel[12*18];
    TH1F* _hZPOCAErrorPanel[12*18];
    TH1F* _hZSensitivityPanel[12*18]; // |(wireDir x trackDir)_Z|, unitless
    TProfile* _hDocaErrorUDOCAPanel[12*18];
    TProfile* _hDocaErrorSquaredUDOCAPanel[12*18];
    TH1F* _hUDOCAEntriesPanel[12*18];
    TProfile* _hDocaErrorSquaredUPosPanel[12*18];
    TH1F* _hUPosEntriesPanel[12*18];

    TH1F* _hNHits      = nullptr;
    TH1F* _hNPanels      = nullptr;
    TH1F* _hNPlanes      = nullptr;
    TH1F* _hmomX      = nullptr;
    TH1F* _hmomY      = nullptr;
    TH1F* _hmomZ      = nullptr;
    TH1F* _hchi2_dof      = nullptr;

    TH1F* _hTrackSigma      = nullptr;


    TH1F* _hqualityScore      = nullptr;
    TH2F* _hchi2_nhits      = nullptr;

    TH1F* _hHitLong      = nullptr;
    TH1F* _hTrackLong      = nullptr;
    TH1F* _hLongError      = nullptr;
    TProfile * _hZPOCAErrorY = nullptr;

    Int_t *nStrawsPerPanel = nullptr;
    Int_t *nPanelsPerPlane = nullptr;
    Int_t *nPlanesPerStation = nullptr;
    Int_t *nStations = nullptr;
    Int_t *nCounter = nullptr;
    TH1F *hWire[18*12];
    TH1F *hULocal[18*12];
    TH1F *hVErrorPanel[18*12];       // local V residual over the full panel
    TH1F *hXErrorPanel[18*12];       // global X residual over the full panel
    TH1F *hYErrorPanel[18*12];       // global Y residual over the full panel
    TProfile *hXErrorZPanelP[18*12]; // <global X residual> vs tshs._uupos
    TProfile *hYErrorZPanelP[18*12]; // <global Y residual> vs tshs._uupos
    TProfile *hVErrorZPanelP[18*12]; // <local V residual> vs tshs._uupos
    TProfile *hdVULocal[18*12];

    // Z residual spatial dependence, grouped by station (18 stations).
    TProfile *hZErrorVsXStation[18];
    TProfile *hZErrorVsYStation[18];
    TProfile2D *hZErrorXYStation[18];

    TH2I* _hPanelID2D      = nullptr;
    TH1I* _hPanelID      = nullptr;
    TH1I* _hPanelIntersected = nullptr;
    TH1I* _hPanelIntersectedWithHit = nullptr;
    TH2F* _hnPanels2D      = nullptr;
    TH2I* _hEvenOddPanel      = nullptr;
    TProfile *hZPOCAErrorWire[18*12];
    TProfile *hZPOCAErrorWireEven[18*12];
    TProfile *hZPOCAErrorWireOdd[18*12];

    TProfile *hZPOCAErrorWireEven_PosZ[18*12];
    TProfile *hZPOCAErrorWireOdd_PosZ[18*12];
    TProfile *hZPOCAErrorWireEven_NegZ[18*12];
    TProfile *hZPOCAErrorWireOdd_NegZ[18*12];

    TProfile *hVPOCAErrorWire[18*12];
    TProfile *hVPOCAErrorWireEven[18*12];
    TProfile *hVPOCAErrorWireOdd[18*12];

    TProfile *hVPOCAErrorWireEven_PosZ[18*12];
    TProfile *hVPOCAErrorWireOdd_PosZ[18*12];
    TProfile *hVPOCAErrorWireEven_NegZ[18*12];
    TProfile *hVPOCAErrorWireOdd_NegZ[18*12];

    TProfile *hVErrorZWireP[96*12*18];
    TProfile *hWErrorZWireP[96*12*18];
    TH2F *hVErrorZWire2D[96*12*18];
    TH2F *hWErrorZWire2D[96*12*18];

    TH2F *hXYPOCA  = nullptr;
    TH2F *hYZPOCA  = nullptr;

};

AlignmentTrackAnalyzer::AlignmentTrackAnalyzer(const Parameters& conf)
  : art::EDAnalyzer(conf),
  _conf(conf()),
  _kalSeedsToken{consumes<KalSeedCollection>(conf().kalSeedsTag())},
  _tmin(conf().tmin()),
  _maxPrint(conf().maxPrint()),
  _maxChi2Dof(conf().maxChi2Dof()),
  _minAbsFitDOCA(conf().minAbsFitDOCA()),
  _minEvenStrawsPerPanel(conf().minEvenStrawsPerPanel()),
  _minOddStrawsPerPanel(conf().minOddStrawsPerPanel()),
  _minGoodHitsPerQualityPanel(conf().minGoodHitsPerQualityPanel()),
  _minGoodHits(conf().minGoodHits()),
  _minQualityPanelIntersections(conf().minQualityPanelIntersections()),
  _maxQualityPanelIntersections(conf().maxQualityPanelIntersections()),
  _minXYSpan(conf().minXYSpan()),
  _minZSpan(conf().minZSpan()){
  }

double AlignmentTrackAnalyzer::computeQualityScore(const std::vector<double>& values, double maxPossibleSpan, double k) {
    size_t n = values.size();
    if (n < 2) return 0.0;

    std::vector<double> sortedValues = values;
    std::sort(sortedValues.begin(), sortedValues.end());

    double minVal = sortedValues.front();
    double maxVal = sortedValues.back();
    double span = maxVal - minVal;

    if (span == 0.0) return 0.0;

    if (maxPossibleSpan <= 0.0) {
      maxPossibleSpan = span;
    }

    // Find max gap
    double maxGap = 0.0;
    for (size_t i = 1; i < n; ++i) {
      double gap = sortedValues[i] - sortedValues[i - 1];
      if (gap > maxGap) maxGap = gap;
    }

    double spanScore = span / maxPossibleSpan;
    if (spanScore > 1.0) spanScore = 1.0;

    double gapScore = 1.0 - (maxGap / span);
    if (gapScore < 0.0) gapScore = 0.0;

    double countScore = 1.0 - std::exp(-k * (n - 1));

    return spanScore * gapScore * countScore;
  }

// Build one trapezoidal active-area proxy for every tracker panel.
// The indexing is panel + 6*plane, matching the existing panel histograms.
void AlignmentTrackAnalyzer::InitializePanelTrapezoids(const Tracker& trk) {
  for (int plane = 0; plane < 36; ++plane) {
    for (int panel = 0; panel < 6; ++panel) {
      const int panelIndex = panel + 6*plane;

      StrawId sid0(plane, panel, 0);
      StrawId sid94(plane, panel, 94);

      Straw const& straw0 = trk.getStraw(sid0);
      Straw const& straw94 = trk.getStraw(sid94);

      auto const pA = straw0.wirePosition(-straw0.halfLength());
      auto const pB = straw0.wirePosition( straw0.halfLength());
      auto const pC = straw94.wirePosition( straw94.halfLength());
      auto const pD = straw94.wirePosition(-straw94.halfLength());

      _panelA[panelIndex] = TVector3(pA.x(), pA.y(), pA.z());
      _panelB[panelIndex] = TVector3(pB.x(), pB.y(), pB.z());
      _panelC[panelIndex] = TVector3(pC.x(), pC.y(), pC.z());
      _panelD[panelIndex] = TVector3(pD.x(), pD.y(), pD.z());

      // Constant-Z approximation for this straw layer.
      _panelZ[panelIndex] = (_panelA[panelIndex].Z() +
                             _panelB[panelIndex].Z() +
                             _panelC[panelIndex].Z() +
                             _panelD[panelIndex].Z()) / 4.0;
    }
  }

  _panelTrapezoidsInitialized = true;
}

// 2D point-in-convex-quadrilateral test using the sign of the cross product
// for each successive edge A->B->C->D->A.  Points on an edge count as inside.
bool AlignmentTrackAnalyzer::PointInsidePanelTrapezoid(double x, double y,
                                                         int panelIndex) const {
  const TVector3& A = _panelA[panelIndex];
  const TVector3& B = _panelB[panelIndex];
  const TVector3& C = _panelC[panelIndex];
  const TVector3& D = _panelD[panelIndex];

  const double d1 = (B.X()-A.X())*(y-A.Y()) - (B.Y()-A.Y())*(x-A.X());
  const double d2 = (C.X()-B.X())*(y-B.Y()) - (C.Y()-B.Y())*(x-B.X());
  const double d3 = (D.X()-C.X())*(y-C.Y()) - (D.Y()-C.Y())*(x-C.X());
  const double d4 = (A.X()-D.X())*(y-D.Y()) - (A.Y()-D.Y())*(x-D.X());

  const bool hasNeg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0) || (d4 < 0.0);
  const bool hasPos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0) || (d4 > 0.0);

  return !(hasNeg && hasPos);
}

TVector3 AlignmentTrackAnalyzer::ComputePOCA2D(
  double DOCA,
  const TVector3& trackDir,
  int leftRight,
  const TVector3& wireStart,
  const TVector3& wireEnd
  ) {
    // 1. Local X along wire
    TVector3 xLocal = (wireEnd - wireStart).Unit();

    // 2. Local Z along global Z
    TVector3 zLocal(0,0,1);

    // 3. Local Y to complete right-handed system
    TVector3 yLocal = zLocal.Cross(xLocal).Unit();

    // 4. Perpendicular DOCA vector
    TVector3 vDOCA = xLocal.Cross(trackDir).Unit() * DOCA * leftRight;

    // 5. Project into local YZ plane
    double A = 0;                  // along wire
    double B = vDOCA.Dot(yLocal);  // local Y
    double C = vDOCA.Dot(zLocal);  // local Z

    return TVector3(A, B, C);
  }

  // Converts a local position (wire frame) to global coordinates.
  // The local (0,0,0) corresponds to the midpoint between wire endpoints.
TVector3 AlignmentTrackAnalyzer::LocalToGlobal(const TVector3 &pLocal,
  const TVector3 &wireEnd1,
  const TVector3 &wireEnd2)
  {
    // 1. Define origin = wire midpoint
    TVector3 origin = 0.5 * (wireEnd1 + wireEnd2);

    // 2. Direction of wire projected onto XY plane
    TVector3 dirXY(wireEnd2.X() - wireEnd1.X(),
    wireEnd2.Y() - wireEnd1.Y(),
    0.0);

    // 3. Build local coordinate axes
    TVector3 xHat;
    if (dirXY.Mag() < 1e-12)
    {
      // wire is vertical; pick arbitrary horizontal direction for xHat
      xHat = TVector3(1.0, 0.0, 0.0);
    }
    else
    {
      xHat = dirXY.Unit();
    }

    TVector3 zHat(0.0, 0.0, 1.0);
    TVector3 yHat = zHat.Cross(xHat).Unit();

    // 4. Convert local point to global coordinates
    TVector3 pGlobal = origin
    + pLocal.X() * xHat
    + pLocal.Y() * yHat
    + pLocal.Z() * zHat;

    return pGlobal;
  }

// Book histograms and profiles.
void AlignmentTrackAnalyzer::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    _hNTracks    = tfs->make<TH1F>("hNTracks", "Number of tracks per event.",                               10,  0.,    10. );
    _hnDOF       = tfs->make<TH1F>("hnDOF",    "Number of degrees of freedom in fit.",                     100,  0.,   100. );
    _hHasCalo    = tfs->make<TH1F>("hHasCalo", "Number of calorimeter hits.",                                2,  0.,     2. );
    _ht0         = tfs->make<TH1F>("ht0",      "Track time at mid-point of Tracker ;(ns)",                 100,  0.,  2000. );
    _hp          = tfs->make<TH1F>("hp",       "Track momentum at mid-point of tracker;( MeV/c)",          100, 70.,   120. );
    _hpErr       = tfs->make<TH1F>("hpErr",    "Error on track momentum at mid-point of tracker;( MeV/c)", 100,  0.,     2. );
    _hnSkip      = tfs->make<TH1F>("hnSkip",   "Cut tree for skipped tracks;( MeV/c)",                       3,  0.,     3. );
    _hTD         = tfs->make<TProfile>("hTD",   "DOCA [mm]: Drift Time",                                     100, -30, 70, 0, 3);
    _hdriftTime         = tfs->make<TH1F>("hdriftTime",   "Drift Time [ns]",                                     100, -30, 70);
    _hUDOCA         = tfs->make<TH1F>("_hUDOCA",   "_hUDOCA [mm]",                                     100, -1, 4);
    _hCDrift         = tfs->make<TH1F>("_hCDrift",   "_hCDrift [mm]",                                     100, -1, 4);
    _hRDrift         = tfs->make<TH1F>("_hRDrift",   "_hRDrift [mm]",                                     100, -1, 4);
    _hNHits         = tfs->make<TH1F>("_hNHits",   "nHits",                                     51, 0, 50);
    _hNPanels         = tfs->make<TH1F>("_hNPanels",   "_hNPanels",                                     13, 0, 12);
    _hNPlanes         = tfs->make<TH1F>("_hNPlanes",   "_hNPlanes",                                     2, 0, 1);
    _hmomX         = tfs->make<TH1F>("_hmomX",   "_hmomX",                                     50, -1, 1);
    _hmomY         = tfs->make<TH1F>("_hmomY",   "_hmomY",                                     50, -1, 1);
    _hmomZ         = tfs->make<TH1F>("_hmomZ",   "_hmomZ",                                     50, -1, 1);

    _hchi2_dof         = tfs->make<TH1F>("_hchi2_dof",   "_hchi2_dof",                                     100, 0, 5);
    _hXYSpan           = tfs->make<TH1F>("hXYSpan",      "Maximum pairwise XY span of quality hits;XY span [mm]", 100, 0., 1600.);
    _hZSpan            = tfs->make<TH1F>("hZSpan",       "Z span of quality hits;Z span [mm]", 100, 0., 4000.);
    _hPanelID         = tfs->make<TH1I>("_hPanelID",   "_hPanelID",                                     12, 0, 12);
    _hPanelIntersected = tfs->make<TH1I>("_hPanelIntersected",
                                         "Geometrically intersected panels;Panel ID;Tracks",
                                         216, -0.5, 215.5);
    _hPanelIntersectedWithHit = tfs->make<TH1I>("_hPanelIntersectedWithHit",
                                                "Geometrically intersected panels with >=1 track hit;Panel ID;Tracks",
                                                216, -0.5, 215.5);
    _hPanelID2D         = tfs->make<TH2I>("_hPanelID2D",   "_hPanelID2D",                                     12, 0, 12, 12, 0, 12);
    _hnPanels2D         = tfs->make<TH2F>("_hnPanels2D",   "_hnPanels2D",                                     50, 0, 49, 50, 0, 49);

    _hchi2_nhits         = tfs->make<TH2F>("_hchi2_nhits",   "_hchi2_nhits",                                     50,0, 10, 50, 0, 50);

    _hEvenOddPanel         = tfs->make<TH2I>("_hEvenOddPanel",   "_hEvenOddPanel",                                     50, 0, 49, 50, 0, 49);
    _hqualityScore         = tfs->make<TH1F>("_hqualityScore",   "_hqualityScore",                                     100, 0, 1);
    _hHitLong         = tfs->make<TH1F>("_hHitLong",   "_hHitLong",                                     100, -1000, 1000);
    _hTrackLong         = tfs->make<TH1F>("_hTrackLong",   "_hTrackLong",                                     100, -1000, 1000);
    _hLongError         = tfs->make<TH1F>("_hLongError",   "_hLongError",                                     100, -500, 500);
    _hFullSpan         = tfs->make<TH1F>("_hFullSpan",   "_hFullSpan",                                     100, 0, 1500);

    
    _hNHitsQuality         = tfs->make<TH1F>("_hNHitsQuality",   "_hNHitsQuality",                                     100, 0, 100);
    _hNPanelsQuality         = tfs->make<TH1F>("_hNPanelsQuality",   "_hNPanelsQuality",                                     50, 0, 50);
    _hTrackPhi         = tfs->make<TH1F>("_hTrackPhi",   "_hTrackPhi",                                     100, -2*3.14,-2*3.14);
    _hTrackTheta         = tfs->make<TH1F>("_hTrackTheta",   "_hTrackTheta",                                     100,  -2*3.14,-2*3.14);

    
    
    
    
    hXYPOCA = tfs->make<TH2F>("hXYPOCA",   "X POCA YPOCA",                      100,  -700, 700, 100, -700, 700);
    hYZPOCA = tfs->make<TH2F>("hYZPOCA",   "Y POCA ZPOCA",                      100,  -700, 700, 100, -100, 100);

    _hDocaErrorRDrift   = tfs->make<TProfile>("_hDocaErrorRDrift",   "Hit DOCA - Track DOCA[um]:Hit RDRIFT", 100, -0.5, 3, -750, 750);
    _hDocaErrorCDrift   = tfs->make<TProfile>("_hDocaErrorCDrift",   "Hit DOCA - Track DOCA[um]:Hit CDRIFT", 100, -0.5, 3, -750, 750);
    _hDocaErrorUDOCA   = tfs->make<TProfile>("_hDocaErrorUDOCA",   "Hit DOCA - Track DOCA[um]:Track UDOCA", 100, -0.5, 3, -750, 750);
    _hTrackSigma = tfs->make<TH1F>("_hTrackSigma",   "_hTrackSigma",                      100,  0, 1000);
    _hDOCAError = tfs->make<TH1F>("hDOCAError",   "Hit DOCA - Track DOCA [um]",                      100,  -2000, 2000);
    _hXPOCAError = tfs->make<TH1F>("_hXPOCAError",   "Hit POCA X - Track POCA X [um]",                      100,  -2000, 2000);
    _hYPOCAError = tfs->make<TH1F>("_hYPOCAError",   "Hit POCA Y - Track POCA Y [um]",                      100,  -2000, 2000);
    _hZPOCAError = tfs->make<TH1F>("_hZPOCAError",   "Hit POCA Z - Track POCA Z [um]",                      100,  -2000, 2000);

    _hDOCAErrorCleaned = tfs->make<TH1F>("_hDOCAErrorCleaned",   "Hit DOCA - Track DOCA [um]",                      100,  -2000, 2000);

    nStrawsPerPanel = tfs->make<Int_t>(98);
    nPanelsPerPlane = tfs->make<Int_t>(6);
    nPlanesPerStation = tfs->make<Int_t>(2);
    nStations = tfs->make<Int_t>(18);
    nCounter = tfs->make<Int_t>(0);

    _hZPOCAErrorY = tfs->make<TProfile>("_hZPOCAErrorY",   "_hZPOCAErrorY",                      10,-700, 700,  -2000, 2000);

    _hXPOCAErrorRatio = tfs->make<TH1F>("_hXPOCAErrorRatio",   "X POCA Error / DOCA Error",                      100,  -1, 1);
    _hYPOCAErrorRatio = tfs->make<TH1F>("_hYPOCAErrorRatio",   "Y POCA Error / DOCA Error",                      100,  -1, 1);
    _hZPOCAErrorRatio = tfs->make<TH1F>("_hZPOCAErrorRatio",   "Z POCA Error / DOCA Error",                      100,  -1, 1);

    char title[600];
    for (Int_t ihist=0;ihist<12*18;ihist++){
      sprintf(title, "hdriftTimePanel%i", ihist);
      _hdriftTimePanel[ihist]= tfs->make<TH1F>(title, title,  100, -30, 70);
      sprintf(title, "hEDepPanel%i", ihist);
      _hEDepPanel[ihist]= tfs->make<TH1F>(title, title, 100, 0, 5);
      _hEDepPanel[ihist]->GetXaxis()->SetTitle("tshs._edep #times 1000");
      sprintf(title, "hDOCAErrorPanel%i", ihist);
      _hDOCAErrorPanel[ihist]= tfs->make<TH1F>(title, title,  100, -2000, 2000);
      sprintf(title, "_hXPOCAErrorPanel%i", ihist);
      _hXPOCAErrorPanel[ihist]= tfs->make<TH1F>(title, title,  100, -2000, 2000);
      sprintf(title, "_hYPOCAErrorPanel%i", ihist);
      _hYPOCAErrorPanel[ihist]= tfs->make<TH1F>(title, title,  100, -2000, 2000);
      sprintf(title, "_hZPOCAErrorPanel%i", ihist);
      _hZPOCAErrorPanel[ihist]= tfs->make<TH1F>(title, title,  100, -700, 700);
      sprintf(title, "hZSensitivityPanel%i", ihist);
      _hZSensitivityPanel[ihist]= tfs->make<TH1F>(title, title, 100, 0., 1.);
      _hZSensitivityPanel[ihist]->GetXaxis()->SetTitle("|(wireDir #times trackDir)_{Z}|");

      sprintf(title, "hDocaErrorUDOCAPanel%i", ihist);
      _hDocaErrorUDOCAPanel[ihist]= tfs->make<TProfile>(title, title, 100, -0.5, 3., -750., 750.);
      _hDocaErrorUDOCAPanel[ihist]->GetXaxis()->SetTitle("Track UDOCA [mm]");
      _hDocaErrorUDOCAPanel[ihist]->GetYaxis()->SetTitle("<Hit DOCA - Track DOCA> [#mu m]");

      sprintf(title, "hDocaErrorSquaredUDOCAPanel%i", ihist);
      _hDocaErrorSquaredUDOCAPanel[ihist]= tfs->make<TProfile>(title, title, 100, -0.5, 3.);
      _hDocaErrorSquaredUDOCAPanel[ihist]->GetXaxis()->SetTitle("Track UDOCA [mm]");
      _hDocaErrorSquaredUDOCAPanel[ihist]->GetYaxis()->SetTitle("<(Hit DOCA - Track DOCA)^{2}> [#mu m^{2}]");

      sprintf(title, "hUDOCAEntriesPanel%i", ihist);
      _hUDOCAEntriesPanel[ihist]= tfs->make<TH1F>(title, title, 100, -0.5, 3.);
      _hUDOCAEntriesPanel[ihist]->GetXaxis()->SetTitle("Track UDOCA [mm]");
      _hUDOCAEntriesPanel[ihist]->GetYaxis()->SetTitle("Entries");

      sprintf(title, "hDocaErrorSquaredUPosPanel%i", ihist);
      _hDocaErrorSquaredUPosPanel[ihist]= tfs->make<TProfile>(title, title, 50, -1000., 1000.);
      _hDocaErrorSquaredUPosPanel[ihist]->GetXaxis()->SetTitle("tshs._uupos [mm]");
      _hDocaErrorSquaredUPosPanel[ihist]->GetYaxis()->SetTitle("<(Hit DOCA - Track DOCA)^{2}> [#mu m^{2}]");

      sprintf(title, "hUPosEntriesPanel%i", ihist);
      _hUPosEntriesPanel[ihist]= tfs->make<TH1F>(title, title, 50, -1000., 1000.);
      _hUPosEntriesPanel[ihist]->GetXaxis()->SetTitle("tshs._uupos [mm]");
      _hUPosEntriesPanel[ihist]->GetYaxis()->SetTitle("Entries");
    }

    for (Int_t ihist=0;ihist<12*18;ihist++){
      sprintf(title, "hWire%i", ihist);
      hWire[ihist]= tfs->make<TH1F>(title, title, (*nStrawsPerPanel), 0., double(*nStrawsPerPanel));

      sprintf(title, "hULocal%i", ihist);
      hULocal[ihist]= tfs->make<TH1F>(title, title, 100, -700, 700);

      sprintf(title, "hdVULocal%i", ihist);
      hdVULocal[ihist]= tfs->make<TProfile>(title, title, 20, -700, 700, -500, 500);

      sprintf(title, "hVErrorPanel%i", ihist);
      hVErrorPanel[ihist]= tfs->make<TH1F>(title, title, 100, -1000, 1000);
      hVErrorPanel[ihist]->GetXaxis()->SetTitle("Hit V - Track V [#mu m]");

      sprintf(title, "hXErrorPanel%i", ihist);
      hXErrorPanel[ihist]= tfs->make<TH1F>(title, title, 100, -1000, 1000);
      hXErrorPanel[ihist]->GetXaxis()->SetTitle("Hit X - Track X [#mu m]");

      sprintf(title, "hYErrorPanel%i", ihist);
      hYErrorPanel[ihist]= tfs->make<TH1F>(title, title, 100, -1000, 1000);
      hYErrorPanel[ihist]->GetXaxis()->SetTitle("Hit Y - Track Y [#mu m]");

      sprintf(title, "hXErrorZPanelP%i", ihist);
      hXErrorZPanelP[ihist]= tfs->make<TProfile>(title, title, 50, -1000, 1000, -1000, 1000);
      hXErrorZPanelP[ihist]->GetXaxis()->SetTitle("tshs._uupos [mm]");
      hXErrorZPanelP[ihist]->GetYaxis()->SetTitle("<Hit X - Track X> [#mu m]");

      sprintf(title, "hYErrorZPanelP%i", ihist);
      hYErrorZPanelP[ihist]= tfs->make<TProfile>(title, title, 50, -1000, 1000, -1000, 1000);
      hYErrorZPanelP[ihist]->GetXaxis()->SetTitle("tshs._uupos [mm]");
      hYErrorZPanelP[ihist]->GetYaxis()->SetTitle("<Hit Y - Track Y> [#mu m]");

      sprintf(title, "hVErrorZPanelP%i", ihist);
      hVErrorZPanelP[ihist]= tfs->make<TProfile>(title, title, 50, -1000, 1000, -1000, 1000);
      hVErrorZPanelP[ihist]->GetXaxis()->SetTitle("tshs._uupos [mm]");
      hVErrorZPanelP[ihist]->GetYaxis()->SetTitle("<Hit V - Track V> [#mu m]");

      sprintf(title, "hZPOCAErrorWire%i", ihist);
      hZPOCAErrorWire[ihist]= tfs->make<TProfile>(title, title, (*nStrawsPerPanel), 0., double(*nStrawsPerPanel), -500, 500);
      sprintf(title, "hZPOCAErrorWireEven%i", ihist);
      hZPOCAErrorWireEven[ihist]= tfs->make<TProfile>(title, title, (*nStrawsPerPanel), 0., double(*nStrawsPerPanel), -500, 500);
      sprintf(title, "hZPOCAErrorWireOdd%i", ihist);
      hZPOCAErrorWireOdd[ihist]= tfs->make<TProfile>(title, title, (*nStrawsPerPanel), 0., double(*nStrawsPerPanel), -500, 500);

      sprintf(title, "hZPOCAErrorWireEven_PosZ%i", ihist);
      hZPOCAErrorWireEven_PosZ[ihist]= tfs->make<TProfile>(title, title, (*nStrawsPerPanel), 0., double(*nStrawsPerPanel), -500, 500);
      sprintf(title, "hZPOCAErrorWireOdd_PosZ%i", ihist);
      hZPOCAErrorWireOdd_PosZ[ihist]= tfs->make<TProfile>(title, title, (*nStrawsPerPanel), 0., double(*nStrawsPerPanel), -500, 500);

      sprintf(title, "hZPOCAErrorWireEven_NegZ%i", ihist);
      hZPOCAErrorWireEven_NegZ[ihist]= tfs->make<TProfile>(title, title, (*nStrawsPerPanel), 0., double(*nStrawsPerPanel), -500, 500);
      sprintf(title, "hZPOCAErrorWireOdd_NegZ%i", ihist);
      hZPOCAErrorWireOdd_NegZ[ihist]= tfs->make<TProfile>(title, title, (*nStrawsPerPanel), 0., double(*nStrawsPerPanel), -500, 500);

      sprintf(title, "hVPOCAErrorWire%i", ihist);
      hVPOCAErrorWire[ihist]= tfs->make<TProfile>(title, title, (*nStrawsPerPanel), 0., double(*nStrawsPerPanel), -500, 500);
      sprintf(title, "hVPOCAErrorWireEven%i", ihist);
      hVPOCAErrorWireEven[ihist]= tfs->make<TProfile>(title, title, (*nStrawsPerPanel), 0., double(*nStrawsPerPanel), -500, 500);
      sprintf(title, "hVPOCAErrorWireOdd%i", ihist);
      hVPOCAErrorWireOdd[ihist]= tfs->make<TProfile>(title, title, (*nStrawsPerPanel), 0., double(*nStrawsPerPanel), -500, 500);

      sprintf(title, "hVPOCAErrorWireEven_PosZ%i", ihist);
      hVPOCAErrorWireEven_PosZ[ihist]= tfs->make<TProfile>(title, title, (*nStrawsPerPanel), 0., double(*nStrawsPerPanel), -500, 500);
      sprintf(title, "hVPOCAErrorWireOdd_PosZ%i", ihist);
      hVPOCAErrorWireOdd_PosZ[ihist]= tfs->make<TProfile>(title, title, (*nStrawsPerPanel), 0., double(*nStrawsPerPanel), -500, 500);

      sprintf(title, "hVPOCAErrorWireEven_NegZ%i", ihist);
      hVPOCAErrorWireEven_NegZ[ihist]= tfs->make<TProfile>(title, title, (*nStrawsPerPanel), 0., double(*nStrawsPerPanel), -500, 500);
      sprintf(title, "hVPOCAErrorWireOdd_NegZ%i", ihist);
      hVPOCAErrorWireOdd_NegZ[ihist]= tfs->make<TProfile>(title, title, (*nStrawsPerPanel), 0., double(*nStrawsPerPanel), -500, 500);

    }

    // Station-level Z residual versus global X/Y.
    // TProfile2D stores the average Z residual in each coarse XY bin.
    for (Int_t istation=0; istation<18; ++istation) {
      sprintf(title, "hZErrorVsXStation%i", istation);
      hZErrorVsXStation[istation] =
        tfs->make<TProfile>(title, title, 50, -1000., 1000., -1000., 1000.);
      hZErrorVsXStation[istation]->GetXaxis()->SetTitle("Track POCA X [mm]");
      hZErrorVsXStation[istation]->GetYaxis()->SetTitle("<Hit global Z - Track global Z> [#mu m]");

      sprintf(title, "hZErrorVsYStation%i", istation);
      hZErrorVsYStation[istation] =
        tfs->make<TProfile>(title, title, 50, -1000., 1000., -1000., 1000.);
      hZErrorVsYStation[istation]->GetXaxis()->SetTitle("Track POCA Y [mm]");
      hZErrorVsYStation[istation]->GetYaxis()->SetTitle("<Hit global Z - Track global Z> [#mu m]");

      sprintf(title, "hZErrorXYStation%i", istation);
      hZErrorXYStation[istation] =
        tfs->make<TProfile2D>(title, title,
                              20, -1000., 1000.,
                              20, -1000., 1000.,
                              -1000., 1000.);
      hZErrorXYStation[istation]->GetXaxis()->SetTitle("Track POCA X [mm]");
      hZErrorXYStation[istation]->GetYaxis()->SetTitle("Track POCA Y [mm]");
      hZErrorXYStation[istation]->GetZaxis()->SetTitle("<Hit global Z - Track global Z> [#mu m]");
    }

    for (Int_t ihist=0;ihist<96*12*18;ihist++){
      char title[600];
      sprintf(title, "hVErrorZWire2D%i", ihist);
      hVErrorZWire2D[ihist]= tfs->make<TH2F>(title, title, 50, -1000, 1000, 20, -500, 500);
      hVErrorZWire2D[ihist]->GetXaxis()->SetTitle("Z [mm]");
      hVErrorZWire2D[ihist]->GetYaxis()->SetTitle("Hit V - Track V [#mu m]");
      sprintf(title, "Wire %i", ihist);
      hVErrorZWire2D[ihist]->SetTitle(title);

      sprintf(title, "hVErrorZWireP%i", ihist);
      hVErrorZWireP[ihist]= tfs->make<TProfile>(title, title, 50, -1000, 1000, -1000, 1000);
      hVErrorZWireP[ihist]->GetXaxis()->SetTitle("Z [mm]");
      hVErrorZWireP[ihist]->GetYaxis()->SetTitle("Hit V - Track V [#mu m]");
      sprintf(title, "Wire %i", ihist);
      hVErrorZWireP[ihist]->SetTitle(title);

      sprintf(title, "hWErrorZWire2D%i", ihist);
      hWErrorZWire2D[ihist]= tfs->make<TH2F>(title, title, 50, -1000, 1000, 20, -500, 500);
      hWErrorZWire2D[ihist]->GetXaxis()->SetTitle("Z [mm]");
      hWErrorZWire2D[ihist]->GetYaxis()->SetTitle("Hit W - Track W [#mu m]");
      sprintf(title, "Wire %i", ihist);
      hWErrorZWire2D[ihist]->SetTitle(title);

      sprintf(title, "hWErrorZWireP%i", ihist);
      hWErrorZWireP[ihist]= tfs->make<TProfile>(title, title, 50, -1000, 1000, -1000, 1000);
      hWErrorZWireP[ihist]->GetXaxis()->SetTitle("Z [mm]");
      hWErrorZWireP[ihist]->GetYaxis()->SetTitle("Hit W - Track W [#mu m]");
      sprintf(title, "Wire %i", ihist);
      hWErrorZWireP[ihist]->SetTitle(title);
    }

  }

// Event-level analysis: loop over tracks and then straw hits.
void AlignmentTrackAnalyzer::analyze( const art::Event& event){
    Tracker const& trk = *GeomHandle<Tracker>();

    // Tracker geometry is fixed during the job, so construct the 216 panel
    // trapezoids only once and reuse them for every event/track.
    if (!_panelTrapezoidsInitialized)
      InitializePanelTrapezoids(trk);

    auto const& kalSeeds = event.getProduct(_kalSeedsToken);
    _hNTracks->Fill( kalSeeds.size() );

    Int_t counter = 0;

    // Loop over tracks in the event
    for ( auto const& ks : kalSeeds ) {
        const std::vector<mu2e::TrkStrawHitSeed> &hits = ks.hits();

      counter++;
      // Require final fit to have converged.
      if ( ! ks.status().hasAllProperties( TrkFitFlag::kalmanConverged ) ){
        _hnSkip->Fill(0.);
        continue;
      }
        
      if ((float)ks.chisquared()/ks.nDOF() > _maxChi2Dof)
        continue;

      
      float momX = 0;
      float momY = 0;
      float momZ = 0;
      float momMag = 0;
      double trackX = 0.0;
      double trackY = 0.0;
      double trackZ = 0.0;

      // Use the first KalIntersection to define a point and local direction
      // for the straight-line panel-intersection approximation.
      std::vector<KalIntersection>  Intersections = ks.intersections();
      if (Intersections.empty())
        continue;
      for(unsigned int i = 0; i < Intersections.size(); i++){
        auto intersection_i = Intersections.at(i);
        momX = intersection_i.momentum3().X();
        momY = intersection_i.momentum3().Y();
        momZ = intersection_i.momentum3().Z();
        trackX = intersection_i.position3().X();
        trackY = intersection_i.position3().Y();
        trackZ = intersection_i.position3().Z();
        break;
      }
      momMag = sqrt(momX*momX + momY*momY + momZ*momZ);
      momX = momX/momMag;
      momY = momY/momMag;
      momZ = momZ/momMag;

      std::vector<int> uniquePanelIDs;
      std::vector<int> uniquePlaneIDs;
      std::vector<int> nGoodPanelHitsIDs;

      std::vector<int> nPanelHitsIDs;
      std::vector<int> evenStrawPanelID;
      std::vector<int> oddStrawPanelID;

      std::vector<int> leftStrawPanelID;
      std::vector<int> rightStrawPanelID;

      // First pass over hits: count unique panels/planes and hit categories.
      for(unsigned int i = 0; i < hits.size(); i++){
        const mu2e::TrkStrawHitSeed &tshs = hits.at(i);
        Int_t uniquePanel = tshs.strawId().getPlane()*(*nPanelsPerPlane) + tshs.strawId().getPanel();
        Int_t unique = 1;
        for(unsigned int j = 0; j < uniquePanelIDs.size(); j++){
          if(uniquePanelIDs[j]==uniquePanel){
            unique=0;
            nPanelHitsIDs[j] = nPanelHitsIDs[j]+1;
            if (abs(tshs.strawHitState())==1 &&  abs(tshs.fitDOCA()) > _minAbsFitDOCA)
            nGoodPanelHitsIDs[j]= nGoodPanelHitsIDs[j]+1;
            if ((Int_t)tshs.strawId().getStraw()%2==0){
              evenStrawPanelID[j] = evenStrawPanelID[j] + 1;
            }
            else{
              oddStrawPanelID[j] = oddStrawPanelID[j] + 1;
            }
            if ((tshs.strawHitState() * (2*tshs.earlyEnd()-1)) ==1){
              rightStrawPanelID[j] = rightStrawPanelID[j]  + 1;
            }
            if ((tshs.strawHitState() * (2*tshs.earlyEnd()-1)) ==-1){
              leftStrawPanelID[j] = leftStrawPanelID[j]  + 1;
            }
          }
        }

        if (unique==1){
          uniquePanelIDs.push_back(uniquePanel);
          nPanelHitsIDs.push_back(1);
          if (abs(tshs.strawHitState())==1 &&  abs(tshs.fitDOCA()) > _minAbsFitDOCA)
          nGoodPanelHitsIDs.push_back(1);
          else
          nGoodPanelHitsIDs.push_back(0);
          if ((Int_t)tshs.strawId().getStraw()%2==0){
            evenStrawPanelID.push_back(1);
            oddStrawPanelID.push_back(0);
          }
          else{
            oddStrawPanelID.push_back(1);
            evenStrawPanelID.push_back(0);
          }
          if ((tshs.strawHitState() * (2*tshs.earlyEnd()-1)) ==1){
            rightStrawPanelID.push_back(1);
            leftStrawPanelID.push_back(0);
          }else if ((tshs.strawHitState() * (2*tshs.earlyEnd()-1)) ==-1){
            leftStrawPanelID.push_back(1);
            rightStrawPanelID.push_back(0);
          }else{
            leftStrawPanelID.push_back(0);
            rightStrawPanelID.push_back(0);
          }

        }

        Int_t uniquePlane = tshs.strawId().getPlane();
        Int_t unique2 = 1;
        for(unsigned int j = 0; j < uniquePlaneIDs.size(); j++){
          if(uniquePlaneIDs[j]==uniquePlane)
          unique2=0;
        }
        if (unique2==1)
        uniquePlaneIDs.push_back(uniquePlane);
      }

      Int_t nGoodHits = 0;
      Int_t nQualityPanelIntersections =0;
      Int_t nHitsOnQualityPanel0 = 0;
      Int_t nHitsOnQualityPanel1 = 0;
      for(unsigned int j = 0; j < uniquePanelIDs.size(); j++){
        nGoodHits = nGoodHits+ nGoodPanelHitsIDs[j];
        if (evenStrawPanelID[j] >= _minEvenStrawsPerPanel &&
            oddStrawPanelID[j] >= _minOddStrawsPerPanel &&
            nGoodPanelHitsIDs[j] >= _minGoodHitsPerQualityPanel){
          nQualityPanelIntersections = nQualityPanelIntersections + 1;
	  if (nHitsOnQualityPanel1==0 && nHitsOnQualityPanel0!=0)
	    nHitsOnQualityPanel1=nGoodPanelHitsIDs[j];
	  if (nHitsOnQualityPanel0==0)
            nHitsOnQualityPanel0=nGoodPanelHitsIDs[j];
	}
      }

      _hPanelID->Fill(uniquePanelIDs[0]);
        if (nGoodHits < _minGoodHits)
            continue;

        if (nQualityPanelIntersections < _minQualityPanelIntersections)
            continue;

        if (_maxQualityPanelIntersections >= 0 &&
            nQualityPanelIntersections > _maxQualityPanelIntersections)
            continue;

        // Compute the largest pairwise transverse separation among quality hits.
        // Only hits passing the straw-hit-state and fitted-DOCA requirements enter
        // this calculation.  Positions are the reconstructed hit POCAs in global XY.
        std::vector<TVector3> qualityHitPocaG;
        qualityHitPocaG.reserve(hits.size());

        for (const auto& tshs : hits) {
          if (abs(tshs.strawHitState()) != 1 ||
              abs(tshs.fitDOCA()) < _minAbsFitDOCA)
            continue;

          Straw const& straw = trk.getStraw(tshs.strawId());
          TVector3 wireEnd0(straw.wirePosition(-straw.halfLength()).x(),
                            straw.wirePosition(-straw.halfLength()).y(),
                            straw.wirePosition(-straw.halfLength()).z());
          TVector3 wireEnd1(straw.wirePosition( straw.halfLength()).x(),
                            straw.wirePosition( straw.halfLength()).y(),
                            straw.wirePosition( straw.halfLength()).z());
          TVector3 trackDir(momX, momY, momZ);

          const double hitDoca = tshs.driftRadius();
          const int leftRight = tshs.strawHitState() * (2*tshs.earlyEnd()-1);

          TVector3 hitPocaL = ComputePOCA2D(hitDoca,
                                            trackDir,
                                            leftRight,
                                            wireEnd0,
                                            wireEnd1);
          hitPocaL.SetX(tshs.wireDist());
          qualityHitPocaG.push_back(LocalToGlobal(hitPocaL, wireEnd0, wireEnd1));
        }

        double xySpan = 0.0;
        for (size_t ihit = 0; ihit < qualityHitPocaG.size(); ++ihit) {
          for (size_t jhit = ihit + 1; jhit < qualityHitPocaG.size(); ++jhit) {
            const double dx = qualityHitPocaG[ihit].X() - qualityHitPocaG[jhit].X();
            const double dy = qualityHitPocaG[ihit].Y() - qualityHitPocaG[jhit].Y();
            xySpan = std::max(xySpan, std::sqrt(dx*dx + dy*dy));
          }
        }

        double zSpan = 0.0;
        if (!qualityHitPocaG.empty()) {
          double minZ = qualityHitPocaG.front().Z();
          double maxZ = minZ;
          for (const auto& poca : qualityHitPocaG) {
            minZ = std::min(minZ, poca.Z());
            maxZ = std::max(maxZ, poca.Z());
          }
          zSpan = maxZ - minZ;
        }


        // Fill span distributions before applying their cuts.
        _hXYSpan->Fill(xySpan);
        _hZSpan->Fill(zSpan);
        if (xySpan < _minXYSpan || zSpan < _minZSpan)
          continue;
        
	_hXYSpan->Fill(xySpan);
        _hZSpan->Fill(zSpan);
         
        
        // Selection passed. Keep the original track-summary printout.
        /*
        std::cout << "Track with: " << xySpan << " XY span, " << zSpan
                  << " Z span and " << nQualityPanelIntersections
                  << " good panel intersections each with "
                  << nHitsOnQualityPanel0 << ", " << nHitsOnQualityPanel1 << std::endl;
         */
        
      _hNHits->Fill(hits.size());
      _hNPanels->Fill(uniquePanelIDs.size());
      _hNPlanes->Fill(uniquePlaneIDs.size());
        
        

      double phi   = std::atan2(momY, momX);
      double theta = std::atan2(std::sqrt(momX*momX + momY*momY), momZ);
        
      _hNHitsQuality->Fill(nGoodHits);
      _hNPanelsQuality->Fill(nQualityPanelIntersections);
      _hTrackPhi->Fill(phi);
      _hTrackTheta->Fill(theta);

    
      _hmomX->Fill(momX);
      _hmomY->Fill(momY);
      _hmomZ->Fill(momZ);
      _hchi2_dof->Fill((float)ks.chisquared()/ks.nDOF());

      _hchi2_nhits->Fill((float)ks.chisquared()/ks.nDOF(), hits.size());

      // ------------------------------------------------------------------
      // First-pass panel efficiency diagnostic.
      //
      // 1) Extrapolate the selected track as a straight line to the constant-Z
      //    plane of each panel and record whether it lands inside that panel's
      //    precomputed trapezoid.
      // 2) In a separate hit loop, record which panels actually have >= 1 hit.
      // 3) Compare the two lists in the three requested categories.
      // ------------------------------------------------------------------
      bool expectedPanel[216] = {false};
      bool hitPanel[216] = {false};
      std::vector<int> expectedPanelIDs;

      // A track with essentially zero Z direction cannot be extrapolated to
      // the panel Z planes using this simple approximation.
      if (std::abs(momZ) > 1.0e-12) {
        for (int panelIndex = 0; panelIndex < 216; ++panelIndex) {
          const double tPanel = (_panelZ[panelIndex] - trackZ) / momZ;
          const double xPanel = trackX + tPanel*momX;
          const double yPanel = trackY + tPanel*momY;

          if (PointInsidePanelTrapezoid(xPanel, yPanel, panelIndex)) {
            expectedPanel[panelIndex] = true;
            expectedPanelIDs.push_back(panelIndex);
          }
        }
      }

      // Deliberately separate from the residual-analysis loop below.  For this
      // first efficiency test, any TrkStrawHitSeed on the panel counts as a hit.
      for (const auto& tshs : hits) {
        const int panelIndex = int(tshs.strawId().getPanel() +
                                   6*tshs.strawId().getPlane());
        if (panelIndex >= 0 && panelIndex < 216)
          hitPanel[panelIndex] = true;
      }

      int nHitExpected = 0;
      int nNoHitExpected = 0;
      int nHitNoExpected = 0;

      for (int panelIndex = 0; panelIndex < 216; ++panelIndex) {
        if (expectedPanel[panelIndex]) {
          _hPanelIntersected->Fill(panelIndex);
          if (hitPanel[panelIndex])
            _hPanelIntersectedWithHit->Fill(panelIndex);
        }

        if ( expectedPanel[panelIndex] &&  hitPanel[panelIndex]) ++nHitExpected;
        if ( expectedPanel[panelIndex] && !hitPanel[panelIndex]) ++nNoHitExpected;
        if (!expectedPanel[panelIndex] &&  hitPanel[panelIndex]) ++nHitNoExpected;
      }

      // Second print line per selected track: panel-efficiency diagnostic.
        /*
      std::cout << "PanelEfficiency"
                << " expected+hit=" << nHitExpected
                << " expected+nohit=" << nNoHitExpected
                << " noexpected+hit=" << nHitNoExpected
                << " expectedPanels=[";
      for (size_t ip = 0; ip < expectedPanelIDs.size(); ++ip) {
        if (ip != 0) std::cout << ",";
        std::cout << expectedPanelIDs[ip];
      }
      std::cout << "]" << std::endl;
         */
      // Second pass over hits: compute POCA residuals and fill diagnostics.
      for(unsigned int i = 0; i < hits.size(); i++){
          
        const mu2e::TrkStrawHitSeed &tshs = hits.at(i);

        // Fill timing diagnostics for every hit on a selected track, before
        // applying the individual hit-state / DOCA quality requirements.
        _hdriftTime->Fill(tshs.time() - tshs._ptoca);
        _hdriftTimePanel[int(tshs.strawId().getPanel() + 6 * tshs.strawId().getPlane())]->Fill(tshs.time() - tshs._ptoca);
        _hTD->Fill(tshs.time() - tshs._ptoca, tshs.driftRadius());
          
        if (abs(tshs.strawHitState())!=1 ||  abs(tshs.fitDOCA()) < _minAbsFitDOCA)
            continue;
          
        Straw const& straw = trk.getStraw(tshs.strawId());

        const int panelIndex = int(tshs.strawId().getPanel() + 6 * tshs.strawId().getPlane());
        //std::cout << "Panel " << panelIndex
        //          << " tshs._edep * 1000 = " << tshs._edep * 1000.0 << std::endl;
        _hEDepPanel[panelIndex]->Fill(tshs._edep * 1000.0);

        // Longitudinal and drift-radius diagnostics for quality-selected hits.
        _hTrackLong->Fill(tshs._uupos);
        _hHitLong->Fill(tshs.wireDist());
        _hLongError->Fill(tshs.wireDist() - tshs._uupos);

        _hUDOCA->Fill(abs(tshs.fitDOCA()));
        _hRDrift->Fill(tshs._rdrift);
        _hCDrift->Fill(tshs._cdrift);
	//std::cout << tshs.fitDOCA() << ", " << sqrt(tshs._udocavar) << std::endl;
        // Compare drift radius to the fitted unsigned DOCA.
	_hTrackSigma->Fill(sqrt(tshs._udocavar)*1000);
	
        const double trackUDOCA = abs(tshs.fitDOCA());
        const double docaResidualUm = (tshs.driftRadius() - trackUDOCA) * 1000.0;
        const double docaResidualSquaredUm2 = docaResidualUm * docaResidualUm;

	_hDocaErrorUDOCA->Fill(trackUDOCA, docaResidualUm);
        _hDocaErrorUDOCAPanel[panelIndex]->Fill(trackUDOCA, docaResidualUm);
        _hDocaErrorSquaredUDOCAPanel[panelIndex]->Fill(trackUDOCA, docaResidualSquaredUm2);
        _hUDOCAEntriesPanel[panelIndex]->Fill(trackUDOCA);
        _hDocaErrorSquaredUPosPanel[panelIndex]->Fill(tshs._uupos, docaResidualSquaredUm2);
        _hUPosEntriesPanel[panelIndex]->Fill(tshs._uupos);
        _hDocaErrorRDrift->Fill(abs(tshs._rdrift), (tshs.driftRadius() - abs(tshs.fitDOCA())) * 1000);
        _hDocaErrorCDrift->Fill(abs(tshs._cdrift), (tshs.driftRadius() - abs(tshs.fitDOCA())) * 1000);

        _hDOCAError->Fill((tshs.driftRadius() - abs(tshs.fitDOCA())) * 1000);
        _hDOCAErrorPanel[int(tshs.strawId().getPanel() + 6 * tshs.strawId().getPlane())]->Fill((tshs.driftRadius() - abs(tshs.fitDOCA())) * 1000);

        TVector3 wireEnd0(straw.wirePosition(-straw.halfLength()).x(), straw.wirePosition(-straw.halfLength()).y(), straw.wirePosition(-straw.halfLength()).z());
        TVector3 wireEnd1(straw.wirePosition(straw.halfLength()).x(), straw.wirePosition(straw.halfLength()).y(), straw.wirePosition(straw.halfLength()).z());
        TVector3 trackDir(momX, momY, momZ);

        double trackDoca = abs(tshs.fitDOCA());
        double hitDoca = tshs.driftRadius();
        int leftRight = tshs.strawHitState() * (2*tshs.earlyEnd()-1); // +1 or -1

        TVector3 trackPocaL =ComputePOCA2D(trackDoca,
        trackDir,
        leftRight,
        wireEnd0,
        wireEnd1);

        TVector3 hitPocaL   = ComputePOCA2D(hitDoca,
        trackDir,
        leftRight,
        wireEnd0,
        wireEnd1);

        hitPocaL.SetX(tshs.wireDist());
        trackPocaL.SetX(tshs.wireDist());

        TVector3 trackPocaG = LocalToGlobal(trackPocaL, wireEnd0, wireEnd1);
        TVector3 hitPocaG = LocalToGlobal(hitPocaL, wireEnd0, wireEnd1);

        hXYPOCA->Fill(trackPocaG.X(), trackPocaG.Y());
        hYZPOCA->Fill(trackPocaG.Y(), trackPocaG.Z());

        _hZPOCAErrorY->Fill(straw.strawPosition(tshs.wireDist()).y(),(hitPocaL.Z() -trackPocaL.Z())*1000);
        _hXPOCAError->Fill((hitPocaL.X() -trackPocaL.X())*1000);
        _hYPOCAError->Fill((hitPocaL.Y() -trackPocaL.Y())*1000);
        _hZPOCAError->Fill((hitPocaL.Z() -trackPocaL.Z())*1000);

        _hXPOCAErrorRatio->Fill((hitPocaL.X() -trackPocaL.X())/abs(tshs.driftRadius() -  abs(tshs.fitDOCA())));
        _hYPOCAErrorRatio->Fill((hitPocaL.Y() -trackPocaL.Y())/abs(tshs.driftRadius() -  abs(tshs.fitDOCA())));
        _hZPOCAErrorRatio->Fill((hitPocaL.Z() -trackPocaL.Z())/abs(tshs.driftRadius() -  abs(tshs.fitDOCA())));
        _hZPOCAError->Fill((hitPocaL.Z() -trackPocaL.Z())*1000);

        hZPOCAErrorWire[int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())]->Fill((Int_t)tshs.strawId().getStraw(), (hitPocaL.Z() -trackPocaL.Z())*1000.);
        if ((Int_t)tshs.strawId().getStraw()%2==0)
        hZPOCAErrorWireEven[int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())]->Fill((Int_t)tshs.strawId().getStraw(), (hitPocaL.Z() -trackPocaL.Z())*1000.);
        else
        hZPOCAErrorWireOdd[int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())]->Fill((Int_t)tshs.strawId().getStraw(), (hitPocaL.Z() -trackPocaL.Z())*1000.);

        if (momZ > 0){
          if ((Int_t)tshs.strawId().getStraw()%2==0)
          hZPOCAErrorWireEven_PosZ[int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())]->Fill((Int_t)tshs.strawId().getStraw(), (hitPocaL.Z() -trackPocaL.Z())*1000.);
          else
          hZPOCAErrorWireOdd_PosZ[int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())]->Fill((Int_t)tshs.strawId().getStraw(), (hitPocaL.Z() -trackPocaL.Z())*1000.);
        }else{
          if ((Int_t)tshs.strawId().getStraw()%2==0)
          hZPOCAErrorWireEven_NegZ[int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())]->Fill((Int_t)tshs.strawId().getStraw(), (hitPocaL.Z() -trackPocaL.Z())*1000.);
          else
          hZPOCAErrorWireOdd_NegZ[int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())]->Fill((Int_t)tshs.strawId().getStraw(), (hitPocaL.Z() -trackPocaL.Z())*1000.);
        }

        hVPOCAErrorWire[int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())]->Fill((Int_t)tshs.strawId().getStraw(), (hitPocaL.Y() -trackPocaL.Y())*1000.);
        if ((Int_t)tshs.strawId().getStraw()%2==0)
        hVPOCAErrorWireEven[int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())]->Fill((Int_t)tshs.strawId().getStraw(), (hitPocaL.Y() -trackPocaL.Y())*1000.);
        else
        hVPOCAErrorWireOdd[int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())]->Fill((Int_t)tshs.strawId().getStraw(), (hitPocaL.Y() -trackPocaL.Y())*1000.);

        if (momZ > 0){
          if ((Int_t)tshs.strawId().getStraw()%2==0)
          hVPOCAErrorWireEven_PosZ[int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())]->Fill((Int_t)tshs.strawId().getStraw(), (hitPocaL.Y() -trackPocaL.Y())*1000.);
          else
          hVPOCAErrorWireOdd_PosZ[int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())]->Fill((Int_t)tshs.strawId().getStraw(), (hitPocaL.Y() -trackPocaL.Y())*1000.);
        }else{
          if ((Int_t)tshs.strawId().getStraw()%2==0)
          hVPOCAErrorWireEven_NegZ[int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())]->Fill((Int_t)tshs.strawId().getStraw(), (hitPocaL.Y() -trackPocaL.Y())*1000.);
          else
          hVPOCAErrorWireOdd_NegZ[int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())]->Fill((Int_t)tshs.strawId().getStraw(), (hitPocaL.Y() -trackPocaL.Y())*1000.);
        }

        //const int panelIndex = int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane());
        const double xErrorGlobal = (hitPocaG.X() - trackPocaG.X()) * 1000.;
        const double yErrorGlobal = (hitPocaG.Y() - trackPocaG.Y()) * 1000.;
        const double vErrorLocal  = (hitPocaL.Y() - trackPocaL.Y()) * 1000.;

        // Existing local-coordinate panel residuals.
        _hXPOCAErrorPanel[panelIndex]->Fill((hitPocaL.X() - trackPocaL.X()) * 1000.);
        _hYPOCAErrorPanel[panelIndex]->Fill(vErrorLocal);
        // Z sensitivity is the magnitude of the global-Z component of the
        // unit DOCA direction: (wireDir x trackDir).Unit().  This directly
        // measures how strongly the DOCA responds to a global-Z displacement.
        const TVector3 wireDir = (wireEnd1 - wireEnd0).Unit();
        const TVector3 docaDir = wireDir.Cross(trackDir).Unit();
        const double zSensitivity = std::abs(docaDir.Z());
        _hZSensitivityPanel[panelIndex]->Fill(zSensitivity);

        // Keep only hits with appreciable global-Z sensitivity.
        if (zSensitivity > 0.4) {
          const double zErrorUm = (hitPocaG.Z() - trackPocaG.Z()) * 1000.;
          _hZPOCAErrorPanel[panelIndex]->Fill(zErrorUm);

          const int stationIndex = int(tshs.strawId().getPlane()) / 2;
          hZErrorVsXStation[stationIndex]->Fill(trackPocaG.X(), zErrorUm);
          hZErrorVsYStation[stationIndex]->Fill(trackPocaG.Y(), zErrorUm);
          hZErrorXYStation[stationIndex]->Fill(trackPocaG.X(), trackPocaG.Y(), zErrorUm);
        }

        // Full-panel residual distributions requested in local V and global X/Y.
        hVErrorPanel[panelIndex]->Fill(vErrorLocal);
        hXErrorPanel[panelIndex]->Fill(xErrorGlobal);
        hYErrorPanel[panelIndex]->Fill(yErrorGlobal);

        // Mean residual versus tshs._uupos for each panel.
        hXErrorZPanelP[panelIndex]->Fill(tshs._uupos, xErrorGlobal);
        hYErrorZPanelP[panelIndex]->Fill(tshs._uupos, yErrorGlobal);
        hVErrorZPanelP[panelIndex]->Fill(tshs._uupos, vErrorLocal);

        hVErrorZWireP[int(tshs.strawId().getStraw())+int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())*96]->Fill(tshs._uupos,
        (hitPocaL.Y() -trackPocaL.Y())*1000);
        hWErrorZWireP[int(tshs.strawId().getStraw())+int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())*96]->Fill(tshs._uupos,
        (hitPocaL.Z() -trackPocaL.Z())*1000);
        hVErrorZWire2D[int(tshs.strawId().getStraw())+int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())*96]->Fill(tshs._uupos,
        (hitPocaL.Y() -trackPocaL.Y())*1000);
        hWErrorZWire2D[int(tshs.strawId().getStraw())+int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())*96]->Fill(tshs._uupos,
        (hitPocaL.Z() -trackPocaL.Z())*1000);

        hWire[int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())]->SetBinContent(int(tshs.strawId().getStraw()),
        hWire[int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())]->GetBinContent(int(tshs.strawId().getStraw()))+1);

        hULocal[int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())]->Fill(trackPocaL.X());

        hdVULocal[int(tshs.strawId().getPanel() + 6*tshs.strawId().getPlane())]->Fill(trackPocaL.X(), (hitPocaL.Y() -trackPocaL.Y())*1000);

        _hDOCAErrorCleaned->Fill((tshs.driftRadius() - abs(tshs.fitDOCA())) * 1000);

      }

    } // end loop over kalSeeds

  } // end analyze

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::AlignmentTrackAnalyzer)
