void show_result(const char* file="PHM_out/test.root") {
  TFile *f = TFile::Open(file);
  gStyle->SetOptStat(0);
  new TBrowser("Browser", "PixelHistoMaker output", 1200, 800);
}
