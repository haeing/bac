void test(){
  
  TFile *beam70_file = new TFile("pattern_file.root", "read");
  TH1I *hit_pattern;
  TH1D *hist = new TH1D("hist","hist",100,-60,60);
  hit_pattern = (TH1I*)beam70_file->Get("hist_bh2_pattern");

  for(int i=0;i<1000;i++){
    double x =-1*( (hit_pattern->GetRandom()-1)/8*112-56);


  hist->Fill(x);
  }
  hist->Draw();
}
