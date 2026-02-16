const int Nboard = 4;
const int Nmppc = 16;

void convert(){
  TFile *file_ped[Nmppc];
  TFile *file_led[Nboard][Nmppc];

  TTree *tree_ped[Nmppc];
  TTree *tree_led[Nboard][Nmppc];

  Double_t ADC_ped[Nmppc][4];
  Double_t ADC_led[Nboard][Nmppc][4];
  
  int ped_run[Nmppc] = {137, 139, 143, 145, 148, 151, 155, 159, 161, 163, 165, 168, 172, 175, 177, 179};

  int led_run[Nboard][Nmppc] = {{138, 142, 144, 147, 150, 154, 158, 160, 162, 164, 167, 169, 174, 176, 178, 180},{138, 142, 144, 147, 150, 154, 158, 160, 162, 164, 167, 169, 174, 176, 178, 180},{138, 142, 144, 147, 150, 153, 158, 160, 162, 164, 167, 171, 174, 176, 178, 180},{138, 142, 144, 147, 150, 153, 158, 160, 162, 164, 167, 169, 174, 176, 178, 180}};

  TString dir = "../KEKAR2023Dec_root";
  for(int j=0;j<Nmppc;j++){
    file_ped[j] = new TFile(Form("%s/kekar_run00%d.root",dir.Data(),ped_run[j]));
    tree_ped[j] = (TTree*)file_ped[j]->Get("tree");
    tree_ped[j]->SetBranchAddress("baca",ADC_ped[j]);
    for(int i=0;i<Nboard;i++){
      file_led[i][j] = TFile::Open(Form("%s/kekar_run00%d.root",dir.Data(),led_run[i][j]), "READ");
      //file_led[i][j] = new TFile(Form("kekar_run00%d.root",led_run[i][j]), "READ");
      tree_led[i][j] = (TTree*)file_led[i][j]->Get("tree");
      tree_led[i][j]->SetBranchAddress("baca",ADC_led[i][j]);
    }
  }
  
  TFile *file = new TFile("BAC_LED_57V.root","recreate");
  TTree *data = new TTree("tree","LED w/ 57 V");

  Double_t ADCp[Nboard][Nmppc];
  Double_t ADCl[Nboard][Nmppc];

  
  data->Branch("ADC_ped",&ADCp,Form("ADC_ped[%d][%d]/D",Nboard,Nmppc));
  data->Branch("ADC_led",&ADCl,Form("ADC_led[%d][%d]/D",Nboard,Nmppc));
  int total = 100000;
  for(int n=0;n<total;n++){
    for(int i=0;i<Nboard;i++){
      for(int j=0;j<Nmppc;j++){
	tree_ped[j]->GetEntry(n);
	tree_led[i][j]->GetEntry(n);

	ADCp[i][j] = ADC_ped[j][i];
	ADCl[i][j] = ADC_led[i][j][i];
      }
    }
    data->Fill();
  }

  data->Write();
  file->Close();
  
				

}
