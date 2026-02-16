const int Nboard = 4;
const int Nhv = 3;

void convert_kek(){
  TFile *file_ped[Nhv];
  TFile *file_led[Nboard][Nhv];

  TTree *tree_ped[Nhv];
  TTree *tree_led[Nboard][Nhv];

  Double_t ADC_ped[Nhv][4];
  Double_t ADC_led[Nboard][Nhv][4];
  
  int ped_run[Nhv] = {243,242,241};

  int led_run[Nboard][Nhv] = {{239,238,240},{239,238,240},{239,238,240},{239,238,240}};

  TString dir = "/Users/ihaein/Work/E72/KEKAR2023Dec/KEKAR2023Dec_root";
  for(int j=0;j<Nhv;j++){
    file_ped[j] = new TFile(Form("%s/kekar_run00%d.root",dir.Data(),ped_run[j]));
    tree_ped[j] = (TTree*)file_ped[j]->Get("tree");
    tree_ped[j]->SetBranchAddress("baca",ADC_ped[j]);
    for(int i=0;i<Nboard;i++){
      file_led[i][j] = TFile::Open(Form("%s/kekar_run00%d.root",dir.Data(),led_run[i][j]), "READ");
      tree_led[i][j] = (TTree*)file_led[i][j]->Get("tree");
      tree_led[i][j]->SetBranchAddress("baca",ADC_led[i][j]);
    }
  }
  
  TFile *file = new TFile("BAC_LED_KEK.root","recreate");
  TNamed *desc = new TNamed(
  "Description",
  "BAC LED calibration run\n"
  "\n"
  "[Board][HV]\n"
  "HV settings (index order):\n"
  "  HV[0] = 56 V\n"
  "  HV[1] = 57 V\n"
  "  HV[2] = 58 V\n"
  "\n"
  "Gate width : 130 ns\n"
  "\n"
  "Channel assignment (0-origin):\n"
  "  All boards use channel = 1"
			    );
  desc->Write();
  TTree *data = new TTree("tree","tree");

  Double_t ADCp[Nboard][Nhv];
  Double_t ADCl[Nboard][Nhv];

  
  data->Branch("ADC_ped",&ADCp,Form("ADC_ped[%d][%d]/D",Nboard,Nhv));
  data->Branch("ADC_led",&ADCl,Form("ADC_led[%d][%d]/D",Nboard,Nhv));
  int total = 100000;
  for(int n=0;n<total;n++){
    for(int i=0;i<Nboard;i++){
      for(int j=0;j<Nhv;j++){
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
