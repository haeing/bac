#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>

struct WStat {
  double sumw  = 0.0;
  double sumwx = 0.0;

  void Fill(double x, double err){
    if (!std::isfinite(x) || !std::isfinite(err)) return;
    if (err <= 0) return;
    const double w = 1.0/(err*err);
    sumw  += w;
    sumwx += w*x;
  }

  bool Has() const { return sumw > 0; }
  double Mean() const { return sumwx/sumw; }
  double Err()  const { return std::sqrt(1.0/sumw); }
};

struct Key {
  int board;
  int ch;
  bool operator<(const Key& o) const {
    if (board != o.board) return board < o.board;
    return ch < o.ch;
  }
};

static std::vector<TString> SplitCSV(const TString& line){
  std::vector<TString> out;
  TObjArray* arr = line.Tokenize(",");
  for (int i=0;i<arr->GetEntriesFast();++i)
    out.push_back(((TObjString*)arr->At(i))->GetString().Strip(TString::kBoth));
  arr->Delete(); delete arr;
  return out;
}

void analysis_t110()
{

  const char* csv = "led_gain/spe_result_update/spe_summary.csv";
  double target_mppc_hv = 58.0;
  double hv_tol = 1e-6;
  std::ifstream fin(csv);
  if(!fin.is_open()){
    std::cerr << "[ERR] cannot open " << csv << std::endl;
    return;
  }


  std::map<Key, WStat> q1stat;
  std::map<Key, WStat> pxtstat;

  std::string line;

  std::getline(fin,line);
  while (std::getline(fin, line)){
    TString t(line.c_str());
    if (t.Strip(TString::kBoth).IsNull()) continue;

    auto tok = SplitCSV(t);
    //if ((int)tok.size() <= 13)continue;
    
    int board = tok[1].Atoi();
    int ch = tok[2].Atoi();
    double mhv = tok[5].Atoi();

    double q1 = tok[9].Atof();
    double q1e = tok[10].Atof();
    double pxt = tok[11].Atof();
    double pxte = tok[12].Atof();

    if(ch == -1)continue;
    Key k{board, ch};
    q1stat[k].Fill(q1, q1e);
    pxtstat[k].Fill(pxt, pxte);
  }
  fin.close();


  std::vector<int>    v_board, v_ch;
  std::vector<double> v_q1, v_q1err, v_pxt, v_pxterr;

  for (const auto& it : q1stat){
    const Key& k = it.first;
    if (!q1stat[k].Has() || !pxtstat[k].Has()) continue;

    v_board.push_back(k.board);
    v_ch.push_back(k.ch);

    v_q1.push_back(q1stat[k].Mean());
    v_q1err.push_back(q1stat[k].Err());

    v_pxt.push_back(pxtstat[k].Mean());
    v_pxterr.push_back(pxtstat[k].Err());
  }

  std::cout << "board  ch   <Q1> ± err     <PXT> ± err\n";
  for (size_t i=0;i<v_board.size();++i){
    std::cout
      << v_board[i] << "  "
      << v_ch[i] << "   "
      << v_q1[i] << " ± " << v_q1err[i] << "   "
      << v_pxt[i] << " ± " << v_pxterr[i] << "\n";
  }
  
}

