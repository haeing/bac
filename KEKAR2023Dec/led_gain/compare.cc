#include <array>
#include <cctype>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

constexpr int NBOARD = 4;
constexpr int NCH    = 16;
constexpr int NHV    = 3;   // 56,57,58

struct ValErr {
  double v  = std::numeric_limits<double>::quiet_NaN();
  double e  = std::numeric_limits<double>::quiet_NaN();
  bool ok   = false;
};

static std::string trim(std::string s) {
  auto issp = [](unsigned char c){ return std::isspace(c); };
  while (!s.empty() && issp((unsigned char)s.front())) s.erase(s.begin());
  while (!s.empty() && issp((unsigned char)s.back()))  s.pop_back();
  return s;
}

static std::vector<std::string> split_csv_line(const std::string& line) {
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) out.push_back(trim(item));
  return out;
}

static bool looks_like_header(const std::string& line) {
  auto cols = split_csv_line(line);
  if (cols.empty()) return true;
  if (cols[0].empty()) return true;
  char c = cols[0][0];
  return !(std::isdigit((unsigned char)c) || c=='-');
}

static int hv_to_index(int hv) {
  if (hv == 56) return 0;
  if (hv == 57) return 1;
  if (hv == 58) return 2;
  return -1;
}

bool LoadBaseline57V(const std::string& csvPath,
                     std::array<std::array<ValErr, NCH>, NBOARD>& baselineQ1,
                     std::array<std::array<ValErr, NCH>, NBOARD>& baselinePxt)
{
  std::ifstream fin(csvPath);
  if (!fin) {
    std::cerr << "[ERROR] Cannot open baseline csv: " << csvPath << "\n";
    return false;
  }

  std::string line;
  bool first = true;

  while (std::getline(fin, line)) {
    line = trim(line);
    if (line.empty()) continue;

    if (first && looks_like_header(line)) { first = false; continue; }
    first = false;

    auto c = split_csv_line(line);
    // 기대: 0:board 1:channel 2:Q1 3:Q1_err 4:pxt 5:pxt_err ...
    if (c.size() < 6) continue;

    int board   = std::stoi(c[0]);
    int channel = std::stoi(c[1]);
    if (board < 0 || board >= NBOARD || channel < 0 || channel >= NCH) continue;

    double Q1     = std::stod(c[2]);
    double Q1_err = std::stod(c[3]);
    double pxt    = std::stod(c[4]);
    double pxt_err= std::stod(c[5]);

    baselineQ1[board][channel]  = {Q1,  Q1_err,  true};
    baselinePxt[board][channel] = {pxt, pxt_err, true};
  }

  return true;
}

// -------------------------
// KEK HV-scan CSV 로드
// (board,channel,HV,Q1,Q1_err,pxt,pxt_err,...)
// channel=1만 사용한다고 가정
// -------------------------
bool LoadKekHVScan(const std::string& csvPath,
                   std::array<std::array<ValErr, NHV>, NBOARD>& kekQ1,
                   std::array<std::array<ValErr, NHV>, NBOARD>& kekPxt,
                   int targetChannel = 1)
{
  std::ifstream fin(csvPath);
  if (!fin) {
    std::cerr << "[ERROR] Cannot open kek csv: " << csvPath << "\n";
    return false;
  }

  std::string line;
  bool first = true;

  while (std::getline(fin, line)) {
    line = trim(line);
    if (line.empty()) continue;

    if (first && looks_like_header(line)) { first = false; continue; }
    first = false;

    auto c = split_csv_line(line);
    // 기대: 0:board 1:channel 2:HV 3:Q1 4:Q1_err 5:pxt 6:pxt_err ...
    if (c.size() < 7) continue;

    int board   = std::stoi(c[0]);
    int channel = std::stoi(c[1]);
    int HV      = std::stoi(c[2]);
    if (board < 0 || board >= NBOARD) continue;
    if (channel != targetChannel) continue;

    int ihv = hv_to_index(HV);
    if (ihv < 0) continue;

    double Q1     = std::stod(c[3]);
    double Q1_err = std::stod(c[4]);
    double pxt    = std::stod(c[5]);
    double pxt_err= std::stod(c[6]);

    kekQ1[board][ihv]  = {Q1,  Q1_err,  true};
    kekPxt[board][ihv] = {pxt, pxt_err, true};
  }

  return true;
}

// -------------------------
// 사용 예시 main
// -------------------------
void compare()
{
  const std::string baselineCsv = "spe_result/Q1_pxt_all.csv";
  const std::string kekCsv      = "spe_result/Q1_pxt_all_kek.csv";

  std::array<std::array<ValErr, NCH>, NBOARD> baselineQ1{}, baselinePxt{};
  std::array<std::array<ValErr, NHV>, NBOARD> kekQ1{}, kekPxt{};

  if (!LoadBaseline57V(baselineCsv, baselineQ1, baselinePxt)) cout<<"error"<<endl;
  if (!LoadKekHVScan(kekCsv, kekQ1, kekPxt, /*targetChannel=*/1)) cout<<"error"<<endl;

  double q1_factor[NBOARD];
  double pxt_factor[NBOARD];

  for(int i=0;i<NBOARD;i++){
    q1_factor[i] = kekQ1[i][2].v / baselineQ1[i][1].v;
    pxt_factor[i] = kekPxt[i][2].v / baselinePxt[i][1].v;
    cout<<"q1 factor : "<<q1_factor[i]<<", pxt factor : "<<pxt_factor[i]<<endl;
  }
  
  
  
  double q1_ave[NBOARD] = {0.};
  double pxt_ave[NBOARD] = {0.};
  for(int i=0;i<NBOARD;i++){
    for(int j=0;j<NCH;j++){
      q1_ave[i] +=baselineQ1[i][j].v;
      pxt_ave[i] +=baselinePxt[i][j].v;
    }
    q1_ave[i] /= (double)NCH;
    pxt_ave[i] /= (double)NCH;
    
    cout<<"Board "<<i<<" : Q1 = "<<q1_ave[i] * q1_factor[i]<<", Pxt = "<<pxt_ave[i] * pxt_factor[i]<<endl;
  }

  
  
  

 }
