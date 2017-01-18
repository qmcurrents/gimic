#if 0

  Convert output of Q-Chem or FermiONs++ to TURBOMOLE
  format

  This Source Code Form is subject to the terms of the
  Mozilla Public License, v. 2.0. If a copy of the MPL
  was not distributed with this file, You can obtain
  one at https://mozilla.org/MPL/2.0/.

  (C) J. Kussmann

#endif

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <sstream>
#include <algorithm>

#define ANG2BOHR 1.88972612456506198632e0

#define die(x) dieLoc(x,__FILE__,__LINE__)

/*  TURBOMOLE ordering...

  D: 200 020 002 110 101 011
  F: 300 030 003 210 201 120 021 102 012 111
  G: 400 040 004 310 301 130 031 103 013 220
     202 022 211 121 112
  H: 500 050 005 410 401 140 041 104 014 320
     302 230 032 203 023 311 131 113 221 212
     122

*/

// from a fermions-perspective...
static int tm_d[6] = {0,3,4,1,5,2};
static int tm_f[10] = {0,3,4,5,9,7,1,6,8,2};
static int tm_g[15] = {0,3,4,9,12,10,5,13,14,7,1,6,11,8,2};

void dieLoc(const char* bye, const char* fff, int ll){
  printf("\nFatal error:\n\n");
  printf("  %s\n\n",bye);
  printf("  in file '%s', line %d.\n\n",fff,ll);
  exit(1);
}

int lxyz2pos (int LX, int LY, int LZ){

  int ltot = LX+LY+LZ;

  int ret = LZ;
  for(int o=0;o<=LY+LZ;o++){
    ret += o;
  }
  return ret;
}

class atom{
  public:
    int ioz;
    double x;
    double y;
    double z;
    int nbf;
    int nbfc;
    int offp;
    int offc;
};

class atomshell{
  public:
    int lqn;
    std::vector<double> linc;
    std::vector<double> expc;
};

class atombasis{
  public:
    int ioz;
    int max_lqn;
    std::vector<atomshell> shells;
};

typedef enum{
  qchem = 0,
  fermions = 1
} prog_type;


void print_mat(int nrow, int ncol, double* scr, const char* name);
void read_output_qchem(std::string& fname, std::string& basfile, std::vector<atom>& theAtoms, std::vector<atombasis>& theBasis, int& nbf, int& nbfc);
void read_output_fermions(std::string& fname, std::string& basfile, std::vector<atom>& theAtoms, std::vector<atombasis>& theBasis, int& nbf, int& nbfc);
void read_density_qchem(std::string& dens_file, std::string& xdens_file, std::vector<atom>& theAtoms, std::vector<atombasis>& theBasis, int nbf, int nbfc, bool s2c, bool opensh);
void read_density_fermions(std::string& scr_dir, std::vector<atom>& theAtoms, std::vector<atombasis>& theBasis, int nbf, int nbfc, bool s2c, bool opensh);
void conv_s2c(prog_type typ, std::vector<double>& matrix, std::vector<atom>& theAtoms, std::vector<atombasis>& theBasis, int nbf, int nbfc);
void nonaxc(prog_type typ, std::vector<double>& matrix, std::vector<atom>& theAtoms, std::vector<atombasis>& theBasis, int nbf, int nbfc);
void qchem_reorder(std::vector<double>& matrix, std::vector<atom>& theAtoms, std::vector<atombasis>& theBasis, int nbf, int nbfc);
void fermions2turbo(std::vector<double>& matrix, std::vector<atom>& theAtoms, std::vector<atombasis>& theBasis, int nbf, int nbfc);
void build_c2p_fermions(double*** c2p_tm, int max_lqn);
void form_nax_fermions(double*** nax, int max_lqn);
void build_c2p_qchem(double*** c2p_tm, int max_lqn);
void form_nax_qchem(double*** nax, int max_lqn);

int get_lqn(const std::string& clqn){

  int ret = -1;
  if (toupper(clqn.c_str()[0]) == 'S')
    ret = 0;
  else if (toupper(clqn.c_str()[0]) == 'P')
    ret = 1;
  else if (toupper(clqn.c_str()[0]) == 'D')
    ret = 2;
  else if (toupper(clqn.c_str()[0]) == 'F')
    ret = 3;
  else if (toupper(clqn.c_str()[0]) == 'G')
    ret = 4;
  else if (toupper(clqn.c_str()[0]) == 'H')
    ret = 5;
  else if (toupper(clqn.c_str()[0]) == 'I')
    ret = 6;
  else if (toupper(clqn.c_str()[0]) == 'J')
    ret = 7;
  else if (toupper(clqn.c_str()[0]) == 'K')
    ret = 8;
  else if (toupper(clqn.c_str()[0]) == 'L')
    ret = 9;
  else{
   printf("lqn: '%d'\n",clqn.c_str()[0]);
   die("lqn too high!");
  }
  return ret;

}

void get_coz(int ioz, std::string& out){

  if (ioz > 111){
    printf("Can't find COZ for >>>%i<<<!\n",ioz);
    die("Error in get_coz...");
  }

  char* a1 = new char[112];
  char* a2 = new char[112];
  strncpy(a1," HHLBBCNOFNNMASPSCAKCSTVCMFCNCZGGASBKRSYZNMTRRPACISSTIXCBLCPNPSEGTDHETYLHTWROIPAHTPBPAR FRATPUNPACBCEFMNLRDSBHM",112);
  strncpy(a2,"  EIE     EAGLI  LR ACI RNEOIUNAESERRBR RBOCUHDGDNNBE ESAAERDMMUDBYORMBUFA ESRTUGLBIOTN RACHA PUMMKFSMDORFBGHST",112);

  out = a1[ioz];
  if (a2[ioz] != ' ') out += a2[ioz];

  delete[] a1;
  delete[] a2;

}

int get_ioz(const std::string& oz){

  char* a1 = new char[113];
  char* a2 = new char[113];
  strncpy(a1," HHLBBCNOFNNMASPSCAKCSTVCMFCNCZGGASBKRSYZNMTRRPACISSTIXCBLCPNPSEGTDHETYLHTWROIPAHTPBPAR FRATPUNPACBCEFMNLRDSBHM",112);
  strncpy(a2,"  EIE     EAGLI  LR ACI RNEOIUNAESERRBR RBOCUHDGDNNBE ESAAERDMMUDBYORMBUFA ESRTUGLBIOTN RACHA PUMMKFSMDORFBGHST",112);

  int ret = -1;

  if (isalpha(oz.c_str()[1])){
    for(int i=1;i<112;i++){
      if (a1[i] == toupper(oz.c_str()[0]) && a2[i] == toupper(oz.c_str()[1])){
        ret = i;
        break;
      }
    }
  }else{
    for(int i=1;i<112;i++){
      if (a1[i] == toupper(oz.c_str()[0]) && !isalpha(a2[i])){
        ret = i;
        break;
      }
    }
  }

  delete[] a1;
  delete[] a2;

  if (ret == -1){
    printf("Can't find OZ for >>>%s<<<!\n",oz.c_str());
    die("Error in get_ioz...");
  }

  return ret;

}

// http://stackoverflow.com/a/25829233
// trim from left
inline std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(0, s.find_first_not_of(t));
    return s;
}

// trim from right
inline std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

// trim from left & right
inline std::string& trim(std::string& s, const char* t = " \t\n\r\f\v")
{
    return ltrim(rtrim(s, t), t);
}

int main(int argc, char* argv[]){

  std::string qcout_name;
  std::string dens_file;
  std::string xdens_file;
  std::string scr_dir; // fermions
  prog_type type = qchem;

  std::string basfile;
  bool do_s2c = false;
  bool opensh = false;

  int ii=1;
  while(ii < argc){
    if (strncmp(argv[ii],"-h",2) == 0){
      printf("USAGE <cmdname> -t <qchem or fermions> -qcout <output-file> -scr <scratch-directory> -s2c (opt.) -openshell (opt.)\n");
      exit(0);
    }else if (strncmp(argv[ii],"-t",2) == 0){
      if (ii == argc-1) die("Illegal arguments!");
      ii++;
      if (strncmp(argv[ii],"qchem",5) == 0) type = qchem;
      else if (strncmp(argv[ii],"fermions",8) == 0) type = fermions;
      else die("Illegal argument for '-t'!");
    }else if (strncmp(argv[ii],"-qcout",6) == 0){
      if (ii == argc-1) die("Illegal arguments!");
      ii++;
      qcout_name = argv[ii];
    }else if (strncmp(argv[ii],"-dens",5) == 0){
      if (ii == argc-1) die("Illegal arguments!");
      ii++;
      dens_file = argv[ii];
    }else if (strncmp(argv[ii],"-xdens",6) == 0){
      if (ii == argc-1) die("Illegal arguments!");
      ii++;
      xdens_file = argv[ii];
    }else if (strncmp(argv[ii],"-bas",4) == 0){
      if (ii == argc-1) die("Illegal arguments!");
      ii++;
      basfile = argv[ii];
    }else if (strncmp(argv[ii],"-scr",4) == 0){
      if (ii == argc-1) die("Illegal arguments!");
      ii++;
      scr_dir = argv[ii];
    }else if (strncmp(argv[ii],"-s2c",4) == 0){
      do_s2c = true;
    }else if (strncmp(argv[ii],"-openshell",10) == 0){
      opensh = true;
    }else{
      die("Illegal arguments!");
    }
    ii++;
  }

  if (qcout_name.size() == 0) die("Need output-file!");
  if (type == qchem){
    if (dens_file.size() == 0)  die("Need density-file!");
    if (xdens_file.size() == 0) die("Need perturbed density-file!");
  }else if (type == fermions){
    if (scr_dir.size() == 0)  die("Need scratch-directory!");
  }
  std::vector<atom> theAtoms;
  std::vector<atombasis> theBasis;

  int nbf  = -1;
  int nbfc = -1;
  if (type == qchem){
    read_output_qchem(qcout_name,basfile,theAtoms,theBasis,nbf,nbfc);
    read_density_qchem(dens_file,xdens_file,theAtoms,theBasis,nbf,nbfc,do_s2c,opensh);
  }else if (type == fermions){
    read_output_fermions(qcout_name,basfile,theAtoms,theBasis,nbf,nbfc);
    read_density_fermions(scr_dir,theAtoms,theBasis,nbf,nbfc,do_s2c,opensh);
  }else{
    die("Illegal output-type!");
  }

}

void read_output_qchem(std::string& fname, std::string& basfile, std::vector<atom>& theAtoms, std::vector<atombasis>& theBasis, int& nbf, int& nbfc){

  theAtoms.clear();
  theBasis.clear();

  std::ifstream input;
  input.open(fname);
  if (!input.good()){
    printf("File: %s\n",fname.c_str());
    die("Can't open output-file...");
  }

  std::string line;

  bool go_on = true;
  bool got_bas = false;
  bool got_coords = false;
  int state = 0;
  atombasis ascr;
  atomshell sscr;
  std::vector<std::string> buff;
  std::string buf;
  int ioz_scr = -1;
  std::set<int> the_oz;
  int np_act = -1;
  int np_count = -1;
  while(go_on){
    getline(input,line);
    trim(line);
    if (input.eof()) go_on = false;
    if (line.length() == 0) continue;
    if (state == 0){
      if (line.compare("Standard Nuclear Orientation (Angstroms)") == 0){
        state = 1;
        got_coords = true;
      }else if (line.compare("$basis") == 0){
        state = 4;
        got_bas = true;
      }else if (line.compare("Atom    Z    L     Exponents        Contraction Coefficients") == 0){
        state = 8;
        got_bas = true;
      }
    }else if (state == 1 || state == 2){ // skip and inc
      state++;
    }else if (state == 3){ // read geo
      if (std::search_n(line.begin(), line.end(), 50, '-') == line.begin()) {
        state = 0;
        if (got_bas) go_on = false;
      }else{
        std::transform(line.begin(), line.end(), line.begin(), tolower);
        // add atom...
        std::stringstream ss(line);
        buff.clear();
        while(ss >> buf) buff.push_back(buf);
        if (buff.size() != 5){
          printf("line: %s\n",line.c_str());
          die("Illegal format!");
        }
        atom toadd;
        toadd.ioz = get_ioz(buff[1]);
        toadd.x   = atof(buff[2].c_str());
        toadd.y   = atof(buff[3].c_str());
        toadd.z   = atof(buff[4].c_str());
        theAtoms.push_back(toadd);
        the_oz.insert(toadd.ioz);
      }
    }else if (state == 4){ // basis: get atom-type
      if (line.compare("$end") == 0){
        state = 0;
        if (got_coords) go_on = false;
        if (ioz_scr > 0){ // store last set...
          if (sscr.linc.size() > 0){
            atomshell sadd;
            sadd.lqn = sscr.lqn;
            for(size_t ii=0;ii<sscr.linc.size();ii++){
              sadd.linc.push_back(sscr.linc[ii]);
              sadd.expc.push_back(sscr.expc[ii]);
            }
            sscr.linc.clear();
            sscr.expc.clear();
            ascr.shells.push_back(sadd);
          }
          atombasis toadd;
          toadd.ioz = ioz_scr;
          for(auto it : ascr.shells){
            atomshell ashell;
            ashell.lqn = it.lqn;
            for(size_t ii=0;ii<it.linc.size();ii++){
              ashell.linc.push_back(it.linc[ii]);
              ashell.expc.push_back(it.expc[ii]);
            }
            toadd.shells.push_back(ashell);
          }
          toadd.max_lqn = toadd.shells[toadd.shells.size()-1].lqn;
          theBasis.push_back(toadd);
          ascr.shells.clear();
        }
        ioz_scr = -1;
        continue;
      }
      std::transform(line.begin(), line.end(), line.begin(), tolower);
      // add atom...
      std::stringstream ss(line);
      buff.clear();
      while(ss >> buf) buff.push_back(buf);
      if (buff.size() != 2){
        printf("line: %s\n",line.c_str());
        die("Illegal format!");
      }
      int oz_act = get_ioz(buff[0]);
      if (ioz_scr > 0){ // store last set...
        if (sscr.linc.size() > 0){
          atomshell sadd;
          sadd.lqn = sscr.lqn;
          for(size_t ii=0;ii<sscr.linc.size();ii++){
            sadd.linc.push_back(sscr.linc[ii]);
            sadd.expc.push_back(sscr.expc[ii]);
          }
          sscr.linc.clear();
          sscr.expc.clear();
          ascr.shells.push_back(sadd);
        }
        atombasis toadd;
        toadd.ioz = ioz_scr;
        for(auto it : ascr.shells){
          atomshell ashell;
          ashell.lqn = it.lqn;
          for(size_t ii=0;ii<it.linc.size();ii++){
            ashell.linc.push_back(it.linc[ii]);
            ashell.expc.push_back(it.expc[ii]);
          }
          toadd.shells.push_back(ashell);
        }
        toadd.max_lqn = toadd.shells[toadd.shells.size()-1].lqn;
        theBasis.push_back(toadd);
        ascr.shells.clear();
      }
      ioz_scr = -1;
      if (the_oz.find(oz_act) != the_oz.end()){
        ioz_scr = oz_act;
        state = 5;
      }
    }else if (state == 5){ // basis: shell-start
      std::transform(line.begin(), line.end(), line.begin(), tolower);
      std::stringstream ss(line);
      buff.clear();
      while(ss >> buf) buff.push_back(buf);
      if (buff.size() != 3) die("Illegal format!");
      if (sscr.linc.size() > 0){
        atomshell sadd;
        sadd.lqn = sscr.lqn;
        for(size_t ii=0;ii<sscr.linc.size();ii++){
          sadd.linc.push_back(sscr.linc[ii]);
          sadd.expc.push_back(sscr.expc[ii]);
        }
        sscr.linc.clear();
        sscr.expc.clear();
        ascr.shells.push_back(sadd);
      }
      sscr.lqn = get_lqn(buff[0]);
      np_act = atoi(buff[1].c_str());
      np_count = 0;
      state = 6;
    }else if (state == 6){ // basis: read shell
      std::transform(line.begin(), line.end(), line.begin(), tolower);
      std::stringstream ss(line);
      buff.clear();
      while(ss >> buf) buff.push_back(buf);
      if (buff.size() != 2) die("Illegal format!");
      np_count++;
      sscr.linc.push_back(atof(buff[1].c_str()));
      sscr.expc.push_back(atof(buff[0].c_str()));
      if (np_count == np_act) state = 7;
    }else if (state == 7){ // basis: shell/atom start
      std::transform(line.begin(), line.end(), line.begin(), tolower);
      // add atom...
      std::stringstream ss(line);
      buff.clear();
      while(ss >> buf) buff.push_back(buf);
      if (buff.size() == 2){
        int oz_act = get_ioz(buff[0]);
        if (ioz_scr > 0){ // store last set...
          if (sscr.linc.size() > 0){
            atomshell sadd;
            sadd.lqn = sscr.lqn;
            for(size_t ii=0;ii<sscr.linc.size();ii++){
              sadd.linc.push_back(sscr.linc[ii]);
              sadd.expc.push_back(sscr.expc[ii]);
            }
            sscr.linc.clear();
            sscr.expc.clear();
            ascr.shells.push_back(sadd);
          }
          atombasis toadd;
          toadd.ioz = ioz_scr;
          for(auto it : ascr.shells){
            atomshell ashell;
            ashell.lqn = it.lqn;
            for(size_t ii=0;ii<it.linc.size();ii++){
              ashell.linc.push_back(it.linc[ii]);
              ashell.expc.push_back(it.expc[ii]);
            }
            toadd.shells.push_back(ashell);
          }
          toadd.max_lqn = toadd.shells[toadd.shells.size()-1].lqn;
          theBasis.push_back(toadd);
          ascr.shells.clear();
        }
        ioz_scr = -1;
        if (the_oz.find(oz_act) != the_oz.end()){
          ioz_scr = oz_act;
          state = 5;
        }else{
          state = 4;
        }
      }else if (buff.size() == 3){
        if (sscr.linc.size() > 0){
          atomshell sadd;
          sadd.lqn = sscr.lqn;
          for(size_t ii=0;ii<sscr.linc.size();ii++){
            sadd.linc.push_back(sscr.linc[ii]);
            sadd.expc.push_back(sscr.expc[ii]);
          }
          sscr.linc.clear();
          sscr.expc.clear();
          ascr.shells.push_back(sadd);
        }
        sscr.lqn = get_lqn(buff[0]);
        np_act = atoi(buff[1].c_str());
        np_count = 0;
        state = 6;
      }else if (line.compare("$end") == 0){
        state = 0;
        if (got_coords) go_on = false;
        if (ioz_scr > 0){ // store last set...
          if (sscr.linc.size() > 0){
            atomshell sadd;
            sadd.lqn = sscr.lqn;
            for(size_t ii=0;ii<sscr.linc.size();ii++){
              sadd.linc.push_back(sscr.linc[ii]);
              sadd.expc.push_back(sscr.expc[ii]);
            }
            sscr.linc.clear();
            sscr.expc.clear();
            ascr.shells.push_back(sadd);
          }
          atombasis toadd;
          toadd.ioz = ioz_scr;
          for(auto it : ascr.shells){
            atomshell ashell;
            ashell.lqn = it.lqn;
            for(size_t ii=0;ii<it.linc.size();ii++){
              ashell.linc.push_back(it.linc[ii]);
              ashell.expc.push_back(it.expc[ii]);
            }
            toadd.shells.push_back(ashell);
          }
          toadd.max_lqn = toadd.shells[toadd.shells.size()-1].lqn;
          theBasis.push_back(toadd);
          ascr.shells.clear();
        }
        ioz_scr = -1;
      }else if (line.compare("****") == 0){
        if (ioz_scr > 0){ // store last set...
          if (sscr.linc.size() > 0){
            atomshell sadd;
            sadd.lqn = sscr.lqn;
            for(size_t ii=0;ii<sscr.linc.size();ii++){
              sadd.linc.push_back(sscr.linc[ii]);
              sadd.expc.push_back(sscr.expc[ii]);
            }
            sscr.linc.clear();
            sscr.expc.clear();
            ascr.shells.push_back(sadd);
          }
          atombasis toadd;
          toadd.ioz = ioz_scr;
          for(auto it : ascr.shells){
            atomshell ashell;
            ashell.lqn = it.lqn;
            for(size_t ii=0;ii<it.linc.size();ii++){
              ashell.linc.push_back(it.linc[ii]);
              ashell.expc.push_back(it.expc[ii]);
            }
            toadd.shells.push_back(ashell);
          }
          toadd.max_lqn = toadd.shells[toadd.shells.size()-1].lqn;
          theBasis.push_back(toadd);
          ascr.shells.clear();
        }
        ioz_scr = -1;
        state = 4;
      }else{
        printf("line: %s\n",line.c_str());
        die("Illegal format!");
      }

    }else if (state == 8){ // basis: inc
      state = 9;
    }else if (state == 9){ // basis: start
      if (std::search_n(line.begin(), line.end(), 50, '-') == line.begin()) {
        state = 0;
        if (got_coords) go_on = false;
        if (ioz_scr > 0){
          if (sscr.linc.size() > 0){
            atomshell sadd;
            sadd.lqn = sscr.lqn;
            for(size_t ii=0;ii<sscr.linc.size();ii++){
              sadd.linc.push_back(sscr.linc[ii]);
              sadd.expc.push_back(sscr.expc[ii]);
            }
            sscr.linc.clear();
            sscr.expc.clear();
            ascr.shells.push_back(sadd);
          }
          atombasis toadd;
          toadd.ioz = ioz_scr;
          for(auto it : ascr.shells){
            atomshell ashell;
            ashell.lqn = it.lqn;
            for(size_t ii=0;ii<it.linc.size();ii++){
              ashell.linc.push_back(it.linc[ii]);
              ashell.expc.push_back(it.expc[ii]);
            }
            toadd.shells.push_back(ashell);
          }
          toadd.max_lqn = toadd.shells[toadd.shells.size()-1].lqn;
          theBasis.push_back(toadd);
          ascr.shells.clear();
        }
      }else{
        std::transform(line.begin(), line.end(), line.begin(), tolower);
        // add atom...
        std::stringstream ss(line);
        buff.clear();
        while(ss >> buf) buff.push_back(buf);
        if (buff.size() == 5){
          if (ioz_scr > 0){
            if (sscr.linc.size() > 0){
              atomshell sadd;
              sadd.lqn = sscr.lqn;
              for(size_t ii=0;ii<sscr.linc.size();ii++){
                sadd.linc.push_back(sscr.linc[ii]);
                sadd.expc.push_back(sscr.expc[ii]);
              }
              sscr.linc.clear();
              sscr.expc.clear();
              ascr.shells.push_back(sadd);
            }
            atombasis toadd;
            toadd.ioz = ioz_scr;
            for(auto it : ascr.shells){
              atomshell ashell;
              ashell.lqn = it.lqn;
              for(size_t ii=0;ii<it.linc.size();ii++){
                ashell.linc.push_back(it.linc[ii]);
                ashell.expc.push_back(it.expc[ii]);
              }
              toadd.shells.push_back(ashell);
            }
            toadd.max_lqn = toadd.shells[toadd.shells.size()-1].lqn;
            theBasis.push_back(toadd);
            ascr.shells.clear();
          }
          ioz_scr = get_ioz(buff[1]);
          sscr.lqn = get_lqn(buff[2]);
          sscr.linc.push_back(atof(buff[4].c_str()));
          sscr.expc.push_back(atof(buff[3].c_str()));
        }else if (buff.size() == 3){
          atomshell sadd;
          sadd.lqn = sscr.lqn;
          for(size_t ii=0;ii<sscr.linc.size();ii++){
            sadd.linc.push_back(sscr.linc[ii]);
            sadd.expc.push_back(sscr.expc[ii]);
          }
          sscr.linc.clear();
          sscr.expc.clear();
          ascr.shells.push_back(sadd);

          sscr.lqn = get_lqn(buff[0]);
          sscr.linc.push_back(atof(buff[2].c_str()));
          sscr.expc.push_back(atof(buff[1].c_str()));
          
        }else if (buff.size() == 2){
          sscr.linc.push_back(atof(buff[1].c_str()));
          sscr.expc.push_back(atof(buff[0].c_str()));
        }
      }
    }
  }

  if (!got_bas){ // read it from file, however, I also could simply...
    die("No basis specification in output-file!");
  }

  // ok, write to 'mol'
  FILE* fout = fopen("mol","w");
  if (fout == NULL) die("Error creating file 'mol'!");
  fprintf(fout,"INTGRL        1    0    1    0    0    0    0    0    0\n");
  fprintf(fout,"Q-CHEM\n");
  fprintf(fout,"              Generated by gimic_conv\n");
  fprintf(fout,"%i    0            0.10E-08              0    0\n",(int)theAtoms.size());
  fprintf(fout,"9999.00      3.00\n");
  nbf  = 0;
  nbfc = 0;
  int ac = 0;
  for(auto it : theAtoms){
    // basis
    size_t which = 0;
    for(size_t ii=0;ii<theBasis.size();ii++){
       if (it.ioz == theBasis[ii].ioz) which = ii;
    }
    fprintf(fout,"%.1f    1%2i",(float)it.ioz,(int)theBasis[which].max_lqn+1);
    int nl = 0;
    int l_act = -1;
    for(size_t ii=0; ii < theBasis[which].shells.size(); ii++){
      if (theBasis[which].shells[ii].lqn != l_act){
        if (l_act >= 0){
          fprintf(fout,"%3i",nl);
        }
        l_act = theBasis[which].shells[ii].lqn;
        nl = 1;
      }else{
        nl++;
      }
    }
    fprintf(fout,"%3i\n",nl);
    std::string coz;
    get_coz(it.ioz,coz);
    fprintf(fout,"%2s 1%20.12f%20.12f%20.12f\n",coz.c_str(),it.x * ANG2BOHR,it.y * ANG2BOHR,it.z * ANG2BOHR);
    // print basis
    int nbf_act  = 0;
    int nbfc_act = 0;
    for(auto it2 : theBasis[which].shells){
      fprintf(fout,"%6i  1\n",(int)it2.linc.size());
      for(size_t ii=0;ii<it2.linc.size();ii++){
        fprintf(fout,"%18.10f %15.10f\n",it2.expc[ii],it2.linc[ii]);
      }
      nbf_act  += 2 * it2.lqn + 1;
      nbfc_act += ((it2.lqn+1)*(it2.lqn+2))/2;
    }
    theAtoms[ac].nbf   = nbf_act;
    theAtoms[ac].nbfc  = nbfc_act;
    theAtoms[ac].offp  = nbf;
    theAtoms[ac].offc  = nbfc;
    nbf     += nbf_act;
    nbfc    += nbfc_act;
    ac++;
  }
  fprintf(fout,"\n");
  fclose(fout);

}


void read_density_qchem(std::string& dens_file, std::string& xdens_file, std::vector<atom>& theAtoms, std::vector<atombasis>& theBasis, int nbf, int nbfc, bool do_s2c, bool opensh){

  // ACES ordering == fermions ordering

  std::vector<double> matrix0(nbfc*nbfc);
  
  FILE* fin = fopen(dens_file.c_str(),"rb");
  if (fin == NULL) die("Could not open density-file!");

  int dim_in = nbfc*nbfc;
  if (do_s2c) dim_in = nbf*nbf;

  if (fread(&matrix0[0],dim_in*8,1,fin) != 1){
    die("Failed to read from density-file!");
  }

  fclose(fin);

  if (do_s2c) conv_s2c(qchem,matrix0,theAtoms,theBasis,nbf,nbfc);
  nonaxc(qchem,matrix0,theAtoms,theBasis,nbf,nbfc);
  qchem_reorder(matrix0,theAtoms,theBasis,nbf,nbfc);


  FILE* fout = fopen("xdens","w");
  if (fout == NULL) die("Could not create density-output 0!");

  for(int elem=0;elem<nbfc*nbfc;elem++)
    fprintf(fout,"%.14E\n",matrix0[elem]);
  fprintf(fout,"\n");

  fin = fopen(xdens_file.c_str(),"rb");
  if (fin == NULL) die("Could not open density-file!");
  size_t off = 0;
  for(int ii=0;ii<3;ii++){

    if (fseek(fin,off*8,SEEK_SET) != 0)
      die("fseek-error!");

    if (fread(&matrix0[0],dim_in*8,1,fin) != 1){
      die("Failed to read from density-file!");
    }
  
    if (do_s2c) conv_s2c(qchem,matrix0,theAtoms,theBasis,nbf,nbfc);
    nonaxc(qchem,matrix0,theAtoms,theBasis,nbf,nbfc);
    qchem_reorder(matrix0,theAtoms,theBasis,nbf,nbfc);
  
  
    for(int elem=0;elem<nbfc*nbfc;elem++)
      fprintf(fout,"%.14E\n",matrix0[elem]);
    fprintf(fout,"\n");
    if (do_s2c) off += nbf*nbf;
    else        off += nbfc*nbfc;
  }
  fclose(fin);
  fclose(fout);



}

void qchem_reorder(std::vector<double>& matrix, std::vector<atom>& theAtoms, std::vector<atombasis>& theBasis, int nbf, int nbfc){

  std::vector<double> mat1(nbfc*nbfc);
  for(int at0=0;at0<nbfc*nbfc;at0++) mat1[at0] = matrix[at0];

  for(int at0=theAtoms.size()-1;at0>=0;at0--){
    size_t which = 0;
    for(size_t ii=0;ii<theBasis.size();ii++){
       if (theAtoms[at0].ioz == theBasis[ii].ioz) which = ii;
    }
    int offc1 = nbfc;
    if (at0 < theAtoms.size()-1){
      offc1 = theAtoms[at0+1].offc;
    }
    for(int at1=theBasis[which].shells.size()-1;at1>=0;at1--){
      int lqn1  = theBasis[which].shells[at1].lqn;
      int dim1c = ((lqn1+1)*(lqn1+2))/2;
      offc1 -= dim1c;
      for(int bt0=theAtoms.size()-1;bt0>=0;bt0--){
        size_t which2 = 0;
        for(size_t ii=0;ii<theBasis.size();ii++){
           if (theAtoms[bt0].ioz == theBasis[ii].ioz) which2 = ii;
        }
        int offc2 = nbfc;
        if (bt0 < theAtoms.size()-1){
          offc2 = theAtoms[bt0+1].offc;
        }
        for(int bt1=theBasis[which2].shells.size()-1;bt1>=0;bt1--){
          int lqn2  = theBasis[which2].shells[bt1].lqn;
          int dim2c = ((lqn2+1)*(lqn2+2))/2;
          offc2 -= dim2c;
          // resort...
          int cc0=0;
          for(int Lz0=0;Lz0<=lqn1;Lz0++){
            for(int Ly0=0;Ly0<=lqn1-Lz0;Ly0++){
              int Lx0 = lqn1-Ly0-Lz0;
              int pos_fermions0 = lxyz2pos(Lx0,Ly0,Lz0);
              int cc1=0;
              for(int Lz1=0;Lz1<=lqn2;Lz1++){
                for(int Ly1=0;Ly1<=lqn2-Lz1;Ly1++){
                  int Lx1 = lqn2-Ly1-Lz1;
                  int pos_fermions1 = lxyz2pos(Lx1,Ly1,Lz1);
                  matrix[cc0 +  offc1 + (cc1+offc2)*nbfc] = mat1[pos_fermions0 +  offc1 + (pos_fermions1+offc2)*nbfc];
                  cc1++;
                }
              }
              cc0++;
            }
          }
        } 
      } 
    }
  }

}

void nonaxc(prog_type typ, std::vector<double>& matrix, std::vector<atom>& theAtoms, std::vector<atombasis>& theBasis, int nbf, int nbfc){

  int max_lqn = 0;
  for(auto it : theBasis){
    for(auto it2 : it.shells) max_lqn = (max_lqn > it2.lqn) ? max_lqn : it2.lqn;
  }

  double** nax;
  if (typ == qchem)
    form_nax_qchem(&nax,max_lqn);
  else if (typ == fermions)
    form_nax_fermions(&nax,max_lqn);
  else
    die("Illegal prog-type!");
  double* mat = &matrix[0];

  for(int at0=theAtoms.size()-1;at0>=0;at0--){
    size_t which = 0;
    for(size_t ii=0;ii<theBasis.size();ii++){
       if (theAtoms[at0].ioz == theBasis[ii].ioz) which = ii;
    }
    int offc1 = nbfc;
    if (at0 < theAtoms.size()-1){
      offc1 = theAtoms[at0+1].offc;
    }
    for(int at1=theBasis[which].shells.size()-1;at1>=0;at1--){
      int lqn1  = theBasis[which].shells[at1].lqn;
      int dim1c = ((lqn1+1)*(lqn1+2))/2;
      offc1 -= dim1c;
      for(int bt0=theAtoms.size()-1;bt0>=0;bt0--){
        size_t which2 = 0;
        for(size_t ii=0;ii<theBasis.size();ii++){
           if (theAtoms[bt0].ioz == theBasis[ii].ioz) which2 = ii;
        }
        int offc2 = nbfc;
        if (bt0 < theAtoms.size()-1){
          offc2 = theAtoms[bt0+1].offc;
        }
        for(int bt1=theBasis[which2].shells.size()-1;bt1>=0;bt1--){
          int lqn2  = theBasis[which2].shells[bt1].lqn;
          int dim2c = ((lqn2+1)*(lqn2+2))/2;
          offc2 -= dim2c;
  
          if (lqn1 > 1){
            for(int J=0;J<dim2c;J++){
              for(int ii=0;ii<dim1c;ii++){
                mat[ii+offc1 + (J+offc2)*nbfc] *= nax[lqn1][ii];
              }
            }
          }
          if (lqn2 > 1){
            for(int J=0;J<dim1c;J++){
              for(int ii=0;ii<dim2c;ii++){
                mat[J+offc1 + (ii+offc2)*nbfc] *= nax[lqn2][ii];
              }
            }
          }
        }
      }
    }
  }

  for (int l=0; l<=max_lqn; l++)
    delete[] nax[l];

  delete[] nax;

}

void conv_s2c(prog_type typ, std::vector<double>& matrix, std::vector<atom>& theAtoms, std::vector<atombasis>& theBasis, int nbf, int nbfc){

  int max_lqn = 0;
  for(auto it : theBasis){
    for(auto it2 : it.shells) max_lqn = (max_lqn > it2.lqn) ? max_lqn : it2.lqn;
  }

  double** c2p;
  if (typ == qchem)
    build_c2p_qchem(&c2p,max_lqn);
  else if (typ == fermions)
    build_c2p_fermions(&c2p,max_lqn);
  else
    die("Illegal prog-type!");

  std::vector<double> scratch(((max_lqn+1)*(max_lqn+2)/2)*nbfc);

  double* mat = &matrix[0];
  // transform rows, i.e. first expand to cartesian size...
  for(int i=(nbf-1);i>0;i--){
    double* cscr = mat + nbfc * i;
    double* pscr = mat + nbf  * i;

    for(int j=(nbf-1);j>=0;j--){
      cscr[j] = pscr[j];
    }

  }

  for(int i=0;i<nbf ;i++){ // COLUMNS

    int coffcol = i*nbfc;
    double* pscr = mat + coffcol;

    for(int k=0;k<nbf ;k++) scratch[k] = pscr[k];

    for(int at0=theAtoms.size()-1;at0>=0;at0--){
      size_t which = 0;
      for(size_t ii=0;ii<theBasis.size();ii++){
         if (theAtoms[at0].ioz == theBasis[ii].ioz) which = ii;
      }
      int offp1 = nbf;
      int offc1 = nbfc;
      if (at0 < theAtoms.size()-1){
        offp1 = theAtoms[at0+1].offp;
        offc1 = theAtoms[at0+1].offc;
      }
      for(int at1=theBasis[which].shells.size()-1;at1>=0;at1--){
        int lqn1  = theBasis[which].shells[at1].lqn;
        int dim1p = 2*lqn1+1;
        int dim1c = ((lqn1+1)*(lqn1+2))/2;
        offc1 -= dim1c;
        offp1 -= dim1p;
    
        double* cscr = mat     + coffcol + offc1;
        double* pscr = &scratch[0] + offp1;
        if (lqn1 > 1){
          for(int l=0;l<dim1c;l++){
            double sum = 0.e0;
            for(int k=0;k<dim1p;k++){
              sum += c2p[lqn1][k+l*dim1p] * pscr[k];
            }
            cscr[l] = sum;
          }
        }else{ // just copy
          for(int k=0;k<dim1c;k++) cscr[k] = pscr[k];
        }
      }
    }
  }
  // transform columns shell by shell...
  for(int at0=theAtoms.size()-1;at0>=0;at0--){
    size_t which = 0;
    for(size_t ii=0;ii<theBasis.size();ii++){
       if (theAtoms[at0].ioz == theBasis[ii].ioz) which = ii;
    }
    int offp1 = nbf;
    int offc1 = nbfc;
    if (at0 < theAtoms.size()-1){
      offp1 = theAtoms[at0+1].offp;
      offc1 = theAtoms[at0+1].offc;
    }
    for(int at1=theBasis[which].shells.size()-1;at1>=0;at1--){
      int lqn1  = theBasis[which].shells[at1].lqn;
      int dim1p = 2*lqn1+1;
      int dim1c = ((lqn1+1)*(lqn1+2))/2;
      offc1 -= dim1c;
      offp1 -= dim1p;
  
      double* cscr = mat + offc1*nbfc;
      double* pscr = mat + offp1*nbfc;
      if (lqn1 > 1){
        for(int i=0;i<nbfc*dim1p;i++) scratch[i] = pscr[i];
        for(int i=0;i<nbfc;i++){ // row for row...
          for(int l=0;l<dim1c;l++){
            double sum = 0.e0;
            for(int k=0;k<dim1p;k++){
              sum += c2p[lqn1][k+l*dim1p] * scratch[i+k*nbfc];
            }
            cscr[i+l*nbfc] = sum;
          }
        }
      }else{ // just copy
        for(int i=(nbfc*dim1c-1);i>=0;i--) cscr[i] = pscr[i];
      }
    }
  }

  for (int l=0; l<=max_lqn; l++)
    delete[] c2p[l];

  delete[] c2p;

}

double xyz2lm_Coeff(int, int, int, int, int);

#define parity(m) ((m)%2 ? -1 : 1)

double factorial(int n)
{

   if (n == 0 || n == 1) return(1.0);
   if (n < 0) return(0.0) ;
   else {
      return ((double) n * factorial(n-1)) ;
      }
}

double combinations(int n, int k)
{
   double comb ;

   if (n == k) return (1.0) ;
   else if (k > n) return(0.0) ;
   else if (k == 0) return(1.0) ;
   comb = factorial(n) / (factorial(k) * factorial(n-k)) ;

   return(comb) ;
}

// IJQC 54, 83 (1995)
double xyz2lm_Coeff(int l, int m, int lx, int ly, int lz, double* bc, int bc_dim, double* fac) {
  int i, j, k, i_max;
  int k_min, k_max;
  int abs_m;
  int comp;
  double pfac, pfac1, sum, sum1;
  
  abs_m = abs(m);
  if ((lx + ly - abs(m))%2)
    return 0.0;
  else
    j = (lx + ly - abs(m))/2;
  
  if (j < 0)
    return 0.0;
  
  comp = (m >= 0) ? 1 : -1;
  i = abs_m-lx;
  if (comp != parity(i))
    return 0.0;
  
  pfac = sqrt(fac[2*lx]*fac[2*ly]*fac[2*lz]*fac[l-abs_m]/(fac[2*l]*fac[l]
      *fac[lx]*fac[ly]*fac[lz]*fac[l+abs_m]));
  pfac /= (1 << l);
  
  if (m < 0)
    pfac *= parity((i-1)/2);
  else
    pfac *= parity(i/2);
  
  i_max = (l-abs_m)/2;
  sum = 0.0;
  for (i=0; i<=i_max; i++) {
    pfac1 = bc[l+bc_dim*i]*bc[i+bc_dim*j];
    if (fabs(pfac1) < 1e-20)
      continue;
    else
      pfac1 *= (parity(i)*fac[2*(l-i)]/fac[l-abs_m-2*i]);
    sum1 = 0.0;
    k_min = ((lx-abs_m)/2 > 0) ? (lx-abs_m)/2 : 0;
    k_max = (j < lx/2) ? j : lx/2;
    for (k=k_min; k<=k_max; k++)
      sum1 += bc[j+bc_dim*k]*bc[abs_m+bc_dim*(lx-2*k)]*parity(k);
    sum += pfac1*sum1;
  }
  
  if (m == 0)
    return pfac*sum;
  else
    return sqrt(2.e0)*pfac*sum;
}

void build_c2p_qchem(double*** c2p_tm, int max_lqn) {

    int i, j, m;
    int ao, l;
    
    *c2p_tm = new double*[max_lqn+1];
    for (l=0; l<=max_lqn; l++)
      (*c2p_tm)[l] = new double[(2*l+1)*((l+1)*(l+2)/2)];
    double* fac = new double[2*max_lqn+1];
    fac[0] = 1.0;
    for (i=1; i<=2*max_lqn; i++)
      fac[i] = i*fac[i-1];
    int bc_dim = max_lqn+1;
    double* bc = new double[(max_lqn+1)*(max_lqn+1)];
    for (i=0; i<(max_lqn+1)*(max_lqn+1); i++) bc[i] = 0.e0;
    for (i=0; i<=max_lqn; i++)
      for (j=0; j<=i; j++)
        bc[i+j*bc_dim] = combinations(i, j);
    
    for (l=0; l<=max_lqn; l++) {
      for (m=-l; m<=l; m++){
        ao = 0;      

          for(int Lz=0;Lz<=l;Lz++){
            for(int Ly=0;Ly<=l-Lz;Ly++){
              int Lx = l-Ly-Lz;
              (*c2p_tm)[l][m+l + ao*(2*l+1)] = xyz2lm_Coeff(l,m,Lx,Ly,Lz,bc,bc_dim,fac);
              ao++;
            }
          }

      }
    }
    
    delete[] fac;
    delete[] bc;
    return;

}

void build_c2p_fermions(double*** c2p_tm, int max_lqn) {
    int i, j, m;
    int ao, l;
    
    *c2p_tm = new double*[max_lqn+1];
    for (l=0; l<=max_lqn; l++)
      (*c2p_tm)[l] = new double[(2*l+1)*((l+1)*(l+2)/2)];
    double* fac = new double[2*max_lqn+1];
    fac[0] = 1.0;
    for (i=1; i<=2*max_lqn; i++)
      fac[i] = i*fac[i-1];
    int bc_dim = max_lqn+1;
    double* bc = new double[(max_lqn+1)*(max_lqn+1)];
    for (i=0; i<(max_lqn+1)*(max_lqn+1); i++) bc[i] = 0.e0;
    for (i=0; i<=max_lqn; i++)
      for (j=0; j<=i; j++)
        bc[i+j*bc_dim] = combinations(i, j);
    
    for (l=0; l<=max_lqn; l++) {
      for (m=-l; m<=l; m++){
        ao = 0;      
        for(int Lx=l;Lx>=0;Lx--){
          int m2 = l-Lx;
          for(int Ly=m2;Ly>=0;Ly--){
            int Lz = l-Lx-Ly;
            (*c2p_tm)[l][m+l + ao*(2*l+1)] = xyz2lm_Coeff(l,m,Lx,Ly,Lz,bc,bc_dim,fac);
            ao++;
          }
        }
      }
    }
    
    delete[] fac;
    delete[] bc;
    return;
}

void form_nax_qchem(double*** nax, int max_lqn){

  double** nax0;
  form_nax_fermions(&nax0,max_lqn);
  *nax = new double*[max_lqn+1];
  for(int ii=0;ii<=max_lqn;ii++) (*nax)[ii] = new double[((ii+1)*(ii+2))/2];
  
  // resort...
  for(int L=0;L<=max_lqn;L++){
    int cc=0;
    for(int Lz=0;Lz<=L;Lz++){
      for(int Ly=0;Ly<=L-Lz;Ly++){
        int Lx = L-Ly-Lz;
        int pos_fermions = lxyz2pos(Lx,Ly,Lz);
        (*nax)[L][cc] = nax0[L][pos_fermions];
        cc++;
      }
    }
    for(int ii=0;ii<(L+1)*(L+2)/2;ii++) (*nax)[L][ii] = sqrt((*nax)[L][ii]);
  }

  for (int l=0; l<=max_lqn; l++)
    delete[] nax0[l];

  delete[] nax0;

}

void form_nax_fermions(double*** nax, int max_lqn){

  *nax = new double*[max_lqn+1];
  for(int ii=0;ii<=max_lqn;ii++) (*nax)[ii] = new double[(ii+1)*(ii+2)/2];
  
  double sq3i = 1.e0/sqrt(3.e0);

  for(int L=0;L<=max_lqn;L++){
    double d1 = 1.e0;
    int k=0;
    for(int Lx=L;Lx>=0;Lx--){
      int    i2 = L - Lx;
      double d2 = d1;
      for(int Ly=i2;Ly>=0;Ly--){
        int Lz = L - Lx - Ly;
        (*nax)[L][k] = d2;
        d2 *= ((double)(2*Ly-1))/((double)(2*Lz+1));
        k++;
      }
      d1 *= (((double)(2*Lx-1))/((double)(2*i2+1)));
    }
    for(int ii=0;ii<(L+1)*(L+2)/2;ii++) (*nax)[L][ii] = sqrt((*nax)[L][ii]);
    if (L > 1) for(int ii=0;ii<(L+1)*(L+2)/2;ii++) (*nax)[L][ii] *= sq3i;
  }
}

void print_mat(int nrow, int ncol, double* scr, const char* name){

  if (nrow == -1){ // vectorized, print columns of 25 elements (ncol == 25)
    int numel = ncol;
    nrow = 25;
    ncol = numel / nrow;
    int lastcol = numel - nrow*ncol;
    if (lastcol > 0) ncol++;

    int n6   = ncol / 6;
    int rest = ncol - 6*n6;

    printf("\n%s --- vectorized (%i elements)\n\n",name,numel);

    int i;
    for(i=0;i<n6;i++){
      printf("           %4i         %4i         %4i         %4i         %4i         %4i\n\n",
             i*6+1,i*6+2,i*6+3,i*6+4,i*6+5,i*6+6);
      for(int j=0;j<nrow;j++){
        printf(" %4i %13.7f%13.7f%13.7f%13.7f%13.7f%13.7f\n",j+1,
               scr[j + i*6*nrow],scr[j + (i*6+1)*nrow],scr[j + (i*6+2)*nrow],scr[j + (i*6+3)*nrow],
               scr[j + (i*6+4)*nrow],scr[j + (i*6+5)*nrow]);
      }
      printf("\n");
    }
    i = ncol - rest;
    if (rest > 0){
      printf("           %4i     ",i+1);
      for(int j=1;j<rest;j++) printf("    %4i     ",i+1+j);
      printf("\n\n");
      for(int j=0;j<nrow;j++){
        printf(" %4i ",j+1);
        if (j < lastcol) for(int k=0;k<rest;k++) printf("%13.7f",scr[j+(i+k)*nrow]);
        else for(int k=0;k<rest-1;k++) printf("%13.7f",scr[j+(i+k)*nrow]);
        printf("\n");
      }
    }
    printf("\n");

  }else{
// 6 cols
    int n6   = ncol / 6;
    int rest = ncol - 6*n6;
 
 
    printf("%s\n",name);
 
    int i;
    for(i=0;i<n6;i++){
    printf("         %4i        %4i        %4i        %4i        %4i         %4i\n",
             i*6+1,i*6+2,i*6+3,i*6+4,i*6+5,i*6+6);
      for(int j=0;j<nrow;j++){
      printf(" %3i %12.7f%12.7f%12.7f%12.7f%12.7f%12.7f\n",j+1,
               scr[j + i*6*nrow],scr[j + (i*6+1)*nrow],scr[j + (i*6+2)*nrow],scr[j + (i*6+3)*nrow],
               scr[j + (i*6+4)*nrow],scr[j + (i*6+5)*nrow]);
      }
    }
    i = ncol - rest;
    if (rest > 0){
      printf("         %4i",i+1);
      for(int j=1;j<rest;j++) printf("        %4i",i+1+j);
      printf("\n");
      for(int j=0;j<nrow;j++){
      printf(" %3i ",j+1);
        for(int k=0;k<rest;k++) printf("%12.7f",scr[j+(i+k)*nrow]);
        printf("\n");
      }
    }
  }
  fflush(stdout);
}

/*



   
    QCMINE




*/
void read_output_fermions(std::string& fname, std::string& basfile, std::vector<atom>& theAtoms, std::vector<atombasis>& theBasis, int& nbf, int& nbfc){

  theAtoms.clear();
  theBasis.clear();

  std::ifstream input;
  input.open(fname);
  if (!input.good()){
    printf("File: %s\n",fname.c_str());
    die("Can't open output-file...");
  }

  std::string line;

  bool go_on = true;
  bool got_bas = false;
  bool got_coords = false;
  int state = 0;
  atombasis ascr;
  atomshell sscr;
  std::vector<std::string> buff;
  std::string buf;
  int ioz_scr = -1;
  std::set<int> the_oz;
  int np_act = -1;
  int np_count = -1;
  std::string bas_name;
  while(go_on){
    getline(input,line);
    if (input.eof()) go_on = false;
    if (line.length() == 0) continue;
    if (state == 0){
      if (line.compare("  Coordinates in Bohr:") == 0){
        state = 1;
        got_coords = true;
      }else if (line.find("  BasisSet:") != std::string::npos){
        std::stringstream ss(line);
        buff.clear();
        while(ss >> buf) buff.push_back(buf);
        if (buff.size() != 2){
          printf("line: %s\n",line.c_str());
          die("Illegal format!");
        }
        bas_name = buff[1];
        if (got_coords) go_on = false;
      }
    }else if (state == 1 || state == 2 || state == 3){ // skip and inc
      state++;
    }else if (state == 4){ // read geo
      if (line.compare(" ----------------------------------------------------------------------------------------") == 0){
        state = 0;
        if (got_bas) go_on = false;
      }else{
        std::transform(line.begin(), line.end(), line.begin(), tolower);
        // add atom...
        std::stringstream ss(line);
        buff.clear();
        while(ss >> buf) buff.push_back(buf);
        if (buff.size() != 7){
          printf("line: %s\n",line.c_str());
          die("Illegal format!");
        }
        atom toadd;
        toadd.ioz = atoi(buff[2].c_str());
        toadd.x   = atof(buff[3].c_str());
        toadd.y   = atof(buff[4].c_str());
        toadd.z   = atof(buff[5].c_str());
        theAtoms.push_back(toadd);
        the_oz.insert(toadd.ioz);
      }
    }
  }

  // read basis
  if (bas_name.size() == 0) die("No basis-name?!?");

  if (std::getenv("QC_BASIS") == NULL) die("Env 'QC_BASIS' is not set!");
  std::string bas_file = std::getenv("QC_BASIS");
  bas_file += "/basis/";
  bas_file += bas_name;

  std::ifstream bas;
  bas.open(bas_file);
  if (!bas.good()){
    printf("File: %s\n",bas_file.c_str());
    die("Can't open basisset-file...");
  }
  go_on = true;
  state = 0;
  while(go_on){
    getline(bas,line);
    if (bas.eof()) go_on = false;
    if (line.length() == 0) continue;
    if (state == 0){
      if (line.compare("****") == 0){
        state = 1;
      }
    }else if (state == 1){ // new atom
      std::transform(line.begin(), line.end(), line.begin(), tolower);
      std::stringstream ss(line);
      buff.clear();
      while(ss >> buf) buff.push_back(buf);
      if (buff.size() != 2){
        state = 0;
      }else{
        int ioz_act = get_ioz(buff[0]);
        bool got_this = false;
        size_t which = 0;
        for(size_t ii=0;ii<theAtoms.size() && !got_this;ii++){
           if (theAtoms[ii].ioz == ioz_act) got_this = true;;
        }
        if (got_this){
          ioz_scr = ioz_act;
          state = 2;
        }else          state = 0;
      }
    }else if (state == 2){ // new shell
      if (line.compare("****") == 0){
        if (sscr.linc.size() > 0){
          atomshell sadd;
          sadd.lqn = sscr.lqn;
          for(size_t ii=0;ii<sscr.linc.size();ii++){
            sadd.linc.push_back(sscr.linc[ii]);
            sadd.expc.push_back(sscr.expc[ii]);
          }
          sscr.linc.clear();
          sscr.expc.clear();
          ascr.shells.push_back(sadd);
        }
        if (ascr.shells.size() > 0){
          atombasis toadd;
          toadd.ioz = ioz_scr;
          for(auto it : ascr.shells){
            atomshell ashell;
            ashell.lqn = it.lqn;
            for(size_t ii=0;ii<it.linc.size();ii++){
              ashell.linc.push_back(it.linc[ii]);
              ashell.expc.push_back(it.expc[ii]);
            }
            toadd.shells.push_back(ashell);
          }
          toadd.max_lqn = toadd.shells[toadd.shells.size()-1].lqn;
          theBasis.push_back(toadd);
          ascr.shells.clear();
        }
        state = 1;
      }else{
        std::transform(line.begin(), line.end(), line.begin(), tolower);
        std::stringstream ss(line);
        buff.clear();
        while(ss >> buf) buff.push_back(buf);
        if (buff.size() == 2){
          sscr.linc.push_back(atof(buff[1].c_str()));
          sscr.expc.push_back(atof(buff[0].c_str()));
        }else if (buff.size() == 3){
          if (sscr.linc.size() > 0){
            atomshell sadd;
            sadd.lqn = sscr.lqn;
            for(size_t ii=0;ii<sscr.linc.size();ii++){
              sadd.linc.push_back(sscr.linc[ii]);
              sadd.expc.push_back(sscr.expc[ii]);
            }
            sscr.linc.clear();
            sscr.expc.clear();
            ascr.shells.push_back(sadd);
          }
          sscr.lqn = get_lqn(buff[0]);
        }else{
          die("Illegal format!");
        }
      }
    }
  }

  

  // ok, write to 'mol'
  FILE* fout = fopen("mol","w");
  if (fout == NULL) die("Error creating file 'mol'!");
  fprintf(fout,"INTGRL        1    0    1    0    0    0    0    0    0\n");
  fprintf(fout,"TURBOMOLE\n");
  fprintf(fout,"              Generated by gimic_conv\n");
  fprintf(fout,"%i    0            0.10E-08              0    0\n",(int)theAtoms.size());
  fprintf(fout,"9999.00      3.00\n");
  nbf  = 0;
  nbfc = 0;
  int ac = 0;
  for(auto it : theAtoms){
    // basis
    size_t which = 0;
    for(size_t ii=0;ii<theBasis.size();ii++){
       if (it.ioz == theBasis[ii].ioz) which = ii;
    }
    fprintf(fout,"%.1f    1%2i",(float)it.ioz,(int)theBasis[which].max_lqn+1);
    int nl = 0;
    int l_act = -1;
    for(size_t ii=0; ii < theBasis[which].shells.size(); ii++){
      if (theBasis[which].shells[ii].lqn != l_act){
        if (l_act >= 0){
          fprintf(fout,"%3i",nl);
        }
        l_act = theBasis[which].shells[ii].lqn;
        nl = 1;
      }else{
        nl++;
      }
    }
    fprintf(fout,"%3i\n",nl);
    std::string coz;
    get_coz(it.ioz,coz);
    fprintf(fout,"%2s 1%20.12f%20.12f%20.12f\n",coz.c_str(),it.x,it.y,it.z);
    int nbf_act  = 0;
    int nbfc_act = 0;
    for(auto it2 : theBasis[which].shells){
      fprintf(fout,"%6i  1\n",(int)it2.linc.size());
      for(size_t ii=0;ii<it2.linc.size();ii++){
        fprintf(fout,"%18.10f %15.10f\n",it2.expc[ii],it2.linc[ii]);
      }
      nbf_act  += 2 * it2.lqn + 1;
      nbfc_act += ((it2.lqn+1)*(it2.lqn+2))/2;
    }
    theAtoms[ac].nbf   = nbf_act;
    theAtoms[ac].nbfc  = nbfc_act;
    theAtoms[ac].offp  = nbf;
    theAtoms[ac].offc  = nbfc;
    nbf     += nbf_act;
    nbfc    += nbfc_act;
    ac++;
  }
  fprintf(fout,"\n");
  fclose(fout);

}

void read_density_fermions(std::string& scr_dir, std::vector<atom>& theAtoms, std::vector<atombasis>& theBasis, int nbf, int nbfc, bool do_s2c, bool opensh){

  // ACES ordering == fermions ordering

  std::vector<double> matrix0(nbfc*nbfc);
  
  std::string toread = scr_dir;
  toread += "/data_38.0";
  FILE* fin = fopen(toread.c_str(),"rb");
  if (fin == NULL) die("Could not open density-file!");

  int dim_in = nbfc*nbfc;
  if (do_s2c) dim_in = nbf*nbf;

  if (fread(&matrix0[0],dim_in*8,1,fin) != 1){
    die("Failed to read from density-file!");
  }

  fclose(fin);

  if (do_s2c) conv_s2c(fermions,matrix0,theAtoms,theBasis,nbf,nbfc);
  nonaxc(fermions,matrix0,theAtoms,theBasis,nbf,nbfc);
  fermions2turbo(matrix0,theAtoms,theBasis,nbf,nbfc);


  FILE* fout = fopen("xdens","w");
  if (fout == NULL) die("Could not create density-output 0!");

  for(int elem=0;elem<nbfc*nbfc;elem++)
    if (fabs(matrix0[elem]) < 1e-12)
      fprintf(fout,"%.14E\n",0.e0);
    else
      fprintf(fout,"%.14E\n",2.e0*matrix0[elem]);
  fprintf(fout,"\n");

  std::string toread0 = scr_dir;
  toread0 += "/GIMIC_PB";
  fin = fopen(toread0.c_str(),"rb");
  if (fin == NULL){
    printf("file: '%s'\n",toread0.c_str());
    die("Could not open xdensity-file!");
  }

  for(int ii=0;ii<3;ii++){

    if (fread(&matrix0[0],dim_in*8,1,fin) != 1){
      die("Failed to read from density-file!");
    }
  
    if (do_s2c) conv_s2c(fermions,matrix0,theAtoms,theBasis,nbf,nbfc);
    nonaxc(fermions,matrix0,theAtoms,theBasis,nbf,nbfc);
    fermions2turbo(matrix0,theAtoms,theBasis,nbf,nbfc);
  
    for(int elem=0;elem<nbfc*nbfc;elem++){
      if (fabs(matrix0[elem]) < 1e-12)
        fprintf(fout,"%.14E\n",0.e0);
      else
        fprintf(fout,"%.14E\n",4e0*matrix0[elem]);
    }
    fprintf(fout,"\n");
  }
  fclose(fin);
  fclose(fout);



}

void fermions2turbo(std::vector<double>& matrix, std::vector<atom>& theAtoms, std::vector<atombasis>& theBasis, int nbf, int nbfc){

  std::vector<double> mat1(nbfc*nbfc);
  //for(int at0=0;at0<nbfc*nbfc;at0++) mat1[at0] = matrix[at0];

  std::vector<std::pair<int,int> > mapgto;

  int max_lqn = 0;
  for(auto it : theBasis){
    for(auto it2 : it.shells) max_lqn = (max_lqn > it2.lqn) ? max_lqn : it2.lqn;
  }

  int cc = 0;
  for(int ll=0;ll<=max_lqn;ll++){
    int offc1 = 0;
    for(int at0=0;at0<theAtoms.size();at0++){
      size_t which = 0;
      for(size_t ii=0;ii<theBasis.size();ii++){
         if (theAtoms[at0].ioz == theBasis[ii].ioz) which = ii;
      }
      for(int at1=0;at1<theBasis[which].shells.size();at1++){
        int lqn1  = theBasis[which].shells[at1].lqn;
        int dim1c = ((lqn1+1)*(lqn1+2))/2;
        if (lqn1 == ll){
          if (ll == 2){
            for(int oo=0;oo<dim1c;oo++){
              mapgto.push_back(std::make_pair(offc1+oo,cc+tm_d[oo]));
            }
            cc += dim1c;
          }else if (ll == 3){
            for(int oo=0;oo<dim1c;oo++){
              mapgto.push_back(std::make_pair(offc1+oo,cc+tm_f[oo]));
            }
            cc += dim1c;
          }else if (ll == 4){
            for(int oo=0;oo<dim1c;oo++){
              mapgto.push_back(std::make_pair(offc1+oo,cc+tm_g[oo]));
            }
            cc += dim1c;
          }else if (ll < 2){
            for(int oo=0;oo<dim1c;oo++){
              mapgto.push_back(std::make_pair(offc1+oo,cc));
              cc++;
            }
          }else{
            die("Only up to g-funcs...");
          }
        }
        offc1 += dim1c;
      }
    }
  }

  for(auto it : mapgto){
    int from = it.first * nbfc;
    int to   = it.second * nbfc;
    for(int oo=0;oo<nbfc;oo++) mat1[oo+to] = matrix[oo+from];
  }
  for(auto it : mapgto){
    int from = it.first;
    int to   = it.second;
    for(int oo=0;oo<nbfc;oo++) matrix[oo*nbfc+to] = mat1[oo*nbfc+from];
  }

}

