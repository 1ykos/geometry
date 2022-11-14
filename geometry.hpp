#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "wmath.hpp"

#include <dlib/matrix.h>

namespace geometry{
  using dlib::abs;
  using dlib::cholesky_decomposition;
  using dlib::diag;
  using dlib::diagm;
  using dlib::dot;
  using dlib::eigenvalue_decomposition;
  using dlib::identity_matrix;
  using dlib::inv;
  using dlib::length;
  using dlib::length_squared;
  using dlib::make_symmetric;
  using dlib::matrix;
  using dlib::matrix_exp;
  using dlib::matrix_op;
  using dlib::normalize;
  using dlib::op_make_symmetric;
  using dlib::round;
  using dlib::set_colm;
  using dlib::tmp;
  using dlib::trace;
  using dlib::zeros_matrix;
  using std::abs;
  using std::array;
  using std::cerr;
  using std::cin;
  using std::cout;
  using std::endl;
  using std::fill;
  using std::fixed;
  using std::floor;
  using std::get;
  using std::getline;
  using std::ifstream;
  using std::isnan;
  using std::istream;
  using std::lower_bound;
  using std::map;
  using std::max_element;
  using std::numeric_limits;
  using std::ofstream;
  using std::pow;
  using std::round;
  using std::setprecision;
  using std::setw;
  using std::sort;
  using std::stod;
  using std::streamsize;
  using std::string;
  using std::stringstream;
  using std::swap;
  using std::tuple;
  using std::unordered_map;
  using std::vector;
  using wmath::long_mul;
  using wmath::clz;
  using wmath::digits;

  constexpr bool VERBOSE = false;

  void inline get_inverse(
      const double& a,const double& b,const double& c,
      const double& d,const double& e,const double& f,
      const double& g,const double& h,const double& i,
      double& j,double& k,double& l,
      double& m,double& n,double& o,
      double& p,double& q,double& r
      ){
    const double A = (e*i-f*h);
    const double B =-(d*i-f*g);
    const double C = (d*h-e*g);
    const double D =-(b*i-c*h);
    const double E = (a*i-c*g);
    const double F =-(a*h-b*g);
    const double G = (b*f-c*e);
    const double H =-(a*f-c*d);
    const double I = (a*e-b*d);
    const double det = a*A+b*B+c*C;
    const double idt = 1.0/det;
    j = A*idt;
    k = D*idt;
    l = G*idt;
    m = B*idt;
    n = E*idt;
    o = H*idt;
    p = C*idt;
    q = F*idt;
    r = I*idt;
  }

  void inline solve3x4(double *m){
    for (size_t i=0;i!=3;++i){
      double max = -1; size_t k;
      for (size_t j=i;j!=3;++j) if (abs(m[j*4+i])>max){max=abs(m[j*4+i]); k=j;}
      if (k!=i) for (size_t l=0;l!=4;++l) swap(m[k*4+l],m[i*4+l]);
      if (m[i*4+i]<0) for (size_t l=i;l!=4;++l) m[i*4+l]/=-m[i*4+l];
      else for (size_t l=i;l!=4;++l)            m[i*4+l]/= m[i*4+i];
      for (size_t j=i+1;j!=3;++j)
        for (size_t l=i;l!=4;++l) m[j*4+l]-=m[j*4+i]/m[i*4+l];
    }
  }

  struct sector2d{
    uint32_t min_fs;
    uint32_t max_fs;
    uint32_t min_ss;
    uint32_t max_ss;
  };

  struct linear_panel_transform{
    uint32_t min_fs,max_fs,min_ss,max_ss;
    double fs_2_x,ss_2_x,xoffset, // detector coordinates
           fs_2_y,ss_2_y,yoffset, // to real world coordinates
           fs_2_z,ss_2_z,zoffset; 
    double x_2_fs,y_2_fs,z_2_fs,  // real world coordinates
           x_2_ss,y_2_ss,z_2_ss,  // to detector panel
           x_2_on,y_2_on,z_2_on;  // coordinates
    double n0,n1,n2; // normal vector of the detector plane
    size_t n_panel; //running panel id
    template<typename T>
    bool inline is_valid(const T&fs, const T&ss) const {
      //cout << min_fs << "<=" << fs << "<=" << max_fs << "&&"
      //     << min_ss << "<=" << ss << "<=" << max_ss << endl;
      return (ss<=max_ss)&&(fs<=max_fs)&&(fs>=min_fs)&&(ss>=min_ss);
    }
    template<typename T>
    void inline operator()(
        const T& fs,
        const T& ss,
        double& x,
        double& y,
        double& z
        ) const {
      if (VERBOSE) cerr << min_fs << " " << fs << " " << max_fs << endl;
      if (VERBOSE) cerr << min_ss << " " << ss << " " << max_ss << endl; 
      if (fs<min_fs||fs>max_fs) return;
      if (ss<min_ss||ss>max_ss) return;
      if (ss<min_ss) return;
      if (fs<min_fs) return;
      if (ss>max_ss) return;
      if (fs>max_fs) return;
      const T _fs=fs-min_fs;
      const T _ss=ss-min_ss;
      if (VERBOSE) cerr << "_fs = " << _fs << endl;
      if (VERBOSE) cerr << "_ss = " << _ss << endl;
      if (VERBOSE){
        cerr << fs_2_x << " " << ss_2_x << endl;
        cerr << fs_2_y << " " << ss_2_y << endl;
      }
      x=_fs*fs_2_x+_ss*ss_2_x+xoffset;
      y=_fs*fs_2_y+_ss*ss_2_y+yoffset;
      //cerr << "# "  << z << " " << zoffset << " " << z+zoffset << endl;
      z=_fs*fs_2_z+_ss*ss_2_z+zoffset+z;
      //z=zoffset;
    }
    void inline operator()(
        double x,double y,double z,
        uint32_t& fs, uint32_t& ss) const {
      const double d = (xoffset*n0+yoffset*n1+zoffset*n2)/(x*n0+y*n1+z*n2);
      x*=d;y*=d;z*=d; // now x,y,z truly is on the detector plane
      fs = floor(x*x_2_fs+y*y_2_fs+z*z_2_fs);
      ss = floor(x*x_2_ss+y*y_2_ss+z*z_2_ss);
      return;
      double on = x*x_2_fs+y*y_2_fs+z*z_2_fs;
      if (abs(on-1)>1e-8) cerr << "die matrix ist kaputt :(" << endl;
    }
    bool inline operator==(const linear_panel_transform& a
        )const{
      return
        (a.min_fs==min_fs)&&
        (a.max_fs==max_fs)&&
        (a.min_ss==min_ss)&&
        (a.max_ss==max_ss);
    }
    bool inline operator>(
        const linear_panel_transform& a
        )const{
      if    (min_fs>a.min_fs) return true;
      if    (min_fs<a.min_fs) return false;
      return min_ss>a.min_ss;
    }
    friend bool inline operator>(
        const linear_panel_transform& a,
        const array<uint32_t,2>& b
        ){
      if    (a.min_fs>b[0]) return true;
      if    (a.min_fs<b[0]) return false;
      return a.min_ss>b[1];
    }
    friend bool inline operator>(
        const array<uint32_t,2>& b,
        const linear_panel_transform& a
        ){
      if    (b[0]>a.min_fs) return true;
      if    (b[0]<a.min_fs) return false;
      return b[1]>a.min_ss;
    }
    bool inline operator<(
        const linear_panel_transform& a
        ){
      if (min_fs<a.min_fs) return true;
      if (min_fs>a.min_fs) return false;
      return min_ss<a.min_ss;
    }
    friend bool inline operator<(
        const linear_panel_transform& a,
        const array<uint32_t,2>& b
        ){
      //cout << "the correct fucking operator is in use (0)" << endl;
      if    (a.min_fs<b[0]) return true;
      if    (a.min_fs>b[0]) return false;
      return a.min_ss<b[1];
    }
    friend bool inline operator<(
        const array<uint32_t,2>& b,
        const linear_panel_transform& a
        ){
      //cout << "the correct fucking operator is in use (1)" << endl;
      if    (b[0]<a.min_fs) return true;
      if    (b[0]>a.min_fs) return false;
      return b[1]<a.min_ss;
    }
  };

  struct crystfel_geometry{
    vector<linear_panel_transform> transforms;
    unordered_map<string,size_t> name$idx;
    template<typename T>
    void inline operator()(
        const T& fs,
        const T& ss,
        const string& name,
        double & x,
        double & y,
        double & z) const {
      transforms.at(name$idx.at(name))(fs,ss,x,y,z);
    }
    template<typename T>
    void inline operator()(
        const T& fs,
        const T& ss,
        double & x,
        double & y,
        double & z) const {
      for (auto it=transforms.begin();it!=transforms.end();++it){
        if (it->is_valid(fs,ss)) (*it)(fs,ss,x,y,z);
        else continue;
        break;
      }
      return; // todo faster version
      const auto it = lower_bound(
          transforms.begin(),
          transforms.end(),
          array<uint32_t,2>{
          static_cast<uint32_t>(floor(fs)),
          static_cast<uint32_t>(floor(ss))});
      if (it==transforms.end()) return;
      (*it)(fs,ss,x,y,z);
    }
    inline bool is_valid(const   double& fs, const   double& ss) const {
      for (auto it=transforms.begin();it!=transforms.end();++it){
        if (it->is_valid(fs,ss)) return true;
      }
      return false; // todo faster version
      const auto it = lower_bound(
          transforms.begin(),
          transforms.end(),
          array<uint32_t,2>{
          static_cast<uint32_t>(floor(fs)),
          static_cast<uint32_t>(floor(ss))})-1;
      if (it==transforms.end()) return false;
      //cout << fs << " " << ss << endl;
      //cout << it->min_fs << " " << it->max_fs << " "
      //     << it->min_ss << " " << it->max_ss << endl;
      return it->is_valid(fs,ss);
      bool isvalid=false;
      for (auto it=transforms.begin();it!=transforms.end();++it){
        isvalid=isvalid||it->is_valid(fs,ss);
      }
      return isvalid;
    }
    vector<sector2d> get_sectors(){
      vector<sector2d> sectors;
      for (auto it=transforms.begin(); it!=transforms.end(); ++it){
        sectors.push_back({it->min_fs,it->max_fs,it->min_ss,it->max_ss});
      }
      return sectors;
    }
    linear_panel_transform& get_transform(const uint32_t& fs,const uint32_t ss){
      for (auto it=transforms.begin();it!=transforms.end();++it){
        if (it->is_valid(fs,ss)) return *it;
      }
      return *transforms.end(); // todo faster version
      return *lower_bound(transforms.begin(),transforms.end(),
                 array<uint32_t,2>{
                 static_cast<uint32_t>(floor(fs)),
                 static_cast<uint32_t>(floor(ss))});
    }
    void add(
        const uint32_t&   min_fs,const uint32_t&   max_fs,
        const uint32_t&   min_ss,const uint32_t&   max_ss,
        const   double&   fs_2_x,const   double&   ss_2_x,
        const   double&   fs_2_y,const   double&   ss_2_y,
        const   double&   fs_2_z,const   double&   ss_2_z,
        const   double& corner_x,const   double& corner_y,
        const   double& corner_z,const   size_t& n_panel,
        const   string& name){
      if (VERBOSE) {
        cerr << "add ["<<max_ss<<"]["<<max_fs<<"]"<< endl;
        cerr << fs_2_x << " " << ss_2_x << " " << corner_x << endl;
        cerr << fs_2_y << " " << ss_2_y << " " << corner_y << endl;
        cerr << fs_2_z << " " << ss_2_z << " " << corner_z << endl;
      }
      double x_2_fs,y_2_fs,z_2_fs,
             x_2_ss,y_2_ss,z_2_ss,
             x_2_on,y_2_on,z_2_on;
      get_inverse(
          fs_2_x,ss_2_x,corner_x,
          fs_2_y,ss_2_y,corner_y,
          fs_2_z,ss_2_z,corner_z,
          x_2_fs,y_2_fs,z_2_fs,
          x_2_ss,y_2_ss,z_2_ss,
          x_2_on,y_2_on,z_2_on);
      name$idx[name]=transforms.size();
      transforms.push_back({min_fs,max_fs,min_ss,max_ss,
                            fs_2_x,ss_2_x,corner_x,
                            fs_2_y,ss_2_y,corner_y,
                            fs_2_z,ss_2_z,corner_z,
                            x_2_fs,y_2_fs,z_2_fs,
                            x_2_ss,y_2_ss,z_2_ss,
                            x_2_on,y_2_on,z_2_on,
                            fs_2_y*ss_2_z-fs_2_z*ss_2_y,
                            fs_2_z*ss_2_x-fs_2_x*ss_2_z,
                            fs_2_x*ss_2_y-fs_2_y*ss_2_x,
                            n_panel});
      //sort(transforms.begin(),transforms.end());
    }
    void add(
        const uint32_t&   min_fs,const uint32_t&   max_fs,
        const uint32_t&   min_ss,const uint32_t&   max_ss,
        const   double&   fs_2_x,const   double&   ss_2_x,
        const   double&   fs_2_y,const   double&   ss_2_y,
        const   double& corner_x,const   double& corner_y,
        const   double& corner_z,const   size_t& n_panel,
        const   string& name){
      add(min_fs,max_fs,min_ss,max_ss,
          fs_2_x,ss_2_x,fs_2_y,ss_2_y,0.0,0.0,
          corner_x,corner_y,corner_z,n_panel,name);
    }
  };

  struct paneldescription{
    uint32_t min_fs = 0;
    uint32_t max_fs = 0;
    uint32_t min_ss = 0;
    uint32_t max_ss = 0;
    double x_2_fs = 0;
    double y_2_fs = 0;
    double z_2_fs = 0;
    double x_2_ss = 0;
    double y_2_ss = 0;
    double z_2_ss = 0;
    double corner_x = 0;
    double corner_y = 0;
    double corner_z = 0;
    double coffset = 0;
  };

  crystfel_geometry const inline get_crystfel_geometry(
      istream& detector_geometry_file,
      double coffset=0
      ){
    unordered_map<string,paneldescription> geom;
    vector<string> inputorder;
    size_t pos;
    double res=1;
    double clen=0;
    for (string line; getline(detector_geometry_file,line);){
      if (line[0]==';') continue;
      if (line.substr(0,3).compare("res")==0){
        pos=line.find("=",3);
        res = stod(line.substr(pos+1));
        if (VERBOSE) cerr << "res = " << res << endl;
        continue;
      }
      if (line.substr(0,4).compare("clen")==0){
        pos=line.find("=",4);
        char * end = nullptr;
        clen = strtod(line.substr(pos+1).c_str(),&end);
        if (end==nullptr) clen = 0;
        if (VERBOSE) cerr << "clen = " << clen << endl;
        continue;
      }
      if (line.substr(0,7).compare("coffset")==0){
        pos=line.find("=",7);
        char * end = nullptr;
        coffset = strtod(line.substr(pos+1).c_str(),&end);
        if (end==nullptr) coffset = 0;
        if (VERBOSE) cerr << "coffset = " << coffset << endl;
        continue;
      }
      pos = line.find('/'); 
      if (pos!=string::npos){
        string name = line.substr(0,pos);
        if (geom.count(name)==0) {
          inputorder.push_back(name);
        }
        if (VERBOSE) cerr << name << endl;
        if (line.substr(pos+1,6).compare("min_fs")==0){
          pos=line.find("=",pos+7);
          geom[name].min_fs=stoi(line.substr(pos+1));
          if (VERBOSE) cerr << "min_fs = " << geom[name].min_fs << endl;
          continue;
        }
        if (line.substr(pos+1,6).compare("min_ss")==0){
          pos=line.find("=",pos+7);
          geom[name].min_ss=stoi(line.substr(pos+1));
          if (VERBOSE) cerr << "min_ss = " << geom[name].min_ss << endl;
          continue;
        }
        if (line.substr(pos+1,6).compare("max_fs")==0){
          pos=line.find("=",pos+7);
          geom[name].max_fs=stoi(line.substr(pos+1));
          if (VERBOSE) cerr << "max_fs = " << geom[name].max_fs << endl;
          continue;
        }
        if (line.substr(pos+1,6).compare("max_ss")==0){
          pos=line.find("=",pos+7);
          geom[name].max_ss=stoi(line.substr(pos+1));
          if (VERBOSE) cerr << "max_ss = " << geom[name].max_ss << endl;
          continue;
        }
        if (line.substr(pos+1,2).compare("fs")==0){
          pos=line.find("=",pos+3);
          if (pos==string::npos){
            cerr << "wrong format?" << endl;
            cerr << line << endl;
            break;
          }
          size_t xpos = line.find("x",pos+1);
          size_t ypos = line.find("y",pos+1);
          if (xpos==string::npos||ypos==string::npos){
            cerr << "wrong format?" << endl;
            cerr << line << endl;
            break;
          }
          //while(xpos<line.size()) if (line[xpos]==' ')++xpos;else break;
          //while(ypos<line.size()) if (line[xpos]==' ')++ypos;else break;
          if (xpos<ypos){
            geom[name].x_2_fs=stod(line.substr(pos+1));
            geom[name].y_2_fs=stod(line.substr(xpos+1));
          }else{
            geom[name].y_2_fs=stod(line.substr(pos+1));
            geom[name].x_2_fs=stod(line.substr(ypos+1));
          }
          //cerr << line << endl;
          if (VERBOSE)
            cerr << "fs = " << geom[name].x_2_fs << "x" << " "
              << geom[name].y_2_fs << "y" << endl;
          continue;
        }
        if (line.substr(pos+1,2).compare("ss")==0){
          pos=line.find("=",pos+3);
          if (pos==string::npos){
            cerr << "wrong format? (ss)" << endl;
            cerr << line << endl;
            break;
          }
          size_t xpos = line.find("x",pos+1);
          size_t ypos = line.find("y",pos+1);
          if (xpos==string::npos||ypos==string::npos){
            cerr << "wrong format?" << endl;
            cerr << line << endl;
            break;
          }
          if (xpos<ypos){
            //cerr << "a" << endl;
            //cerr << line.substr(pos) << endl;
            geom[name].x_2_ss=stod(line.substr(pos+1));
            //cerr << line.substr(xpos) << endl;
            geom[name].y_2_ss=stod(line.substr(xpos+1));
          }else{
            geom[name].y_2_ss=stod(line.substr(pos+1));
            geom[name].x_2_ss=stod(line.substr(ypos+1));
          }
          //cerr << line << endl;
          if (VERBOSE)
            cerr << "ss = " << geom[name].x_2_ss << "x" << " " << 
              geom[name].y_2_ss << "y" << endl;
          continue;
        }
        if (line.substr(pos+1,8).compare("corner_x")==0){
          pos=line.find("=",pos+9);
          if (pos==string::npos){
            cerr << "wrong format? (corner_x)" << endl;
            cerr << line << endl;
            break;
          }
          geom[name].corner_x=stod(line.substr(pos+1));
          if (VERBOSE) cerr << "corner_x = " << geom[name].corner_x << endl;
          continue;
        }
        if (line.substr(pos+1,8).compare("corner_y")==0){
          pos=line.find("=",pos+8);
          if (pos==string::npos){
            cerr << "wrong format? (corner_y)" << endl;
            cerr << line << endl;
            break;
          }
          geom[name].corner_y=stod(line.substr(pos+1));
          if (VERBOSE) cerr << "corner_y = " << geom[name].corner_y << endl;
          continue;
        }
        if (line.substr(pos+1,7).compare("coffset")==0){
          pos=line.find("=",pos+8); 
          if (pos==string::npos){
            cerr << "wrong format? (coffset)" << endl;
            cerr << line << endl;
            break;
          }
          //cerr << "found coffset" << endl;
          //cerr << line.substr(pos+1) << endl;
          //cerr << stod(line.substr(pos+1)) << endl;
          geom[name].coffset=stod(line.substr(pos+1));
          if (VERBOSE) cerr << "coffset = " << geom[name].coffset << endl;
          continue;
        }
      }
    }
    coffset+=clen*res;
    crystfel_geometry f;
    size_t n_panel=0;
    for (auto it=inputorder.begin();it!=inputorder.end();++it){
      const string& name = *it;
      if (geom.count(name)==0) continue;
      const paneldescription& p = geom.at(*it);
      if (VERBOSE){
        cerr << name << " " << n_panel << endl;
        cerr << p.x_2_fs << " " << p.y_2_fs << endl;
        cerr << p.x_2_ss << " " << p.y_2_ss << endl;
        //cerr << " inverted= " << endl;
        //cerr << it->second.fs_2_x << " " << it->second.ss_2_x << endl;
        //cerr << it->second.fs_2_y << " " << it->second.ss_2_y << endl;
      }
      const double rdet = 1.0/(p.x_2_fs*p.y_2_ss-p.y_2_fs*p.x_2_ss);
      const double fs_2_x =  rdet*p.y_2_ss;
      const double ss_2_x = -rdet*p.y_2_fs;
      const double fs_2_y = -rdet*p.x_2_ss;
      const double ss_2_y =  rdet*p.x_2_fs;
      f.add(
          p.min_fs,p.max_fs,
          p.min_ss,p.max_ss,
          fs_2_x,ss_2_x,fs_2_y,ss_2_y,
          p.corner_x,p.corner_y,
          p.coffset+coffset,
          n_panel++,name);
      if (VERBOSE) {
        cerr << p.coffset << " " << res << " " << coffset << endl;
        cerr << p.min_fs << " " << p.max_fs << " "
             << p.min_ss << " " << p.max_ss << endl;
      }
    }
    //sort(f.transforms.begin(),f.transforms.end());
    return f;
  }

  struct panel{
    matrix<double,3,2> D = zeros_matrix<double>(3,3);
    matrix<double,3,1> o = zeros_matrix<double>(3,1);
    size_t nfs,nss,d;
    tuple<size_t,size_t> const inline operator()(const size_t& n) const {
      const size_t i = n-d;
      return {i%nfs,i/nfs};
    }
    size_t const inline operator()(const size_t& fs, const size_t& ss) const {
      return d+nfs*ss+fs;
    }
    matrix<double,3,1> const inline operator()(
        const matrix<double,2,1>& fsss) const {
      return D*fsss+o;
    }
    matrix<double,2,1> const inline operator()(
        matrix<double,3,1>& w) const {
      const matrix<double,3,3> M =
        inv(matrix<double,3,3>{
            D(0,0),D(0,1),w(0),
            D(1,0),D(1,1),w(1),
            D(2,0),D(2,1),w(2)
          });
      const matrix<double,3,1> fssse = -M*o;
      w*=-fssse(2);
      return matrix<double,2,1>{fssse(0),fssse(1)};
    }
    bool inline isvalid(const size_t j) const {
      if (j<d)         return false;
      if (j<d+nfs*nss) return true;
      return false;
    }
    bool inline isvalid(const size_t& fs,const size_t& ss) const {
      if (fs<nfs&&ss<nss) return true;
      return false;
    }
    bool inline isvalid(const matrix<double,2,1>& fsss) const {
      return isvalid(size_t(fsss(0)),size_t(fsss(1)));
    }
  };

  struct geometry{
    vector<panel> panels;
    size_t num_pixels=0;
    std::tuple<matrix<double,2,1>,const panel&> const inline map_to_fsss(
        matrix<double,3,1>& w
        ) const {
      for (auto it=panels.begin();it!=panels.end();++it){
        const matrix<double,2,1> fsss = (*it)(w);
        if ( !it->isvalid(fsss) ) continue;
        return {fsss,*it};
      }
      throw std::out_of_range("fs ss out of range for all panels");
    }
    matrix<double,3,1> const inline map_to_xyz(
        const panel& p,
        const matrix<double,2,1>& fsss) const {
      return p(fsss);
    }
    matrix<double,3,1> const inline map_to_xyz(
        const panel&  p,
        const size_t& i, // pointer to start of panel array
        const size_t& j) const {
      const double fs  = (j-i)%p.nfs;
      const double ss  = (j-i)/p.nfs;
      return p(matrix<double,2,1>{fs,ss});
    }
    size_t inline interpol(
        const size_t& ok,
        const size_t& olo,
        const size_t& ohi,
        const size_t& lo,
        const size_t& hi
        ) const {
      auto lm          = long_mul(size_t(ok-olo),hi-lo);
      const size_t n   = clz(get<0>(lm));
      const size_t m   = digits<size_t>()-n;
      const size_t den = (ohi-olo)>>m;
      const size_t nom = (get<0>(lm)<<n)+(get<1>(lm)>>m);
      return lo+nom/den;
    }
    inline const panel&
    get_panel(
        const size_t& j
        ) const {
      size_t lo = 0, hi = panels.size()-1;
      while(true){
        //cerr << "in scary get_panel loop" << endl;
        //cerr << lo << " " << hi << endl;
        if ((hi<lo)||(hi>panels.size()-1))
          throw std::out_of_range("index ouf of range, there is no such panel");
        if (lo==hi) return panels[lo];
        size_t mi = interpol(
            j,
            panels[lo].d,
            panels[hi].d+panels[hi].nfs*panels[hi].nss,
            lo,hi);
        if (hi-lo==1) {
          if ((panels[lo].d<=j)
            &&(j<(panels[lo].d+panels[lo].nfs*panels[lo].nss)))
            return panels[lo];
          if ((panels[hi].d<=j)
            &&(j<(panels[hi].d+panels[hi].nfs*panels[hi].nss)))
            return panels[hi];
          throw std::out_of_range("index ouf of range, there is no such panel");
        }
        if (hi-lo>1) {
          if (mi==lo) ++mi;
          if (mi==hi) --mi;
        }
        if (panels[mi].d>j) hi=mi;
        else if (panels[mi].d+panels[mi].nfs*panels[mi].nss<=j) lo=mi;
        else return panels[mi];
      }
    }
    std::tuple<matrix<double,3,1>,const panel&> 
    const inline map_to_xyz(
        const size_t& j
        ) const {
      for (auto it = panels.begin();it!=panels.end();++it){
        if (j<it->d){ 
          auto fsss = (*it)(j);
          return {
            (*it)(matrix<double,2,1>{
                double(get<0>(fsss)),
                double(get<1>(fsss))
                }),
            *it
          };
        }
      }
      throw std::out_of_range("index ouf of range, there is no such panel"); 
    }
    size_t const inline operator()(
        const size_t& fs,
        const size_t& ss,
        const size_t& n) const {
      return panels[n](fs,ss);
    }
  };

  bool const inline read_geometry(istream& file,geometry& geom){
    do{
      panel p;
      if (!(file >> p.nfs   )) return false;
      if (!(file >> p.nss   )) return false;
      if (!(file >> p.D(0,0))) return false;
      if (!(file >> p.D(0,1))) return false;
      if (!(file >> p.o(0  ))) return false;
      if (!(file >> p.D(1,0))) return false;
      if (!(file >> p.D(1,1))) return false;
      if (!(file >> p.o(1  ))) return false;
      if (!(file >> p.D(2,0))) return false;
      if (!(file >> p.D(2,1))) return false;
      if (!(file >> p.o(2  ))) return false;
      p.d = geom.num_pixels;
      geom.panels.push_back(p);
      if (file.eof()) return true;
      file.ignore(numeric_limits<streamsize>::max(),'\n');
      if (file.eof()) return true;
      if (file.peek()=='#') return true;
      geom.num_pixels+=p.nfs*p.nss;
    }while(true);
    geom.panels.shrink_to_fit();
    return true;
  }
  
  bool const inline read_geometry_bin(istream& file,geometry& geom){
    uint64_t n;
    file.read(reinterpret_cast<char*>(&n),8);
    for (size_t i=0;i!=n;++i){
      struct panel p;
      uint64_t nfs,nss;
      if(!file.read(reinterpret_cast<char*>(&nfs     ),8)) return false;
      if(!file.read(reinterpret_cast<char*>(&nss     ),8)) return false;
      if(!file.read(reinterpret_cast<char*>(&p.D(0,0)),8)) return false;
      if(!file.read(reinterpret_cast<char*>(&p.D(0,1)),8)) return false;
      if(!file.read(reinterpret_cast<char*>(&p.o(0  )),8)) return false;
      if(!file.read(reinterpret_cast<char*>(&p.D(1,0)),8)) return false;
      if(!file.read(reinterpret_cast<char*>(&p.D(1,1)),8)) return false;
      if(!file.read(reinterpret_cast<char*>(&p.o(1  )),8)) return false;
      if(!file.read(reinterpret_cast<char*>(&p.D(2,0)),8)) return false;
      if(!file.read(reinterpret_cast<char*>(&p.D(2,1)),8)) return false;
      if(!file.read(reinterpret_cast<char*>(&p.o(2  )),8)) return false;
      p.nfs=nfs;
      p.nss=nss;
      p.d = geom.num_pixels;
      geom.panels.push_back(p);
      geom.num_pixels+=p.nfs*p.nss;
    }
    geom.panels.shrink_to_fit();
    return true;
  }

  template<class ifstream,class ofstream>
  const inline void geom_ascii2bin(ifstream& in,ofstream& out){
    vector<tuple<size_t,size_t,array<double,9>>> buffer;
    while (in&&(!in.eof())){
      tuple<size_t,size_t,array<double,9>> entry;
      if (!(in >> get<0>(entry))) break;
      if (!(in >> get<1>(entry))) break;
      for (size_t i=0;i!=9;++i) if (!(in >> (get<2>(entry))[i])) break;
      buffer.push_back(entry);
    }
    uint64_t n = buffer.size();
    out.write(reinterpret_cast<char*>(&n),8);
    for (auto it=buffer.begin();it!=buffer.end();++it){
      out.write(reinterpret_cast<char*>(&(get<0>(*it))),8);
      out.write(reinterpret_cast<char*>(&(get<1>(*it))),8);
      for (size_t i=0;i!=9;++i)
        out.write(reinterpret_cast<char*>(&((get<2>(*it))[i])),8);
    }
  }
}
#endif // GEOMETRY_H
