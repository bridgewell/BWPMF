#ifndef __BWPMF_H__
#define __BWPMF_H__

#ifdef NOISY_DDEBUG
#define NOISY_DEBUG
#endif

#include <string>
#include <vector>
#include <unordered_map>

typedef std::unordered_map<std::string, size_t> Dictionary;
typedef std::pair<size_t, int> HostnameCount;
typedef std::vector<HostnameCount> UserData;

struct History {
  
  std::vector<UserData> userdata;
  
  size_t cookie_size;
  
  size_t hostname_size;
  
};

typedef std::pair<size_t, size_t> UserHostname;
typedef std::pair<UserHostname, int> UserHostnameCount;

struct Prior {
  double a1, a2, b2, c1, c2, d2;
  Prior(double _a1, double _a2, double _b2, double _c1, double _c2, double _d2) : 
    a1(_a1), a2(_a2), b2(_b2), c1(_c1), c2(_c2), d2(_d2) { }
  Prior(const Prior& src) : a1(src.a1), a2(src.a2), b2(src.b2), c1(src.c1), c2(src.c2), d2(src.d2) { }
};

struct Param {
  std::vector<double> shp1, rte1;
  double shp2, rte2;
  Param() { }
  Param(int k) : shp1(k, 0.0), rte1(k, 0.0) { }
};

struct Model {
  
  int K;
  Prior prior;

  std::vector< Param > user_param;
  std::vector< Param > item_param;

  Model(const Prior& _prior, int _k, size_t user_size, size_t item_size);
  
  ~Model() { }

};

#endif // __BWPMF_H__