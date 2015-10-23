#ifndef __BWPMF_H__
#define __BWPMF_H__

#ifdef NOISY_DDEBUG
#define NOISY_DEBUG
#endif

#include <string>
#include <vector>
#include <unordered_map>
#include "list_of_list.h"

typedef float DTYPE;

// encoding data
typedef std::unordered_map<std::string, size_t> Dictionary;

extern Dictionary cookie_dict, hostname_dict;

// history
struct ItemCount {
  
  size_t item;
  
  int count;
  
  ItemCount() : item(0), count(0) { }
  
  ItemCount(size_t _item, int _count) : item(_item), count(_count) { }
  
  ItemCount(const ItemCount& src) : ItemCount() {
    this->operator=(src);
  }
  
  void operator=(const ItemCount& src) {
    item = src.item;
    count = src.count;
  }
  
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar & item;
    ar & count;
  }
};

struct History {
  
  size_t user_size, item_size;

  ListOfList<ItemCount> data;

  History() : user_size(0), item_size(0), data() { }
    
  History(const std::vector<std::vector<ItemCount> >& src, const size_t _item_size)
    : user_size(src.size()), item_size(_item_size), data(src)
    { }
  
  ~History() { }
  
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar & user_size;
    ar & item_size;
    ar & data;
  }

};
  
struct Prior {
  double a1, a2, b2, c1, c2, d2;
  Prior() { }
  Prior(double _a1, double _a2, double _b2, double _c1, double _c2, double _d2) : 
    a1(_a1), a2(_a2), b2(_b2), c1(_c1), c2(_c2), d2(_d2) { }
  Prior(const Prior& src) : a1(src.a1), a2(src.a2), b2(src.b2), c1(src.c1), c2(src.c2), d2(src.d2) { }
};

struct Param {
  
  static int K;
  
  static size_t current_param_count;
  
  static void set_K(int k) {
    if (k != K) {
      if (current_param_count != 0) throw std::logic_error("K should be set before initializing Param");
    }
    K = k;
  }
  
  DTYPE *shp1, *rte1;
  
  DTYPE shp2, rte2;
  
  Param() : shp1(new DTYPE[K]), rte1(new DTYPE[K]), shp2(0.0), rte2(0.0) {
    std::fill_n(shp1, K, 0.0);
    std::fill_n(rte1, K, 0.0);
    current_param_count++;
  }
  
  Param(const Param& src) : Param() {
    this->operator=(src);
  }
  
  void operator=(const Param& src) {
    std::copy(src.shp1, src.shp1 + K, shp1);
    std::copy(src.rte1, src.rte1 + K, rte1);
    shp2 = src.shp2;
    rte2 = src.rte2;
  }
  
  ~Param() {
    delete [] rte1;
    delete [] shp1;
  }
  
};

struct Model {
  
  int K;
  Prior prior;
  size_t user_size, item_size;
  Param *user_param, *item_param;

  Model();
  
  Model(const Prior& _prior, int _k, size_t user_size, size_t item_size);
  
  Model(const Model& m);
  
  void operator=(const Model& m);
  
  ~Model();

};

#endif // __BWPMF_H__