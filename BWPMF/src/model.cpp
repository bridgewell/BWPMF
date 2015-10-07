#include <Rcpp.h>
#include "bwpmf.h"
#include "rcpp_serialization.h"
#include "omp.h"

Model::Model(const Prior& _prior, int _k, size_t user_size, size_t item_size)
  : K(_k), user_param(user_size, Param(K)), item_param(item_size, Param(K)), prior(_prior)
  {
    const Prior& prior(this->prior);
#pragma omp parallel default(none) shared(_k) shared(prior)
    {
      int i;
      unsigned int seed = omp_get_thread_num() + (int) time(NULL);
#pragma omp for
      for(size_t i = 0;i < user_param.size();i++) {
        Param& param(user_param[i]);
        param.shp2 = prior.a2 + _k * prior.a1;
        param.rte2 = prior.a2 / prior.b2* (0.9 + rand_r(&seed) * 0.2 / RAND_MAX);
        for(int k = 0;k < _k;k++) {
          param.shp1[k] = prior.a1 * (0.9 + rand_r(&seed) * 0.2 / RAND_MAX);
          // param.rte1[k] = param.shp2 / param.rte2* (0.9 + rand_r(&seed) * 0.2 / RAND_MAX);
          param.rte1[k] = prior.b2 * (0.9 + rand_r(&seed) * 0.2 / RAND_MAX);
        }
        param.shp1.reserve(_k);
        param.rte1.reserve(_k);
      }
#pragma omp for
      for(size_t item = 0;item < item_param.size();item++) {
        Param& param(item_param[item]);
        param.shp2 = prior.c2 + _k * prior.c1;
        param.rte2 = prior.c2 / prior.d2 * (0.9 + rand_r(&seed) * 0.2 / RAND_MAX);
        for(int k = 0;k < _k;k++) {
          param.shp1[k] = prior.c1 * (0.9 + rand_r(&seed) * 0.2 / RAND_MAX);
          // param.rte1[k] = param.shp2 / param.rte2 * (0.9 + rand_r(&seed) * 0.2 / RAND_MAX);
          param.rte1[k] = prior.d2 * (0.9 + rand_r(&seed) * 0.2 / RAND_MAX);
        }
        param.shp1.reserve(_k);
        param.rte1.reserve(_k);
      }
    }
    user_param.reserve(user_param.size());
    item_param.reserve(item_param.size());
  }

namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive& ar, Prior& p, const unsigned int version) {
  ar & p.a1;
  ar & p.a2;
  ar & p.b2;
  ar & p.c1;
  ar & p.c2;
  ar & p.d2;
}

template<class Archive>
void serialize(Archive& ar, Param& param, const unsigned int version) {
  ar & param.rte1;
  ar & param.shp1;
  ar & param.rte2;
  ar & param.shp2;
}

template<class Archive>
void serialize(Archive& ar, Model& m, const unsigned int version) {
  ar & m.K;
  ar & m.prior;
  ar & m.user_param;
  ar & m.item_param;
}

}
}

RCPP_EXPOSED_CLASS(Prior)
RCPP_EXPOSED_CLASS(Param)
RCPP_EXPOSED_CLASS(Model)

using namespace Rcpp;

double user_size(Model* m) {
  return m->user_param.size();
}

SEXP user_param(Model* m, double i) {
  return wrap(m->user_param[i]);
}

double item_size(Model* m) {
  return m->item_param.size();
}

SEXP item_param(Model* m, double i) {
  return wrap(m->item_param[i]);
}

void prior_show(Prior* p) {
  Rprintf("Prior: a1:%f a2:%f b2:%f c1:%f c2:%f d2:%f\n", 
          p->a1, p->a2, p->b2, p->c1, p->c2, p->d2);
}

void param_show(Param* p) {
  Rprintf("Param: \n\tshp1: ");
  for(auto d : p->shp1) {
    Rprintf("%f ", d);
  }
  Rprintf("\n\tshp2: %f\n\trte1: ", p->shp2);
  for(auto d : p->rte1) {
    Rprintf("%f ", d);
  }
  Rprintf("\n\trte2: %f\n", p->rte2);
}

void model_serialize(Model* m, const std::string& path) {
  std::ofstream ofs(path.c_str());
  boost::iostreams::filtering_stream<boost::iostreams::output> f;
  f.push(boost::iostreams::gzip_compressor());
  f.push(ofs);
  boost::archive::binary_oarchive oa(f);
  oa << *m;
}

void model_deserialize(Model* m, const std::string& path) {
  std::ifstream ifs(path.c_str());
  boost::iostreams::filtering_stream<boost::iostreams::input> f;
  f.push(boost::iostreams::gzip_decompressor());
  f.push(ifs);
  boost::archive::binary_iarchive ia(f);
  ia >> *m;
}

RCPP_MODULE(model) {

  class_<Prior>("Prior")
    .constructor<double,double,double,double,double,double>()
    .field_readonly("a1", &Prior::a1)
    .field_readonly("a2", &Prior::a2)
    .field_readonly("b2", &Prior::b2)
    .field_readonly("c1", &Prior::c1)
    .field_readonly("c2", &Prior::c2)
    .field_readonly("d2", &Prior::d2)
    .method("show", &prior_show)
  ;
  
  class_<Param>("Param")
    .field_readonly("shp1", &Param::shp1)
    .field_readonly("rte1", &Param::rte1)
    .field_readonly("shp2", &Param::shp2)
    .field_readonly("rte2", &Param::rte2)
    .method("show", &param_show)
  ;
  
  class_<Model>("Model")
    .constructor<Prior,int,size_t,size_t>()
    .field_readonly("K", &Model::K)
    .field("prior", &Model::prior)
    .method("user_size", &user_size)
    .method("user_param", &user_param)
    .method("item_size", &item_size)
    .method("item_param", &item_param)
    .method("serialize", &model_serialize)
    .method("deserialize", &model_deserialize)
  ;

}