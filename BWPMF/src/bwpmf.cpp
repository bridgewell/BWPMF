#include "stdafx.h"

size_t Param::current_param_count = 0;

int Param::K = 0;

Model::Model() 
  : K(0), prior(), user_size(0), item_size(0), user_param(NULL), item_param(NULL)
  { }

Model::Model(const Model& m) 
  : K(m.K), prior(m.prior), user_size(m.user_size), item_size(m.item_size),
    user_param(new Param[m.user_size]), item_param(new Param[m.item_size])
  {
    std::copy(m.user_param, m.user_param + m.user_size, user_param);
    std::copy(m.item_param, m.item_param + m.item_size, item_param);
  }

void Model::operator=(const Model& m) {
  delete [] user_param;
  delete [] item_param;
  K = m.K;
  prior = m.prior;
  user_size = m.user_size;
  item_size = m.item_size;
  user_param = new Param[user_size];
  item_param = new Param[item_size];
  std::copy(m.user_param, m.user_param + m.user_size, user_param);
  std::copy(m.item_param, m.item_param + m.item_size, item_param);
}

Model::Model(const Prior& _prior, int _k, size_t _user_size, size_t _item_size)
  : K(_k), prior(_prior), user_size(_user_size), item_size(_item_size),
    user_param(new Param[user_size]), item_param(new Param[item_size])
  {
    Param::set_K(_k);
#pragma omp parallel
    {
      unsigned int seed = omp_get_thread_num() + (int) time(NULL);
#pragma omp for
      for(size_t i = 0;i < user_size;i++) {
        Param& param(user_param[i]);
        param.shp2 = prior.a2 + _k * prior.a1;
        param.rte2 = prior.a2 / prior.b2* (0.9 + rand_r(&seed) * 0.2 / RAND_MAX);
        for(int k = 0;k < _k;k++) {
          param.shp1[k] = prior.a1 * (0.9 + rand_r(&seed) * 0.2 / RAND_MAX);
          // param.rte1[k] = param.shp2 / param.rte2* (0.9 + rand_r(&seed) * 0.2 / RAND_MAX);
          param.rte1[k] = prior.b2 * (0.9 + rand_r(&seed) * 0.2 / RAND_MAX);
        }
      }
#pragma omp for
      for(size_t item = 0;item < item_size;item++) {
        Param& param(item_param[item]);
        param.shp2 = prior.c2 + _k * prior.c1;
        param.rte2 = prior.c2 / prior.d2 * (0.9 + rand_r(&seed) * 0.2 / RAND_MAX);
        for(int k = 0;k < _k;k++) {
          param.shp1[k] = prior.c1 * (0.9 + rand_r(&seed) * 0.2 / RAND_MAX);
          // param.rte1[k] = param.shp2 / param.rte2 * (0.9 + rand_r(&seed) * 0.2 / RAND_MAX);
          param.rte1[k] = prior.d2 * (0.9 + rand_r(&seed) * 0.2 / RAND_MAX);
        }
      }
    }
  }

Model::~Model() {
  delete [] item_param;
  delete [] user_param;
}

BOOST_SERIALIZATION_SPLIT_FREE(Model)

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
  for(int k = 0;k < Param::K;k++) {
    ar & param.rte1[k];
    ar & param.shp1[k];
  }
  ar & param.rte2;
  ar & param.shp2;
}

template<class Archive>
void save(Archive& ar, const Model& m, const unsigned int version) {
  ar & m.K;
  ar & m.prior;
  ar & m.user_size;
  for(size_t user = 0;user < m.user_size;user++) {
    ar & m.user_param[user];
  }
  ar & m.item_size;
  for(size_t item = 0;item < m.item_size;item++) {
    ar & m.item_param[item];
  }
}

template<class Archive>
void load(Archive& ar, Model& m, const unsigned int version) {
  ar & m.K;
  Param::set_K(m.K);
  ar & m.prior;
  ar & m.user_size;
  delete [] m.user_param;
  m.user_param = new Param[m.user_size];
  for(size_t user = 0;user < m.user_size;user++) {
    ar & m.user_param[user];
  }
  ar & m.item_size;
  delete [] m.item_param;
  m.item_param = new Param[m.item_size];
  for(size_t item = 0;item < m.item_size;item++) {
    ar & m.item_param[item];
  }
}

}
}

RCPP_EXPOSED_CLASS(Prior)
RCPP_EXPOSED_CLASS(Param)
RCPP_EXPOSED_CLASS(Model)

using namespace Rcpp;

double user_size(Model* m) {
  return m->user_size;
}

SEXP user_param(Model* m, double i) {
  return wrap(m->user_param[(size_t )i]);
}

double item_size(Model* m) {
  return m->item_size;
}

SEXP item_param(Model* m, double i) {
  return wrap(m->item_param[(size_t) i]);
}

void prior_show(Prior* p) {
  Rprintf("Prior: a1:%f a2:%f b2:%f c1:%f c2:%f d2:%f\n", 
          p->a1, p->a2, p->b2, p->c1, p->c2, p->d2);
}

void param_show(Param* p) {
  Rprintf("Param: \n\tshp1: ");
  for(int k = 0;k < Param::K;k++) {
    Rprintf("%f ", p->shp1[k]);
  }
  Rprintf("\n\tshp2: %f\n\trte1: ", p->shp2);
  for(int k = 0;k < Param::K;k++) {
    Rprintf("%f ", p->rte1[k]);
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

SEXP param_shp1(Param* param) {
  NumericVector retval(Param::K);
  for(int k = 0;k < Param::K;k++) {
    retval[k] = param->shp1[k];
  }
  return retval;
}

SEXP param_rte1(Param* param) {
  NumericVector retval(Param::K);
  for(int k = 0;k < Param::K;k++) {
    retval[k] = param->rte1[k];
  }
  return retval;
}

//[[Rcpp::export]]
void set_K(int K) {
  Param::set_K(K);
}

NumericMatrix model_export_user(Model* pmodel) {
  NumericMatrix retval(pmodel->user_size, pmodel->K);
#pragma omp parallel for
  for(size_t user = 0;user < pmodel->user_size;user++) {
    Param& user_param(pmodel->user_param[user]);
    for(int k = 0;k < pmodel->K;k++) {
      retval(user, k) = user_param.shp1[k] / user_param.rte1[k];
    }
  }
  return retval;
}

NumericMatrix model_export_item(Model* pmodel) {
  NumericMatrix retval(pmodel->item_size, pmodel->K);
#pragma omp parallel for
  for(size_t item = 0;item < pmodel->item_size;item++) {
    Param& item_param(pmodel->item_param[item]);
    for(int k = 0;k < pmodel->K;k++) {
      retval(item, k) = item_param.shp1[k] / item_param.rte1[k];
    }
  }
  return retval;
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
    .method("shp1", &param_shp1)
    .method("rte1", &param_rte1)
    .field_readonly("shp2", &Param::shp2)
    .field_readonly("rte2", &Param::rte2)
    .method("show", &param_show)
  ;
  
  class_<Model>("Model")
    .constructor<Prior,int,size_t,size_t>()
    .constructor<Model>()
    .field_readonly("K", &Model::K)
    .field("prior", &Model::prior)
    .method("user_size", &user_size)
    .method("user_param", &user_param)
    .method("item_size", &item_size)
    .method("item_param", &item_param)
    .method("serialize", &model_serialize)
    .method("deserialize", &model_deserialize)
    .method("export_user", &model_export_user)
    .method("export_item", &model_export_item)
  ;

}