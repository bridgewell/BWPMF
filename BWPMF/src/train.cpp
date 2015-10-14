#include "stdafx.h"

using namespace Rcpp;

RCPP_EXPOSED_CLASS(Model)

// typedef std::vector<double> Phi;
// typedef std::vector< std::pair<size_t, Phi> > UserDataPhi;
// typedef std::vector<UserDataPhi> UserItemPhi;
// static UserItemPhi phi;
// 
// typedef std::vector< std::vector< std::pair<size_t, size_t> > > ItemInvertedIndex;
// 
// // [[Rcpp::export]]
// SEXP compute_inverted_index(SEXP Rmodel, SEXP Rhistory) {
//   Model* pmodel(as<Model*>(Rmodel));
//   Model& model(*pmodel);
//   XPtr<History> phistory(Rhistory);
//   const History& history(*phistory);
//   XPtr<ItemInvertedIndex> retval(new ItemInvertedIndex(model.item_param.size()));
//   ItemInvertedIndex& item_inverted_index(*retval);
//   for(size_t user = 0;user < history.userdata.size();user++) {
//     const UserData& user_data(history.userdata[user]);
//     for(size_t j = 0;j < user_data.size();j++) {
//       const auto& element(user_data[j]);
//       const size_t item = element.first;
//       item_inverted_index[item].push_back(std::make_pair(user, j));
//     }
//   }
// #pragma omp parallel for
//   for(size_t item = 0;item < model.item_param.size();item++) {
//     auto& v(item_inverted_index[item]);
//     v.reserve(v.size());
//   }
//   return retval;
// }
// 

//[[Rcpp::export]]
SEXP init_phi(SEXP Rmodel, SEXP Rhistory) {
  Model* pmodel(as<Model*>(Rmodel));
  Model& model(*pmodel);
  XPtr<History> phistory(Rhistory);
  History& history(*phistory);
  return XPtr<PhiList>(new PhiList(history.data.get_index(), history.data.get_index_size(), true));
}

//[[Rcpp::export]]
void train_once(SEXP Rmodel, SEXP Rhistory, SEXP Rtesting_history, SEXP Rphi, Function logger) {
  Model *pmodel(as<Model*>(Rmodel));
  Model& model(*pmodel);
  XPtr<History> phistory(Rhistory), ptesting_history(Rtesting_history);
  History &history(*phistory), &testing_history(*ptesting_history);
  if (model.user_size != history.user_size) throw std::invalid_argument("user_size is inconsistent");
  XPtr<PhiList> pphi_list(Rphi);
  PhiList& phi_list(*pphi_list);
  if (phi_list.get_index_size() != model.user_size) throw std::invalid_argument("index_size of phi_list is inconsistent");
  const int K(Param::K);
  static std::vector<double> user_sum, item_sum;
  user_sum.resize(K);
  user_sum.shrink_to_fit();
  item_sum.resize(K);
  item_sum.shrink_to_fit();
#pragma omp parallel
  {
    std::vector<double> local_item_sum(K, 0.0), local_user_sum(K, 0.0);
#pragma omp single
    logger(Rf_mkString("Calculating phi..."));
#pragma omp for
    for(size_t user = 0;user < history.user_size;user++) {
#ifdef NOISY_DEBUG
      Rcout << __FILE__ << "(" << __LINE__ << ")" << std::endl;
#endif
      // Phi *pphi_start = phi_list(user), *pphi_end = phi_list(user + 1);
      auto pphi_range = phi_list.range(user);
#ifdef NOISY_DEBUG
      Rcout << __FILE__ << "(" << __LINE__ << ")" << std::endl;
#endif
      const ItemCount *item_count = history.data(user);
#ifdef NOISY_DEBUG
      if (history.data.size(user) != phi_list.size(user)) throw std::logic_error(
        boost::str(boost::format("Inconsistent history size(%1%) and phi size(%2%)") % history.data.size(user) % phi_list.size(user))
        );
#endif
      for(Phi *pphi = pphi_range.first; pphi != pphi_range.second;pphi++) {
        size_t item = item_count->item;
#ifdef NOISY_DEBUG
        Rprintf("user: %zu item: %zu \n", user, item);
#endif
        Phi& phi(*pphi);
        Param& user_param(model.user_param[user]), item_param(model.item_param[item]);
        for(int k = 0;k < K;k++) {
          phi.data[k] = exp(Rf_digamma(user_param.shp1[k]) - log(user_param.rte1[k]) + Rf_digamma(item_param.shp1[k]) - log(item_param.rte1[k]));
        }
        double denom = std::accumulate(phi.data, phi.data + K, 0.0);
        std::transform(phi.data, phi.data + K, phi.data, [&denom](const double input) {
          return input / denom;
        });
        item_count++;
      }
    }

#pragma omp single
    logger(Rf_mkString("Updating user parameters..."));
    
    std::fill(local_item_sum.begin(), local_item_sum.end(), 0.0);
#pragma omp for
    for(size_t item = 0;item < model.item_size;item++) {
      Param& item_param(model.item_param[item]);
      for(int k = 0;k < K;k++) {
        local_item_sum[k] += item_param.shp1[k] / item_param.rte1[k];
      }
    }
    for(int k = 0;k < K;k++) {
#pragma omp atomic
      item_sum[k] += local_item_sum[k];
    }
#pragma omp for
    for(size_t user = 0;user < history.user_size;user++) {
      Param& user_param(model.user_param[user]);
      std::fill(user_param.shp1, user_param.shp1 + K, model.prior.a1);
      std::transform(item_sum.begin(), item_sum.end(), user_param.rte1, [&user_param](const double input) {
        return input + user_param.shp2 / user_param.rte2;
      });
      const auto range = history.data.range(user);
      // const ItemCount *start = history.data(user), *end = history.data(user + 1);
      const Phi* pphi = phi_list(user);
      for(const ItemCount *pitem_count = range.first; pitem_count != range.second;pitem_count++) {
        const size_t item = pitem_count->item;
        const int y = pitem_count->count;
        for(int k = 0;k < K;k++) {
          user_param.shp1[k] += y * pphi->data[k];
        }
        pphi++;
      }
      user_param.rte2 = model.prior.a2 / model.prior.b2;
      for(int k = 0;k < K;k++) {
        user_param.rte2 += user_param.shp1[k] / user_param.rte1[k];
      }
    }
    
#pragma omp single
    logger(Rf_mkString("Updating item parameters..."));
    
    std::fill(local_user_sum.begin(), local_user_sum.end(), 0.0);
#pragma omp for
    for(size_t user = 0;user < model.user_size;user++) {
      Param& user_param(model.user_param[user]);
      for(int k = 0;k < K;k++) {
        local_user_sum[k] += user_param.shp1[k] / user_param.rte1[k];
      }
    }
    for(int k = 0;k < K;k++) {
#pragma omp atomic
      user_sum[k] += local_user_sum[k];
    }
#pragma omp for
    for(size_t item = 0;item < model.item_size;item++) {
      Param& item_param(model.item_param[item]);
      std::fill(item_param.shp1, item_param.shp1 + K, model.prior.c1);
      std::transform(user_sum.begin(), user_sum.end(), item_param.rte1, [&item_param](const double input) {
        return input + item_param.shp2 / item_param.rte2;
      });
    }
#pragma omp for
    for(size_t user = 0;user < history.user_size;user++) {
      auto range = history.data.range(user);
      // const ItemCount *start = history.data(user), *end = history.data(user + 1);
      const Phi* pphi = phi_list(user);
      for(const ItemCount *pitem_count = range.first; pitem_count != range.second;pitem_count++) {
        const size_t item = pitem_count->item;
        Param& item_param(model.item_param[item]);
        const int y = pitem_count->count;
        for(int k = 0;k < K;k++) {
          double tmp = y * pphi->data[k];
#pragma omp atomic
          item_param.shp1[k] += tmp;
        }
        pphi++;
      }
    }
#pragma omp for
    for(size_t item = 0;item < model.item_size;item++) {
      Param& item_param(model.item_param[item]);
      item_param.rte2 = model.prior.c2 / model.prior.d2;
      for(int k = 0;k < K;k++) {
        item_param.rte2 += item_param.shp1[k] / item_param.rte1[k];
      }
    }
    
  } // #pragma omp parallel
}

//[[Rcpp::export]]
double pmf_logloss(SEXP Rmodel, SEXP Rhistory) {
  Model* pmodel(as<Model*>(Rmodel));
  Model& model(*pmodel);
  XPtr<History> phistory(Rhistory);
  History& history(*phistory);
  static std::vector<double> user_sum, item_sum;
  user_sum.resize(model.K, 0.0);
  std::fill(user_sum.begin(), user_sum.end(), 0.0);
  item_sum.resize(model.K, 0.0);
  std::fill(item_sum.begin(), item_sum.end(), 0.0);
  double retval = 0.0;
#pragma omp parallel
  {
  std::vector<double> local_user_sum(model.K, 0.0), local_item_sum(model.K, 0.0);
  double local_retval = 0.0;
    // y log(lambda)
#pragma omp for
    for(size_t user = 0;user < history.user_size;user++) {
      const Param& user_param(model.user_param[user]);
      auto range = history.data.range(user);
      // const ItemCount *start = history.data(user), *end = history.data(user + 1);
      for(const ItemCount *item_count = range.first; item_count != range.second;item_count++) {
        const size_t item = item_count->item;
        const int y = item_count->count;
        const Param& item_param(model.item_param[item]);
        double lambda = 0.0;
        for(int k = 0;k < model.K;k++) {
          double user_score = user_param.shp1[k] / user_param.rte1[k];
          double item_score = item_param.shp1[k] / item_param.rte1[k];
          lambda += user_score * item_score;
        }
        local_retval += y * log(lambda);
      }
    }
    // sum(theat_{u,k})
#pragma omp for
    for(size_t user = 0;user < model.user_size;user++) {
      const auto& param(model.user_param[user]);
      for(int k = 0;k < model.K;k++) {
        local_user_sum[k] += param.shp1[k] / param.rte1[k];
      }
    }
    // sum(beta_{i,k})
#pragma omp for
    for(size_t item = 0;item < model.item_size;item++) {
      const auto& param(model.item_param[item]);
      for(int k = 0;k < model.K;k++) {
        local_item_sum[k] += param.shp1[k] / param.rte1[k];
      }
    }
    
    for(int k = 0;k < model.K;k++) {
#pragma omp atomic
      user_sum[k] += local_user_sum[k];
#pragma omp atomic
      item_sum[k] += local_item_sum[k];
    }
#pragma omp atomic
    retval += local_retval;
  }
  for(int k = 0;k < model.K;k++) {
    retval -= user_sum[k] * item_sum[k];
  }
  return retval;
}
// 
// //[[Rcpp::export]]
// double pmf_mae(SEXP Rmodel, SEXP Rhistory) {
//   Model* pmodel(as<Model*>(Rmodel));
//   const Model& model(*pmodel);
//   XPtr<History> phistory(Rhistory);
//   const History& history(*phistory);
//   std::vector<double> user_sum(model.K, 0.0), item_sum(model.K, 0.0);
//   double retval = 0.0, adjustment = 0.0;
// #pragma omp parallel
//   {
//   std::vector<double> local_user_sum(model.K, 0.0), local_item_sum(model.K, 0.0);
//   double local_retval = 0.0, local_adjustment = 0.0;
//     // y log(lambda)
// #pragma omp for
//     for(size_t user = 0;user < history.userdata.size();user++) {
//       const auto& user_data(history.userdata[user]);
//       const auto& user_param(model.user_param[user]);
//       for(const auto& element : user_data) {
//         const size_t item = element.first;
//         const auto& item_param(model.item_param[item]);
//         const int y = element.second;
//         double lambda = 0.0;
//         for(int k = 0;k < model.K;k++) {
//           double user_score = user_param.shp1[k] / user_param.rte1[k];
//           double item_score = item_param.shp1[k] / item_param.rte1[k];
//           lambda += user_score * item_score;
//         }
//         local_retval += abs(lambda - y);
//         local_adjustment += lambda;
//       }
//     }
//     // sum(theat_{u,k})
// #pragma omp for
//     for(size_t user = 0;user < model.user_param.size();user++) {
//       const auto& param(model.user_param[user]);
//       for(int k = 0;k < model.K;k++) {
//         local_user_sum[k] += param.shp1[k] / param.rte1[k];
//       }
//     }
//     // sum(beta_{i,k})
// #pragma omp for
//     for(size_t item = 0;item < model.item_param.size();item++) {
//       const auto& param(model.item_param[item]);
//       for(int k = 0;k < model.K;k++) {
//         local_item_sum[k] += param.shp1[k] / param.rte1[k];
//       }
//     }
//     
//     for(int k = 0;k < model.K;k++) {
// #pragma omp atomic
//       user_sum[k] += local_user_sum[k];
// #pragma omp atomic
//       item_sum[k] += local_item_sum[k];
//     }
// #pragma omp atomic
//     retval += local_retval;
// #pragma omp atomic
//     adjustment += local_adjustment;
//   }
//   for(int k = 0;k < model.K;k++) {
//     retval += user_sum[k] * item_sum[k];
//   }
//   return (retval - adjustment) / (((double)model.item_param.size()) * model.user_param.size());
// }
