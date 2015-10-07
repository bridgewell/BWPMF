#include <Rcpp.h>
#include "bwpmf.h"

using namespace Rcpp;

RCPP_EXPOSED_CLASS(Model)

typedef std::vector<double> Phi;
typedef std::vector< std::pair<size_t, Phi> > UserDataPhi;
typedef std::vector<UserDataPhi> UserItemPhi;
static UserItemPhi phi;

typedef std::vector< std::vector< std::pair<size_t, size_t> > > ItemInvertedIndex;

// [[Rcpp::export]]
SEXP compute_inverted_index(SEXP Rmodel, SEXP Rhistory) {
  Model* pmodel(as<Model*>(Rmodel));
  Model& model(*pmodel);
  XPtr<History> phistory(Rhistory);
  const History& history(*phistory);
  XPtr<ItemInvertedIndex> retval(new ItemInvertedIndex(model.item_param.size()));
  ItemInvertedIndex& item_inverted_index(*retval);
  for(size_t user = 0;user < history.userdata.size();user++) {
    const UserData& user_data(history.userdata[user]);
    for(size_t j = 0;j < user_data.size();j++) {
      const auto& element(user_data[j]);
      const size_t item = element.first;
      item_inverted_index[item].push_back(std::make_pair(user, j));
    }
  }
#pragma omp parallel for
  for(size_t item = 0;item < model.item_param.size();item++) {
    auto& v(item_inverted_index[item]);
    v.reserve(v.size());
  }
  return retval;
}

//[[Rcpp::export]]
void init_phi(SEXP Rmodel, SEXP Rtraining_history) {
  Model* pmodel(as<Model*>(Rmodel));
  Model& model(*pmodel);
  XPtr<History> ptraining_history(Rtraining_history);
  History& training_history(*ptraining_history);

  phi.resize(model.user_param.size());
#pragma omp parallel for
  for(size_t i = 0;i < training_history.userdata.size();i++) {
    phi[i].resize(training_history.userdata[i].size());
    phi[i].reserve(training_history.userdata[i].size());
    for(size_t j = 0;j < training_history.userdata[i].size();j++) {
      HostnameCount& user_data(training_history.userdata[i][j]);
      std::pair<size_t, Phi>& user_data_phi(phi[i][j]);
      user_data_phi.first = user_data.first;
      user_data_phi.second.resize(model.K);
    }
  }  
}

//[[Rcpp::export]]
void train_once(SEXP Rmodel, SEXP Rtraining_history, SEXP Rtesting_history, SEXP Ritem_inverted_index) {
  Model* pmodel(as<Model*>(Rmodel));
  Model& model(*pmodel);
  XPtr<History> ptraining_history(Rtraining_history), ptesting_history(Rtesting_history);
  History& training_history(*ptraining_history), testing_history(*ptesting_history);
  XPtr<ItemInvertedIndex> pitem_inverted_index(Ritem_inverted_index);
  ItemInvertedIndex& item_inverted_index(*pitem_inverted_index);

  if (training_history.userdata.size() != testing_history.userdata.size()) throw std::invalid_argument("Inconsistent training and testing data");
  if (training_history.userdata.size() != model.user_param.size()) throw std::invalid_argument("Inconsistent training and cookie dictionary");
  // calculate phi
#ifdef NOISY_DEBUG
  Rprintf("Calculating phi...\n");
#endif
#pragma omp parallel for
  for(size_t i = 0;i < training_history.userdata.size();i++) {
    for(size_t j = 0;j < training_history.userdata[i].size();j++) {
      HostnameCount& user_data(training_history.userdata[i][j]);
      std::pair<size_t, Phi>& user_data_phi(phi[i][j]);
      for(size_t k = 0; k < model.K; k++) {
        const Param& user_param(model.user_param[i]), item_param(model.item_param[user_data.first]);
        user_data_phi.second[k] = exp(Rf_digamma(user_param.shp1[k]) - log(user_param.rte1[k]) + Rf_digamma(item_param.shp1[k]) - log(item_param.rte1[k]));
      }
      double denom = std::accumulate(user_data_phi.second.begin(), user_data_phi.second.end(), 0.0);
      std::transform(user_data_phi.second.begin(), user_data_phi.second.end(), user_data_phi.second.begin(), [&denom](double input) {
        return input / denom;
      });
#ifdef NOISY_DDEBUG
      Rprintf("sum of phi: %f\n", std::accumulate(user_data_phi.second.begin(), user_data_phi.second.end(), 0.0));
#endif
    }
  }
  // update user
  {
#ifdef NOISY_DEBUG
    Rprintf("Updating user...\n");
#endif
    const double a1 = model.prior.a1, a2 = model.prior.a2, b2 = model.prior.b2;
    static std::vector<double> item_sum;
    item_sum.resize(model.K);
    std::fill(item_sum.begin(), item_sum.end(), 0.0);
#pragma omp parallel
    {
      std::vector<double> local_item_sum(model.K, 0.0);
#pragma omp for
      for(size_t item = 0;item < model.item_param.size();item++) {
        const auto& item_param(model.item_param[item]);
        for(int k = 0;k < model.K;k++) {
          local_item_sum[k] += item_param.shp1[k] / item_param.rte1[k];
        }
      }
      for(int k = 0;k < model.K;k++) {
#pragma omp atomic
        item_sum[k] += local_item_sum[k];
      }
    }
    
#pragma omp parallel for
    for(size_t i = 0;i < training_history.userdata.size();i++) {
      const auto& user_data(training_history.userdata[i]);
      auto& user_param(model.user_param[i]);
      for(int k = 0;k < model.K;k++) {
        user_param.shp1[k] = a1;
        user_param.rte1[k] = user_param.shp2 / user_param.rte2 + item_sum[k];
      }
      for(size_t j = 0;j < user_data.size();j++) {
        size_t item = user_data[j].first;
        const auto& item_param(model.item_param[item]);
        int y = user_data[j].second;
        for(int k = 0;k < model.K;k++) {
          user_param.shp1[k] += y * phi[i][j].second[k];
        }
      }
      user_param.rte2 = a2 / b2;
      for(int k = 0;k < model.K;k++) {
        user_param.rte2 += user_param.shp1[k] / user_param.rte1[k];
      }
    }
  }
  // update item
  {
#ifdef NOISY_DEBUG
    Rprintf("Updating item...\n");
#endif
    double c1 = model.prior.c1, c2 = model.prior.c2, d2 = model.prior.d2;
    static std::vector<double> user_sum;
    user_sum.resize(model.K);
    std::fill(user_sum.begin(), user_sum.end(), 0.0);
#pragma omp parallel
    {
      std::vector<double> local_user_sum(model.K, 0.0);
#pragma omp for
      for (size_t user = 0;user < model.user_param.size();user++) {
        const auto& user_param(model.user_param[user]);
        for(int k = 0;k < model.K;k++) {
          local_user_sum[k] += user_param.shp1[k] / user_param.rte1[k];
        }
      }
      for(int k = 0;k < model.K;k++) {
#pragma omp atomic
        user_sum[k] += local_user_sum[k];
      }
    }
#pragma omp parallel for
    for(size_t item = 0;item < item_inverted_index.size();item++) {
      const auto& index(item_inverted_index[item]);
      auto& item_param(model.item_param[item]);
      for(int k = 0;k < model.K;k++) {
        item_param.shp1[k] = c1;
        item_param.rte1[k] = item_param.shp2 / item_param.rte2 + user_sum[k];
      }
      for(const auto& position : index) {
        const size_t user = position.first;
        const size_t user_data_index = position.second;
        const int y = training_history.userdata[user][user_data_index].second;
#ifdef NOISY_DEBUG
        if (training_history[user][user_data_index].first != item) throw std::logic_error("The inverted index is invalid");
        if (phi[user][user_data_index].first != item) throw std::logic_error("The inverted index and phi are inconsistent");
#endif
        const auto& user_param(model.user_param[user]);
        for(int k = 0;k < model.K;k++) {
          item_param.shp1[k] += y * phi[user][user_data_index].second[k];
        }
      }
      item_param.rte2 = c2 / d2;
      for(int k = 0;k < model.K;k++) {
        item_param.rte2 += item_param.shp1[k] / item_param.rte1[k];
      }
    }
  }
}

//[[Rcpp::export]]
double pmf_logloss(SEXP Rmodel, SEXP Rhistory) {
  Model* pmodel(as<Model*>(Rmodel));
  const Model& model(*pmodel);
  XPtr<History> phistory(Rhistory);
  const History& history(*phistory);
  std::vector<double> user_sum(model.K, 0.0), item_sum(model.K, 0.0);
  double retval = 0.0;
#pragma omp parallel
  {
  std::vector<double> local_user_sum(model.K, 0.0), local_item_sum(model.K, 0.0);
  double local_retval = 0.0;
    // y log(lambda)
#pragma omp for
    for(size_t user = 0;user < history.userdata.size();user++) {
      const auto& user_data(history.userdata[user]);
      const auto& user_param(model.user_param[user]);
      for(const auto& element : user_data) {
        const size_t item = element.first;
        const int y = element.second;
        const auto& item_param(model.item_param[item]);
        double lambda = 0.0;
        for(int k = 0;k < model.K;k++) {
          double user_score = user_param.shp1[k] / user_param.rte1[k];
          double item_score = item_param.shp1[k] / item_param.rte1[k];
          lambda += user_score * item_score;
#ifdef NOISY_DEBUG
          Rprintf("(user: %zu item: %zu k: %d ) user_score: %f item_score: %f\n", user, item, k, user_score, item_score);
#endif
        }
        local_retval += y * log(lambda);
      }
    }
    // sum(theat_{u,k})
#pragma omp for
    for(size_t user = 0;user < model.user_param.size();user++) {
      const auto& param(model.user_param[user]);
      for(int k = 0;k < model.K;k++) {
        local_user_sum[k] += param.shp1[k] / param.rte1[k];
      }
    }
    // sum(beta_{i,k})
#pragma omp for
    for(size_t item = 0;item < model.item_param.size();item++) {
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

//[[Rcpp::export]]
double pmf_mae(SEXP Rmodel, SEXP Rhistory) {
  Model* pmodel(as<Model*>(Rmodel));
  const Model& model(*pmodel);
  XPtr<History> phistory(Rhistory);
  const History& history(*phistory);
  std::vector<double> user_sum(model.K, 0.0), item_sum(model.K, 0.0);
  double retval = 0.0, adjustment = 0.0;
#pragma omp parallel
  {
  std::vector<double> local_user_sum(model.K, 0.0), local_item_sum(model.K, 0.0);
  double local_retval = 0.0, local_adjustment = 0.0;
    // y log(lambda)
#pragma omp for
    for(size_t user = 0;user < history.userdata.size();user++) {
      const auto& user_data(history.userdata[user]);
      const auto& user_param(model.user_param[user]);
      for(const auto& element : user_data) {
        const size_t item = element.first;
        const auto& item_param(model.item_param[item]);
        const int y = element.second;
        double lambda = 0.0;
        for(int k = 0;k < model.K;k++) {
          double user_score = user_param.shp1[k] / user_param.rte1[k];
          double item_score = item_param.shp1[k] / item_param.rte1[k];
          lambda += user_score * item_score;
        }
        local_retval += abs(lambda - y);
        local_adjustment += lambda;
      }
    }
    // sum(theat_{u,k})
#pragma omp for
    for(size_t user = 0;user < model.user_param.size();user++) {
      const auto& param(model.user_param[user]);
      for(int k = 0;k < model.K;k++) {
        local_user_sum[k] += param.shp1[k] / param.rte1[k];
      }
    }
    // sum(beta_{i,k})
#pragma omp for
    for(size_t item = 0;item < model.item_param.size();item++) {
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
#pragma omp atomic
    adjustment += local_adjustment;
  }
  for(int k = 0;k < model.K;k++) {
    retval += user_sum[k] * item_sum[k];
  }
  return (retval - adjustment) / (((double)model.item_param.size()) * model.user_param.size());
}
