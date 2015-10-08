#include <boost/algorithm/string.hpp>
#include <boost/progress.hpp>
#include <Rcpp.h>
#include "bwpmf.h"
#include "rcpp_serialization.h"

namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive& ar, History& history, const unsigned int version) {
  ar & history.cookie_size;
  ar & history.hostname_size;
  ar & history.userdata;
}

}
}

using namespace Rcpp;

Dictionary cookie_dict, hostname_dict;

//[[Rcpp::export]]
SEXP serialize_cookie() {
  return rcpp_serialize(cookie_dict, true, true);
}

//[[Rcpp::export]]
void deserialize_cookie(RawVector src) {
  rcpp_deserialize(cookie_dict, src, true, true);
}

//[[Rcpp::export]]
SEXP serialize_hostname() {
  return rcpp_serialize(hostname_dict, true, true);
}

//[[Rcpp::export]]
void deserialize_hostname(RawVector src) {
  rcpp_deserialize(hostname_dict, src, true, true);
}

void encode(const std::string& key, Dictionary& dict) {
  if (dict.find(key) == dict.end()) {
    size_t value = dict.size();
    dict.insert(std::make_pair(key, value));
  }
}

void encode_cookie(const std::string& cookie) {
  encode(cookie, cookie_dict);
}

void encode_hostname(const std::string& hostname) {
  encode(hostname, hostname_dict);
}

size_t query(const std::string& key, const Dictionary& dict) {
  auto itor = dict.find(key);
  if (itor == dict.end()) return NA_INTEGER;
  return itor->second;
}

SEXP query_vector(const CharacterVector& src, const Dictionary& dict) {
  NumericVector retval(src.size());
  for(int i = 0;i < src.size();i++) {
    retval[i] = query(CHAR(src[i]), dict);
  }
  return retval;
}

//[[Rcpp::export]]
SEXP query_cookie(CharacterVector cookie) {
  return query_vector(cookie, cookie_dict);
}

//[[Rcpp::export]]
SEXP query_hostname(CharacterVector hostname) {
  return query_vector(hostname, hostname_dict);
}

//[[Rcpp::export]]
SEXP count_cookie() {
  return wrap(cookie_dict.size());
}

//[[Rcpp::export]]
SEXP count_hostname() {
  return wrap(hostname_dict.size());
}

//[[Rcpp::export]]
void clean_cookie() {
  cookie_dict.clear();
  cookie_dict.reserve(0);
}

//[[Rcpp::export]]
void clean_hostname() {
  hostname_dict.clear();
  hostname_dict.reserve(0);
}

const std::string& first_comp(const std::string& src) {
  static std::string buf;
  auto delim = src.find((char) 3);
  // Rprintf("%zd\n", delim);
  if (delim == std::string::npos) {
    buf.resize(0);
  } else {
    buf.assign(src, 0, delim);
  }
  return buf;
}

const std::pair<std::string, int>& two_comp(const std::string& src) {
  static std::pair<std::string, int> retval;
  static std::string buf;
  auto delim = src.find((char) 3);
  buf.assign(src, delim + 1, src.size() - delim - 1);
  if (delim == std::string::npos) {
    retval.first.resize(0);
  } else {
    retval.first.assign(src, 0, delim);
    retval.second = std::atoi(buf.c_str());
  }
  return retval;
}

//[[Rcpp::export]]
void encode(const std::string& path, size_t user_visit_lower_bound = 0, double progress = 0) {
  std::ifstream input(path.c_str());
  static std::vector<std::string> buf1, buf2;
  std::shared_ptr<boost::progress_display> pb(NULL);
  if (progress > 0) pb.reset(new boost::progress_display(progress));
  for(std::string str ; std::getline(input, str);) {
    if (progress > 0) pb->operator++();
    boost::split(buf1, str, boost::is_any_of("\1"));
    if (buf1.size() < 2) throw std::logic_error("invalid data");
    std::string& cookie(buf1[0]);
    boost::split(buf2, buf1[1], boost::is_any_of("\2"));
    if (buf2.size() < 1) throw std::logic_error("invalid user data");
    if (buf2.size() < user_visit_lower_bound) continue;
    encode_cookie(cookie);
    for(const std::string& s : buf2) {
      const std::string& s2(first_comp(s));
      if (s2.size() > 0) {
        encode_hostname(s2);
      }
    }
  }
}


//[[Rcpp::export]]
SEXP encode_data(const std::string& path, double progress = 0) {
  std::ifstream input(path.c_str());
  static std::vector<std::string> buf1, buf2;
  std::shared_ptr<boost::progress_display> pb(NULL);
  if (progress > 0) pb.reset(new boost::progress_display(progress));
  XPtr<History> retval(new History());
  History& history(*retval);
  history.userdata.resize(cookie_dict.size());
  history.cookie_size = cookie_dict.size();
  history.hostname_size = hostname_dict.size();
  for(std::string str ; std::getline(input, str);) {
    if (progress > 0) pb->operator++();
    boost::split(buf1, str, boost::is_any_of("\1"));
    if (buf1.size() < 2) throw std::logic_error("invalid data");
    std::string& cookie(buf1[0]);
    auto itor_cookie = cookie_dict.find(cookie);
    if (itor_cookie == cookie_dict.end()) continue;
    UserData& user_data(history.userdata[itor_cookie->second]);
    if (user_data.size() > 0) throw std::logic_error("Duplicated cookie");
    boost::split(buf2, buf1[1], boost::is_any_of("\2"));
    if (buf2.size() < 1) throw std::logic_error("invalid user data");
//     std::vector<double> hostname_id;
//     std::vector<int> count;
    for(const std::string& s : buf2) {
      const std::pair<std::string, int>& comp(two_comp(s));
      if (comp.first.size() > 0) {
        auto itor_hostname = hostname_dict.find(comp.first);
        if (itor_hostname == hostname_dict.end()) throw std::logic_error("Unknown hostname");
//         hostname_id.push_back(itor_hostname->second);
//         count.push_back(comp.second);
        user_data.push_back(HostnameCount(itor_hostname->second, comp.second));
      }
    }
  }
  return retval;
}

//[[Rcpp::export]]
SEXP serialize_history(SEXP Rhistory) {
  XPtr<History> phistory(Rhistory);
  History& history(*phistory);
  return rcpp_serialize(history, true, true);
}

//[[Rcpp::export]]
SEXP deserialize_history(RawVector src) {
  XPtr<History> retval(new History());
  rcpp_deserialize(*retval, src, true, true);
  return retval;
}

//[[Rcpp::export]]
void print_history(SEXP Rhistory) {
  XPtr<History> phistory(Rhistory);
  const History& history(*phistory);
  for(size_t user = 0;user < history.userdata.size();user++) {
    const auto& user_data(history.userdata[user]);
    if (user_data.size() == 0) continue;
    Rprintf("user(%zu): ", user);
    for(const auto& element : user_data) {
      Rprintf("(item: %zu, count: %d) ", element.first, element.second);
    }
    Rprintf("\n");
  }
}


//[[Rcpp::export]]
NumericVector check_history(SEXP Rhistory) {
  XPtr<History> phistory(Rhistory);
  History& history(*phistory);
  size_t error = 0, count = 0;
#pragma omp parallel for reduction( + : error, count )
  for(size_t i = 0;i < history.userdata.size();i++) {
    UserData& user_data(history.userdata[i]);
    std::set<int> check_tmp;
    int local_error = 0, local_count = 0;
    for(const HostnameCount& x : user_data) {
      if (check_tmp.find(x.first) != check_tmp.end()) local_error += 1;
      check_tmp.insert(x.first);
      local_count += x.second;
    }
    error += local_error;
    count += local_count;
  }
  if (error > 0) throw std::logic_error("Duplicated hostname in an user data");
  return wrap(count);
}

//[[Rcpp::export]]
SEXP count_non_zero_of_history(SEXP Rhistory) {
  XPtr<History> phistory(Rhistory);
  History& history(*phistory);
  size_t non_zero_count = 0;
#pragma omp parallel for reduction( + : non_zero_count)   
  for(size_t i = 0;i < history.userdata.size();i++) {
    non_zero_count += history.userdata[i].size();
  }
  return wrap(non_zero_count);
}

//[[Rcpp::export]]
size_t count_cookie_history(SEXP Rhistory) {
  XPtr<History> phistory(Rhistory);
  History& history(*phistory);
  return history.cookie_size;
}

//[[Rcpp::export]]
size_t count_hostname_history(SEXP Rhistory) {
  XPtr<History> phistory(Rhistory);
  History& history(*phistory);
  return history.hostname_size;
}

//[[Rcpp::export]]
SEXP extract_history(SEXP Rhistory, NumericVector id) {
  XPtr<History> phistory(Rhistory);
  History& history(*phistory);
  std::sort(id.begin(), id.end());
  std::vector<size_t> group_size(history.userdata.size(), 0);
  // count size
#pragma omp parallel for
  for(size_t i = 0;i < history.userdata.size();i++) {
    group_size[i] = history.userdata[i].size();
  }
#ifdef NOISY_DEBUG
  Function("print")(wrap(group_size));
#endif
  // cumulative sum
  std::vector<size_t> cumulative_group_size(group_size);
  for(size_t i = 1;i < group_size.size();i++) {
    cumulative_group_size[i] = cumulative_group_size[i-1] + group_size[i];
  }
#ifdef NOISY_DEBUG
  Function("print")(wrap(cumulative_group_size));
#endif
  // calculate group of id
  std::vector<size_t> id_group(id.size(), 0);
#pragma omp parallel for
  for(size_t i = 0;i < id.size();i++) {
    auto j = std::upper_bound(cumulative_group_size.begin(), cumulative_group_size.end(), id[i] - 1);
    auto group_idx = j - cumulative_group_size.begin();
    id_group[i] = group_idx;
  }
  // pick out these samples
  XPtr<History> retval(new History);
  retval->cookie_size = history.cookie_size;
  retval->hostname_size = history.hostname_size;
  retval->userdata.resize(history.cookie_size);
  {
    size_t i = id.size();
    while(i > 0) {
      i--;
#ifdef NOISY_DEBUG
      Rprintf("i: %zu\n", i);
      Rprintf("id[i]: %zu\n", (size_t) id[i]);
      Rprintf("id_group[i]: %zu\n", id_group[i]);
      Rprintf("cumulative_group_size[id_group[i] - 1]): %zu\n", cumulative_group_size[id_group[i] - 1]);
#endif
      size_t inner_group_id = (id_group[i] > 0 ? id[i] - cumulative_group_size[id_group[i] - 1] : id[i]);
#ifdef NOISY_DEBUG
      Rprintf("inner_group_id: %zu\n", inner_group_id);
#endif
      if (inner_group_id > history.userdata[id_group[i]].size()) throw std::logic_error("Invalid size ( > size )");
      if (inner_group_id == 0) throw std::logic_error("Invalid size ( == 0 )");
      auto target = history.userdata[id_group[i]].begin() + (inner_group_id - 1);
      retval->userdata[id_group[i]].push_back(*target);
      history.userdata[id_group[i]].erase(target);
    }
  }
#pragma omp parallel for
  for(size_t user = 0;user < retval->cookie_size;user++) {
    auto& user_data(retval->userdata[user]);
    user_data.reserve(user_data.size());
  }
  return retval;
}