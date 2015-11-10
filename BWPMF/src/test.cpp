#include "stdafx.h"
#include <iterator>

std::ostream& operator<<(std::ostream& os, const ListOfList<int>& dt) {
  dt.print(os);
  return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& dt) {
  os << "\n";
  std::ostream_iterator<T> iot(os, ", ");
  std::copy(dt.begin(), dt.end(), iot);
  os << std::endl;
  return os;
}

using namespace Rcpp;

template<typename T>
void report_error(const T& a, const std::string& file = "", int line = 0) {
  Rcpp::Rcout << a;
  std::string msg(file);
  msg.append(" ");
  msg.append(std::to_string(line));
  throw std::logic_error(msg);
}

//[[Rcpp::export]]
void test_list_of_list() {
  
  typedef ListOfList<int> LOLI;
  std::vector< std::vector<int> > a;
  a.push_back({1, 2, 3});
  a.push_back({1, 2});
  a.push_back({2, 3, 4});
  {
    LOLI a2(a);
    {
      Rcout << "Testing constructor" << std::endl;
      if (a2.get_index_size() != 3) report_error(a2, __FILE__, __LINE__);
      if (a2.get_total_size() != 8) report_error(a2, __FILE__, __LINE__);
      if (a2.size(0) != 3 | a2.size(1) != 2 | a2.size(2) != 3) report_error(a2, __FILE__, __LINE__);
      const int *p2 = a2.operator()(0);
      if (p2[0] != 1 | p2[1] != 2 | p2[2] != 3 | p2[7] != 4) report_error(a2, __FILE__, __LINE__);
    }
    {
      Rcout << "Testing for_each" << std::endl;
      std::vector<size_t> retval(3, 0);
      for(size_t i = 0;i < 3;i++) {
        a2(i, [&retval, &i](const int& j) {
          Rcout << j << std::endl;
          retval[i] += j;
        });
      }
      if (retval[0] != 6 | retval[1] != 3 | retval[2] != 9) report_error(retval, __FILE__, __LINE__);
    }
    {
      Rcout << "Testing clean" << std::endl;
      a2.clean([](const int& j) {
        Rcout << j << std::endl;
        return j == 1;
      });
      if (a2.get_index_size() != 3 | a2.get_total_size() != 2) report_error(a2, __FILE__, __LINE__);
      if (a2.size(0) != 1 | a2.size(1) != 1 | a2.size(2) != 0) report_error(a2, __FILE__, __LINE__);
      const int* start = a2(0);
      for (size_t i = 0;i < a2.get_total_size();i++) {
        if (start[i] != 1) report_error(a2, __FILE__, __LINE__);
      }
    }
  }
  {
    LOLI a2(a);
    Rcout << "Testing clean (2)" << std::endl;
    a2.clean([](const int& j) {
      // Rcout << j << std::endl;
      return j < 3;
    });
    if (a2.get_index_size() != 3 | a2.get_total_size() != 5) report_error(a2, __FILE__, __LINE__);
    if (a2.size(0) != 2 | a2.size(1) != 2 | a2.size(2) != 1) report_error(a2, __FILE__, __LINE__);
    const int* start = a2(0);
    for (size_t i = 0;i < a2.get_total_size();i++) {
      if (start[i] > 2) report_error(a2, __FILE__, __LINE__);
    }
  }  
  
  
}

//[[Rcpp::export]]
void test_phi(SEXP Rphi, SEXP Rhistory) {
  {
    RObject phi(Rphi);
    const std::string storage(as<std::string>(phi.attr("storage")));
    if (storage.compare("memory") != 0) return;
  }
  PhiList &phi_list(*XPtr<PhiList>(Rphi));
  History &history(*XPtr<History>(Rhistory));
  if (phi_list.get_index_size() != history.data.get_index_size()) throw std::logic_error(
    boost::str(boost::format("index_size of phi (%1%) and history (%2%) are inconsistent") % phi_list.get_index_size() % history.data.get_index_size())
  );
  const size_t *phi_index = phi_list.get_index(), *history_index = history.data.get_index();
  for(size_t i = 0;i < phi_list.get_index_size();i++) {
    if (phi_index[i] != history_index[i]) throw std::logic_error(boost::str(
      boost::format("index value of phi (%1%) and history (%2%) are inconsistent ant position %3%)") % phi_index[i] % history_index[i] % i
    ));
  }
}

template<typename T>
void print_list_of_list_index(const ListOfList<T>& src) {
  std::ostream_iterator<size_t> it(Rcpp::Rcout, ", ");
  std::copy(src.get_index(), src.get_index() + src.get_index_size(), it);
  Rcpp::Rcout << std::endl;
}

//[[Rcpp::export]]
void print_phi_index(SEXP Rphi) {
  {
    RObject phi(Rphi);
    if (as<std::string>(phi.attr("storage")).compare("memory") != 0) return;
  }
  PhiList &phi_list(*XPtr<PhiList>(Rphi));
  return print_list_of_list_index(phi_list);
}

//[[Rcpp::export]]
void print_history_index(SEXP Rhistory) {
  History &history(*XPtr<History>(Rhistory));
  return print_list_of_list_index(history.data);
}

NumericMatrix dump_phi_memory(SEXP Rphi) {
  XPtr<PhiList> pphi_list(Rphi);
  const PhiList& phi_list(*pphi_list);
  size_t total_size = phi_list.get_total_size();
  NumericMatrix retval(total_size, Param::K);
  size_t j = 0;
  for(size_t user = 0;user < phi_list.get_index_size();user++) {
    const auto& user_range(phi_list.range(user));
    for(auto p = user_range.first;p != user_range.second;p++) {
      const Phi& phi(*p);
      for(int k = 0;k < Param::K;k++) {
        retval(j, k) = phi.data[k];
      }
      j++;
    }
  }
  return retval;
}

NumericMatrix dump_phi_disk(SEXP Rphi, SEXP Rhistory, int K) {
  XPtr<History> phistory(Rhistory);
  History &history(*phistory);
  XPtr<pPhiOnDiskVec> pphi_disk_vec(Rphi);
  pPhiOnDiskVec& phi_disk_vec(*pphi_disk_vec);
  size_t total_size = history.data.get_total_size();
  std::vector< std::vector< std::vector<double> > > retval_buffer;
  retval_buffer.resize(history.user_size);
  bool is_valid = true;
#pragma omp parallel
  {
#pragma omp master 
    {
      Rprintf("Checking threads...\n");
      if (as<int>(pphi_disk_vec.attr("threads")) != omp_get_num_threads()) {
        is_valid = false;
      }
    }
  }
  if (!is_valid) throw std::runtime_error("The threads of phi and openmp are inconsistent!");
#pragma omp parallel
  {
    size_t thread_id = omp_get_thread_num();
    PhiOnDisk& phi_disk(*pphi_disk_vec->operator[](thread_id).get());
    {
      auto read_flag(phi_disk.get_read_flag());
#pragma omp for
      for(size_t user = 0;user < history.user_size;user++) {
        auto& retval_element(retval_buffer[user]);
        auto item_range = history.data.range(user);
        for(const ItemCount *pitem_count = item_range.first; pitem_count != item_range.second;pitem_count++) {
          // size_t item = pitem_count->item;
          std::vector<double> element(K, 0);
          const Phi& phi(phi_disk.get_read_target());
          for(int k = 0;k < K;k++) {
            element[k] = phi.data[k];
          }
          retval_element.push_back(element);
        }
      }
    }
  }
  NumericMatrix retval(total_size, K);
  size_t i = 0;
  for(const auto& retval_element : retval_buffer) {
    for(const auto& element : retval_element) {
      for(int k = 0;k < K;k++) {
        retval(i, k) = element[k];
      }
      i++;
    }
  }
  return retval;
}

//[[Rcpp::export]]
NumericMatrix dump_phi(SEXP Rphi, SEXP Rhistory = R_NilValue, SEXP K = R_NilValue) {
  RObject phi(Rphi);
  const std::string storage(as<std::string>(phi.attr("storage")));
  if (storage.compare("memory") == 0) {
    return dump_phi_memory(Rphi);
  } else if (storage.compare("disk") == 0) {
    return dump_phi_disk(Rphi, Rhistory, as<int>(K));
  } else {
    throw std::invalid_argument("Unknown storage type");
  }
}

//[[Rcpp::export]]
SEXP test_phi_on_disk(const std::string& path, NumericMatrix value) {
  Param::set_K(value.ncol());
  PhiOnDisk phi_disk(path, 10);
  // write
  {
    auto write_flag(phi_disk.get_write_flag());
    for(size_t i = 0;i < value.nrow();i++) {
      Phi& phi(phi_disk.get_write_target());
      for(int k = 0;k < value.ncol();k++) {
        phi.data[k] = (float) value(i,k);
      }
    }
  }
  if (value.nrow() != phi_disk.get_total_size()) {
    Rcout << "value.nrow(): " << value.nrow() << " phi_disk.get_total_size(): " << phi_disk.get_total_size() << std::endl;
    throw std::logic_error("total size is incorrect");
  }
  NumericMatrix retval1(value.nrow(), value.ncol());
  // read
  {
    auto read_flag(phi_disk.get_read_flag());
    for(size_t i = 0;i < value.nrow();i++) {
      const Phi& phi(phi_disk.get_read_target());
      for(int k = 0;k < value.ncol();k++) {
        retval1(i, k) = phi.data[k];
      }
    }
  }
  // write again
  {
    auto write_flag(phi_disk.get_write_flag());
    for(size_t i = 0;i < value.nrow();i++) {
      Phi& phi(phi_disk.get_write_target());
      for(int k = 0;k < value.ncol();k++) {
        phi.data[k] = (float) value(i,k) + 1;
      }
    }
  }
  if (value.nrow() != phi_disk.get_total_size()) throw std::logic_error("total size is incorrect");
  NumericMatrix retval2(value.nrow(), value.ncol());
  // read again
  // read
  {
    auto read_flag(phi_disk.get_read_flag());
    for(size_t i = 0;i < value.nrow();i++) {
      const Phi& phi(phi_disk.get_read_target());
      for(int k = 0;k < value.ncol();k++) {
        retval2(i, k) = phi.data[k];
      }
    }
  }
  return List::create(Named("retval1") = retval1, Named("retval2") = retval2);
}