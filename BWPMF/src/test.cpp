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
  PhiList &phi_list(*XPtr<PhiList>(Rphi));
  return print_list_of_list_index(phi_list);
}

//[[Rcpp::export]]
void print_history_index(SEXP Rhistory) {
  History &history(*XPtr<History>(Rhistory));
  return print_list_of_list_index(history.data);
}