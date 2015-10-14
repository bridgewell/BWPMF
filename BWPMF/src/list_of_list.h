#ifndef __LIST_OF_LIST_H__
#define __LIST_OF_LIST_H__

#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <boost/serialization/split_member.hpp>
#ifdef NOISY_DEBUG
#include <Rcpp.h>
#endif // NOISY_DEBUG

template<typename T>
class ListOfList {
  
  size_t total_size;
  size_t index_size;
  size_t *index;
  T *data;
  
  ListOfList(const ListOfList&);
  void operator=(const ListOfList&);
  
  ListOfList(size_t _total_size, size_t _index_size) 
    : total_size(_total_size), index_size(_index_size), 
      index(new size_t[_index_size + 1]), data(new T[_total_size])
  { }

public:
  
  ListOfList() : total_size(0), index_size(0), index(nullptr), data(nullptr) 
  { }
  
  ListOfList(const std::vector< std::vector<T> >& src) 
    : ListOfList(std::accumulate(src.begin(), src.end(), 0, [](const size_t& retval, const std::vector<T>& i) {
        return retval + i.size();
    }), src.size())
  {
    index[0] = 0;
    size_t counter = 0;
    for(size_t i = 0;i < index_size;i++) {
      index[i + 1] = index[i] + src[i].size();
      for(const auto& element : src[i]) {
        data[counter++] = element;
      }
    }
  }
  
  ListOfList(const std::vector<size_t>& _size)
    : ListOfList(std::accumulate(_size.begin(), _size.end(), 0), _size.size())
  {
    index[0] = 0;
    for(size_t i = 0;i < _size.size();i++) {
      index[i + 1] = index[i] + _size[i];
    }
  }
  
  ListOfList(const size_t* _size, size_t _index_size, bool diff = false)
    : ListOfList((diff ? _size[_index_size] : std::accumulate(_size, _size + _index_size, 0)), _index_size)
  {
    if (diff) {
      std::copy(_size, _size + _index_size + 1, index);
    } else {
      index[0] = 0;
      for(size_t i = 0;i < _index_size;i++) {
        index[i + 1] = index[i]  + _size[i];
      }
    }
  }
  
  ~ListOfList() {
    delete [] index;
    delete [] data;
  }

  T* operator()(size_t i) {
#ifdef CHECK_BOUNDARY
    if (i >= index_size) throw std::invalid_argument("i exceeds the index_size");
#endif
    return data + index[i];
  }
  
  std::pair<T*,T*> range(size_t i) {
#ifdef CHECK_BOUNDARY
    if (i >= index_size) throw std::invalid_argument("i exceeds the index_size");
#endif
    return std::make_pair(data + index[i], data + index[i+1]);
  }
  
  T& operator()(size_t i, size_t j) {
#ifdef CHECK_BOUNDARY
    if (i >= index_size) throw std::invalid_argument("i exceeds the index_size");
    if (index[i] + j >= index[i + 1]) throw std::invalid_argument("j exceeds the size of i-th list");
#endif
    return data[index[i] + j];
  }
  
  template<class UnaryFunction>
  void operator()(size_t i, UnaryFunction f) {
#ifdef CHECK_BOUNDARY
    if (i >= index_size) throw std::invalid_argument("i exceeds the index_size");
#endif
    T 
      *begin = data + index[i],
      *end = data + index[i + 1];
    std::for_each(begin, end, f);
  }
  
  const size_t get_total_size() const {
    return total_size;
  }
  
  const size_t get_index_size() const {
    return index_size;
  }
  
  const size_t* get_index() const {
    return index;
  }
  
  const size_t size(size_t i) const {
#ifdef CHECK_BOUNDARY
    if (i >= index_size) throw std::invalid_argument("i exceeds the index_size");
#endif
    return index[i+1] - index[i];
  }
  
  template<class UnaryOperator>
  void clean(UnaryOperator f) {
    size_t total_adj = 0;
    size_t k = 0;
    size_t index_begin = index[0];
    for(size_t i = 0;i < index_size;i++) {
      for(size_t j = index_begin;j < index[i+1];j++) {
        if (!f(data[j])) total_adj += 1;
        if (j > k) {
          data[k] = data[j];
        }
        if (f(data[k])) k++; 
      }
#ifdef NOISY_DEBUG
      Rprintf("%uz - %uz \n", index[i+1], total_adj);
#endif
      index_begin = index[i + 1];
      index[i + 1] -= total_adj;
    }
    total_size -= total_adj;
  }
  
  void print(std::ostream& os) const {
    os << "total_size: " << total_size << "\n";
    os << "index_size: " << index_size << "\n";
    os << "index:" << "\n\t";
    std::ostream_iterator<size_t> out_it1 ( os, "\n\t");
    std::copy(index, index + index_size + 1, out_it1);
    os << "\n";
    os << "data:" << "\n\t";
    std::ostream_iterator<T> out_it2 ( os, "\n\t" );
    std::copy(data, data + total_size, out_it2);
    os << std::endl;
  }
  
private:
  friend class boost::serialization::access;
  
  template<class Archive>
  void save(Archive &ar, const unsigned int version) const {
    ar & total_size;
    ar & index_size;
    for(size_t i = 0;i < index_size + 1;i++) {
      ar & index[i];
    }
    for(size_t i = 0;i < total_size;i++) {
      ar & data[i];
    }
  }
  
  template<class Archive>
  void load(Archive &ar, const unsigned int version) {
    delete [] index;
    delete [] data;
    ar & total_size;
    data = new T[total_size];
    ar & index_size;
    index = new size_t[index_size + 1];
    for(size_t i = 0;i < index_size + 1;i++) {
      ar & index[i];
    }
    for(size_t i = 0;i < total_size;i++) {
      ar & data[i];
    }
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()

};

#endif //__LIST_OF_LIST_H__