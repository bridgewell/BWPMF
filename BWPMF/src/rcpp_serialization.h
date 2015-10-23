#ifndef __RCPP_SERIALIZATION_H__
#define __RCPP_SERIALIZATION_H__

#include <iostream>
#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/null.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#include <Rcpp.h>

template<typename T>
std::string serialize(const T& m, bool is_binary, bool is_gzip) {
  std::stringstream os;
  {
    boost::iostreams::filtering_stream<boost::iostreams::output> f;
    if (is_gzip) f.push(boost::iostreams::gzip_compressor());
    f.push(os);
    if (is_binary) {
      boost::archive::binary_oarchive oa(f);
      oa << m;
    } else {
      boost::archive::text_oarchive oa(f);
      oa << m;
    }
  }
  return os.str();
}

template<typename T>
void deserialize(const std::string& path, T& target) {
  std::ifstream is(path.c_str());
  boost::iostreams::filtering_stream<boost::iostreams::input> f;
  f.push(boost::iostreams::gzip_decompressor());
  f.push(is);
  boost::archive::binary_iarchive ia(f);
  ia >> target;
}

template <typename T>
SEXP serialize(SEXP Rpath, const T& target) {
  const std::string path(Rcpp::as<std::string>(Rpath));
  std::ofstream os(path.c_str());
  boost::iostreams::filtering_stream<boost::iostreams::output> f;
  f.push(boost::iostreams::gzip_compressor());
  f.push(os);
  boost::archive::binary_oarchive oa(f);
  oa << target;
  return R_NilValue;
}

template<typename T>
SEXP rcpp_serialize(const T& m, bool is_binary, bool is_gzip) {
  const std::string src(serialize<T>(m, is_binary, is_gzip));
  Rcpp::RawVector retval(src.size());
  memcpy(&retval[0], src.c_str(), src.size());
  return retval;  
}

template<typename T>
void rcpp_deserialize(T& m, Rcpp::RawVector src, bool is_binary, bool is_gzip) {
  typedef boost::iostreams::basic_array_source<char> Device;
  boost::iostreams::stream<Device> stream((char *) &src[0], src.size());
  boost::iostreams::filtering_stream<boost::iostreams::input> f;
  if (is_gzip) f.push(boost::iostreams::gzip_decompressor());
  f.push(stream);
  if (is_binary) {
    boost::archive::binary_iarchive ia(f);
    ia >> m;
  } else {
    boost::archive::text_iarchive ia(f);
    ia >> m;
  }  
}

#endif //__RCPP_SERIALIZATION_H__