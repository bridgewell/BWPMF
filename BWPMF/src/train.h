#ifndef __TRAIN_H__
#define __TRAIN_H__

#include <memory>
#include <boost/serialization/split_member.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <Rcpp.h>
#include "bwpmf.h"

struct Phi {
  
  DTYPE *data;
  
  Phi() : data(new DTYPE[Param::K]) 
  { }
  
  ~Phi() { delete [] data; }

  template<class Archive>
  void serialize(Archive &ar, const unsigned int version) const {
    for(size_t i = 0;i < Param::K;i++) {
      ar & data[i];
    }
  }

};

class PhiOnDisk {
  
  enum Mode {
    read,
    write
  } mode;
  
  size_t buffer_size;
  
  Phi *buffer;
  
  std::string path;
  
  size_t current_position;
  
  size_t total_size;
  
  size_t read_size;
  
  std::shared_ptr<std::ifstream> ifs;
  
  std::shared_ptr<boost::archive::binary_iarchive> iar;
  
  std::shared_ptr<std::ofstream> ofs;
  
  std::shared_ptr<boost::archive::binary_oarchive> oar;
  
  void flush() {
    switch (mode) {
    case Mode::read: {
#ifdef NOISY_DDEBUG
      Rcpp::Rcout << __FILE__ << "(" << __LINE__ << ")" << std::endl;
#endif
      size_t i = 0;
      while(i < buffer_size & i + read_size < total_size) {
        Phi& phi(buffer[i]);
#ifdef NOISY_DDEBUG
        Rprintf("%zu(total: %zu)\n", i, total_size);
#endif
        *iar >> phi;
        i += 1;
      }
      read_size += i;
#ifdef NOISY_DDEBUG
      Rcpp::Rcout << __FILE__ << "(" << __LINE__ << ") read_size: " << read_size <<  std::endl;
#endif
      return;
    }
    case Mode::write: {
#ifdef NOISY_DDEBUG
      Rcpp::Rcout << __FILE__ << "(" << __LINE__ << ")" << std::endl;
#endif
      for(size_t i = 0;i < current_position;i++) {
        Phi& phi(buffer[i]);
        *oar << phi;
        total_size++;
      }
#ifdef NOISY_DDEBUG
      Rcpp::Rcout << __FILE__ << "(" << __LINE__ << ") total_size: " << total_size << std::endl;
#endif
      return;
    }
    }
  }
  
  void reset() {
    current_position = 0;
  }
  
  void check() {
    if (current_position == buffer_size) {
      flush();
      reset();
    }
  }
  
  void start_write() {
    mode = Mode::write;
    current_position = 0;
    total_size = 0;
    ofs.reset(new std::ofstream(path.c_str()));
    oar.reset(new boost::archive::binary_oarchive(*ofs));
    reset();
  }

  void end_write() {
    flush();
    reset();
    oar.reset();
    ofs->flush();
    ofs->close();
    ofs.reset();
  }

  void start_read() {
    mode = Mode::read;
    ifs.reset(new std::ifstream(path.c_str()));
    iar.reset(new boost::archive::binary_iarchive(*ifs));
    read_size = 0;
    flush();
    reset();
  }
  
  void end_read() {
    reset();
    iar.reset();
    ifs->close();
    ifs.reset();
  }

public:
  PhiOnDisk(const std::string _path, size_t _buffer_size = 10000) 
    : buffer_size(_buffer_size), buffer(new Phi[buffer_size]), path(_path),
      current_position(0), iar(NULL), oar(NULL), total_size(0), read_size(0)
  {  
#ifdef NOISY_DEBUG
    Rcpp::Rcout << "PhiOnDisk with K: " << Param::K << std::endl;
#endif
  }
  
  ~PhiOnDisk() {
    delete [] buffer;
  }
  
  struct WriteFlag {
    PhiOnDisk& p;
    WriteFlag(PhiOnDisk& _p) : p(_p) {
      p.start_write();
    }
    ~WriteFlag() {
      p.end_write();
    }
  };
  
  WriteFlag get_write_flag() {
    return WriteFlag(*this);
  }
  
  struct ReadFlag {
    PhiOnDisk& p;
    ReadFlag(PhiOnDisk& _p) : p(_p) {
      p.start_read();
    }
    ~ReadFlag() {
      p.end_read();
    }
  };
  
  ReadFlag get_read_flag() {
    return ReadFlag(*this);
  }
  
  Phi& get_write_target() {
    Phi& retval(buffer[current_position++]);
    check();
    return retval;
  }
  

  const Phi& get_read_target() {
    const Phi& retval(buffer[current_position++]);
    check();
    return retval;
  }
  
  const size_t get_total_size() const {
    return total_size;
  }
  
};

typedef ListOfList<Phi> PhiList;

#endif // __TRAIN_H__