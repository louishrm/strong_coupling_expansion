#pragma once

#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>
#include <cerrno>
#include <cstring>
#include <stdexcept>
#include <string>

namespace sc_expansion {

  // Append one row to a CSV file, writing the header first if the file is empty.
  // Safe under concurrent appenders: takes an exclusive flock for the duration of the
  // size check + write, so the header is written exactly once even when many MPI jobs
  // finish simultaneously.
  inline void append_csv_row(std::string const &path, std::string const &header, std::string const &row) {
    int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_APPEND, 0644);
    if (fd < 0) { throw std::runtime_error("csv append: open failed for " + path + ": " + std::strerror(errno)); }

    if (::flock(fd, LOCK_EX) != 0) {
      ::close(fd);
      throw std::runtime_error("csv append: flock failed for " + path + ": " + std::strerror(errno));
    }

    off_t sz = ::lseek(fd, 0, SEEK_END);
    std::string buf;
    if (sz == 0) buf = header + "\n";
    buf += row + "\n";

    ssize_t n = ::write(fd, buf.data(), buf.size());
    int werr  = errno;

    ::flock(fd, LOCK_UN);
    ::close(fd);

    if (n < 0 || (size_t)n != buf.size()) { throw std::runtime_error("csv append: short write to " + path + ": " + std::strerror(werr)); }
  }

} // namespace sc_expansion
