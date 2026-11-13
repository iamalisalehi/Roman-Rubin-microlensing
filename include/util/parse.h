// A buffered whitespace-token number reader for the big text inputs (the CMD star lists).
#ifndef ROMAN_UTIL_PARSE_H
#define ROMAN_UTIL_PARSE_H

#include <cstddef>
#include <cstdio>
#include <vector>

// Reads doubles exactly as `std::ifstream >> double` does, but through a bounded buffer and
// std::from_chars. The fast path takes plain decimal tokens; anything else (nan, inf, a leading '+',
// hex, out of range, a token istream would split) is handed to an istringstream, so the result is
// today's code by construction.
class TextReader {
public:
    explicit TextReader(const char* path, size_t bufBytes = size_t(1) << 20);
    ~TextReader();
    TextReader(const TextReader&) = delete;
    TextReader& operator=(const TextReader&) = delete;

    bool open(const char* path);    // switch to another file (the buffer is reused)
    bool isOpen() const { return fp_ != nullptr; }
    void skipLine();            // consume through the next '\n' (std::getline)
    bool next(double& out);     // false at end of file or on a token `>>` would reject (out is still set as `>>` sets it)

private:
    bool more();                // compact to the front and read more; false when nothing was read
    std::FILE* fp_;
    std::vector<char> buf_;
    size_t pos_ = 0, len_ = 0;
};

// True if the fast from_chars path is compiled in (the compiler's libstdc++ has floating-point from_chars);
// false means every token is read through an istringstream: identical values, slower.
bool textReaderFastPath();

#endif // ROMAN_UTIL_PARSE_H
