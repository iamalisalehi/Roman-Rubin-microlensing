#include "util/parse.h"

#include <charconv>
#include <cstring>
#include <sstream>
#include <string>

// istream's skipws in the C locale.
static inline bool isWs(char c)
{
    return c == ' ' or c == '\t' or c == '\n' or c == '\v' or c == '\f' or c == '\r';
}

TextReader::TextReader(const char* path, size_t bufBytes)
    : fp_(nullptr), buf_(bufBytes < 4 ? 4 : bufBytes)
{
    open(path);
}

TextReader::~TextReader()
{
    if (fp_) std::fclose(fp_);
}

bool TextReader::open(const char* path)
{
    if (fp_) std::fclose(fp_);
    fp_ = std::fopen(path, "rb");
    pos_ = len_ = 0;
    return fp_ != nullptr;
}

bool TextReader::more()
{
    if (!fp_) return false;
    if (pos_ > 0) {                       // keep the unread tail, drop what was consumed
        std::memmove(buf_.data(), buf_.data() + pos_, len_ - pos_);
        len_ -= pos_;
        pos_ = 0;
    }
    if (len_ == buf_.size()) buf_.resize(2 * buf_.size());   // one token longer than the buffer
    const size_t got = std::fread(buf_.data() + len_, 1, buf_.size() - len_, fp_);
    len_ += got;
    return got > 0;
}

void TextReader::skipLine()
{
    for (;;) {
        const void* nl = pos_ < len_ ? std::memchr(buf_.data() + pos_, '\n', len_ - pos_) : nullptr;
        if (nl) {
            pos_ = static_cast<const char*>(nl) - buf_.data() + 1;
            return;
        }
        pos_ = len_;
        if (!more()) return;
    }
}

bool TextReader::next(double& out)
{
    for (;;) {                            // skip whitespace
        while (pos_ < len_ and isWs(buf_[pos_])) ++pos_;
        if (pos_ < len_) break;
        if (!more()) return false;
    }
    size_t e = pos_;                      // the token: the maximal run of non-whitespace
    for (;;) {
        while (e < len_ and !isWs(buf_[e])) ++e;
        if (e < len_) break;
        const size_t tokLen = e - pos_;
        if (!more()) break;               // end of file ends the token
        e = pos_ + tokLen;                // more() moved the token to the front
    }
    const char* p = buf_.data() + pos_;
    const char* end = buf_.data() + e;

    // Fast path: a plain decimal token that from_chars consumes whole.
    const char c0 = *p;
    bool plain = (c0 >= '0' and c0 <= '9') or c0 == '-' or c0 == '.';
    for (const char* q = p; plain and q < end; ++q)
        plain = (*q >= '0' and *q <= '9') or *q == '.' or *q == 'e' or *q == 'E' or *q == '+' or *q == '-';
    if (plain) {
        double v;
        const auto r = std::from_chars(p, end, v, std::chars_format::general);
        if (r.ec == std::errc{} and r.ptr == end) {
            out = v;
            pos_ = e;
            return true;
        }
    }

    // Slow path: what `>>` does, and it leaves the rest of the token unread, as `>>` does.
    std::istringstream ss(std::string(p, end));
    double x = 0.0;
    const bool ok = static_cast<bool>(ss >> x);
    out = x;
    if (!ok) {
        pos_ = e;
        return false;
    }
    pos_ += ss.eof() ? size_t(e - pos_) : static_cast<size_t>(ss.tellg());
    return true;
}
