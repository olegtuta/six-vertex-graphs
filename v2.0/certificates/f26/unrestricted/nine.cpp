// Exact order-nine certificate: a lift of 169 even-root Gram blocks plus a
// new root-one block. Two sampling modes (cached partitions / fresh 5-sets).
#include "wide_integer.h"
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
using Big = WideInteger;
int positions9[9][9]{};
Big integer(const string &s) {
  Big v = 0;
  int start = s.size() && s[0] == '-' ? 1 : 0;
  for (int i = start; i < (int)s.size(); i++) {
    if (s[i] < '0' || s[i] > '9')
      throw runtime_error("Bad integer");
    v = v * 10 + s[i] - '0';
  }
  return start ? -v : v;
}
string decimal(Big x) {
  if (!x)
    return "0";
  bool negative = x < 0;
  if (negative)
    x = -x;
  string s;
  while (x) {
    s += char('0' + static_cast<int>(x % 10));
    x /= 10;
  }
  if (negative)
    s += '-';
  reverse(s.begin(), s.end());
  return s;
}
int flagcode(uint64_t mask, int root, unsigned others) {
  int v[5];
  v[0] = root;
  int n = 1;
  for (int i = 0; i < 9; i++)
    if ((others >> i) & 1)
      v[n++] = i;
  if (n != 5)
    throw runtime_error("Not four extension vertices");
  int code = 0, b = 0;
  for (int i = 0; i < 5; i++)
    for (int j = i + 1; j < 5; j++, b++)
      if ((mask >> positions9[v[i]][v[j]]) & 1)
        code |= 1 << b;
  return code;
}
int main(int argc, char **argv) {
  try {
    if (argc != 6)
      throw runtime_error(
          "root-factor eight-values nine-prefix mode output required");
    int bit = 0;
    for (int i = 0; i < 9; i++)
      for (int j = i + 1; j < 9; j++)
        positions9[i][j] = positions9[j][i] = bit++;
    string mode = argv[4];
    if (mode != "cached" && mode != "independent")
      throw runtime_error("Unknown partition mode");
    ifstream in(argv[1]);
    int d;
    int64_t D;
    in >> d >> D;
    if (d != 90 || D < 1 || D > (1ll << 42))
      throw runtime_error("Invalid root Gram resource");
    array<int, 1024> classifier;
    for (int &c : classifier) {
      in >> c;
      if (c < 0 || c >= 90)
        throw runtime_error("Invalid root classifier");
    }
    vector<int64_t> F(8100);
    for (auto &v : F) {
      in >> v;
      if (v <= -(1ll << 47) || v >= (1ll << 47))
        throw runtime_error("Root factor exceeds integer resource");
    }
    if (!in)
      throw runtime_error("Incomplete root factor");
    array<Big, 8100> Q;
    for (int i = 0; i < 90; i++)
      for (int j = 0; j < 90; j++) {
        Big v = 0;
        for (int k = 0; k < 90; k++)
          v += Big(F[i * 90 + k]) * F[j * 90 + k];
        Q[i * 90 + j] = v;
      }
    ifstream old(argv[2]);
    int n8;
    string den;
    old >> n8 >> den;
    if (n8 != 12346 || integer(den) != Big(40320) * D * D)
      throw runtime_error("Eight/nine Gram denominators differ");
    vector<Big> values8(n8);
    for (auto &v : values8) {
      string s;
      old >> s;
      v = integer(s);
    }
    if (!old)
      throw runtime_error("Incomplete eight-values");
    string prefix = argv[3];
    ifstream hosts(prefix + ".u64", ios::binary),
        sub(prefix + ".sub8.u16", ios::binary),
        ptr(prefix + ".root1.ptr.u64", ios::binary),
        ix(prefix + ".root1.indices.u16", ios::binary),
        counts(prefix + ".root1.counts.u16", ios::binary);
    uint64_t previous = 0;
    if (mode == "cached") {
      ptr.read((char *)&previous, 8);
      if (previous)
        throw runtime_error("Bad sparse pointer origin");
    }
    uint64_t mask;
    int nh = 0, maxhost = -1;
    Big maximum = 0;
    ofstream progress(string(argv[5]) + ".progress.jsonl");
    while (hosts.read((char *)&mask, 8)) {
      if (mask >= (1ull << 36))
        throw runtime_error("Invalid nine host");
      Big value = 0;
      for (int i = 0; i < 9; i++) {
        uint16_t id;
        sub.read((char *)&id, 2);
        if (!sub || id >= 12346)
          throw runtime_error("Invalid eight-subhost");
        value += values8[id];
      }
      Big root = 0;
      if (mode == "cached") {
        uint64_t next;
        ptr.read((char *)&next, 8);
        if (!ptr || next < previous || next - previous > 4095)
          throw runtime_error("Invalid root sparse pointer");
        int total = 0;
        for (uint64_t p = previous; p < next; p++) {
          uint16_t id, count;
          ix.read((char *)&id, 2);
          counts.read((char *)&count, 2);
          if (!ix || !counts || id >= 4095 || count > 630)
            throw runtime_error("Invalid root coefficient");
          int j = 0;
          while ((j + 1) * (j + 2) / 2 <= id)
            j++;
          int i = id - j * (j + 1) / 2;
          root += Big(count) * Q[i * 90 + j];
          total += count;
        }
        if (total != 630)
          throw runtime_error("Wrong root normalisation");
        previous = next;
      } else {
        int total = 0;
        // Independent enumeration by 5-subsets and their selected root, rather
        // than by roots and complementary 4-subset partition tables.
        for (unsigned five = 0; five < 512; five++)
          if (popcount(five) == 5)
            for (int r = 0; r < 9; r++)
              if ((five >> r) & 1) {
                unsigned left = five ^ (1u << r), right = 511u ^ five;
                int i = classifier[flagcode(mask, r, left)],
                    j = classifier[flagcode(mask, r, right)];
                root += Q[i * 90 + j];
                total++;
              }
        if (total != 630)
          throw runtime_error("Independent root partition count failed");
      }
      value += Big(576) * root;
      if (maxhost < 0 || value > maximum) {
        maximum = value;
        maxhost = nh;
      }
      nh++;
      if (nh % 2000 == 0) {
        progress << "{\"verified_nine_hosts\":" << nh << ",\"mode\":\"" << mode
                 << "\"}\n";
        progress.flush();
      }
    }
    if (nh != 274668)
      throw runtime_error("Incomplete nine-host resource");
    ofstream out(argv[5]);
    out << "{\"status\":\"exact_nine_gram_upper_recomputed\",\"mode\":\""
        << mode << "\",\"host_classes\":" << nh << ",\"numerator\":\""
        << decimal(maximum) << "\",\"denominator\":\""
        << decimal(Big(362880) * D * D) << "\",\"maximizing_host\":" << maxhost
        << ",\"root_one_gram_psd_exact\":true,\"problem_closed\":false}\n";
    return 0;
  } catch (exception &e) {
    cerr << e.what() << "\n";
    return 1;
  }
}
