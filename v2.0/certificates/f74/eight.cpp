#include <bit>
// A second sampling algorithm: average over all 8! ordered host vertices.
// Recounts the target, proves host orbit coverage, constructs Gram PSD exactly.
#include "../f26/unrestricted/wide_integer.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
using Big = WideInteger;
string decimal(Big x) {
  if (!x)
    return "0";
  bool neg = x < 0;
  if (neg)
    x = -x;
  string s;
  while (x) {
    s += static_cast<char>('0' + static_cast<int>(x % 10));
    x /= 10;
  }
  if (neg)
    s += '-';
  reverse(s.begin(), s.end());
  return s;
}
struct Block {
  int roots, dimension;
  vector<int> index;
  vector<Big> q;
};
int packed(const bool adj[8][8], const int *vertices, int n) {
  int m = 0, b = 0;
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++, b++)
      if (adj[vertices[i]][vertices[j]])
        m |= 1 << b;
  return m;
}
int extension(const bool adj[8][8], const array<int, 8> &p, int k, int offset,
              int u) {
  int result = 0, b = 0;
  for (int i = 0; i < k; i++)
    for (int j = 0; j < u; j++, b++)
      if (adj[p[i]][p[offset + j]])
        result |= 1 << b;
  for (int i = 0; i < u; i++)
    for (int j = i + 1; j < u; j++, b++)
      if (adj[p[offset + i]][p[offset + j]])
        result |= 1 << b;
  return result;
}
Big integer(string s) {
  Big v = 0;
  if (s.empty())
    throw runtime_error("missing integer");
  for (char c : s) {
    if (c < '0' || c > '9')
      throw runtime_error("invalid bound");
    v = v * 10 + c - '0';
  }
  return v;
}
int main(int argc, char **argv) {
  try {
    if (argc != 5 && argc != 6)
      throw runtime_error("input, expected numerator, result, progress [exact "
                          "host rows] required");
    ifstream in(argv[1]);
    int blocks, hosts;
    int64_t d;
    in >> blocks >> hosts >> d;
    if (blocks != 169 || hosts != 12346 || d <= 0 || d > (1ll << 42))
      throw runtime_error("invalid dimensions");
    vector<Block> bs;
    map<int, vector<int>> type;
    for (int k : {2, 4, 6})
      type[k] = vector<int>(1 << (k * (k - 1) / 2), -1);
    for (int b = 0; b < blocks; b++) {
      int n, k, root, size, width;
      in >> n >> k >> root >> size >> width;
      if (!type.count(k) || 2 * n - k != 8 || size < 1 || size > 272 ||
          width < 1 || width > 9 || root < 0 || root >= (int)type[k].size() ||
          type[k][root] != -1)
        throw runtime_error("invalid block");
      type[k][root] = b;
      Block block;
      block.roots = k;
      block.dimension = size;
      block.index.resize(1 << width);
      for (int &i : block.index) {
        in >> i;
        if (i < 0 || i >= size)
          throw runtime_error("flag index outside matrix");
      }
      vector<int64_t> f(size * size);
      for (auto &v : f) {
        in >> v;
        if (v <= -(1ll << 47) || v >= (1ll << 47))
          throw runtime_error("factor resource exceeded");
      }
      block.q.resize(size * size);
      for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++) {
          Big dot = 0;
          for (int c = 0; c < size; c++)
            dot += Big(f[i * size + c]) * f[j * size + c];
          block.q[i * size + j] = dot;
        }
      bs.push_back(move(block));
    }
    bool motif[8][8]{};
    for (auto e :
         vector<pair<int, int>>{{0, 1}, {0, 3}, {0, 4}, {1, 2}, {1, 3}, {3, 5}, {4, 5}})
      motif[e.first][e.second] = motif[e.second][e.first] = true;
    array<int, 6> p6{0, 1, 2, 3, 4, 5};
    set<int> motif_orbit;
    do {
      motif_orbit.insert(packed(motif, p6.data(), 6));
    } while (next_permutation(p6.begin(), p6.end()));
    Big bound = integer(argv[2]), observed = 0;
    int maximizer = -1;
    set<int> canonical_graphs;
    uint64_t covered = 0;
    ofstream progress(argv[4]);
    vector<string> exact_rows;
    for (int h = 0; h < hosts; h++) {
      uint32_t graph;
      int declared_count;
      in >> graph >> declared_count;
      if (!in || graph >= (1u << 28))
        throw runtime_error("invalid host");
      bool adj[8][8]{};
      array<int, 8> degree{};
      int b = 0;
      for (int i = 0; i < 8; i++)
        for (int j = i + 1; j < 8; j++, b++)
          if ((graph >> b) & 1) {
            adj[i][j] = adj[j][i] = true;
            degree[i]++;
            degree[j]++;
          }
      int count = 0;
      for (int omit1 = 0; omit1 < 8; omit1++)
        for (int omit2 = omit1 + 1; omit2 < 8; omit2++) {
          int vertices[6], j = 0;
          for (int v = 0; v < 8; v++)
            if (v != omit1 && v != omit2)
              vertices[j++] = v;
          if (motif_orbit.count(packed(adj, vertices, 6)))
            count++;
        }
      if (count != declared_count)
        throw runtime_error("independent target recount differs");
      Big value = Big(count) * 1440 * d * d;
      array<int, 8> p{0, 1, 2, 3, 4, 5, 6, 7};
      int canonical = numeric_limits<int>::max(), aut = 0;
      do {
        bool sorted_degrees = true;
        for (int i = 1; i < 8; i++)
          if (degree[p[i - 1]] > degree[p[i]]) {
            sorted_degrees = false;
            break;
          }
        if (sorted_degrees) {
          int image = packed(adj, p.data(), 8);
          if (image < canonical) {
            canonical = image;
            aut = 1;
          } else if (image == canonical)
            aut++;
        }
        for (int k : {2, 4, 6}) {
          int root = packed(adj, p.data(), k), block_id = type[k][root];
          if (block_id < 0)
            continue;
          const Block &block = bs[block_id];
          int u = (8 - k) / 2;
          int left = block.index[extension(adj, p, k, k, u)],
              right = block.index[extension(adj, p, k, k + u, u)];
          value += block.q[left * block.dimension + right];
        }
      } while (next_permutation(p.begin(), p.end()));
      if (aut < 1 || 40320 % aut || !canonical_graphs.insert(canonical).second)
        throw runtime_error("host orbit certificate inconsistent");
      covered += 40320 / aut;
      if (value > bound)
        throw runtime_error(
            "an exact host inequality exceeds proposed upper bound at host " +
            to_string(h));
      if (argc == 6)
        exact_rows.push_back(decimal(value));
      if (maximizer < 0 || value > observed) {
        observed = value;
        maximizer = h;
      }
      if (h % 100 == 0) {
        progress << "{\"independent_hosts_verified\":" << h + 1
                 << ",\"total\":" << hosts << "}\n";
        progress.flush();
      }
    }
    if (covered != (1ull << 28) || canonical_graphs.size() != 12346 ||
        observed != bound)
      throw runtime_error("incomplete host coverage or unequal computed bound");
    ofstream out(argv[3]);
    out << "{\"status\":\"independently_verified_integer_gram_upper\",\"target_"
           "recounted\":true,\"all_host_orbits_complete\":true,\"gram_psd_by_"
           "construction\":true,\"host_classes\":12346,\"labelled_graphs_"
           "covered\":"
        << covered << ",\"numerator\":\"" << decimal(observed)
        << "\",\"denominator\":\"" << decimal(Big(40320) * d * d)
        << "\",\"maximizing_host_index\":" << maximizer
        << "}\n";
    if (argc == 6) {
      ofstream rows(argv[5]);
      rows << "{\"status\":\"independent_exact_eight_host_rows\","
              "\"denominator\":\""
           << decimal(Big(40320) * d * d) << "\",\"host_numerators\":[";
      for (int h = 0; h < hosts; h++) {
        if (h)
          rows << ",";
        rows << "\"" << exact_rows[h] << "\"";
      }
      rows << "]}\n";
      rows.flush();
      if (!rows)
        throw runtime_error("Exact host row output failed");
    }
    return 0;
  } catch (const exception &e) {
    cerr << e.what() << "\n";
    return 1;
  }
}
