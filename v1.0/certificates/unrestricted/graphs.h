#include <algorithm>
#include <bit>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
using U = uint64_t;
struct G {
  int n;
  U mask;
  vector<unsigned> adj;
};
G graph(int n, U mask) {
  G g{n, mask, vector<unsigned>(n)};
  int b = 0;
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++, b++)
      if ((mask >> b) & 1) {
        g.adj[i] |= 1u << j;
        g.adj[j] |= 1u << i;
      }
  return g;
}
G induced(const G &g, const vector<int> &v) {
  U m = 0;
  int b = 0;
  for (int i = 0; i < (int)v.size(); i++)
    for (int j = i + 1; j < (int)v.size(); j++, b++)
      if ((g.adj[v[i]] >> v[j]) & 1)
        m |= 1ull << b;
  return graph((int)v.size(), m);
}
vector<int> vertices(unsigned s) {
  vector<int> r;
  for (int i = 0; s; s >>= 1, i++)
    if (s & 1)
      r.push_back(i);
  return r;
}
U canonical(const G &g) {
  vector<int> v(g.n), d(g.n);
  iota(v.begin(), v.end(), 0);
  for (int i = 0; i < g.n; i++)
    d[i] = popcount(g.adj[i]);
  sort(v.begin(), v.end(),
       [&](int a, int b) { return make_pair(d[a], a) < make_pair(d[b], b); });
  vector<pair<int, int>> groups;
  for (int l = 0; l < g.n;) {
    int r = l + 1;
    while (r < g.n && d[v[r]] == d[v[l]])
      r++;
    groups.push_back({l, r});
    l = r;
  }
  U best = ~0ull;
  function<void(int)> visit = [&](int t) {
    if (t == (int)groups.size()) {
      best = min(best, induced(g, v).mask);
      return;
    }
    auto [l, r] = groups[t];
    sort(v.begin() + l, v.begin() + r);
    do {
      visit(t + 1);
    } while (next_permutation(v.begin() + l, v.begin() + r));
  };
  visit(0);
  return best;
}
