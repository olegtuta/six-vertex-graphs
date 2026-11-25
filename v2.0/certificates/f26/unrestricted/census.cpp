// Independent target recognition and complete order-nine orbit coverage.
#include "graphs.h"
#include <array>
pair<U, int> orbit(const G &g) {
  vector<int> v(g.n), d(g.n);
  iota(v.begin(), v.end(), 0);
  for (int i = 0; i < g.n; i++)
    d[i] = popcount(g.adj[i]);
  sort(v.begin(), v.end(),
       [&](int a, int b) { return make_pair(d[a], a) < make_pair(d[b], b); });
  vector<pair<int, int>> groups;
  for (int l = 0; l < g.n;) {
    int r = l + 1;
    while (r < g.n && d[v[l]] == d[v[r]])
      r++;
    groups.push_back({l, r});
    l = r;
  }
  U best = ~0ull;
  int multiplicity = 0;
  function<void(int)> visit = [&](int t) {
    if (t == (int)groups.size()) {
      U m = induced(g, v).mask;
      if (m < best) {
        best = m;
        multiplicity = 1;
      } else if (m == best)
        multiplicity++;
      return;
    }
    auto [l, r] = groups[t];
    sort(v.begin() + l, v.begin() + r);
    do {
      visit(t + 1);
    } while (next_permutation(v.begin() + l, v.begin() + r));
  };
  visit(0);
  return {best, multiplicity};
}
int main(int argc, char **argv) {
  try {
    if (argc != 4)
      throw runtime_error("nine-prefix eight-hosts output-json required");
    string prefix = argv[1];
    ifstream input(prefix + ".u64", ios::binary),
        targetfile(prefix + ".target.u8", ios::binary),
        subfile(prefix + ".sub8.u16", ios::binary);
    ifstream eights(argv[2]);
    map<U, int> id8;
    for (int h = 0; h < 12346; h++) {
      U m;
      eights >> m;
      if (!eights || !id8.emplace(canonical(graph(8, m)), h).second)
        throw runtime_error("Bad eight-host list");
    }
    bool motif[6][6]{};
    for (auto [a, b] :
         vector<pair<int, int>>{{0, 1}, {0, 2}, {0, 5}, {1, 3}, {4, 5}})
      motif[a][b] = motif[b][a] = true;
    array<bool, 32768> is_f{};
    array<int, 6> p{0, 1, 2, 3, 4, 5};
    do {
      int bit = 0, m = 0;
      for (int i = 0; i < 6; i++)
        for (int j = i + 1; j < 6; j++, bit++)
          if (motif[p[i]][p[j]])
            m |= 1 << bit;
      is_f[m] = true;
    } while (next_permutation(p.begin(), p.end()));
    set<U> classes;
    uint64_t total = 0;
    int rows = 0;
    U mask;
    ofstream progress(string(argv[3]) + ".progress.jsonl");
    while (input.read((char *)&mask, 8)) {
      if (mask >= (1ull << 36))
        throw runtime_error("Invalid nine-mask");
      G g = graph(9, mask);
      auto [c, aut] = orbit(g);
      if (aut < 1 || 362880 % aut || !classes.insert(c).second)
        throw runtime_error("Repeated orbit or invalid automorphism group");
      total += 362880 / aut;
      uint8_t observed;
      targetfile.read((char *)&observed, 1);
      if (!targetfile)
        throw runtime_error("Missing target");
      int exact = 0;
      for (unsigned s = 1; s < 512; s++)
        if (popcount(s) == 6)
          exact += is_f[induced(g, vertices(s)).mask];
      if (exact != observed)
        throw runtime_error("Independent F26 target count mismatch");
      for (int omit = 0; omit < 9; omit++) {
        uint16_t declared;
        subfile.read((char *)&declared, 2);
        if (!subfile)
          throw runtime_error("Missing subhost");
        vector<int> v;
        for (int i = 0; i < 9; i++)
          if (i != omit)
            v.push_back(i);
        if (id8.at(canonical(induced(g, v))) != declared)
          throw runtime_error("Induced eight-host mismatch");
      }
      rows++;
      if (rows % 1000 == 0) {
        progress << "{\"audited_nine_hosts\":" << rows
                 << ",\"labelled_orbits_covered\":" << total << "}\n";
        progress.flush();
      }
    }
    if (rows != 274668 || total != (1ull << 36))
      throw runtime_error("Nine-host orbit cover incomplete");
    ofstream out(argv[3]);
    out << "{\"status\":\"complete_nine_host_orbit_coverage_and_targets_"
           "verified\",\"host_classes\":"
        << rows << ",\"labelled_graphs_covered\":" << total
        << ",\"independent_F26_recognition\":true,\"all_induced_eight_hosts_"
           "checked\":true,\"problem_closed\":false}\n";
    return 0;
  } catch (exception &e) {
    cerr << e.what() << "\n";
    return 1;
  }
}
