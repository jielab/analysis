/*
 * Tiny native lookup regression: compile with the Kraken source used by metag.
 * No production data or existing index is read or changed.
 *
 * g++ -std=c++11 -O2 -fopenmp -I KRAKEN_SRC -I f \
 *   f/test_union_memory.cc KRAKEN_SRC/taxonomy.cc \
 *   KRAKEN_SRC/mmap_file.cc KRAKEN_SRC/compact_hash.cc -o /tmp/test_union_memory
 * mkdir /tmp/union-memory-test-NEW
 * /tmp/test_union_memory /tmp/union-memory-test-NEW
 */
#include "metag_union.h"
#include <cassert>

using namespace kraken2;

template<typename Cell>
void check_width(const std::string &folder, unsigned bits, const Taxonomy &taxonomy) {
  const auto a = taxonomy.GetInternalID(3);
  const auto b = taxonomy.GetInternalID(4);
  const auto c = taxonomy.GetInternalID(5);
  const auto genus = taxonomy.GetInternalID(2);
  const auto root = taxonomy.GetInternalID(1);
  assert(a && b && c && genus && root);
  std::vector<std::string> paths;
  for (unsigned t = 0; t < 5; ++t) {
    const auto path = folder + "/hash" + std::to_string(bits) + "." + std::to_string(t) + ".k2d";
    // Use a realistic taxonomy value width supported by both cell types.
    CompactHashTable<Cell> table(1009, bits - 16, 16);
    auto add = [&table](hkey_t key, hvalue_t value) {
      hvalue_t old = 0;
      assert(table.CompareAndSet(key, value, &old));
    };
    for (unsigned i = 0; i < 200; ++i)
      add(1000 + t * 500 + i, i % 3 == 0 ? a : i % 3 == 1 ? b : c);
    if (t == 0) { add(17, a); add(18, a); add(19, a); add(20, a); }
    if (t == 1) { add(18, b); add(20, a); }
    if (t == 2) add(19, c);
    table.WriteTable(path.c_str());
    paths.push_back(path);
  }
  metag_union::UnionHashTable ram(paths, taxonomy, false);
  metag_union::UnionHashTable mapped(paths, taxonomy, true);
  assert(ram.Get(17) == a && mapped.Get(17) == a);
  assert(ram.Get(18) == genus && mapped.Get(18) == genus);
  assert(ram.Get(19) == root && mapped.Get(19) == root);
  assert(ram.Get(20) == a && mapped.Get(20) == a);
  assert(ram.Get(21) == 0 && mapped.Get(21) == 0);
  std::vector<hkey_t> keys;
  for (unsigned i = 0; i < 10000; ++i) keys.push_back(i);
  std::vector<hvalue_t> expected(keys.size()), actual(keys.size());
  ram.GetBatch(keys.data(), expected.data(), keys.size());
  mapped.GetBatch(keys.data(), actual.data(), keys.size());
  assert(expected == actual);
  for (size_t i = 0; i < keys.size(); ++i) {
    assert(expected[i] == ram.Get(keys[i]));
    assert(actual[i] == mapped.Get(keys[i]));
  }
  ram.GetBatch(nullptr, nullptr, 0);
  mapped.GetBatch(nullptr, nullptr, 0);
  // Concurrent lookup must remain identical: every worker owns its scratch.
  #pragma omp parallel for num_threads(4)
  for (int run = 0; run < 20; ++run) {
    std::vector<hvalue_t> batch(keys.size());
    mapped.GetBatch(keys.data(), batch.data(), keys.size());
    assert(batch == expected);
  }
  std::cout << "PASS " << bits << "-bit cells, five categories, 10000 keys; "
            << "RAM=mmap, scalar=batch, known species/genus/root/unmatched, 4 threads" << std::endl;
}

int main(int argc, char **argv) {
  if (argc != 2) return 2;
  const std::string folder(argv[1]);
  const auto nodes = folder + "/nodes.dmp", names = folder + "/names.dmp";
  {
    std::ofstream n(nodes), s(names);
    if (!n || !s) return 2;
    n << "1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tgenus\t|\n"
      << "3\t|\t2\t|\tspecies\t|\n4\t|\t2\t|\tspecies\t|\n5\t|\t1\t|\tspecies\t|\n";
    for (int i = 1; i <= 5; ++i)
      s << i << "\t|\tTest " << i << "\t|\t\t|\tscientific name\t|\n";
  }
  const auto taxfile = folder + "/taxo.k2d";
  NCBITaxonomy ncbi(nodes, names);
  for (uint64_t i = 1; i <= 5; ++i) ncbi.MarkNode(i);
  ncbi.ConvertToKrakenTaxonomy(taxfile.c_str());
  Taxonomy taxonomy(taxfile, true);
  taxonomy.GenerateExternalToInternalIDMap();
  check_width<CompactHashCell>(folder, 32, taxonomy);
  check_width<CompactHashCell40>(folder, 40, taxonomy);
  return 0;
}
