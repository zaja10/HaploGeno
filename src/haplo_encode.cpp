#include <Rcpp.h>
#include <map>
#include <vector>
using namespace Rcpp;

//' Encode Haplotypes (Optimized)
//'
//' Maps unique rows of a genotype matrix to integers without string conversion.
//'
//' @param mat A numeric matrix of genotypes (0, 1, 2).
//' @return Integer vector of haplotype IDs.
//' @export
// [[Rcpp::export]]
 IntegerVector encode_block_fast(NumericMatrix mat) {
   int n = mat.nrow();
   int m = mat.ncol();
   
   // Use vector<int> as key to avoid string conversion overhead
   std::map<std::vector<int>, int> haplo_map;
   IntegerVector ids(n);
   
   int current_id = 1;
   
   for (int i = 0; i < n; i++) {
     // Construct vector key for current row
     std::vector<int> row_key(m);
     for (int j = 0; j < m; j++) {
       row_key[j] = (int)mat(i, j);
     }
     
     // Lookup or Insert
     // Note: lower_bound is slightly more efficient than find + insert
     std::map<std::vector<int>, int>::iterator it = haplo_map.lower_bound(row_key);
     if (it == haplo_map.end() || haplo_map.key_comp()(row_key, it->first)) {
       haplo_map.insert(it, std::make_pair(row_key, current_id));
       ids[i] = current_id;
       current_id++;
     } else {
       ids[i] = it->second;
     }
   }
   
   return ids;
 }