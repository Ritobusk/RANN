import "knn-driver"
import "knn-bruteforce"


def main [m] [n] [d] (Tval: i32) (k: i64) (test_set: [m][d]f32) (queries: [n][d]f32) =
  let RANN_knns  = RANN Tval k test_set queries
  let brute_knns = bruteNNs k test_set queries

  let correct_nns =
        let cor_nns_tmp = 
          map2 (\R_knn b_knn ->
                    let isInTrue  = map (\neighbour -> any (\trueNs -> neighbour == trueNs) b_knn) R_knn
                    let isInTrueInt = map (\iit -> if iit then 1i64 else 0i64) isInTrue
                    in reduce (+) 0i64 isInTrueInt
                ) RANN_knns brute_knns
        
        in reduce (+) 0i64 cor_nns_tmp
  let accuracy =  (f64.i64 correct_nns) / (f64.i64 (k * n))

  in (accuracy)
