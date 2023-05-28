import "knn-driver"
import "knn-bruteforce"


def main [m] [n] [d] (Tval: i32) (k: i64) (test_set: [m][d]f32) (queries: [n][d]f32) =
  let queries' = (queries[:100])
  let RANN_knns  = testRANN Tval k test_set queries'
  let brute_knns = bruteNNs k test_set queries'

  let correct_nns =
        let correct_nns_tmp = 
          map2 (\R_knn b_knn ->
                    let isInTrueNN  = map (\neighbour -> any (\trueNs -> neighbour == trueNs) b_knn) R_knn
                    let isInTrueNNInt = map (\iit -> if iit then 1i64 else 0i64) isInTrueNN
                    in reduce (+) 0i64 isInTrueNNInt
                ) RANN_knns brute_knns
        
        in reduce (+) 0i64 correct_nns_tmp
  let accuracy =  (f64.i64 correct_nns) / (f64.i64 (k * 100))

  in (accuracy)
