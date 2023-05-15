import "mass"
import "thetaPOT"
import "treeProcessIrr"
import "kdTreeIrregularRankK"

def RANN [m] [n] [d] (Tval: i32) (k: i64) (test_set: [m][d]f32) (queries: [n][d]f32) =
  -- Step 1: shift points
  let shifted_points = shiftPoints test_set

  -- Setup for loop 
  let init_knns = replicate n (replicate k (-1i32, f32.inf))
  let height =  ( log2Int (m / 256))

  -- Step 2-6 The loop:
  let new_knns =
    loop curr_nns = init_knns for t < Tval do
      -- Step 2 Perform the pseudo random orthogonal transformation on the test set and quiery set
      let M1 = i64.i32 <| log2Int d
      let transformed_test_set = pseudoRandomOrthogonalTransformation M1 t test_set
      let transformed_queries  = (pseudoRandomOrthogonalTransformation M1 t queries)
      
      -- Step 3 Build the kd-tree
      let (leaves, indir, median_dims, median_vals, shp_arr) = buildKdTree height transformed_test_set

      -- Step 4 & 5 Search the tree and find new candidates
      in searchForKnns transformed_queries curr_nns
                       leaves indir median_dims median_vals shp_arr 
                       height

  -- Step 7 omitted here
  let (k_inds, _) =  unzip <| map (\i_knn -> unzip i_knn) new_knns
  in k_inds


def superRANN [m] [n] [d] (Tval: i32) (k: i64) (test_set: [m][d]f32) (queries: [n][d]f32) =
  -- Step 1: shift points
  let shifted_points = shiftPoints test_set

  -- Setup for loop 
  let init_knns_q = replicate n (replicate k (-1i32, f32.inf))
  let init_knns_t = replicate m (replicate k (-1i32, f32.inf))
  let height =  ( log2Int (m / 256))

  -- Step 2-6 The loop:
  let (new_knns_q, new_knns_t) =
    loop (curr_knns_q, curr_knns_t)= (init_knns_q, init_knns_t) for t < Tval do
      -- Step 2 Perform the pseudo random orthogonal transformation on the test set and quiery set
      let M1 = i64.i32 <| log2Int d
      let transformed_test_set = pseudoRandomOrthogonalTransformation M1 t test_set
      let transformed_queries  = (pseudoRandomOrthogonalTransformation M1 t queries)
      
      -- Step 3 Build the kd-tree
      let (leaves, indir, median_dims, median_vals, shp_arr) = buildKdTree height transformed_test_set

      -- Step 4 & 5 Search the tree and find new candidates
      let curr_knns_q' = searchForKnns transformed_queries curr_knns_q
                          leaves indir median_dims median_vals shp_arr 
                          height
      let curr_knns_t' = searchForKnns transformed_test_set curr_knns_t
                          leaves indir median_dims median_vals shp_arr 
                          height
      in (curr_knns_q', curr_knns_t')

  -- Step 7 perform depth one search "supercharging" on the found knns of queries
  let (knn_inds_q, _) =  unzip <| map (\i_knn -> unzip i_knn) new_knns_q
  let (knn_inds_t, _) =  unzip <| map (\i_knn -> unzip i_knn) new_knns_t
  let supercharging = 
    let ksqr = k**2
    let ksqr_inds = map (\knn_ind_q -> (flatten <| 
                            (map (\ind -> knn_inds_t[ind]) knn_ind_q)) :> [ksqr]i32 ) knn_inds_q
    
    let ksqr_points = map (\inds -> map (\ind -> (ind, test_set[ind])) inds) ksqr_inds
    in map2 (\refs ind -> bruteForce queries[ind] new_knns_q[ind] refs ) ksqr_points (iota n)

  let (super_knn_inds, _) =  unzip <| map (\i_knn -> unzip i_knn) supercharging
  in (super_knn_inds)


def main [m] [n] [d] (Tval: i32) (k: i64) (test_set: [m][d]f32) (queries: [n][d]f32) =
  superRANN Tval k test_set queries