import "mass"
import "thetaPOT"
import "treeProcessIrr"
import "kdTreeIrregularRankK"

def main [m] [n] [d] (Tval: i32) (k: i64) (height: i32) (test_set: [m][d]f32) (queries: [n][d]f32) =
  -- Step 1: shift points
  let shifted_points = shiftPoints test_set

  -- Setup for loop (maybe it should be a loop instead of map)
  let init_knns = replicate n (replicate k (-1i32, f32.inf))

  -- Step 2-5 (and 6) Loop:
  let new_knns =
    loop curr_nns = init_knns for t < Tval do
      -- Step 2 Perform the pseudo random orthogonal transformation on the test set and quiery set
      let M1 = i64.i32 <| log2Int d
      let transformed_test_set = pseudoRandomOrthogonalTransformation M1 t test_set
      let transformed_queries  = (pseudoRandomOrthogonalTransformation M1 t queries)
      
      -- Step 3 Build the kd-tree
      let (leaves, indir, median_dims , median_vals, shp_arr) = buildKdTree height test_set

      -- Step 4 & 5 Search the tree and find new candidates
      in  searchForKnns transformed_queries curr_nns
                        leaves indir median_dims median_vals shp_arr 
                        height

  -- Step 7 perform depth one search "supercharging" on the found knns of quiries
  let (k_inds, k_dists) =  unzip <| map (\i_knn -> unzip i_knn) new_knns
  in (k_inds[0], k_dists[0])