import "mass"
import "thetaPOT"
import "treeProcess"

def main [m] [n] [d] (Tval: i64) (k: i64) (defppl: i32) (test_set: [m][d]f32) (queries: [n][d]f32) =
    -- Step 1: shift points
    let shifted_points = shiftPoints test_set

    -- Setup for loop (maybe it should be a loop instead of map)
    let init_knns_q = replicate n (replicate k (-1i32, f32.inf))

    -- Step 2-5 (and 6) Loop:
    let new_knns =
        loop curr_nn_set_after_t_iterations = init_knns_q for t < Tval do
            -- Perform the pseudo random orthogonal transformation on the test set and quiery set
            let t = i32.i64 (t_val+1)
            let transformed_test_set = pseudoRandomOrthogonalTransformation t test_set
            let transformed_queries  = pseudoRandomOrthogonalTransformation t queries

            -- Find knn for this iteration

    -- Step 7 perform depth one search "supercharging" on the found knns of quiries

        