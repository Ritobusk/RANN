-- ==
-- compiled input { 
--   [2.0f32, -1.0f32, -1.0f32, 2.0f32]
-- } 

import "lib/github.com/diku-dk/cpprandom/random"
import "lib/github.com/diku-dk/cpprandom/shuffle"

local module rng_engine = pcg32
local module shuffle = mk_shuffle rng_engine

let Calculate_Pj [d] (point: [d]f32) (permutation : [d]i64) : [d]f32 =
    map (\per -> point[per]) permutation

let main [n][d]
         (points : [n][d]f32)  =
    
    let M1 =  d / 2
    let M2 = M1

    let rng_origin = rng_engine.rng_from_seed [1]
	let rngs_M1M2 = rng_engine.split_rng (M1+M2) rng_origin
    
    let dim_arr = iota d

    let (rngs_P, permutations_P) = unzip <| map (\r -> shuffle.shuffle r dim_arr) rngs_M1M2
    let points_Pj = map (\vec -> Calculate_Pj vec permutations_P[0]) points 
    in points_Pj