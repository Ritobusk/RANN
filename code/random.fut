-- ==
-- compiled input { 
--   [2.0f32, -1.0f32, -1.0f32, 2.0f32]
-- } 

import "lib/github.com/diku-dk/cpprandom/random"
import "lib/github.com/diku-dk/cpprandom/shuffle"

local module rng_engine = pcg32
local module shuffle = mk_shuffle rng_engine

let main [n]
         (point : [n]f32)  =
	let a = rng_engine.rand
    in a