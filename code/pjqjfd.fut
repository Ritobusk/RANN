
import "lib/github.com/diku-dk/cpprandom/random"
import "lib/github.com/diku-dk/cpprandom/shuffle"
import "lib/github.com/diku-dk/complex/complex"

local module real = { open f64 }
local module rng_engine = pcg32
local module shuffle = mk_shuffle rng_engine
module c64 = mk_complex f32

-- Permutation function
let calculate_Pj [d] (point: [d]f32) (permutation : [d]i64) : [d]f32 = 
    map (\i -> point[i]) permutation


let call_Qjk [d] [d1] (point: [d]f32) (rand_numbers: [d1]f32) (i : i64) : [d]f32 =
    let Qjk = map (\k ->
                    if      i == k then ((f32.cos rand_numbers[i]) * point[i]  + (f32.sin rand_numbers[i]) * point[i+1])
                    else if i+1 == k then  ((-f32.sin rand_numbers[i]) * point[i] + (f32.cos rand_numbers[i]) * point[i+1] )
                    else point[k] 
                ) (iota d)
    in Qjk
 
let calculate_Qj [d] [d1] (point: [d]f32) (rand_numbers: [d1]f32) : [d]f32 =
    loop acc = call_Qjk point rand_numbers 0 for i < d1 -1  do
        call_Qjk acc rand_numbers (i+1)

let Z [d] (point: [d]f32) (d2 : i64) : [d2](f32,f32)  =
    map (\i -> c64.mk (point[i*2]) (point[(i*2) + 1])) (iota d2)

let Z_inv [d2] (point: [d2](f32,f32)) : []f32 =
    map (\i -> if i%2 == 0  then c64.re point[i/2]
                            else c64.im point[i/2] ) (iota (2*d2))

let Fd [d] (point: [d]f32)  =
    let d2 = d / 2
    let bim = Z point d2
    in Z_inv bim


let main [n][d]
         (points : [n][d]f32)  =
    
    let M1 =  d / 2
    let M2 = M1

    let points' = points

    -- RNGs for the P and Q functions
    let rng_origin_P = rng_engine.rng_from_seed [1]
	let rngs_M1M2 = rng_engine.split_rng (M1+M2) rng_origin_P
    let rng_origin_Q = rng_engine.rng_from_seed [2]
	let rngs_dM1M2 = rng_engine.split_rng ((M1+M2) * (d-1)) rng_origin_Q
    
    -- (d-1)*(m1+m2) random numbers for Q normalized and then multiplied by 2 pi
    let (rngs_Q, rand_numbers_Q)     = unzip <| map (\r -> rng_engine.rand r) rngs_dM1M2 
    let rand_numbers_Q_normalized2pi = map (\num -> (f32.f64 (((f64.u32 num) / (f64.u32 4294967295)) * (2*real.pi) ))) rand_numbers_Q
    --rand_numbers_Q_normalized2pi[0:d-1:1] the first d-1 random numbers [i*(d-1):(i+1)(d-1):1] when iterating over m1 and m2


    let (rngs_P, permutations_P) = unzip <| map (\r -> shuffle.shuffle r (iota d)) rngs_M1M2
    let points_Pj = map (\vec -> calculate_Pj vec permutations_P[0]) points' 
    let points_Qj = map (\vec -> calculate_Qj vec rand_numbers_Q_normalized2pi[0:d-1:1] ) points'
    let points_Fd =  Fd [1.0f32, 2.0f32, 3.0f32, 4.0f32] 


    in points_Fd
-- [[-0.924159f32, -1.109819f32, 0.607398f32, -3.863729f32], [-4.771064f32, -1.707392f32, -0.522218f32, -7.393312f32]]
    --3.104016f32, 4.048981f32, 3.253206f32,