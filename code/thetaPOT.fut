
import "lib/github.com/diku-dk/cpprandom/random"
import "lib/github.com/diku-dk/cpprandom/shuffle"
import "lib/github.com/diku-dk/complex/complex"

local module real = { open f64 }
local module rng_engine = pcg32
local module shuffle = mk_shuffle rng_engine
module c32 = mk_complex f32

-- Permutation function
let calculate_Pj [d] (point: [d]f32) (permutation : [d]i64) : [d]f32 = 
    map (\i -> point[i]) permutation

-- Rotation functions
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

-- Fourier transform functions
let mat_vec_complex [d2] (mat : [d2][d2](f32,f32)) (vec : [d2](f32,f32))  =
    map (\row -> reduce (c32.+) (c32.mk 0.0 0.0) <| map2 (\m v -> m c32.* v) row vec ) mat

let Z [d] (point: [d]f32) (d2 : i64) : [d2](f32,f32)  =
    map (\i -> c32.mk (point[i*2]) (point[(i*2) + 1])) (iota d2)

let Z_inv [d2] [d] (TZ: [d2](f32,f32))  (point : [d]f32) : [d]f32 =
    map (\i ->  if i%2==0 && i == d-1 then point[d-1]
                else if i%2 == 0  then c32.re TZ[i/2]
                else c32.im TZ[i/2] 
               
        ) (iota (d))

let Fd [d] (point: [d]f32) (d2 : i64) : [d]f32  =
    let pZ = Z point d2
    let one_over_d2sqr = c32.mk_re (1.0 / ( (f32.sqrt (f32.i64 d2)))) 
    let T = map (\k ->  
                    map (\l ->
                        let exponent_im = c32.exp (c32.mk_im 
                                (-(2.0 * f32.pi * (f32.i64 k) * (f32.i64 l)) / (f32.i64 d2)))
                        in one_over_d2sqr c32.* exponent_im
                        ) (iota d2)  
                ) (iota d2)
    let pTZ = mat_vec_complex T pZ
    in Z_inv pTZ point

let Theta [d] (point :  [d]f32) (permutations : [][]i64) (random_numbers : []f32) (m1 : i64) (m2 : i64) =
    let m2_PQ = loop acc = point for i < m2 do
                    let pointP = calculate_Pj acc permutations[i]
                    let ind_1 = i * (d-1)
                    let ind_2 = (i+1)*(d-1) 
                    in calculate_Qj pointP random_numbers[ind_1:ind_2:1]
    let d2 = d / 2
    let Fd_m2_PQ = Fd m2_PQ d2
    in loop acc = Fd_m2_PQ for i < m1 do
        let pointP = calculate_Pj acc permutations[i + m2]
        let ind_1 = (i * (d-1)) + m2 * (d-1)
        let ind_2 = (i+1) * (d-1) + m2 * (d-1)
        in calculate_Qj pointP random_numbers[ind_1:ind_2:1]


   let main [n][d]
         (points : [n][d]f32)  =
    
    let M1 =  d / 2
    let M2 = M1

    let points' = points

    -- RNGs for the P and Q functions
    let rng_origin_P = rng_engine.rng_from_seed [1]
	let rngs_M1M2 = rng_engine.split_rng (M1+M2) rng_origin_P
    let (rngs_P, permutations_P) = unzip <| map (\r -> shuffle.shuffle r (iota d)) rngs_M1M2

    let rng_origin_Q = rng_engine.rng_from_seed [2]
	let rngs_dM1M2 = rng_engine.split_rng ((M1+M2) * (d-1)) rng_origin_Q
    
    -- (d-1)*(m1+m2) random numbers for Q normalized and then multiplied by 2 pi
    let (rngs_Q, rand_numbers_Q)     = unzip <| map (\r -> rng_engine.rand r) rngs_dM1M2 
    let rand_numbers_Q_normalized2pi = map (\num -> (f32.f64 (((f64.u32 num) / (f64.u32 4294967295)) * (2*real.pi) ))) rand_numbers_Q
    --rand_numbers_Q_normalized2pi[0:d-1:1] the first d-1 random numbers [i*(d-1):(i+1)(d-1):1] when iterating over m1 and m2


    --let points_Pj = map (\vec -> calculate_Pj vec permutations_P[0]) points' 
    --let points_Qj = map (\vec -> calculate_Qj vec rand_numbers_Q_normalized2pi[0:d-1:1] ) points'
    --let points_Fd =  Fd [1.0f32, 2.0f32, 3.0f32, 4.0f32] 


    in map (\p ->  Theta p  permutations_P rand_numbers_Q_normalized2pi M1 M2) points'
-- [[-0.924159f32, -1.109819f32, 0.607398f32, -3.863729f32], [-4.771064f32, -1.707392f32, -0.522218f32, -7.393312f32]]
    --3.104016f32, 4.048981f32, 3.253206f32,