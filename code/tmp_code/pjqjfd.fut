
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
let calculate_Qj [d] [d1] (point: *[d]f32) (rand_numbers: [d1]f32) : *[d]f32 =
    loop acc = point for i < d1  do
        let tmp = acc[i]
        let acc[i]   = ((f32.cos rand_numbers[i]) * acc[i]  + (f32.sin rand_numbers[i]) * acc[i+1])
        let acc[i+1] = ((-f32.sin rand_numbers[i]) * tmp + (f32.cos rand_numbers[i]) * acc[i+1])
        in  acc
 
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


--let Z_inv [d2] (point: [d2](f32,f32)) : []f32 =
--    map (\i -> if i%2 == 0  then c32.re point[i/2]
--                            else c32.im point[i/2] ) (iota (2*d2))

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


--let Fd [d] (point: [d]f32)  =
--    let d2 = d / 2
--    let pZ = Z point d2
--    let one_over_d2sqr = c32.mk_re (1.0 / ( (f32.sqrt (f32.i64 d2)))) 
--    let T = map (\k ->  
--                    map (\l ->
--                        let exponent_im = c32.exp (c32.mk_im 
--                                (-(2.0 * f32.pi * (f32.i64 k) * (f32.i64 l)) / (f32.i64 d2)))
--                        in one_over_d2sqr c32.* exponent_im
--                        ) (iota d2)  
--                ) (iota d2)
--    let pTZ = mat_vec_complex T pZ
--
--    in Z_inv pTZ


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
    --let points_Qj = map (\vec -> calculate_Qj vec rand_numbers_Q_normalized2pi[0:d-1:1] ) points'
    let points_Fd =  map (\vec -> Fd vec (d/2)) points  


    in points_Fd