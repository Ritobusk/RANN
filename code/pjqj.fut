import "lib/github.com/diku-dk/cpprandom/random"
import "lib/github.com/diku-dk/cpprandom/shuffle"

local module real = { open f64 }
local module rng_engine = pcg32
local module shuffle = mk_shuffle rng_engine

-- Permutation function
let calculate_Pj [d] (point: [d]f32) (permutation : [d]i64) : [d]f32 = 
    map (\i -> point[i]) permutation



---- I do not think inplace updates will work... ----
let modify_Qjk [d]  (point: *[d]f32) (rand_number: f32) (k : i32) : *[d]f32 =
    let a =  (f32.cos rand_number) * point[k]  + (f32.sin rand_number) * point[k+1]
    in point with [k]  = a
    

let modify_Qjk1 [d] (point: *[d]f32) (rand_number: f32) (k : i32) : *[d]f32 =
        let b =  (-f32.sin rand_number) * point[k] + (f32.cos rand_number) * point[k+1] 
        in  point with [k+1] = b



---- FAILED Q ATTEMPTS ----


let calculate_Qj_tmp [d] [d1] (point: [d]f32) (rand_numbers: [d1]f32) : [d]f32 =
    map (\k ->
            if k == 0          then ((f32.cos rand_numbers[k]) * point[k]  + (f32.sin rand_numbers[k]) * point[k+1])
            else if k == (d-1) then ((-f32.sin rand_numbers[k-1]) * point[k-1] + (f32.cos rand_numbers[k-1]) * point[k]) 
            else let tmp =  ((-f32.sin rand_numbers[k-1]) * point[k-1] + (f32.cos rand_numbers[k-1]) * point[k] )
                    in ((f32.cos rand_numbers[k]) * tmp  + (f32.sin rand_numbers[k]) * point[k+1] )
        ) (iota d)

let calculate_Qj2_tmp [d] [d1] (point: [d]f32) (rand_numbers: [d1]f32) (k: i64) (acc : f32)  =

    if k == 0   then ((f32.cos rand_numbers[k]) * point[k]  + (f32.sin rand_numbers[k]) * point[k+1])
    else if k == 1   then  
                    let tmp =  ((-f32.sin rand_numbers[k-1]) * point[k-1] + (f32.cos rand_numbers[k-1]) * point[k] )
                    in ((f32.cos rand_numbers[k]) * tmp  + (f32.sin rand_numbers[k]) * point[k+1])
                    
    else if k == (d-1) then ((-f32.sin rand_numbers[k-1]) * point[k-1] + (f32.cos rand_numbers[k-1]) * point[k]) 
    else let tmp =  ((-f32.sin rand_numbers[k-1]) * acc + (f32.cos rand_numbers[k-1]) * point[k] )
        let acc_val = ((f32.cos rand_numbers[k]) * tmp  + (f32.sin rand_numbers[k]) * point[k+1] )
        in acc_val
--         
--let calculate_Qj2 [d] [d1] (point: [d]f32) (rand_numbers: [d1]f32) : [d]f32 =
--    let arr = map (\i -> f32.i64 i) (iota d)
--    let a =
--        loop acc = point[0] for x in (iota d1) do
--            let arr[x] = f32.i64 (x+2)
--            in x
--    in arr

--let calculate_Qj3 [d] [d1] (point: [d]f32) (rand_numbers: [d1]f32) =
--    let i = 0
--    let x = loop a = point[0] for i < 3 do
--                let b = calculate_Qj2 point rand_numbers i a
--                in concat a b
--    in x
---- FAILED Q ATTEMPTS END ----

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
    in points_Qj
-- [[-0.924159f32, -1.109819f32, 0.607398f32, -3.863729f32], [-4.771064f32, -1.707392f32, -0.522218f32, -7.393312f32]]
    --3.104016f32, 4.048981f32, 3.253206f32,