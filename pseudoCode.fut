-- reorders the values of a vector/point randomly
let Calculate_Pj(Vector: [d]f32) : [d]f32 =
    -- We generate permutations π1;...;π_{M1+M2} of the numbers {1,...,d} uniformly
    -- at random and independent of each other. Using the shuffle module in Futhark
    -- https://futhark-lang.org/pkgs/github.com/diku-dk/cpprandom/1.3.1/doc/lib/github.com/diku-dk/cpprandom/shuffle.html

    let rng                 = Make an "rng_engine" with min=0 and max=d-1
    let shuffler            = Make a "shuffle"x 
    let (rng_state, shuffled_array) = shuffler.shuffle (rng.seed) (Vector)
    in shuffled_array

-- Does the planar rotation as exlpained in equation 6 in the paper.
-- Then it inserts the rotation at index i and i+1 while keeping the other 
-- coordinates intact.
let Calculate_Qjks (Vector: [d]f32) (random_number: f32) (i: i32.64) : [d]f32 =
    let rotation_ks = [[f32.cos random_number, f32.sin random_number], 
                       [-f32.sin random_number, f32.cos random_number]] * [[Vector[i]],[Vector[i + 1]]]

    in Qjk = map (\k ->
                    if      i == k then rotation_ks[0]
                    else if i+1 == k+1 then rotation_ks[1]
                    else Vector[k] 
                ) map i32.i64 (iota d) 
    -- !!I am using map here instead of just changing the array at the kth index directly. Is this the best method ? !!

-- Calls Calculate_Qjks as shown in equation 7 from the paper. ↓
let Calculate_Qj (Vector: [d]f32) (random_numbers: [d]f32) : [d]f32 =
    -- Use scan to do the composition of the "Q_{j,k}"s? This looks wrong... (Vector is given to Calculate_Qjks)
    scan (\i -> 
            Calculate_Qjks Vector random_numbers[i] i) 0 (iota d-1)    


-- We need complex numbers for this:
-- https://futhark-lang.org/examples/complex-numbers.html
module c32 = mk_complex f32
let Calculate_Zinverse (Z: [d2]c32) : [d]f32 =
    map (\i ->
            if i%2 == 0 then Z[i/2].re
                        else Z[i/2].im
            ) iota (2*d2)

-- Calculate F^(d) 
let Calculate_Fd (Vector: [d]f32) : [d]f32 =
    -- Is d odd? maybe not nessasary
    let d2 =  if (d%2) then d / 2 else (d-1) / 2
    -- Construct Z ↓
    let d2index = iota d2
    let Z = map (\i -> c32.mk (Vector[i*2]) (Vector[(i*2)+1]) ) d2index
    -- Construct T
    let T = map(\k -> 
                    map(l -> 
                        c32.mk_re (1.0f / sqrt(d2)) c32.*                                     -- <-- Maybe this should be complex too? (With no imaginary part)
                        (c32.exp(  (c32.mk (0) (-(2*pi*(k-1)*(l-1))/d2)  ))) 
                        ) d2index ) d2index
    -- Now we can calculate eq. 10 ↓
    let T' = T * Z -- matix vector mul
    let Fd = Calculate_Zinverse T'
    -- If d was odd
    if (!(d%2)) then Fd :: Vector[d-1]
    in Fd

let Calculate_theta (Vector : [d]f32) (pseudorandom_numbers : [d]f32) (m2: i32) (m1:i32) : [d]f32 = 
    -- !! I dont know how function composition should work with scan !! --
    let PQm2 = scan (\i -> 
                        Calculate_Qj (Calculate_Pj Vector) pseudorandom_numbers[m1 + i]) 0 (iota m2) 
    let FdPQm2 = Calculate_Fd PQm2
    in scan (\i -> 
                    Calculate_Qj (Calculate_Pj Vector) pseudorandom_numbers[i]) 0 (iota m1) 
    

--------------------------------------------------
--------------------------------------------------

let Main (points : [N][d]f32) =
    -- construct (d−1)·(m1+m2) independent pseudorandom numbers for Calculate_Qj
    -- https://futhark-lang.org/pkgs/github.com/diku-dk/cpprandom/1.3.1/doc/lib/github.com/diku-dk/cpprandom/random.html
    let m1 = (log d) / 2
    let m2 = (log d) / 2
    let seed_amount = (d-1) * m1 + m2
    let rng         = Make an "rng_engine" with min=0 and max=2π
    let rng_seeds   = rng.split(seed_amount)
    let (rng_states, pseudorandom_numbers) =  map (\x -> x.rand) rng_seeds
    -- Use some fancy slicing maybe for the "pseudorandom_numbers"

    let newPoints = map (\point ->
                            Calculate_theta point pseudorandom_numbers m2 m1
                            ) points