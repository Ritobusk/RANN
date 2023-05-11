def iota32 n = (0..1..<i32.i64 n) :> [n]i32

def sumSqrsSeq [d] (xs: [d]f32) (ys: [d]f32) : f32 =
    loop (res) = (0.0f32) for (x,y) in (zip xs ys) do
        let z = x-y in res + z*z

def bruteForce2 [m][d][k] (query: [d]f32) 
                         (knns0: [k](i32,f32))
                         (refs: [m](i32,[d]f32))
                       : [k](i32,f32) =
    loop (knns) = (copy knns0)
      for i < i32.i64 m do
        let dist = sumSqrsSeq query (refs[i].1) in
        if dist >= knns[k-1].1 then knns -- early exit
        else let ref_ind = refs[i].0 in
              let (_, _, knns') =
                loop (dist, ref_ind, knns) for j < k do
                  let cur_nn = knns[j] in
                  if dist >= cur_nn.1
                       then (dist, ref_ind, knns)
                  else let tmp_ind = cur_nn.0
                       let knns[j] = (ref_ind, dist)
                       let ref_ind = tmp_ind
                       in  (cur_nn.1, ref_ind, knns)
              in  knns'

def bruteNNs [m] [n] [d] (k: i64) (test_set: [m][d]f32) (queries: [n][d]f32) =
  -- Setup for loop 
  let init_knns = replicate n (replicate k (-1i32, f32.inf))

  let refs_ind = iota32 m 
  let refs = zip refs_ind test_set

  let f query curr_nn = 
         bruteForce2 query curr_nn refs

  let new_knns = map2 f queries init_knns
    --loop cur_nns =  init_knns for i < n do
    --  let new_knn = bruteForce2 queries[i] cur_nns[i] refs
    --  let cur_nns' = cur_nns with [i] = copy new_knn
    --  in cur_nns' 


  let (k_inds, _) =  unzip <| map (\i_knn -> unzip i_knn) new_knns
  in k_inds

def main [m] [n] [d] (k: i64) (test_set: [m][d]f32) (queries: [n][d]f32) =
  bruteNNs k test_set queries