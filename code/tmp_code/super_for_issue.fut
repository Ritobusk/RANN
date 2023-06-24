def bruteForce [m][d][k] (query: [d]f32) 
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
                  if cur_nn.0 == ref_ind 
                       then (f32.inf, -1, knns) -- already encountered
                  else if dist >= cur_nn.1
                       then (dist, ref_ind, knns)
                  else let tmp_ind = cur_nn.0
                       let knns[j] = (ref_ind, dist)
                       let ref_ind = tmp_ind
                       in  (cur_nn.1, ref_ind, knns)
              in  knns'


def main [m] [n] [d] (k: i64) (test_knns: [m](i32, [d]f32)) (queries_knn: [n](i32,[d]f32)) =
  let (knn_inds_q, _) =  unzip <| map (\i_knn -> unzip i_knn) queries_knn
  let (knn_inds_t, _) =  unzip <| map (\i_knn -> unzip i_knn) test_knns
  let supercharging = 
    --let ksqr = k**2
    let (super_knn_inds_seq) =
      loop (curr_knns_q) = (new_knns_q) for i < k do
        let k_inds = map (\knn_ind_q -> map (\q -> knn_inds_t[knn_ind_q[i],q]) (iota k)) knn_inds_q
        let k_points = --map (\inds -> map (\ind -> 
                       --                 let test_point = map (\k_ind -> test_set[ind, k_ind]) (iota d)
                       --                 in (ind, test_point)
                       --                 ) inds) k_inds
                       map (\inds -> map (\ind -> (ind, test_set[ind])) inds) k_inds
        in  map2 (\refs ind -> bruteForce queries[ind] curr_knns_q[ind] refs ) k_points (iota n)
    in super_knn_inds_seq
  let (super_knn_inds, _) =  unzip <| map (\i_knn -> unzip i_knn) supercharging
  in (super_knn_inds)
