import "kdTreeRankK"


def sumSqrsSeq [d] (xs: [d]f32) (ys: [d]f32) : f32 =
    loop (res) = (0.0f32) for (x,y) in (zip xs ys) do
        let z = x-y in res + z*z

def bruteForce [m][d][k] (query: [d]f32) 
                         (knns0: [k](i32,f32))
                         (refs: [m](i32,[d]f32))
                       : [k](i32,f32) =
    --if query[0] == f32.lowest then copy knns else
    loop (knns) = (copy knns0)
      --for (i,refpt) in (zip (iota m) refs) do
      --  let dist = f32.sqrt <| sumSqrsSeq query refpt in
      for i < i32.i64 m do
        let dist = f32.sqrt <| sumSqrsSeq query (refs[i].1) in
        if dist > knns[k-1].1 then knns -- early exit
        --else if dist == 0.0 then knns  -- Causes compile error?? 0.23.1
        else let ref_ind = refs[i].0 in
             let (_, _, knns') =
               loop (dist, ref_ind, knns) for j < k do
                 let cur_nn = knns[j].1  in
                 if dist >= cur_nn
                 then (dist, ref_ind, knns)
                 else let tmp_ind = knns[j].0
                      let knns[j] = (ref_ind, dist)
                      let ref_ind = tmp_ind
                      in  (cur_nn, ref_ind, knns)
             in  knns'

-- give the median dimensions and values for each internal k-d tree node,
--   finds the leaf to which the query naturally belongs to.
def findLeaf [q][d] (median_dims: [q]i32) (median_vals: [q]f32)
                    (height: i32) (query: [d]f32) =
  let leaf =
    loop (node_index) = (0)
      while !( node_index >= ((1 << (height+1)) - 1)) do
        if query[median_dims[node_index]] < median_vals[node_index]
            then
                (node_index+1)*2-1
            else
                (node_index+1)*2
  let leaf_val = leaf - i32.i64 q
  in i64.i32 leaf_val
      --let qleafs = map (\l -> findLeaf median_dims median_vals height l) input

def findAllPaths [m] (path : [m]i32) (leaf_num : i32) =
    --let init_leafs_V = replicate m leaf_num
    map (\i ->let contact = path[i]
              let offset = if (contact == 1) then -(1 <<(m - i - 1))
                                          else (1 <<(m - i - 1))
              in leaf_num + (i32.i64 offset)
              ) (iota m)

let main [m] [d] (k: i64) (defppl: i32) (input: [m][d]f32) =
    let init_knns = replicate m (replicate k (-1i32, f32.inf))
    --- Build tree
    let (height, num_inner_nodes, _, m') = computeTreeShape (i32.i64 m) defppl
    let m'64 = i64.i32 m'
    let init2_knns = init_knns ++ (replicate (m'64 - m) (replicate k (-1i32, 0.0))) :> [m'64][k](i32, f32)
    let (leafs, indir, median_dims, median_vals) =
            mkKDtree height (i64.i32 num_inner_nodes) (m'64) input
    --- The initial leaf and its path
    let leaf_numbers = map (\i -> i/2) (iota32 m'64)
    let height1 = i64.i32 (height+1)
    --- Maybe dont reverse for readability 
    --- since findAllPaths have to reverse aswell then
    let path_arrs = map (\l -> map (\p -> (l / (2**p)) % 2) (reverse (iota32 height1))) leaf_numbers
    --- The leafs that have contact = 1 for each leaf.
    let defppl64 = i64.i32 defppl
    let Vis = map2 (\pa lnum -> (findAllPaths pa lnum) ) 
                      path_arrs[::defppl64] leaf_numbers[::defppl64]
    let leafs_with_ind = zip indir leafs
    let leafs2d        = unflatten (m'64 / defppl64) defppl64 leafs_with_ind
    let Vi_vals4d = map (\vi -> map (\ind -> leafs2d[ind]) vi) Vis
    --- flatten the values so each Vi is just an array of points and not leaves
    let Vi_leafnum_m1  = defppl64 * height1 
    let Vi_vals3d      = map (\vi -> flatten vi :> [Vi_leafnum_m1](i32,[d]f32) ) Vi_vals4d
    let Vi_for_queries = map (\i -> Vi_vals3d[i]) leaf_numbers :> [m'64][](i32, [d]f32)
    --let Vim = map2 (\i vi -> leafs2d[i] ++ vi) leaf_numbers Vi_for_queries  
    let knns2 =  map3 (\query vi knn0 -> 
                        bruteForce query knn0 vi 
                    ) leafs Vi_for_queries init2_knns
    let knn3  = map3 (\query i knn0 -> 
                        bruteForce query knn0 leafs2d[i] 
                    ) leafs leaf_numbers knns2
    let (knn_inds, knn_dists) = unzip <| map (\knn_tup -> unzip knn_tup) knn3
    in  (leafs, indir, median_dims, median_vals, Vis, knn_inds, knn_dists)