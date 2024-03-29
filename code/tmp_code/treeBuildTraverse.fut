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
    map (\i ->let contact = path[i]
              let offset = if (contact == 1) then -(1 <<( i)) -- -(1 <<(m - i - 1)) if reversed
                                          else (1 <<(i))
              in leaf_num + (i32.i64 offset)
              ) (iota m)

let main [m] [d] (k: i64) (defppl: i32) (input: [m][d]f32) =
    let init_knns = replicate m (replicate k (-1i32, f32.inf))
    --- Build tree
    let (height, num_inner_nodes, _, m') = computeTreeShape (i32.i64 m) defppl
    let m'64 = i64.i32 m'
    let defppl64 = i64.i32 defppl
    let init2_knns = init_knns ++ (replicate (m'64 - m) (replicate k (-1i32, f32.inf))) :> [m'64][k](i32, f32)
    let (leafs, indir, median_dims, median_vals) =
            mkKDtree height (i64.i32 num_inner_nodes) (m'64) input

    -- 2. Find the leaf to which each query "naturally" belongs
    --    If your set of querries is named `querries` this is
    --    achieved with a map:
    --      `let leaf_inds = map (findLeaf kdtree) querries`

    -- 3. sort `(zip querries leaf_inds)` in increasing order of
    --      their leaf_ind (second element)

    -- 4. for each querry compute its knns from its own leaf

    --- The initial leaf and its path
    let leaf_numbers = map (\i -> i/defppl) (iota32 m'64) -- per point in leaf
    let height1     = i64.i32 (height+1)
    let leafs_in_Vi = height1 + 1
    let Vi_leafnum  = defppl64 * leafs_in_Vi 
    let path_arrs = map (\l -> map (\p -> (l / (2**p)) % 2) ((iota32 height1))) leaf_numbers
    

    -- 5. have a loop which goes from [0 .. height - 1] which refines
    --    the nearest neighbors
    --
    -- loop (curr_nn_set) = (init2_knns) for i < height do
    --     let new_leaves = map (reverseBit i) leaf_numbers
    --     let better_nn_set = map3 bruteForce querries curr_nn_set new_leaves 
    --     in  better_nn_set
    --

    --- The leafs that have contact <= 1 for each leaf.
    let Vis = map2 (\pa lnum -> let cont1 =(findAllPaths pa lnum)
                                  in [leaf_numbers[lnum]] ++ cont1 :> [leafs_in_Vi]i32) 
                      path_arrs[::defppl64] leaf_numbers[::defppl64]
    let leafs_with_ind = zip indir leafs
    let leafs2d        = unflatten (m'64 / defppl64) defppl64 leafs_with_ind
    let Vi_vals4d = map (\vi -> map (\ind -> leafs2d[ind]) vi) Vis --try using flatten in here!

    --- flatten the values so each Vi is just an array of points and not leaves
    let Vi_vals3d      = map (\vi -> flatten vi :> [Vi_leafnum](i32,[d]f32) ) Vi_vals4d
    let Vi_for_queries = map (\i -> Vi_vals3d[i]) leaf_numbers :> [m'64][](i32, [d]f32)
    
    --- Find knns in Vi for each point i!
    let knns2 =  map3 (\query vi knn0 -> 
                        bruteForce query knn0 vi 
                    ) leafs Vi_for_queries init2_knns

    --- Now I need to scatter knns2 vals with the indir indicies 

    let (knn_inds, knn_dists) = unzip <| map (\knn_tup -> unzip knn_tup) knns2
    in  (leafs, indir, median_dims, median_vals, Vis, knn_inds, knn_dists)
