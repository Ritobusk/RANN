import "kdTreeRankK"
import "lib/github.com/diku-dk/sorts/radix_sort"

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

def findAllPaths [m] (path : [m]i32) (leaf_num : i32) =
    map (\i ->let contact = path[i]
              let offset = if (contact == 1) then -(1 <<( i)) -- -(1 <<(m - i - 1)) if reversed
                                          else (1 <<(i))
              in leaf_num + (i32.i64 offset)
              ) (iota m)

let main [m] [n] [d] (k: i64) (defppl: i32) (input: [m][d]f32) (queries: [n][d]f32) =
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
    let queries_init_leafs = map (\l -> findLeaf median_dims median_vals height l) queries

    -- 3. sort `(zip querries leaf_inds)` in increasing order of
    --      their leaf_ind (second element)
    let (sorted_query, sorted_query_ind) = 
      zip queries queries_init_leafs |>
      (radix_sort_float_by_key (\(_,l) -> l) i64.num_bits i64.get_bit) |> unzip
    
    -- 4. for each querry compute its knns from its own leaf
    

    -- 5. have a loop which goes from [0 .. height - 1] which refines
    --    the nearest neighbors
    --
    -- loop (curr_nn_set) = (init2_knns) for i < height do
    --     let new_leaves = map (reverseBit i) leaf_numbers
    --     let better_nn_set = map3 bruteForce querries curr_nn_set new_leaves 
    --     in  better_nn_set
    --

    in  (leafs, indir, median_dims, median_vals, queries_init_leafs, sorted_query)
