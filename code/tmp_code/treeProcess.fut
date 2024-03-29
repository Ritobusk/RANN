-- Randomized Approximate K-D Trees
-- ==
-- compiled input @ dataProcess.in

-- output @ LocVolCalib-data/small.out

import "kdTreeRankK"
import "lib/github.com/diku-dk/sorts/radix_sort"

def sumSqrsSeq [d] (xs: [d]f32) (ys: [d]f32) : f32 =
    loop (res) = (0.0f32) for (x,y) in (zip xs ys) do
        let z = x-y in res + z*z

-- Bruteforce knn that returns 
def bruteForce [m][d][k] (query: [d]f32) 
                         (knns0: [k](i32,f32))
                         (refs: [m](i32,[d]f32))
                       : [k](i32,f32) =
    loop (knns) = (copy knns0)
      for i < i32.i64 m do
        let dist = sumSqrsSeq query (refs[i].1) in
        if dist >= knns[k-1].1 then knns -- early exit !!SHOULD BE >=
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

-- Given the median dimensions and values for each internal k-d tree node,
--   finds the leaf to which the query naturally belongs to.
def findLeaf [q][d] (median_dims: [q]i32) (median_vals: [q]f32)
                    (height: i32) (query: [d]f32) =
  let leaf =
    loop (node_index) = (0)
      while !( node_index >= ((1 << (height+1)) - 1)) do -- Maybe remove ! and replace >= with < 
        if query[median_dims[node_index]] < median_vals[node_index]
            then
                (node_index+1)*2-1
            else
                (node_index+1)*2
  let leaf_val = leaf - i32.i64 q
  in i64.i32 leaf_val

def reverseBit (num: i64) (bit: i64) : i64 =
    let bit_val  = 1 << bit
    let bit_rm = (num / bit_val) % 2
    let rev_val = if bit_rm == 1 then -bit_val
                                  else  bit_val
    in num + rev_val
    
def log2Int (n : i64) : i32 =
  let (_, res) =
    loop (n, r) = (n, 0i32)
      while n > 1 do
        (n >> 1, r+1)
  in res 

-- ToDos:
-- 1. Pass the desired `height` of the kd-tree as program parameter instead of `defppl`
-- 2. Fix the construction of the kd-tree to not split in the middle of a set of points
--      having the same value as the median.
--    Representation of kd-tree:
--      - You know that the tree is still fully balanced, but the number of points in
--        each leaf does not have to be the same.
--      - Therefore, when you are building the kd-tree, you need to keep track of the
--        shape of each level of the tree, i.e., level `i` has `2^i` nodes, hence it
--        corresponds to an arrray of `2^i` subarrays. The sum of the shape array has
--        to equal the number of reference points.
--      - on the rightside, you will probably pad with some empty arrays.
-- 3. adjust the rest of the code, i.e., of this file, to work with the new
--    kd-tree representation.
--
def main [m] [n] [d] (k: i64) (defppl: i32) (input: [m][d]f32) (queries: [n][d]f32) =
    let init_knns = replicate n (replicate k (-1i32, f32.inf))
    --- Build tree (height is "0-indexed")
    let (height, num_inner_nodes, _, m') = computeTreeShape (i32.i64 m) defppl
    let m'64 = i64.i32 m'
    let defppl64 = i64.i32 defppl
    let (leafs, indir, median_dims, median_vals) =
            mkKDtree height (i64.i32 num_inner_nodes) (m'64) input
    let num_leafs = length leafs

    -- 2. Find the leaf to which each query "naturally" belongs
    --    If your set of querries is named `querries` this is
    --    achieved with a map:
    --      `let leaf_inds = map (findLeaf kdtree) querries`
    let queries_init_leafs = map (findLeaf median_dims median_vals height) queries

    -- 3. sort `(zip querries leaf_inds)` in increasing order of
    --      their leaf_ind (second element)
    let (sorted_query_with_ind, sorted_query_leaf) = 
      let q_with_ind = zip queries (iota n) 
      let q_ind_with_leaf = zip q_with_ind queries_init_leafs
      in (radix_sort_int_by_key (\(_,l) -> l) (log2Int num_leafs) i64.get_bit q_ind_with_leaf) |> unzip
    
    let (sorted_query, sorted_query_ind) = unzip sorted_query_with_ind 

    -- 4. for each querry compute its knns from its own leaf
    let leafs_with_ind = zip indir leafs
    let leafs2d = unflatten (m'64 / defppl64) defppl64 leafs_with_ind

    -- For flat: use values in sorted_query_leaf to index into scanned shape array. Use that to "slice" the leafs correctly
    let knn_nat_leaf = map3 (\q q_ind l_ind -> bruteForce q init_knns[q_ind] leafs2d[l_ind] ) 
                                                              sorted_query sorted_query_ind sorted_query_leaf

    -- 5. have a loop which goes from [0 .. height - 1] which refines
    --    the nearest neighbors
    --
    -- loop (curr_nn_set) = (init2_knns) for i < height do
    --     let new_leaves = map (reverseBit i) leaf_numbers
    --     let better_nn_set = map3 bruteForce querries curr_nn_set new_leaves 
    --     in  better_nn_set

    let new_knns_sorted =
      loop (curr_nn_set) = (knn_nat_leaf) for i < (i64.i32 height + 1) do
          let new_leaves = map (\l_num -> reverseBit l_num i) sorted_query_leaf
          let better_nn_set = map3 (\q q_knn l_ind -> bruteForce q q_knn leafs2d[l_ind])
                                                        sorted_query curr_nn_set new_leaves
          in better_nn_set

    let new_knns = scatter (init_knns) sorted_query_ind new_knns_sorted -- instead of copy just replicate with 0 and the right size
    let (new_knn_ind, new_knn_dists) = unzip new_knns[0]

    in  (leafs, indir, median_vals, sorted_query, 
                      sorted_query_ind, new_knn_ind, new_knn_dists)
