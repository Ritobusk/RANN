import "kdTreeRankK"

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



let main [m] [d] (defppl: i32) (input: [m][d]f32) =
    let (height, num_inner_nodes, _, m') = computeTreeShape (i32.i64 m) defppl
    let (leafs, indir, median_dims, median_vals) =
            mkKDtree height (i64.i32 num_inner_nodes) (i64.i32 m') input
    let leaf_numbers = map (\i -> i/2) (iota32 m)
    --let qleafs = map (\l -> findLeaf median_dims median_vals height l) input
    let height1 = i64.i32 (height+1)
    let path_arrs = map (\l -> map (\p -> (l / (2**p)) % 2) (reverse (iota32 (height1)))) leaf_numbers --:> [m][height1]i64
    in  (leafs, indir, median_dims, median_vals, leaf_numbers , path_arrs)