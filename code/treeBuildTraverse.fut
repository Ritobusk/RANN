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
  --let pad_query = query ++ (replicate  (q-d) query[d-1]) :> [q]f32
  --let path_arr  = map2 (\pq mv -> if pq < mv then 0 else 1) pad_query median_vals
  in ((leaf - i32.i64 q), path_arr)

let main [m] [d] (defppl: i32) (input: [m][d]f32) =
    let (height, num_inner_nodes, _, m') = computeTreeShape (i32.i64 m) defppl
    let (leafs, indir, median_dims, median_vals) =
        mkKDtree height (i64.i32 num_inner_nodes) (i64.i32 m') input
    let (qleafs, path_arrs) = unzip (map (\l -> findLeaf median_dims median_vals height l) input)
    in  (leafs, indir, median_dims, median_vals, 0, qleafs, path_arrs)