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
      --let qleafs = map (\l -> findLeaf median_dims median_vals height l) input

def findAllPaths [m] (path : [m]i32) (leaf_num : i32) =
    --let init_leafs_V = replicate m leaf_num
    map (\i ->let contact = path[i]
              let offset = if (contact == 1) then -(1 <<(m - i - 1))
                                          else (1 <<(m - i - 1))
              in leaf_num + (i32.i64 offset)
              ) (iota m)

let main [m] [d] (defppl: i32) (input: [m][d]f32) =
    --- Build tree
    let (height, num_inner_nodes, _, m') = computeTreeShape (i32.i64 m) defppl
    let (leafs, indir, median_dims, median_vals) =
            mkKDtree height (i64.i32 num_inner_nodes) (i64.i32 m') input
    --- The initial leaf and its path
    let leaf_numbers = map (\i -> i/2) (iota32 m)
    let height1 = i64.i32 (height+1)
    --- Maybe dont reverse for readability 
    --- since findAllPaths have to reverse aswell then
    let path_arrs = map (\l -> map (\p -> (l / (2**p)) % 2) (reverse (iota32 height1))) leaf_numbers
    --- The leafs that have contact = 1 for each leaf.
    let defppl64 = i64.i32 defppl
    let Vis = map2 (\pa lnum -> (findAllPaths pa lnum) ) 
                      path_arrs[::defppl64] leaf_numbers[::defppl64]
    let leafs2d = unflatten (m / 2) defppl64 leafs
    let Vi_vals4d = map (\vi -> map (\ind -> leafs2d[ind]) vi) Vis
    --- flatten the values so each Vi is just an array of points and not leaves
    let Vi_leafnum_m1 = defppl64 * height1
    let Vi_vals3d = map (\vi -> flatten vi :> [Vi_leafnum_m1][d]f32 ) Vi_vals4d
    in  (leafs, indir, median_dims, median_vals, leaf_numbers , path_arrs, Vis, Vi_vals3d)