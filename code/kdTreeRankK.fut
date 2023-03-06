-- ==
-- entry: main
--
-- compiled input @ valid-data/kdtree-ppl-32-m-2097152.in
-- output @ valid-data/kdtree-ppl-32-m-2097152.out

import "lib/github.com/diku-dk/sorts/radix_sort"
import "util"

def iota32 n = (0..1..<i32.i64 n) :> [n]i32

local def closestLog2 (p: i32) : i32 =
    if p<=1 then 0
    else let (_,res) = loop (q,r) = (p,0) 
                       while q > 1 do
                            (q >> 1, r+1)
         let err_down = p - (1 << res)
         let err_upwd = (1 << (res+1)) - p
         in  if err_down <= err_upwd
             then res else res+1

-- m: the number of reference points
-- defppl: the default number of points per leaf
-- result: (height of tree without leaves, number of points per leaf)
def computeTreeShape (m: i32) (defppl: i32) : (i32, i32, i32, i32) =
    let def_num_leaves = (m + defppl - 1) / defppl
    let hp1 = closestLog2 def_num_leaves in
    if hp1 <= 0 then (-1, 0, m, m)
    else let h = hp1 - 1
         let num_leaves = 1 << (h+1)
         let ppl = (m + num_leaves - 1) / num_leaves
         in  (h, num_leaves-1, ppl, num_leaves*ppl)

local def updateBounds [n] [d2] (level: i32) (median_dims: [n]i32) (median_vals: [n]f32)
                     (node_ind: i32) (lubs_cur: *[d2]f32) : *[d2]f32=
    let d = d2 / 2
    let ancestor = 0
    let (_, res) =
      loop (ancestor,lubs_cur) for i < level do
        let k = level - i - 1
        let ancestor_child = compute_Kth_ancestor k node_ind
        let anc_dim = median_dims[ancestor]
        let lub_ind = if  (ancestor_child & 1) == 0
                      then anc_dim
                      else i32.i64 d+anc_dim
                      -- if right node, then update lower bound
        let anc_med = median_vals[ancestor]
        let lubs_cur[lub_ind] = if !(f32.isinf anc_med) then anc_med else lubs_cur[lub_ind]
        in  (ancestor_child, lubs_cur)
    in  res

local def findClosestMed [n] (cur_dim: i32) (median_dims: [n]i32) (node_ind: i32) : i32 =
    let cur_node = node_ind
    let res_ind  = -1i32
    let (_, res) =
        loop (cur_node, res_ind) 
          while (cur_node != 0) && (res_ind == (-1i32)) do
            let parent = getParent cur_node
            let res_ind = if median_dims[parent] == cur_dim then parent else -1
            in  (parent, res_ind)
    in  res

-- height: the height of the tree excluding leaves
-- q: the number of internal tree nodes (i.e., without leaves)
-- ppl: number of reference points per leaf
-- m' : the number of reference points hold in the tree
--    The following invariants hold: #leaves = q + 1 AND m' = (q+1) * ppl
-- input:  the d-dimensional array of reference points from which the tree is constructed
-- result: a tuple of five arrays
--         1. the reordered points (per leaf)
--         2. the indirect array that holds the original indices of each point
--         3. the index of the dimension that is split
--         4. the median value of the split dimension
-- !NOT NEEDED        5. the closest ancestor node index that splits the same dimension (or -1 if none)!
def mkKDtree [m] [d] (height: i32) (q: i64) (m' : i64)
                     (input: [m][d]f32) :
           (*[m'][d]f32, *[m']i32, *[q]i32, *[q]f32) =

         let num_pads = m' - m
         let input' = input ++ (replicate num_pads (replicate d f32.inf)) :> [m'][d]f32
         let indir  = iota32 m'
         
         let median_vals = replicate q 0.0f32
         let median_dims = replicate q (-1i32)
         let ( indir' : *[m']i32
             , median_dims': *[q]i32
             , median_vals': *[q]f32
             ) =
           loop ( indir  : *[m']i32
                , median_dims: *[q]i32
                , median_vals: *[q]f32 )
             for lev < (height+1) do
               let nodes_this_lvl = 1 << i64.i32 lev
               let pts_per_node_at_lev = m' / nodes_this_lvl
               let indir2d = unflatten nodes_this_lvl pts_per_node_at_lev indir

               -- dimensions to be split for each node at this level is equal to the level
               let med_dims =  if lev >= (i32.i64 d) then replicate nodes_this_lvl (i32.i64 (d-1))
                               else replicate nodes_this_lvl lev
                
               -- sort the chosen dimension for each node
               let chosen_columns = map2 (\indir_chunk dim ->
                                            map (\ind -> input'[ind, dim]
                                                ) indir_chunk
                                         ) indir2d med_dims

               ---! USE rank-k-flat-opt --- 
               let med_vals_rank = computeMedianWithRankK chosen_columns
               -- Call rank-k to get med_vals
               -- Somehow create indir2d'. I need to have the local indicies of sorting chosen_columns
                --  where each node is split by the median. Then I can give these indicies to "let indir2d'"
               -- DONE?
               let initial_local_ind = replicate nodes_this_lvl (iota32 pts_per_node_at_lev) 
               let (sorted_dim_2d, sort_inds_2d) =
                    map2 zip chosen_columns (replicate nodes_this_lvl (iota32 pts_per_node_at_lev))
                    |> map (radix_sort_float_by_key (\(l,_) -> l) f32.num_bits f32.get_bit)
                    |> map unzip |> unzip
               ---! END ---
               ---! I NEED TO CREATE A NEW sort_inds_2d AFTER I HAVE THE MEDIANS ---
               let med_vals = map  (\sorted_dim -> 
                                        let mi = pts_per_node_at_lev/2
                                        in (sorted_dim[mi] + sorted_dim[mi-1])/2
                                   ) sorted_dim_2d
               ---! I NEED TO CREATE A NEW sort_inds_2d ! ---
               let indir2d' = map2(\ indir_chunk sort_inds ->
                                        map (\ind -> indir_chunk[ind]) sort_inds
                                  ) indir2d sort_inds_2d
            
               -- scatter the values of this level in the global result arrays
               let this_lev_inds = map (+ (nodes_this_lvl-1)) (iota nodes_this_lvl)
               let median_dims' = scatter median_dims this_lev_inds med_dims
               let median_vals' = scatter median_vals this_lev_inds med_vals
               let indir'' = flatten indir2d' :> *[m']i32

               in  (indir'', median_dims', median_vals')

         let input'' = map (\ ind -> map (\k -> input'[ind, k]) (iota32 d) ) indir' :> *[m'][d]f32
         in  (input'', indir', median_dims', median_vals')


def main0 (m: i32) (defppl: i32) =
    computeTreeShape m defppl

def main [m] [d] (defppl: i32) (input: [m][d]f32) =
    let (height, num_inner_nodes, _, m') = computeTreeShape (i32.i64 m) defppl
    let (leafs, indir, median_dims, median_vals) =
        mkKDtree height (i64.i32 num_inner_nodes) (i64.i32 m') input
    let r = i64.i32 (m' / 64)
    in  (leafs, indir, median_dims, median_vals)