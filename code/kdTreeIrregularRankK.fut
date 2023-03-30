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
def computeTreeShape (m: i32) (defppl: i32) : (i32, i32, i32) =
    let def_num_leaves = (m + defppl - 1) / defppl
    let hp1 = closestLog2 def_num_leaves in
    if hp1 <= 0 then (-1, 0, m)
    else let h = hp1 - 1
         let num_leaves = 1 << (h+1)
         let ppl = (m + num_leaves - 1) / num_leaves
         in  (h, num_leaves-1, num_leaves*ppl)



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
def mkKDtree [m] [d] (height: i32) (q: i64) (m' : i64)
                     (input: [m][d]f32) :
           (*[m'][d]f32, *[m']i32, *[q]i32, *[q]f32, *[q]i32) =

    let num_pads = m' - m
    let input' = input ++ (replicate num_pads (replicate d f32.inf)) :> [m'][d]f32
    let indir  = iota32 m'
    
    let median_vals = replicate q 0.0f32
    let median_dims = replicate q (-1i32)
    let shape_arr   = replicate q (i32.i64 m')
    let ( indir' : *[m']i32
        , median_dims': *[q]i32
        , median_vals': *[q]f32
        , shape_arr':   *[q]i32 
        ) =
        loop ( indir  : *[m']i32
                , median_dims: *[q]i32
                , median_vals: *[q]f32
                , shape_arr:   *[q]i32)
                for lev < 1 do

            let nodes_this_lvl = 1 << i64.i32 lev
            let average_pts_per_node_at_lvl = m' / nodes_this_lvl

            let start_shp = if lev == 0 then 0 else nodes_this_lvl
            let end_shp   = start_shp + nodes_this_lvl

            let shp_this_lvl       = shape_arr[start_shp:end_shp] :> [nodes_this_lvl]i32
            let scan_shp_this_lvl  = scan (+) 0 shp_this_lvl  


            -- Dimension to be split for each node at this level 
            let med_dim = lev % (i32.i64 d) 

            -- Values of chosen dimension   
            let chosen_column = map (\ind -> input'[ind, med_dim]) indir 

            -- Calculate the median
            let shp_flag_arr = idxs_to_flags shp_this_lvl :> [m']bool
            let sums_of_chosen_vals = segmented_scan (+) (0.0f32) shp_flag_arr chosen_column
            let means = map (\ind -> sums_of_chosen_vals[ind-1] / 
                                                    (f32.i64 average_pts_per_node_at_lvl)) scan_shp_this_lvl
            let ks  = map (\s -> s/2) shp_this_lvl
            let II1 = mkII1 shp_this_lvl :> [m']i32
            let medians_this_lvl = rankSearchBatch means ks shp_this_lvl II1 (copy chosen_column) 
                -- Got an error when not using copy. Don't know why.

            --- Calculate the mask for partition3L
            let medians_for_each_elem = segmented_replicate shp_this_lvl medians_this_lvl :> [m']f32
            let mask_arr = map2 (\p_val m_val -> p_val < m_val) chosen_column medians_for_each_elem

            --- Partition to split each node
            -- let (indir'', new_split) = partition3L mask_arr shp_flag_arr scan_shp_this_lvl <| (zip shp_this_lvl indir)

            --- For the shape I just need to 'weave' isT with a 
            ---   map2 (\shp_val T_val -> shp_val - T_val ) 
            --- let tmp = map2 (\shp_val T_val -> shp_val - (i32.i64 T_val) ) shp new_split
            --- let new_shape_ind = iota (nodes_this_lvl << 1)
            --- let new_shape = map (\ind -> if (ind % 2) == 0 then new_split[ind/2] else tmp[ind/2]) new_shape_ind


            let this_lev_inds = map (+ (nodes_this_lvl-1)) (iota nodes_this_lvl)
            let median_dims' = scatter median_dims this_lev_inds (replicate nodes_this_lvl med_dim)
            let median_vals' = scatter median_vals this_lev_inds medians_this_lvl
            in  (indir, median_dims', median_vals', shape_arr)

    let input'' = map (\ ind -> map (\k -> input'[ind, k]) (iota32 d) ) indir' :> *[m'][d]f32
    in  (input'', indir', median_dims', median_vals', shape_arr')


def main0 (m: i32) (defppl: i32) =
    computeTreeShape m defppl

def main [m] [d] (defppl: i32) (input: [m][d]f32) =
    let (height, num_inner_nodes, m') = computeTreeShape (i32.i64 m) defppl
    let (leafs, indir, median_dims, median_vals, shp_arr) =
        mkKDtree height (i64.i32 num_inner_nodes) (i64.i32 m') input
    in  (leafs, indir, median_dims, median_vals, shp_arr)
