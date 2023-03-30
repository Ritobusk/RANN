def log2 x = (loop (y,c) = (x,0i32) while y > 1i32 do (y >> 1, c+1)).1

def iota32 n = (0..1..<i32.i64 n) :> [n]i32

def partition2Ind [n] (cs: [n]bool) : ([n]i32, i32) =
    let tfs = map (\f -> if f then 1 else 0) cs
    let isT = scan (+) 0 tfs
    let ffs = map (\f -> if f then 0 else 1) cs
    let isF0 = scan (+) 0 ffs

    let i = isT[n-1]
    let isF = map (+ i) isF0 
    let inds = map3 (\ c iT iF ->
                        if c then iT-1 else iF-1
                    ) cs isT isF
    in (inds, i)

def sumSqrs [d] (xs: [d]f32) (ys: [d]f32) : f32 =
    map2 (\x y -> let z = x-y in z*z) xs ys |> reduce (+) 0.0f32

def sumSqrsSeq [d] (xs: [d]f32) (ys: [d]f32) : f32 =
    loop (res) = (0.0f32) for (x,y) in (zip xs ys) do
        let z = x-y in res + z*z

def gather1D 't [m] (arr1D: [m]t) (inds: [m]i32) : *[m]t =
    map (\ind -> arr1D[ind] ) inds

def gather2D 't [m][d] (arr2D: [m][d]t) (inds: [m]i32) : *[m][d]t =
    map (\ind -> map (\j -> arr2D[ind,j]) (iota d) ) inds

def scatter2D [m][k][n] 't (arr2D: *[m][k]t) (qinds: [n]i32) (vals2D: [n][k]t) : *[m][k]t =
  let nk = n*k
  let flat_qinds = map (\i -> let (d,r) = (i32.i64 i / i32.i64 k,
                                           i32.i64 i % i32.i64 k)
                              in i64.i32 (qinds[d]*i32.i64 k + r)
                       ) (iota nk)
  let res1D = scatter (flatten arr2D) flat_qinds ((flatten vals2D) :> [nk]t) 
  in  unflatten m k res1D 


def getParent (node_index: i32) = (node_index-1) / 2

def isLeaf (h: i32) (node_index: i32) =
    node_index >= ((1 << (h+1)) - 1)

-- the k'th ancestor of `node_ind` can be computed with
-- the formula: `(node_ind + 1 - (2^k)) / (2^k)`, for example
-- the parent           (k==1): `(node_ind - 1) / 2`
-- the grandparent      (k==2): `(node_ind - 3) / 4`
-- the grandgrandparent (k==3): `(node_ind - 7) / 8`
def compute_Kth_ancestor (k: i32) (node_ind: i32) =
    let tpk = 1 << k
    in  (node_ind + 1 - tpk) / tpk

def findNodeLevel (node: i32) : i32 =
 ( loop (lev, idx) = (0i32, node)
     while idx > 0i32 do
       (lev+1i32, getParent idx) ).0

-- given a tree `node1` at level `lev` and another tree leaf `leaf`,
-- this function computes the closest common ancestor of `node1` and `node2`
-- `h` is the height of the binary tree (without leaves).
def findClosestCommonAncestor (h: i32) (lev: i32) (node1: i32) (leaf: i32) : i32 =
    let node2 = compute_Kth_ancestor (h+1-lev) leaf
    let (res,_) =
      loop (node1, node2) while node1 != node2 do
        (getParent node1, getParent node2)
    in  res


def imap  as f = map f as
def imap2 as bs f = map2 f as bs
def imap3 as bs cs f = map3 f as bs cs
def imap4 as bs cs ds f = map4 f as bs cs ds
def ifilter as p = filter p as

def ones [q] 't (_xs: [q]t) = replicate q 1i32


--- Took these functions from 'Parallel Programming in Futhark' chapter 8
def segmented_scan 't [n] (g:t->t->t) (ne: t) (flags: [n]bool) (vals: [n]t): [n]t =
  let pairs = scan ( \ (v1,f1) (v2,f2) ->
                       let f = f1 || f2
                       let v = if f2 then v2 else g v1 v2
                       in (v,f) ) (ne,false) (zip vals flags)
  let (res,_) = unzip pairs
  in res

def replicated_iota [n] (reps:[n]i32) : []i32 =
  let s1 = scan (+) 0 reps
  let s2 = map (\i -> if i==0 then 0 else (i64.i32 s1[i-1])) (iota n)
  let tmp = scatter (replicate (i64.i32(reduce (+) 0 reps)) 0) (s2) (iota32 n) -- replace reduce with last elem of s1
  let flags = map (>0) tmp
  in segmented_scan (+) 0 flags tmp

---! changed this to not just work on int32. Will it cause problems?
def segmented_replicate [n] 't (reps:[n]i32) (vs:[n]t) : []t =
  let idxs = replicated_iota reps
  in map (\i -> vs[i]) idxs


def idxs_to_flags [n] (is : [n]i32) : []bool =
  let vs = segmented_replicate is (iota32 n)
  let m = length vs
  in map2 (!=) (vs :> [m]i32) ([0] ++ vs[:m-1] :> [m]i32)
---


-- Write this function in flattened form
def partition3L 't [n] [p]
        (mask : [n]bool) -- mask[i] == True => associated predicate holds on elem i
        (shp_flag_arr: [n]bool)
        (scan_shp: [p]i32)
        (shp : [p]i32, flat_arr : [n]t) -- representation of an irregular array of array
      : ([n]t, [p]i32) = -- result: the flat array reorganized & splitting point of each segment
  -- replace the dummy implementation
  let tfs = map (\f -> if f then 1 else 0) mask
  let isT = segmented_scan (+) 0 shp_flag_arr tfs
  let Ts_per_segment = map (\ind -> isT[ind-1]) scan_shp :> [p]i32
  let T_indicies = replicated_iota shp :> [n]i32
  let ffs = map (\f -> if f then 0 else 1) mask 

  let shp_flag_alt = scatter (copy shp_flag_arr) [0] [true] 
  ---! ffs_2 Doens not work since the first elem of shp_flag_arr is false
  let ffs_2 = map3  (\f_val is_new_segment segment_Tval ->
                      if is_new_segment then  (f_val + Ts_per_segment[segment_Tval])
                      else f_val
                    ) ffs shp_flag_alt T_indicies

  let isF  = segmented_scan (+) 0 shp_flag_arr ffs_2
  let inds = map3 (\c iT iF -> if c then i64.i32(iT -1)
                                       else i64.i32(iF -1)
                  ) mask isT isF
  let r =  scatter (replicate n flat_arr[0]) inds flat_arr

  in (r, Ts_per_segment) 

-- ==
-- compiled input {
--  [false,true,false,true,false,true,false,true,false,true,false,true]
--  [3,4,5]
--  [1,2,3,4,5,6,7,8,9,10,11,12]  
-- }
-- output {
-- [2,1,3,4,6,5,7,8,10,12,9,11]
-- [1,2,3] 
--}
let main [m] [n] (mask: [n]bool) (shp: [m]i32) (f_arr : [n]i32) =
  let scan_shp = scan (+) 0 shp
  let shp_flag_arr = (idxs_to_flags shp) :> [n]bool
  let (new_flat_arr, splits) =  partition3L mask shp_flag_arr scan_shp (shp, f_arr)
  in (new_flat_arr, splits)


def partition3 [ n ] 't 
              ( p : ( t -> bool )) 
              ( arr : [ n ] t ) : ([ n ]t , i64 ) = 
    let cs = map p arr 
    let tfs = map (\ f -> if f then 1 else 0) cs
    let isT = scan (+) 0 tfs
    let i = isT[n-1] 
    let ffs = map (\f -> if f then 0 else 1) cs 
    let isF = map (+i) <| scan (+) 0 ffs 
    let inds = map3 (\c iT iF -> if c then iT -1
                                       else iF -1
                    ) cs isT isF
    let r = scatter (replicate n arr[0]) inds arr
    in (r , i )


let mkFlagArray 't [m] 
            (aoa_shp: [m]i32) (zero: t)
            (aoa_val: [m]t  ) : []t =
  let shp_scn = scan (+) 0 aoa_shp
  let aoa_len = shp_scn[m-1]
  let shp_ind = imap2 aoa_shp (indices aoa_shp)
                      (\ s i ->
                         if s==0 then -1i64
                         else if i==0 then 0i64
                         else i64.i32 shp_scn[i-1]
                      )
  let flags = scatter (replicate (i64.i32 aoa_len) zero)
                      shp_ind aoa_val
  in flags

let mkII1 [m] (shp: [m]i32) : *[]i32 =
    let flags = mkFlagArray shp 0i8 (replicate m 1i8)
    in  map i32.i8 flags
     |> scan (+) 0i32

-- meds: hopefully a decent estimate of the median values for each partition
-- ks:   the k-th smallest element to be searched for each partition (starting from 1)
-- shp, II1, A:  the rep of the iregular array: shape, II1-helper (plus 1) and flat data
-- More precisely, if we have m segments II1 will have the same length as the flat A,
--   and each element will indicate the segment (number plus one) in which the current
--   element of A rezides.
-- E.g., assuming shp = [3,5,7], then the length of A and II1 is 3+5+7=15,
--       if we want the median, ks = [2, 3, 4], and
--       II1 = [1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3].
--
def rankSearchBatch [m][n] (meds: [m]f32) (ks: [m]i32)
                           (shp: [m]i32) (II1: *[n]i32) (A: *[n]f32) : [m]f32 =
  let II1_bak = replicate n 0i32
  let A_bak = replicate n 0f32
  let res = replicate m 0f32
  let q = 0i64

  let (_, _shp, _,_,_,_,_, res) =
    loop (ks : [m]i32, shp : [m]i32, II1, II1_bak, A, A_bak, q, res)
      while (length A > 0) do
        -- compute helpers based on shape
        let shp_sc = scan (+) 0 shp

        -- F(let pivot = last A)
        let pivots =
            imap2 shp_sc (indices shp_sc)
              (\ off i -> if q == 0i64
                          then meds[i]
                          else if off == 0
                               then 0f32
                               else A[off - 1]
              ) |> opaque

        -- compute lt_len and eq_len by means of histograms:
        let h_inds =
            imap2 II1 A
                  (\ sgmindp1 a ->
                    let sgmind = sgmindp1 - 1
                    let pivot  = pivots[sgmind]
                    let h_ind  = sgmind << 1
                    in  i64.i32 <|
                          if a < pivot then h_ind
                          else if pivot == a then h_ind + 1
                          else -1i32
                  )
        let h_vals = ones A
        let lens = reduce_by_index (replicate (2*m) 0i32) (+) 0i32 h_inds h_vals

        --
        let (shp', kinds, ks') =
          imap2 ks (indices ks)
            (\ k i ->
                if k < 0 then (0, 3i8, -1) -- already processed
                else let lt_len = lens[i << 1] in
                     if k < lt_len then (lt_len, 0i8, k)
                     else let eq_len = lens[ (i << 1) + 1]
                          let lteq_len = lt_len + eq_len in
                          if k < lteq_len then (0, 1i8, -1)
                          else (shp[i] - lteq_len, 2i8, k - lteq_len)
            )
          |> unzip3

        -- write the subarrays that have finished
        let (scat_inds, scat_vals) =
            imap2 (indices kinds) kinds
                  (\ i knd ->
                    if knd == 1i8
                    then (i, pivots[i])
                    else (-1, 0.0)
                  )
            |> unzip
        let res' = scatter res scat_inds scat_vals

        -- use a filter to extract elements
        let keepElem sgmindp1 a =
                let sgmind = sgmindp1 - 1
                let pivot = pivots[sgmind]
                let kind  =  kinds[sgmind] in
                if (a < pivot && kind == 0) then true
                else if (a > pivot && kind == 2) then true
                else false

        let conds = map2 keepElem II1 A |> opaque -- strange fusion with duplicating computation

        let tmp_inds = map i32.bool conds
                    |> scan (+) 0i32
        let tot_len = i64.i32 (last tmp_inds)
        let scat_inds = imap2 conds tmp_inds
              (\ c ind -> if c then i64.i32 (ind-1) else -1i64)
        let A'   = scatter A_bak scat_inds A
        let II1' = scatter II1_bak scat_inds II1
        let II1''= II1'[:tot_len]
        let A''  = A'[:tot_len]
        
        in  (ks', shp', II1'', II1, A'', A, q+1, res')
  in res


def computeMedianWithRankK [m][n] (ass: [m][n]f32) =
    let mins = map (reduce_comm f32.min f32.highest) ass
    let maxs = map (reduce_comm f32.max f32.lowest ) ass
    let means = map2 (+) mins maxs |> map (/ 2) |> opaque

    let k  = i32.i64 (n / 2)
    let ks = replicate m k
    let N  = m * n
    let shp = replicate m (i32.i64 n)

    -- II1 should in principle be computed with (mkII1 shp)
    let II1 = (tabulate_2d m n (\i _j -> i32.i64 i + 1) |> flatten) :> *[N]i32
    let A = copy (flatten ass) :> *[N]f32
    let medians = rankSearchBatch means ks shp II1 A
    in medians
