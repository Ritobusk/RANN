def main (k: i64) (d:i64) (test_size:i64) (q_size:i64) :
                      (i64, [](i32, []f32), [](i32, []f32))=
  let test_knns = map (\ind -> 
                          let vals  = map (f32.i64) (iota d)
                          let ind32 = i32.i64 ind
                          in (ind32, vals)) (iota test_size)
  let q_knns = map (\ind -> 
                          let vals  = map (f32.i64) (iota d) --(0..1..<f32.i64 d) :> [d]f32
                          let ind32 = i32.i64 ind
                          in (ind32, vals)) (iota q_size)
  in (k, test_knns, q_knns)

