-- ==
-- compiled input { 
--   [2.0f32, -1.0f32, -1.0f32, 2.0f32]
-- } 
-- output { [4.0f32, 1.0f32, 1.0f32, 4.0f32, 9.0f32] }

let main [n]
         (point : [n]f32)  : [n]f32 =
	map (\x -> x**2) point