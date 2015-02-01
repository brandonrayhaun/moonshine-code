import Data.List

integerVectors :: Int -> Int -> [[Int]]
integerVectors 0 0 = [[]]
integerVectors _ 0 = []
integerVectors k n = 
	[x:t | x <- [-k..k], t <- integerVectors (k - abs x) (n-1)]
	
dot :: [Int] -> [Int] -> Int
dot iv cl = sum $ zipWith (*) iv cl

searchDecomposition :: [Int] -> Int -> [[Int]] -> [[Int]]
searchDecomposition cl s = filter $ \v -> v `dot` cl == s

searchMDecomposition :: [[Int]] -> [Int] -> [[Int]] -> [[Int]]
searchMDecomposition cls ss = filter $ \v -> map (dot v) cls == ss

s :: Int
s = 708938752

cl :: [[Int]]
cl = [[1, 248, 4123, 27000, 27000, 30628, 30875, 61256, 85995, 85995, 147250, 767637, 767637, 779247, 779247, 957125, 1707264, 1707264, 2450240, 2572752, 3376737, 4096000, 4096000, 4123000, 4881384, 4936750, 6669000, 6669000, 6669000, 6669000, 10822875, 11577384, 16539120, 18154500, 21326760, 21326760, 28861000, 30507008, 40199250, 44330496, 51684750, 72925515, 76271625, 77376000, 81153009, 91171899, 111321000, 190373976], 
      [1, 14, 64, -27, -27, 91, 104, 182, 0, 0, 181, 0, 0, -189, -189, 650, 0, 0, 832, 624, 819, 64, 64, 118, 729, 637, -351, -351, -378, -378, 924, 351, 0, -273, 0, 0, 1078, 896, -78, 168, 0, 0, 729, 1560, -729, 0, -1728, 0]]

mutualNub :: [Int] -> [Int] -> ([Int],[Int])
mutualNub (x:[]) (y:[]) = ([x], [y])
mutualNub (x:xs) (y:ys) = if (x == head xs)
    then rs
    else (x:(fst rs), y:(snd rs))
    	where rs = mutualNub xs ys
    	
rD :: [[Int]]
rD = [fst $ mutualNub (cl!!0) (cl!!1), snd $ mutualNub (cl!!0) (cl!!1)]