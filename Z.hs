{-# LANGUAGE BangPatterns #-}
module Main where
import Data.Complex
import Math.NumberTheory.GCD
import Math.NumberTheory.Moduli
import Data.Maybe
import Text.Printf
import Data.List

type I = Int

ep :: I -> Complex Double
ep d = if mod d 4 == 1 
	then 1
	else 0 :+ 1
	
dodd :: I -> Complex Double
dodd c = if mod c 2 == 0 then 0 else 1

klooster :: I -> I -> I -> I -> Complex Double
klooster lamb m n c = go 0 1 where
	c' = fromIntegral c
	go !acc !i
		| coprime i c  =
			let x = fromIntegral (m * (fromInteger $ fromJust (invertMod (fromIntegral i) (fromIntegral c))) + n * i) / c'
			in
			go (acc + (fromIntegral (jacobi c i) * (ep i)^(2*lamb + 1) * exp (0 :+ (2 * pi * x)))) (i + 1)
		| i == c       = acc
		| otherwise    = go acc (i + 1)

bess :: Double -> Complex Double
bess x = ((sqrt (2/(pi*x)))*(sinh x)):+0

b0h :: I -> I -> I -> I -> Complex Double
b0h level m 0 prec = go 0 (4*level) where
	go !acc !c
		| c > prec*level*4 = acc
		| otherwise         = go (acc + ((1 + dodd (div c 4))*(klooster 0 (-m) 0 c))*(((fromIntegral c)**(-1.5)):+0)) (c + 4*level)
b0h level m n prec = go 0 (4*level) where
	go !acc !c
		| c > prec*level*4 = acc
		| otherwise         = go (acc + ((1 + dodd (div c 4))*(klooster 0 (-m) n c)*((fromIntegral c)**(-1))*(bess (4*pi*(sqrt $ fromIntegral (m*n))*((fromIntegral c)**(-1)))))) (c + 4*level)
		
b0 :: I -> I -> I -> I -> Complex Double
b0 level m 0 prec = ((pi**(1.5))*2*(sqrt $ fromIntegral m)*(1:+(-1))*((0.5*sqrt(pi))**(-1)))*(b0h level m 0 prec)
b0 level m n prec = (pi*(sqrt 2)*((fromIntegral n)**(-0.25))*((fromIntegral m)**(0.25))*(1:+(-1)))*(b0h level m n prec)

levelPrint :: I -> I -> [(I, Complex Double)]
levelPrint prec level = zip coeffs $ map (\n -> b0 level 3 n prec) coeffs
	where coeffs = [0,1,4,5,8,9,12,13,16,17,20]
	
printComplex :: Complex Double -> String
printComplex c = if imagPart c < 0.0000001 
	then printf "%.4f" $ realPart c
	else printf "%.4f" (realPart c) ++ "+" ++ show (imagPart c) ++ "I"

representC :: (I, Complex Double) -> String
representC (i, a) = if i > 0 
	then "[" ++ printComplex a ++ "]" ++ "q^{" ++ show i ++ "}" 
	else printComplex a

showQ :: I -> [(I, Complex Double)] -> String
showQ level expansion = "Z_{" ++ show level ++ "} = " ++ (intercalate " + " $ map representC expansion)

doubleMap :: (a -> b -> c) -> [a] -> [b] -> [c]
doubleMap f [] [] = []
doubleMap f _ [] = []
doubleMap f [] _ = []
doubleMap f (x:xs) (y:ys) = (f x y):(doubleMap f xs ys)

main = putStr $ unlines $ doubleMap showQ [11,16,22,26,32,34] $ map (levelPrint 2800) [11,16,22,26,32,34]