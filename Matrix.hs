module Matrix where

import Data.List hiding (transpose)
import Debug.Trace

{-- Модуль с вычислительными методами --}

--Определитель Лапласа
detLaplace :: [[Integer]] -> Integer
detLaplace [[a]] = a
detLaplace a     = sum $ map (\x -> (head $ drop (x - 1) $ last a) * (detLaplace [take (x - 1) y ++ drop x y| y <- init a]) * ((-1)^(x + n))) [1..n] where n = length a

--Все возможные перестановки
permute    :: [a] -> [[a]]
permute [a]      = [[a]]
permute (x:xs)   = concat $ map (\p -> [take p y ++ [x] ++ drop p y|y <- permute xs]) [0..n] where n = length xs

--Четность перестановки
parity     :: [Int] -> Int
parity [a]       = 1
parity (x:xs)    = parity xs * product (map (\p -> signum $ p - x) xs)

--Определитель пересчетом перестановок
detPerm    :: [[Int]] -> Int
detPerm   a      = sum $ map (\p -> parity p * (product $ zipWith (\c d -> head $ drop (c - 1) d) p a)) $ permute [1..n] where n = length a

--Транспонирование
transpose  :: [[a]] -> [[a]]
trasnpose [[]]   = [[]]
transpose []     = []
transpose ([]:_) = []
transpose a      = (map head a):(transpose $ map tail a)

--Квантор всеобщности для списков
forAll     :: (a -> Bool) -> [a] -> Bool
forAll _ []      = True
forAll f (x:xs)  = f x && forAll f xs

--Квантор существования для списков
exists     :: (a -> Bool) -> [a] -> Bool 
exists _ []      = False
exists f (x:xs)  = f x || exists f xs

--Определитель Гаусса
detGauss     :: (RealFloat a) => [[a]] -> a
detGauss [[a]]   = a
detGauss a@((x:_):_) | x == 0    = if exists (0/=) $ head $ transpose a then (*) (-1) $ detGauss $ head [(head $ drop y a):(take y a ++ drop (y + 1) a)|y <- [1..], head (head $ drop y a) /= 0] else 0 
                     | otherwise = (*) x $ detGauss $ map (\k -> tail $ zipWith (\q p -> p - ((q * head k)/x)) (head a) k) $ tail a
--Перемножение матриц
mtxMult      :: (RealFloat a) => [[a]] -> [[a]] -> [[a]]
mtxMult a b      = map (\x -> map (\y -> sum $ zipWith (*) x y) $ transpose b) a

--Домножение матрицы на константу
mtxMultC     :: (RealFloat a) => a -> [[a]] -> [[a]]
mtxMultC a b     = map (map (a*)) b 

--Сумма матриц
mtxAdd       :: (RealFloat a) => [[a]] -> [[a]] -> [[a]]
mtxAdd a b       = zipWith (zipWith (+)) a b

--Разность матриц
mtxSub       :: (RealFloat a) => [[a]] -> [[a]] -> [[a]]
mtxSub a b       = zipWith (zipWith (-)) a b

--Возведение матрицы в степень
mtxPow       :: (RealFloat a) => [[a]] -> Int -> [[a]]
mtxPow a 0       = mtxOne $ length a
mtxPow a n       = mtxMult a $ mtxPow a $ n - 1

--Евклидова норма вектора
vcrNorm2     :: (RealFloat a) => [a] -> a
vcrNorm2 a       = sqrt $ sum $ map (**2) a

--Евклидова норма матрицы
mtxNorm2     :: (RealFloat a) => [[a]] -> a
mtxNorm2 a       = vcrNorm2 $ concat a

--Нулевые матрицы
mtxZero      :: (Num a) => Int -> [[a]]
mtxZero n        = [[0|_<-k]|_<-k] where k = [1..n]

--Единичные матрицы
mtxOne       :: (RealFloat a) => Int -> [[a]]
mtxOne 1         = [[1]]
mtxOne n         = ([1] ++ [0|_<-[1..(n - 1)]]):(map ([0] ++) $ mtxOne $ n - 1)

--Метод вычисления определителя по-умолчанию
det              = detGauss 

--Преобразование вектора в квадратную матрицу
mtx          :: [a] -> [[a]]
mtx v            = [take n $ drop (n * x) v|x <- [0..(n - 1)]] where n = floor $ sqrt $ fromIntegral $ length v

--Отображение матрицы
mtxMap       :: (a -> b) -> [[a]] -> [[b]]
mtxMap f a       = map (map f) a

--Максимум по вектору
vcrMax       :: Ord a => [a] -> a
vcrMax a         = head [x| x <- a, forAll (x>=) a]

--Минимум по вектору
vcrMin       :: Ord a => [a] -> a
vcrMin a         = head [x| x <- a, forAll (x<=) a]

--Максимум по матрице
mtxMax       :: Ord a => [[a]] -> a
mtxMax a         = vcrMax $ concat a

--Минимум по матрице
mtxMin       :: Ord a => [[a]] -> a
mtxMin a         = vcrMin $ concat a

--Сумма набора матриц
mtxSum       :: (RealFloat a) =>  [[[a]]] -> [[a]]
mtxSum [a]       = a
mtxSum (x:xs)    = mtxAdd x $ mtxSum xs

--Произведение набора матриц
mtxProd      :: (RealFloat a) =>  [[[a]]] -> [[a]]  
mtxProd [a]      = a
mtxProd (x:xs)   = mtxMult x $ mtxProd xs

--Многочлен из матриц
mtxPoly      :: (RealFloat a) =>  [a] -> [[a]] -> [[a]]
mtxPoly [c] a    = mtxMultC c $ mtxOne $ length a
mtxPoly c a      = mtxAdd (mtxPoly [head c] a) $ mtxMult a $ mtxPoly (tail c) a

--Союзная матрица
adj          :: [[Double]] -> [[Double]]
adj a            = [[let t q p = take (q - 1) p ++ drop q p in (*) ((^) (-1) $ x + y) $ det $ t y $ map (t x) a| y <- [1..n]] | x <- [1..n]] where n = length a

--Обращение через союзные матрицы
invAdj       :: [[Double]] -> [[Double]]
invAdj a           = mtxMap (/ det a) $ adj a

--Обращение Шульца
invSchultz   :: (RealFloat a) => Int -> a -> [[a]] -> [[a]]
invSchultz n err a = schultz n err a (mtxMultC ((/) 1 $ mtxNorm2 $ mtxMult a t) t) where 
    t = transpose a
    schultz m err a u = let psi = mtxSub (mtxOne (length a)) $ mtxMult a u in if mtxNorm2 psi < err then u else schultz m err a (mtxMult u $ mtxPoly [1|_<-[0..m]] psi) 
--Обращение Жордана-Гаусса
invGaussJordan :: (RealFloat a) => [[a]] -> [[a]]
invGaussJordan a0  = gauss $ jordan (a0, mtxOne $ length a0) where 
    takeLast n xs = foldl (const . tail) xs (drop n xs)
    dropLast n xs = foldl (const . init) xs [1..n]
    gauss ([[1]], b) = b
    gauss (a, b)     = gauss (map init $ init a, (zipWith (\x -> zipWith (\y z -> z - y*x) current) (init $ map last a) $ take (length a - 1) b) ++ (drop (length a - 1) b)) where current = last $ take (length a) b
    jordan  ([[a]], b) = ([[1]], init b ++ [map (/a) $ last b])
    jordan  (a@((x:_):_), b) | x == 0    = jordan $ head [let swap m = (take (length m - length a) m) ++ [last $ dropLast y m] ++ (drop (length m - length a) $ dropLast (y + 1) m ++ takeLast y m) in (swap a, swap b)|y <- [0..], head (last $ dropLast y a) /= 0]
                             | otherwise = ((head a'):(map (0:) $ fst next), snd next) where 
        (a', b')   = ((map (/x) $ head a):(tail a), dropLast n b ++ [map (/x) $ head $ takeLast n b] ++ (takeLast (n - 1) b))
        (a'', b'') = (map (\k -> tail $ zipWith (\q p -> p - q * head k) (head a') k) $ tail a', dropLast (n - 1) b' ++ (zipWith (\k l -> zipWith (\q p -> p - q * head l) (head $ takeLast n b') k) (takeLast (n - 1) b') $ tail a'))
        n          = length a
        next       = jordan (a'', b'')

--Метод вычисления обратной матрицы по-умолчанию
inv                = invGaussJordan

--Умножение матрицы на вектор
opApply      :: (RealFloat a) => [[a]] -> [a] -> [a]
opApply o v      = head $ transpose $ mtxMult o $ transpose [v] 

--Решение методом Крамера
cramer       :: [[Double]] -> [Double] -> [Double]
cramer a0 v      = [(/) (det (take (y - 1) a ++ [v] ++ drop y a)) $ det a|y <- [1..(length a)]] where a = transpose a0 


--Маска матрицы
mtxMask      :: (Num a) => [[Bool]] -> [[a]] -> [[a]]
mtxMask m a = zipWith (zipWith (\x y -> if x then y else 0)) m a

--Служебная функция для упрощения рекурсивных методов
--Возвращает подматрицу данной матрицы, исключая первые столбец и строку
subMtx    :: [[a]] -> [[a]]
subMtx a = map tail $ tail a


--Служебная функция подставляющая подматрицу
mtxMerge     :: [[a]] -> [[a]] -> [[a]]
mtxMerge a b = head a : (zipWith (\x y -> head x : y) (tail a) b)

--Разложение матрицы на диагональную и триугольные
mtxDecomp    :: (Num a) => [[a]] -> ([[a]], [[a]], [[a]])
mtxDecomp [] = ([], [], [])
mtxDecomp a  = (d, l, u) where
    (sd, sl, su) = mtxDecomp $ subMtx a
    zm           = mtxZero $ length a
    d            = (((head (head a) : tail (head zm))) : tail zm) `mtxMerge` sd
    l            = (head zm : tail a) `mtxMerge` sl
    u            =  ((0 : tail (head a)) : tail zm) `mtxMerge` su

--Невязка
residual     :: [[Double]] -> [Double] -> [Double] -> Double
residual a x b = vcrNorm2 $ zipWith (-) (opApply a x) b 

--Максимальная невязка
maxres = 1.0e50

--Верхняя релаксация 
sccOvRl      :: Double -> Double -> [Double] -> [[Double]] -> [Double] -> (Int, [Double])
sccOvRl w eps x0 a b 
    | res > maxres = (-1, map (const (0/0)) b)
    | res > eps    = let (i, x) = sccOvRl w eps current a b in if i /= (-1) then (1 + i, x) else (i, x)
    | otherwise    = (0, current)    
    where (d, l, u) = mtxDecomp a
          pr1       = inv $ d `mtxAdd` mtxMultC w l
          rhs       = opApply pr1 $ zipWith (-) b $ opApply a x0
          delta     = map (w*) rhs
          current   = zipWith (+) x0 delta
          res       = residual a current b

--Решение методом Гаусса
solveGauss  :: [[Double]] -> [Double] -> [Double]
solveGauss a f = reverse $ svTrg $ mkTrg a f where
    mkTrg [] _ = ([], [])
    mkTrg a  f = (mtxMerge newa nga, head newf : ngf) where
        af           = zipWith (\m v -> m ++ [v]) a f
        newaf        = af0 : (map newafn $ tail af) where 
            af0       = head af
            newafn m  = zipWith (-) m $ map ((head m  / head af0)*) af0
        (newa, newf) = (map init newaf, map last newaf)
        (nga, ngf)   = mkTrg (subMtx newa) (tail newf)  
    svTrg ([], _) = []
    svTrg (a,  f) = x : svTrg ((map init $ init a), newf) where
        x    = last f / (last $ last a)
        newf = zipWith (\m v -> v - (last m * x)) (init a) (init f)

--Применение перестановки (или их композиция)
permApply :: [Int] -> [a] -> [a]
permApply p s = map (\x -> s !! (x - 1)) p

--Обращение перестановки
permInvert :: [Int] -> [Int]
permInvert x = map snd $ sort $ zip x [1..(length x)]

--Транспозиция
permSwap :: Int -> Int -> Int -> [Int]
permSwap l n0 m0
    | n0 == m0  = [1..l] 
    | otherwise = (take (n - 1) s) ++ [s !! (m - 1)] ++ (take (m - n - 1) $ drop n s) ++ [s !! (n - 1)] ++ (drop m s) where
        (n, m) = (min n0 m0, max n0 m0)
        s      = [1..l]


--Решение методом Гаусса с выбором главного элемента
solveGaussLE  :: [[Double]] -> [Double] -> [Double]
solveGaussLE a f = permApply (permInvert prm2) $ reverse $ svTrg (na, nf) where
    maxPerms a = head [(permSwap l 1 x, permSwap l 1 y) | x <- [1..l], y <- [1..l], abs (a !! (x - 1) !! (y - 1)) == mtxMax (mtxMap abs a)] where
        l = length a
    mkTrg [] _   = ([], [], [], [])
    mkTrg ao fo  =  (rearTopRow ngp2 $ mtxMerge newa nga, head newf : ngf, newp1, newp2) where
        rearTopRow p m = ((head $ head m) : (permApply p $ tail $ head m)) : tail m
        (p1, p2)       = maxPerms ao
        f1             = permApply p1 fo
        a1             = transpose $ permApply p2 $ transpose $ permApply p1 ao 
        af             = zipWith (\m v -> m ++ [v]) a1 f1
        newaf          = af0 : (map newafn $ tail af) where 
            af0       = head af
            newafn m  = zipWith (-) m $ map ((head m  / head af0)*) af0
        (newa, newf) = (map init newaf, map last newaf)
        (nga, ngf, ngp1, ngp2) = mkTrg (subMtx newa) (tail newf)
        (newp1, newp2)         = let comp x y = head x : (permApply y $ tail x) in (comp p1 ngp1, comp p2 ngp2) 
    svTrg ([], _) = []
    svTrg (a,  f) = x : svTrg ((map init $ init a), newf) where
        x    = last f / (last $ last a)
        newf = zipWith (\m v -> v - (last m * x)) (init a) (init f)
    (na, nf, prm1, prm2) = mkTrg a f

--1-норма вектора
vcrNorm1 :: [Double] -> Double
vcrNorm1 v = sum $ map abs v

--Норма, используемая при вычислении ЧО
mtxNormInf :: [[Double]] -> Double
mtxNormInf m = maximum $ map vcrNorm1 $ transpose m


--Число обусловленности
condNumber :: [[Double]] -> Double
condNumber a = (mtxNormInf a) * (mtxNormInf $ inv a)


--Красивый вывод матрицы
mtxDisp      :: (Show a) => [[a]] -> String
mtxDisp []           = ""
mtxDisp ([]:xs)      = "\n" ++ mtxDisp xs
mtxDisp ((x:xs):xxs) = (show x ++ " ") ++ mtxDisp (xs:xxs)  
