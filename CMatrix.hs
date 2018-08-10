{-# LANGUAGE ForeignFunctionInterface #-}

module CMatrix where
{-- Модуль с межъязыковым интерфейсом  --}

import FFI
import Matrix
import Foreign.C.Types
import Foreign.Ptr
import Debug.Trace

detGauss_hs :: Ptr (Ptr CDouble) -> CInt -> CDouble
detGauss_hs mtx0 n0 = coerce $ detGauss $ mtx where
                     mtx = readMatrix n0 n0 mtx0

solveGauss_hs :: Ptr (Ptr CDouble) -> Ptr CDouble -> CInt -> Ptr CDouble
solveGauss_hs mtx0 v0 n0 = writeVector $ (\x -> trace ("Residual: "  ++ show (residual mtx x v)) x) $ solveGauss mtx v where
    v   = readVector n0 v0
    mtx = readMatrix n0 n0 mtx0

solveGaussLE_hs :: Ptr (Ptr CDouble) -> Ptr CDouble -> CInt -> Ptr CDouble
solveGaussLE_hs mtx0 v0 n0 = writeVector $ (\x -> trace ("Residual: "  ++ show (residual mtx x v)) x) $  solveGaussLE mtx v where
    v   = readVector n0 v0
    mtx = readMatrix n0 n0 mtx0


sccOvRl_hs :: Ptr (Ptr CDouble) -> Ptr CDouble -> CInt -> Double -> Double -> Ptr CDouble
sccOvRl_hs mtx0 v0 n0 w eps = writeVector $ (\x -> trace (if fst x == (-1) then "Diverges." else if isNaN $ head $ snd x then "Computational failure." else "Iterations: " ++ show (fst x)) $ snd x) $ sccOvRl (coerce w) (coerce eps) (map (const 0) v) mtx v  where
    v   = readVector n0 v0
    mtx = readMatrix n0 n0 mtx0


inv_hs :: Ptr (Ptr CDouble) -> CInt -> Ptr (Ptr CDouble)
inv_hs mtx0 n0 = writeMatrix $ inv mtx where
    mtx = readMatrix n0 n0 mtx0


condNumber_hs :: Ptr (Ptr CDouble) -> CInt -> CDouble
condNumber_hs mtx0 n0 = coerce $ condNumber mtx where
    mtx = readMatrix n0 n0 mtx0


foreign export ccall detGauss_hs :: Ptr (Ptr CDouble) -> CInt -> CDouble
foreign export ccall solveGauss_hs :: Ptr (Ptr CDouble) -> Ptr CDouble -> CInt -> Ptr CDouble
foreign export ccall solveGaussLE_hs :: Ptr (Ptr CDouble) -> Ptr CDouble -> CInt -> Ptr CDouble
foreign export ccall sccOvRl_hs :: Ptr (Ptr CDouble) -> Ptr CDouble -> CInt -> Double -> Double -> Ptr CDouble
foreign export ccall inv_hs :: Ptr (Ptr CDouble) -> CInt -> Ptr (Ptr CDouble)
foreign export ccall condNumber_hs :: Ptr (Ptr CDouble) -> CInt -> CDouble


