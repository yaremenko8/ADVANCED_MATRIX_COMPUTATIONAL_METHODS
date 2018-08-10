{-# LANGUAGE ForeignFunctionInterface #-}

module FFI where
{-- Модуль со средствами пострения межъязыковых интерфейсов --}

import Foreign.C.Types
import Foreign.Marshal.Array
import Foreign.Marshal.Unsafe
import Foreign.Ptr
import Foreign.Storable
import Data.Typeable
import Unsafe.Coerce

--Преобразование типов. Используется для преобразования
--типов одного языка в типы второго.
coerce :: a -> b
coerce = unsafeCoerce

--Преобразовать массив заданной длинны в список
readVector :: Storable a => CInt -> Ptr a -> [b]
readVector len0 arr = map coerce rawlist where 
                    len     = fromIntegral len0
                    rawlist = unsafeLocalState $ peekArray len arr  

--Преобразовать двумерный массив заданных размеров в двумерный список
readMatrix :: Storable a => CInt -> CInt -> Ptr (Ptr a) -> [[b]]
readMatrix n0 m0 arr = map (readVector n0) ptrlist where
                    m       = fromIntegral m0
                    ptrlist = unsafeLocalState $ peekArray m arr

--Преобразовать вектор в массив
writeVector :: (Storable a, Storable b) => [a] -> Ptr b
writeVector v = unsafeLocalState $ newArray (map coerce v)

--Преобразовать двумерный список в двумерный массив
writeMatrix :: (Storable a, Storable b) => [[a]] -> Ptr (Ptr b)
writeMatrix m = unsafeLocalState $ newArray (map writeVector m)

