-- | See Hoffman, Gelman (2011) The No U-Turn Sampler: Adaptively Setting Path
--   Lengths in Hamiltonian Monte Carlo.

{-# OPTIONS_GHC -Wall -fno-warn-type-defaults #-}

import Control.Monad
import Control.Monad.Loops
import Control.Monad.Primitive
import Data.Hashable
import Data.HashSet (HashSet)
import qualified Data.HashSet as HashSet
import System.Random.MWC
import System.Random.MWC.Distributions

hmc :: (Enum a, Eq a, Ord a, Num a, PrimMonad m )
    => ([Double] -> Double)   -- ^ Log target function
    -> ([Double] -> [Double]) -- ^ Gradient of log target
    -> [Double]               -- ^ Parameters
    -> a                      -- ^ Epochs to run the chain
    -> a                      -- ^ Number of discretizing steps
    -> Double                 -- ^ Step size
    -> Gen (PrimState m)      -- ^ PRNG
    -> m [[Double]]           -- ^ Chain
hmc lTarget glTarget t n ndisc e g = unfoldrM kernel (n, (t, []))
  where 
    kernel (m, (p, _)) = do
      (p1, r1) <- hmcKernel lTarget glTarget p ndisc e g
      return $ if   m <= 0
               then Nothing
               else Just (p1, (pred m, (p1, r1)))
                   
hmcKernel :: (Enum a, Eq a, Ord a, Num a, PrimMonad m)
          => ([Double] -> Double)   -- ^ Log target function
          -> ([Double] -> [Double]) -- ^ Gradient of log target
          -> [Double]               -- ^ Parameters
          -> a                      -- ^ Number of discretizing steps
          -> Double                 -- ^ Step size
          -> Gen (PrimState m)      -- ^ PRNG
          -> m ([Double], [Double]) -- ^ m (End params, end momenta)
hmcKernel lTarget glTarget t0 ndisc e g = do
  r0 <- replicateM (length t0) (normal 0 1 g)
  z  <- uniform g
  let (t1, r1) = leapfrog glTarget t0 r0 ndisc e 
      a        = min 1 $ hmcAcceptanceRatio lTarget t0 t1 r0 r1
      final | a > z     = (t1, map negate r1)
            | otherwise = (t0, r0)
  return final

-- note that this is not the greatest buildTree we could use
buildTree
  :: (Enum a, Eq a, Floating c, Integral t, Num a, RealFrac c, Hashable c)
  => ([c] -> c)   -- ^ Log target
  -> ([c] -> [c]) -- ^ Gradient
  -> [c]          -- ^ Position
  -> [c]          -- ^ Momentum
  -> c            -- ^ Slice variable
  -> c            -- ^ Direction (-1, +1)
  -> a            -- ^ Depth
  -> c            -- ^ Step size
  -> ([c], [c], [c], [c], HashSet ([c], [c]), t)
buildTree lTarget glTarget = go 
  where 
    go t r u v 0 e = 
      let (t1, r1) = leapfrog glTarget t r 1 (v * e)
          c | u <= auxilliaryTarget lTarget t1 r1 = HashSet.singleton (t1, r1)
            | otherwise                           = HashSet.empty
          s | u < exp 1000 * auxilliaryTarget lTarget t1 r1 = 1
            | otherwise                                     = 0
      in  (t1, r1, t1, r1, c, s)

    go t r u v j e = 
      let (tn, rn, tp, rp, c0, s0)     = go t r u v (pred j) e
          (tnn, rnn, tpp, rpp, c1, s1) = if   roundTo 6 v == -1
                                         then go tn rn u v (pred j) e
                                         else go tp rp u v (pred j) e

          s2 = s0 * s1 * indicator ((tpp .- tnn) `innerProduct` rnn >= 0)
                       * indicator ((tpp .- tnn) `innerProduct` rpp >= 0)

          c2 = c0 `HashSet.union` c1

      in  (tnn, rnn, tpp, rpp, c2, s2) 

leapfrog :: (Enum a, Eq a, Ord a, Fractional c, Num a)
         => ([c] -> [c]) -- ^ Gradient of log target function
         -> [c]          -- ^ List of parameters to target
         -> [c]          -- ^ Momentum variables
         -> a            -- ^ Number of discretizing steps
         -> c            -- ^ Step size
         -> ([c], [c])   -- ^ (End parameters, end momenta)
leapfrog glTarget t0 r0 ndisc e | ndisc < 0 = (t0, r0)
                                          | otherwise = go t0 r0 ndisc
  where go t r 0 = (t, r)
        go t r n = let rm = zipWith (+) r  (map (* (0.5 * e)) (glTarget t))
                       tt = zipWith (+) t  (map (* e) rm)
                       rt = zipWith (+) rm (map (* (0.5 * e)) (glTarget t))
                   in  go tt rt (pred n)

-- | Acceptance ratio for a proposed move.  t0/r0 denote the present state of 
--   the parameters and auxilliary variables, and t1/r1 denote the proposed 
--   state.
hmcAcceptanceRatio :: Floating a => (t -> a) -> t -> t -> [a] -> [a] -> a
hmcAcceptanceRatio lTarget t0 t1 r0 r1 = auxilliaryTarget lTarget t1 r1
                                       / auxilliaryTarget lTarget t0 r0

-- | Augment a log target with some auxilliary variables.
auxilliaryTarget :: Floating a => (t -> a) -> t -> [a] -> a
auxilliaryTarget lTarget t r = exp (lTarget t - 0.5 * innerProduct r r)

findReasonableEpsilon :: PrimMonad m 
                      => ([Double] -> Double) 
                      -> ([Double] -> [Double]) 
                      -> [Double] 
                      -> Gen (PrimState m) 
                      -> m Double
findReasonableEpsilon lTarget glTarget t0 g = do
  r0 <- replicateM (length t0) (normal 0 1 g)
  let (t1, r1) = leapfrog glTarget t0 r0 1 1.0

      a = 2 * indicator (hmcAcceptanceRatio lTarget t0 t1 r0 r1 > 0.5) - 1

      go e t r | (hmcAcceptanceRatio lTarget t0 t r0 r) ^ a > 2 ^ (-a) = 
                   let en       = 2 ^ a * e
                       (tn, rn) = leapfrog glTarget t r 1 e
                   in  go en tn rn 
               | otherwise = e

  return $ go 1.0 t1 r1




-- Utilities ------------------------------------------------------------------

innerProduct :: Num a => [a] -> [a] -> a
innerProduct xs ys = sum $ zipWith (*) xs ys

(.-) :: Num a => [a] -> [a] -> [a]
xs .- ys = zipWith (-) xs ys

indicator :: Integral a => Bool -> a
indicator True  = 1
indicator False = 0

-- | Round to a specified number of digits.
roundTo :: RealFrac a => Int -> a -> a
roundTo n f = fromIntegral (round $ f * (10 ^ n) :: Int) / (10.0 ^^ n)

