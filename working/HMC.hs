-- | See Hoffman, Gelman (2011) The No U-Turn Sampler: Adaptively Setting Path
--   Lengths in Hamiltonian Monte Carlo.

{-# OPTIONS_GHC -fno-warn-type-defaults #-}
{-# LANGUAGE ScopedTypeVariables #-}

import Control.Monad
import Control.Monad.Loops
import Control.Monad.Primitive
import System.Random.MWC
import System.Random.MWC.Distributions

-- TODO what am i
dMax :: Num t => t
dMax = 1000

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

      a = 2 * indicate (hmcAcceptanceRatio lTarget t0 t1 r0 r1 > 0.5) - 1

      go e t r | (hmcAcceptanceRatio lTarget t0 t r0 r) ^ a > 2 ^ (-a) = 
                   let en       = 2 ^ a * e
                       (tn, rn) = leapfrog glTarget t r 1 e
                   in  go en tn rn 
               | otherwise = e

  return $ go 1.0 t1 r1

-- problem
buildTree :: (Enum a,  Eq a, Floating t, Fractional c,  Integral c, Integral d
             , Num a, Num e, RealFrac d,   RealFrac t, PrimMonad m , Variate c)
  => ([t] -> t)
  -> ([t] -> [t])
  -> Gen (PrimState m)
  -> [t]
  -> [t]
  -> t
  -> t
  -> a
  -> t
  -> t1
  -> [t]
  -> m ([t], [t], [t], [t], [t], c, d, t, e)
buildTree lTarget glTarget g = go 
  where
    go t r u v 0 e _ r0 = return $
      let (t1, r1) = leapfrog glTarget t r 1 (v * e)
          n        = indicate (u <= auxilliaryTarget lTarget t1 r1)
          s        = indicate (u <  exp dMax * auxilliaryTarget lTarget t1 r1)
          m        = min 1 (hmcAcceptanceRatio lTarget t1 r1 r0 r0)
      in  (t1, r1, t1, r1, t1, n, s, m, 1)

    go t r u v j e t0 r0 = do
      z <- uniform g
      (tn, rn, tp, rp, t1, n1, s1, a1, na1) <- go t r u v (pred j) e t0 r0

      if   s1 == 1
      then do
        (tnn, rnn, tpp, rpp, t2, n2, s2, a2, na2) <-
          if   v == -1
          then go tn rn u v (pred j) e t0 r0
          else go tp rp u v (pred j) e t0 r0
  
        let p  = n2 / (n1 + n2)
  
            t3 | p > z     = t2
               | otherwise = t1
  
            a3  = a1  + a2
            na3 = na1 + na2
  
            s3  = s2 * indicate ((tpp .- tnn) `innerProduct` rnn >= 0)
                     * indicate ((tpp .- tnn) `innerProduct` rpp >= 0) 
  
            n3  = n1 + n2
        return (tnn, rnn, tpp, rpp, t3, n3, s3, a3, na3)
      else return (tn, rn, tp, rp, t1, n1, s1, a1, na1)

innerNutsKernel :: (PrimMonad m, Variate b) 
                => ([a] -> Double) 
                -> t 
                -> [a]
                -> c
                -> Gen (PrimState m)
                -> m b
innerNutsKernel lTarget glTarget t e g = do
  r0 <- replicateM (length t) (normal 0 1 g)
  u  <- uniformR (0, auxilliaryTarget lTarget t r0) g

  let go (tn, tp, rn, rp, j, tm, n, s) a b gen = do
        vj <- symmetricCategorical [-1, 1] gen
        z  <- uniform gen 

        return z
--  let go (tn, tp, rn, rp, j, tm, n, s) aOrig naOrig g
--        | s == 1 = do
--            vj <- symmetricCategorical [-1, 1] g
--            z  <- uniform g
--
--            (tnn, rnn, tpp, rpp, t1, n1, s, a, na) <- 
--              buildTree lTarget glTarget g tn rn u vj j e t r0 -- FIXME
--
--            return $ (t1, a, na)

  go (t, t, r0, r0, 0, t, 1, 1) 0 0 g

  -- let go (tn, tp, rn, rp, j, tm, n, s) aOrig naOrig g
  --       | s == 1 = do
  --           vj <- symmetricCategorical [-1, 1] g
  --           z  <- uniform g

  --           (tnn, rnn, tpp, rpp, t1, n1, s1, a, na) <-
  --             if   vj == -1
  --             then buildTree lTarget glTarget g tn rn u vj j e t r0
  --             else buildTree lTarget glTarget g tp rp u vj j e t r0

  --           let t2 | s1 == 1 && (min 1 (fromIntegral n1 / fromIntegral n :: Double) > z) = tnn
  --                  | otherwise                       = t

  --               n2 = n + n1
  --               s2 = s1 * indicate ((tpp .- tnn) `innerProduct` rnn >= 0)
  --                       * indicate ((tpp .- tnn) `innerProduct` rpp >= 0) 
  --               j1 = succ j

  --           go (tnn, rnn, tpp, rpp, j1, t2, n2, s2) a na g

  --       | otherwise = return (tm, aOrig, naOrig)

  -- return $ go (t, t, r0, r0, 0, t, 1, 1) 0 0 g





-- Utilities ------------------------------------------------------------------

innerProduct :: Num a => [a] -> [a] -> a
innerProduct xs ys = sum $ zipWith (*) xs ys

(.-) :: Num a => [a] -> [a] -> [a]
xs .- ys = zipWith (-) xs ys

indicate :: Integral a => Bool -> a
indicate True  = 1
indicate False = 0

-- | Round to a specified number of digits.
roundTo :: RealFrac a => Int -> a -> a
roundTo n f = fromIntegral (round $ f * (10 ^ n) :: Int) / (10.0 ^^ n)

symmetricCategorical :: PrimMonad m => [a] -> Gen (PrimState m) -> m a
symmetricCategorical [] _ = error "symmetricCategorical: no candidates"
symmetricCategorical zs g = do
  z <- uniform g
  return $ zs !! truncate (z * fromIntegral (length zs) :: Double)

