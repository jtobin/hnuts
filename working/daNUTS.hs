-- | See Hoffman, Gelman (2011) The No U-Turn Sampler: Adaptively Setting Path
--   Lengths in Hamiltonian Monte Carlo.

import Control.Monad
import Control.Monad.Primitive
import System.Random.MWC
import System.Random.MWC.Distributions
import Statistics.Distribution.Normal

type Parameters = [Double]
type Density    = Parameters -> Double
type Gradient   = Parameters -> Parameters
type Particle   = (Parameters, Parameters)

leapfrogIntegrator :: Int -> Gradient -> Particle -> Double -> Particle
leapfrogIntegrator n glTarget particle e = go particle n
  where go state ndisc 
          | ndisc <= 0 = state
          | otherwise  = go (leapfrog glTarget state e) (pred n)

leapfrog :: Gradient -> Particle -> Double -> Particle
leapfrog glTarget (t, r) e = (tf, rf)
  where rm = zipWith (+) r  ((e / 2) .* glTarget t)
        tf = zipWith (+) t  (e .* rm)
        rf = zipWith (+) rm ((e / 2) .* glTarget tf)

findReasonableEpsilon :: PrimMonad m 
                      => Density
                      -> Gradient
                      -> Parameters
                      -> Gen (PrimState m) 
                      -> m Double
findReasonableEpsilon lTarget glTarget t0 g = do
  r0 <- replicateM (length t0) (normal 0 1 g)
  let (t1, r1) = leapfrog glTarget (t0, r0) 1.0
      a        = 2 * indicate (acceptanceRatio lTarget t0 t1 r0 r1 > 0.5) - 1

      go e t r | (acceptanceRatio lTarget t0 t r0 r) ^ a > 2 ^ (-a) = 
                   let (tn, rn) = leapfrog glTarget (t, r) e
                   in  go (2 ^ a * e) tn rn 
               | otherwise = e

  return $ go 1.0 t1 r1

-- this is the dual averaging buildTree
buildTree 
  :: PrimMonad m 
  => Density
  -> Gradient
  -> Gen (PrimState m)
  -> Parameters
  -> Parameters
  -> Double
  -> Double
  -> Int
  -> Double
  -> Parameters
  -> Parameters
  -> m ([Double], [Double], [Double], [Double], [Double], Int, Int, Double, Int)
buildTree lTarget glTarget g = go
  where 
    go t r u v 0 e _ r0 = return $
      let (t1, r1) = leapfrog glTarget (t, r) (v * e)
          n        = indicate (u <= auxilliaryTarget lTarget t1 r1)
          s        = indicate (u <  exp 1000 * auxilliaryTarget lTarget t1 r1)
          m        = min 1 (acceptanceRatio lTarget t1 r1 r0 r0)
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

        let p   = fromIntegral n2 / fromIntegral (n1 + n2)
            n3  = n1 + n2
            t3  | p > (z :: Double) = t2
                | otherwise         = t1
            a3  = a1 + a2
            na3 = na1 + na2
            s3  = s2 * indicate ((tpp .- tnn) `innerProduct` rnn >= 0)
                     * indicate ((tpp .- tnn) `innerProduct` rpp >= 0)

        return (tnn, rnn, tpp, rpp, t3, n3, s3, a3, na3)
      else return (tn, rn, tp, rp, t1, n1, s1, a1, na1)

relaxingNuts = undefined

-- better idea: wrap this dual averaging scheme around the actual nuts
--              kernel itself.  in fact you'd like to just be able to loosely
--              add dual-averaging to any procedure.
--
-- adaptingNutsKenel lTarget glTarget t m g = do
--   e0 <- findReasonableEpsilon lTarget glTarget t g
-- 
--   let mu      = log (10 * e)
--       epsBar0 = 0
--       h0Bar   = 0
--       gamma   = 0.05
--       delta   = 0.45 -- target mean acceptance probability
--       tau0    = 10
--       kappa   = 0.75
-- 
--       go hBar eNext logEpsBar tToReturn n
--         | n <= 0    = return (tToReturn, logEpsBar, 
-- 
--         | otherwise = do
--             (t0, a, na) <- innerNutsKernel lTarget glTarget t e g
--             let hBarNext = (1 - 1 / (m - n + tau0)) * hBar
--                          + (1 / (m - n + tau0)) * (delta - a) 
-- 
--                 logEpsNext = mu - ((sqrt (m - n)) / gamma) * hmBar
-- 
--                 logEpsBarNext = (m - n) ^ (-kappa) * logEpsNext
--                               + (1 - (m - n) ^ (-kappa)) * logEpsBar
-- 
--             go hBarNext logEpsBarNext t0 (pred n)




innerNutsKernel 
    :: PrimMonad m 
    => Density
    -> Gradient
    -> Parameters
    -> Double
    -> Gen (PrimState m)
    -> m (Parameters, Double, Int)
innerNutsKernel lTarget glTarget t e g = do
  r0 <- replicateM (length t) (normal 0 1 g)
  u  <- uniformR (0, auxilliaryTarget lTarget t r0) g

  let go (tn, tp, rn, rp, j, tm, n, s) aOrig naOrig g
        | s == 1 = do
            vj <- symmetricCategorical [-1, 1] g
            z  <- uniform g

            (tnn, rnn, tpp, rpp, t1, n1, s1, a, na) <-
              if   vj == -1
              then buildTree lTarget glTarget g tn rn u vj j e t r0
              else buildTree lTarget glTarget g tp rp u vj j e t r0

            let t2 | s1 == 1 && min 1 (fi n1 / fi n :: Double) > z = tnn
                   | otherwise = t

                n2 = n + n1
                s2 = s1 * indicate ((tpp .- tnn) `innerProduct` rnn >= 0)
                        * indicate ((tpp .- tnn) `innerProduct` rpp >= 0) 
                j1 = succ j

            go (tnn, rnn, tpp, rpp, j1, t2, n2, s2) a na g

        | otherwise = return (tm, aOrig, naOrig)

  go (t, t, r0, r0, 0, t, 1, 1) 0 0 g

auxilliaryTarget :: Floating a => (t -> a) -> t -> [a] -> a
auxilliaryTarget lTarget t r = exp (lTarget t - 0.5 * innerProduct r r)

acceptanceRatio :: Floating a => (t -> a) -> t -> t -> [a] -> [a] -> a
acceptanceRatio lTarget t0 t1 r0 r1 = auxilliaryTarget lTarget t1 r1
                                    / auxilliaryTarget lTarget t0 r0

innerProduct :: Num a => [a] -> [a] -> a
innerProduct xs ys = sum $ zipWith (*) xs ys

(.*) :: Num b => b -> [b] -> [b]
z .* xs = map (* z) xs

(.-) :: Num a => [a] -> [a] -> [a]
xs .- ys = zipWith (-) xs ys

indicate :: Integral a => Bool -> a
indicate True  = 1
indicate False = 0

symmetricCategorical :: PrimMonad m => [a] -> Gen (PrimState m) -> m a
symmetricCategorical [] _ = error "symmetricCategorical: no candidates"
symmetricCategorical zs g = do
  z <- uniform g
  return $ zs !! truncate (z * fromIntegral (length zs) :: Double)

fi :: (Integral a, Num b) => a -> b
fi = fromIntegral

-- Testing

f :: Density
f _ = log $ 1 / 10

g :: Gradient
g xs = replicate (length xs) 0

