-- | See Hoffman, Gelman (2011) The No U-Turn Sampler: Adaptively Setting Path
--   Lengths in Hamiltonian Monte Carlo.

module Numeric.MCMC.NUTS where

import Control.Monad
import Control.Monad.Loops
import Control.Monad.Primitive
import System.Random.MWC
import System.Random.MWC.Distributions hiding (gamma)
import Statistics.Distribution.Normal

type Parameters = [Double] 
type Density    = Parameters -> Double
type Gradient   = Parameters -> Parameters
type Particle   = (Parameters, Parameters)

newtype BuildTree = BuildTree { 
    getBuildTree :: ([Double], [Double], [Double], [Double], [Double], Int, Int)
  }

instance Show BuildTree where
  show (BuildTree (tm, rm, tp, rp, t', n, s)) = 
       "\n" ++ "tm: " ++ show tm 
    ++ "\n" ++ "tp: " ++ show tp
    ++ "\n" ++ "t': " ++ show t'
    ++ "\n" ++ "n : " ++ show n
    ++ "\n" ++ "s : " ++ show s

data DualAveragingParameters = DualAveragingParameters {
    mAdapt :: Int
  , delta  :: Double
  , mu     :: Double
  , gamma  :: Double
  , tau0   :: Double
  , kappa  :: Double
  } deriving Show

-- | The NUTS sampler.
nuts 
  :: PrimMonad m
  => Density
  -> Gradient
  -> Int
  -> Double
  -> Parameters
  -> Gen (PrimState m)
  -> m [Parameters]
nuts lTarget glTarget n e t g = go t 0 []
  where go position j acc
          | j >= n    = return acc
          | otherwise = do
              nextPosition <- nutsKernel lTarget glTarget e position g
              go nextPosition (succ j) (nextPosition : acc)

-- | The NUTS sampler with dual averaging.
nutsDualAveraging
  :: PrimMonad m
  => Density
  -> Gradient
  -> Int
  -> Int
  -> Parameters
  -> Gen (PrimState m)
  -> m [Parameters]
nutsDualAveraging lTarget glTarget n nAdapt t g = do
    e0 <- findReasonableEpsilon lTarget glTarget t g
    let daParams = basicDualAveragingParameters e0 nAdapt
    chain <- unfoldrM (kernel daParams) (1, e0, 1, 0, t)
    return $ drop nAdapt chain
  where
    kernel params (m, e, eAvg, h, t0) = do
      (eNext, eAvgNext, hNext, tNext) <- 
        nutsKernelDualAvg lTarget glTarget e eAvg h m params t0 g
      return $ if   m > n + nAdapt
               then Nothing
               else Just (t0, (succ m, eNext, eAvgNext, hNext, tNext))

-- | Default DA parameters, given a base step size and burn in period.
basicDualAveragingParameters :: Double -> Int -> DualAveragingParameters
basicDualAveragingParameters step burnInPeriod = DualAveragingParameters {
    mu     = log (10 * step)
  , delta  = 0.5
  , mAdapt = burnInPeriod
  , gamma  = 0.05
  , tau0   = 10
  , kappa  = 0.75
  }

-- | A single iteration of dual-averaging NUTS.
nutsKernelDualAvg 
  :: PrimMonad m 
  => Density
  -> Gradient
  -> Double
  -> Double
  -> Double
  -> Int
  -> DualAveragingParameters
  -> [Double]
  -> Gen (PrimState m)
  -> m (Double, Double, Double, Parameters)
nutsKernelDualAvg lTarget glTarget e eAvg h m daParams t g = do
  r0 <- replicateM (length t) (normal 0 1 g)
  z0 <- exponential 1 g
  let logu = log (auxilliaryTarget lTarget t r0) - z0

  let go (tn, tp, rn, rp, tm, j, n, s, a, na) g
        | s == 1 = do
            vj <- symmetricCategorical [-1, 1] g
            z  <- uniform g

            (tnn, rnn, tpp, rpp, t1, n1, s1, a1, na1) <-
              if   vj == -1
              then do
                (tnn', rnn', _, _, t1', n1', s1', a1', na1') <- 
                  buildTreeDualAvg lTarget glTarget g tn rn logu vj j e t r0
                return (tnn', rnn', tp, rp, t1', n1', s1', a1', na1')
              else do
                (_, _, tpp', rpp', t1', n1', s1', a1', na1') <- 
                  buildTreeDualAvg lTarget glTarget g tp rp logu vj j e t r0
                return (tn, rn, tpp', rpp', t1', n1', s1', a1', na1')

            let accept = s1 == 1 && (min 1 (fi n1 / fi n :: Double)) > z 

                n2 = n + n1
                s2 = s1 * stopCriterion tnn tpp rnn rpp
                j1 = succ j
                t2 | accept    = t1
                   | otherwise = tm

            go (tnn, tpp, rnn, rpp, t2, j1, n2, s2, a1, na1) g

        | otherwise = return (tm, a, na)

  (nextPosition, alpha, nalpha) <- go (t, t, r0, r0, t, 0, 1, 1, 0, 0) g
  
  let (hNext, eNext, eAvgNext) =
          if   m <= mAdapt daParams
          then (hm, exp logEm, exp logEbarM)
          else (h, eAvg, eAvg)
        where
          eta = 1 / (fromIntegral m + tau0 daParams)
          hm  = (1 - eta) * h 
              + eta * (delta daParams - alpha / fromIntegral nalpha)

          zeta = fromIntegral m ** (- (kappa daParams))

          logEm    = mu daParams - sqrt (fromIntegral m) / gamma daParams * hm
          logEbarM = (1 - zeta) * log eAvg + zeta * logEm

  return (eNext, eAvgNext, hNext, nextPosition)

-- | A single iteration of NUTS.
nutsKernel 
  :: PrimMonad m 
  => Density
  -> Gradient
  -> Double
  -> Parameters
  -> Gen (PrimState m)
  -> m Parameters
nutsKernel lTarget glTarget e t g = do
  r0   <- replicateM (length t) (normal 0 1 g)
  z0   <- exponential 1 g
  let logu = log (auxilliaryTarget lTarget t r0) - z0

  let go (tn, tp, rn, rp, tm, j, n, s) g
        | s == 1 = do
            vj <- symmetricCategorical [-1, 1] g
            z  <- uniform g

            (tnn, rnn, tpp, rpp, t1, n1, s1) <- 
              if   vj == -1
              then do
                (tnn', rnn', _, _, t1', n1', s1') <- 
                  buildTree lTarget glTarget g tn rn logu vj j e
                return (tnn', rnn', tp, rp, t1', n1', s1')
              else do
                (_, _, tpp', rpp', t1', n1', s1') <- 
                  buildTree lTarget glTarget g tp rp logu vj j e
                return (tn, rn, tpp', rpp', t1', n1', s1')

            let accept = s1 == 1 && (min 1 (fi n1 / fi n :: Double)) > z

                n2 = n + n1
                s2 = s1 * stopCriterion tnn tpp rnn rpp
                j1 = succ j
                t2 | accept    = t1
                   | otherwise = tm

            go (tnn, tpp, rnn, rpp, t2, j1, n2, s2) g

        | otherwise = return tm

  go (t, t, r0, r0, t, 0, 1, 1) g

-- | Build the 'tree' of candidate states.
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
  -> m ([Double], [Double], [Double], [Double], [Double], Int, Int)
buildTree lTarget glTarget g t r logu v 0 e = do
  let (t0, r0) = leapfrog glTarget (t, r) (v * e)
      joint    = log $ auxilliaryTarget lTarget t0 r0
      n        = indicate (logu < joint)
      s        = indicate (logu - 1000 < joint)
  return (t0, r0, t0, r0, t0, n, s)

buildTree lTarget glTarget g t r logu v j e = do
  z <- uniform g
  (tn, rn, tp, rp, t0, n0, s0) <- 
    buildTree lTarget glTarget g t r logu v (pred j) e

  if   s0 == 1
  then do
    (tnn, rnn, tpp, rpp, t1, n1, s1) <- 
      if   v == -1
      then do
        (tnn', rnn', _, _, t1', n1', s1') <- 
          buildTree lTarget glTarget g tn rn logu v (pred j) e
        return (tnn', rnn', tp, rp, t1', n1', s1')
      else do
        (_, _, tpp', rpp', t1', n1', s1') <- 
          buildTree lTarget glTarget g tp rp logu v (pred j) e
        return (tn, rn, tpp', rpp', t1', n1', s1')

    let accept = (fi n1 / max (fi (n0 + n1)) 1) > (z :: Double)
        n2     = n0 + n1
        s2     = s0 * s1 * stopCriterion tnn tpp rnn rpp
        t2     | accept    = t1
               | otherwise = t0 

    return (tnn, rnn, tpp, rpp, t2, n2, s2)
  else return (tn, rn, tp, rp, t0, n0, s0)

-- | Determine whether or not to stop doubling the tree of candidate states.
stopCriterion :: (Integral a, Num b, Ord b) => [b] -> [b] -> [b] -> [b] -> a
stopCriterion tn tp rn rp = 
      indicate (positionDifference `innerProduct` rn >= 0)
    * indicate (positionDifference `innerProduct` rp >= 0)
  where
    positionDifference = tp .- tn

-- | Build the tree of candidate states under dual averaging.
buildTreeDualAvg
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
buildTreeDualAvg lTarget glTarget g t r logu v 0 e t0 r0 = do
  let (t1, r1) = leapfrog glTarget (t, r) (v * e)
      joint    = log $ auxilliaryTarget lTarget t1 r1
      n        = indicate (logu < joint)
      s        = indicate (logu - 1000 <  joint)
      a        = min 1 (acceptanceRatio lTarget t0 t1 r0 r1)
  return (t1, r1, t1, r1, t1, n, s, a, 1)
      
buildTreeDualAvg lTarget glTarget g t r logu v j e t0 r0 = do
  z <- uniform g
  (tn, rn, tp, rp, t1, n1, s1, a1, na1) <- 
    buildTreeDualAvg lTarget glTarget g t r logu v (pred j) e t0 r0

  if   s1 == 1
  then do
    (tnn, rnn, tpp, rpp, t2, n2, s2, a2, na2) <-
      if   v == -1
      then do 
        (tnn', rnn', _, _, t1', n1', s1', a1', na1') <- 
          buildTreeDualAvg lTarget glTarget g tn rn logu v (pred j) e t0 r0
        return (tnn', rnn', tp, rp, t1', n1', s1', a1', na1')
      else do
        (_, _, tpp', rpp', t1', n1', s1', a1', na1') <-
          buildTreeDualAvg lTarget glTarget g tp rp logu v (pred j) e t0 r0
        return (tn, rn, tpp', rpp', t1', n1', s1', a1', na1')

    let p      = fromIntegral n2 / max (fromIntegral (n1 + n2)) 1
        accept = p > (z :: Double)
        n3     = n1 + n2
        a3     = a1 + a2
        na3    = na1 + na2
        s3     = s1 * s2 * stopCriterion tnn tpp rnn rpp

        t3  | accept    = t2
            | otherwise = t1

    return (tnn, rnn, tpp, rpp, t3, n3, s3, a3, na3)
  else return (tn, rn, tp, rp, t1, n1, s1, a1, na1)

-- | Heuristic for initializing step size.
findReasonableEpsilon 
  :: PrimMonad m 
  => Density
  -> Gradient
  -> Parameters
  -> Gen (PrimState m) 
  -> m Double
findReasonableEpsilon lTarget glTarget t0 g = do
  r0 <- replicateM (length t0) (normal 0 1 g)
  let (t1, r1) = leapfrog glTarget (t0, r0) 1.0
      a        = 2 * indicate (acceptanceRatio lTarget t0 t1 r0 r1 > 0.5) - 1

      go j e t r 
        | j <= 0 = e -- no need to shrink this excessively
        | (acceptanceRatio lTarget t0 t r0 r) ^^ a > 2 ^^ (-a) = 
            let (tn, rn) = leapfrog glTarget (t, r) e
            in  go (pred j) (2 ^^ a * e) tn rn 
        | otherwise = e

  return $ go 10 1.0 t1 r1

-- | Simulate a single step of Hamiltonian dynamics.
leapfrog :: Gradient -> Particle -> Double -> Particle
leapfrog glTarget (t, r) e = (tf, rf)
  where 
    rm = adjustMomentum glTarget e t r
    tf = adjustPosition e rm t
    rf = adjustMomentum glTarget e tf rm

-- | Adjust momentum.
adjustMomentum :: Fractional c => (t -> [c]) -> c -> t -> [c] -> [c]
adjustMomentum glTarget e t r = r .+ ((e / 2) .* glTarget t)

-- | Adjust position.
adjustPosition :: Num c => c -> [c] -> [c] -> [c]
adjustPosition e r t = t .+ (e .* r)

-- | The MH acceptance ratio for a given proposal.
acceptanceRatio :: Floating a => (t -> a) -> t -> t -> [a] -> [a] -> a
acceptanceRatio lTarget t0 t1 r0 r1 = auxilliaryTarget lTarget t1 r1
                                    / auxilliaryTarget lTarget t0 r0

-- | The negative potential. 
auxilliaryTarget :: Floating a => (t -> a) -> t -> [a] -> a
auxilliaryTarget lTarget t r = exp (lTarget t - 0.5 * innerProduct r r)

-- | Simple inner product.
innerProduct :: Num a => [a] -> [a] -> a
innerProduct xs ys = sum $ zipWith (*) xs ys

-- | Vectorized multiplication.
(.*) :: Num b => b -> [b] -> [b]
z .* xs = map (* z) xs

-- | Vectorized subtraction.
(.-) :: Num a => [a] -> [a] -> [a]
xs .- ys = zipWith (-) xs ys

-- | Vectorized addition.
(.+) :: Num a => [a] -> [a] -> [a]
xs .+ ys = zipWith (+) xs ys

-- | Indicator function.
indicate :: Integral a => Bool -> a
indicate True  = 1
indicate False = 0

-- | A symmetric categorical (discrete uniform) distribution.
symmetricCategorical :: PrimMonad m => [a] -> Gen (PrimState m) -> m a
symmetricCategorical [] _ = error "symmetricCategorical: no candidates"
symmetricCategorical zs g = do
  j <- uniformR (0, length zs - 1) g
  return $ zs !! j

-- | Alias for fromIntegral.
fi :: (Integral a, Num b) => a -> b
fi = fromIntegral

