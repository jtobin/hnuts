-- | See Hoffman, Gelman (2011) The No U-Turn Sampler: Adaptively Setting Path
--   Lengths in Hamiltonian Monte Carlo.

module Numeric.MCMC.NUTS where

import Control.Monad
import Control.Monad.Loops
import Control.Monad.Primitive
import System.Random.MWC
import System.Random.MWC.Distributions
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
    -- ++ "\n" ++ "rm: " ++ show rm
    ++ "\n" ++ "tp: " ++ show tp
    -- ++ "\n" ++ "rp: " ++ show rp
    ++ "\n" ++ "t': " ++ show t'
    ++ "\n" ++ "n : " ++ show n
    ++ "\n" ++ "s : " ++ show s

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
nuts lTarget glTarget m e t g = unfoldrM (kernel e) (m, t)
  where
    kernel eps (n, t0) = do
      t1 <- nutsKernel lTarget glTarget eps t0 g
      return $ if   n <= 0
               then Nothing
               else Just (t0, (pred n, t1))

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
  let logu = auxilliaryTarget lTarget t r0 - z0

  let go (tn, tp, rn, rp, j, tm, n, s) g
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

            let t2 | s1 == 1 
                  && (fi n1 / fi n :: Double) > z = t1
                   | otherwise                    = tm

                n2 = n + n1
                s2 = s1 * indicate ((tpp .- tnn) `innerProduct` rnn >= 0)
                        * indicate ((tpp .- tnn) `innerProduct` rpp >= 0) 
                j1 = succ j

            go (tnn, rnn, tpp, rpp, j1, t2, n2, s2) g

        | otherwise = return tm

  go (t, t, r0, r0, 0, t, 1, 1) g

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
  let (t0, r0)   = leapfrog glTarget (t, r) (v * e)
      lAuxTarget = log $ auxilliaryTarget lTarget t0 r0
      n          = indicate (logu < lAuxTarget)
      s          = indicate (logu - 1000 < lAuxTarget)
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
        t2     | accept    = t1
               | otherwise = t0 
        s2     = s1 * indicate ((tpp .- tnn) `innerProduct` rnn >= 0)
                    * indicate ((tpp .- tnn) `innerProduct` rpp >= 0)

    return (tnn, rnn, tpp, rpp, t2, n2, s2)
  else return (tn, rn, tp, rp, t0, n0, s0)

-- | Simulate a single step of Hamiltonian dynamics.
leapfrog :: Gradient -> Particle -> Double -> Particle
leapfrog glTarget (t, r) e = (tf, rf)
  where 
    rm = adjustMomentum glTarget e t r
    tf = adjustPosition e rm t
    rf = adjustMomentum glTarget e tf rm

-- | Adjust momentum.
adjustMomentum :: Fractional c => (t -> [c]) -> c -> t -> [c] -> [c]
adjustMomentum glTarget e t r = zipWith (+) r ((e / 2) .* glTarget t)

-- | Adjust position.
adjustPosition :: Num c => c -> [c] -> [c] -> [c]
adjustPosition e r t = zipWith (+) t (e .* r)

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

-- | Indicator function.
indicate :: Integral a => Bool -> a
indicate True  = 1
indicate False = 0

-- | A symmetric categorical (discrete uniform) distribution.
symmetricCategorical :: PrimMonad m => [a] -> Gen (PrimState m) -> m a
symmetricCategorical [] _ = error "symmetricCategorical: no candidates"
symmetricCategorical zs g = do
  z <- uniform g
  return $ zs !! truncate (z * fromIntegral (length zs) :: Double)

-- | Alias for fromIntegral.
fi :: (Integral a, Num b) => a -> b
fi = fromIntegral

