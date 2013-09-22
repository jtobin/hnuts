-- | See Hoffman, Gelman (2011) The No U-Turn Sampler: Adaptively Setting Path
--   Lengths in Hamiltonian Monte Carlo.

module Numeric.MCMC.NUTS where

import Control.Monad
import Control.Monad.Primitive
import System.Random.MWC -- FIXME change to Prob monad
import System.Random.MWC.Distributions
import Statistics.Distribution.Normal

-- FIXME change to probably api
type Parameters = [Double] 
type Density    = Parameters -> Double
type Gradient   = Parameters -> Parameters
type Particle   = (Parameters, Parameters)

-- FIXME must be streaming
nuts :: PrimMonad m
     => Density
     -> Gradient
     -> Int
     -> Parameters
     -> Gen (PrimState m)
     -> m Parameters
nuts lTarget glTarget m t g = do
  e <- findReasonableEpsilon lTarget glTarget t g
  let go 0 t0 = return t0
      go n t0 = nutsKernel lTarget glTarget e t0 g >>= go (pred n)

  go m t

nutsKernel :: PrimMonad m 
           => Density
           -> Gradient
           -> Double
           -> Parameters
           -> Gen (PrimState m)
           -> m Parameters
nutsKernel lTarget glTarget e t g = do
  r0 <- replicateM (length t) (normal 0 1 g)
  u  <- uniformR (0, auxilliaryTarget lTarget t r0) g

  let go (tn, tp, rn, rp, j, tm, n, s) g
        | s == 1 = do
            vj <- symmetricCategorical [-1, 1] g
            z  <- uniform g

            (tnn, rnn, tpp, rpp, t1, n1, s1) <-
              if   vj == -1
              then do
                (tnn', rnn', _, _, t1', n1', s1') <- 
                  buildTree lTarget glTarget g tn rn u vj j e
                return (tnn', rnn', tp, rp, t1', n1', s1')
              else do
                (_, _, tpp', rpp', t1', n1', s1') <- 
                  buildTree lTarget glTarget g tp rp u vj j e
                return (tn, rn, tpp', rpp', t1', n1', s1')

            let t2 | s1 == 1 && min 1 (fi n1 / fi n :: Double) > z = tnn
                   | otherwise = t

                n2 = n + n1
                s2 = s1 * indicate ((tpp .- tnn) `innerProduct` rnn >= 0)
                        * indicate ((tpp .- tnn) `innerProduct` rpp >= 0) 
                j1 = succ j

            go (tnn, rnn, tpp, rpp, j1, t2, n2, s2) g

        | otherwise = return tm

  go (t, t, r0, r0, 0, t, 1, 1) g

newtype BuildTree = BuildTree { 
    getBuildTree :: ([Double], [Double], [Double], [Double], [Double], Int, Int)
  }

instance Show BuildTree where
  show (BuildTree (tm, rm, tp, rp, t', n, s)) = 
       "\n" ++ "tm: " ++ show tm 
    ++ "\n" ++ "rm: " ++ show rm
    ++ "\n" ++ "tp: " ++ show tp
    ++ "\n" ++ "rp: " ++ show rp
    ++ "\n" ++ "t': " ++ show t'
    ++ "\n" ++ "n : " ++ show n
    ++ "\n" ++ "s : " ++ show s

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
buildTree lTarget glTarget g = go
  where
    go t r u v 0 e = return $
      let (t0, r0) = leapfrog glTarget (t, r) (v * e)
          auxTgt   = auxilliaryTarget lTarget t0 r0
          n        = indicate (u <= auxTgt)
          s        = indicate (auxTgt > log u - 1000)
      in  (t0, r0, t0, r0, t0, n, s)

    go t r u v j e = do
      z <- uniform g
      (tn, rn, tp, rp, t0, n0, s0) <- go t r u v (pred j) e

      if   s0 == 1
      then do
        (tnn, rnn, tpp, rpp, t1, n1, s1) <- 
          if   v == -1
          then do
            (tnn', rnn', _, _, t1', n1', s1') <- go tn rn u v (pred j) e
            return (tnn', rnn', tp, rp, t1', n1', s1')
          else do
            (   _,    _, tpp', rpp', t1', n1', s1') <- go tp rp u v (pred j) e
            return (tn, rn, tpp', rpp', t1', n1', s1')

        let p  = fromIntegral n1 / fromIntegral (n0 + n1)
            n2 = n0 + n1
            t2 | p > (z :: Double) = t1
               | otherwise         = t0 
            s2 = s1 * indicate ((tpp .- tnn) `innerProduct` rnn >= 0)
                    * indicate ((tpp .- tnn) `innerProduct` rpp >= 0)

        return (tnn, rnn, tpp, rpp, t2, n2, s2)
      else return (tn, rn, tp, rp, t0, n0, s0)

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

      go e t r | (acceptanceRatio lTarget t0 t r0 r) ^^ a > 2 ^^ (-a) = 
                   let (tn, rn) = leapfrog glTarget (t, r) e
                   in  go (2 ^^ a * e) tn rn 
               | otherwise = e

  return $ go 1.0 t1 r1

leapfrogIntegrator :: Int -> Gradient -> Particle -> Double -> Particle
leapfrogIntegrator n glTarget particle e = go particle n
  where go state ndisc 
          | ndisc <= 0 = state
          | otherwise  = go (leapfrog glTarget state e) (pred n)

leapfrog :: Gradient -> Particle -> Double -> Particle
leapfrog glTarget (t, r) e = (tf, rf)
  where rm = adjustMomentum glTarget e t r
        tf = adjustPosition e rm t
        rf = adjustMomentum glTarget e tf rm

adjustMomentum :: Fractional c => (t -> [c]) -> c -> t -> [c] -> [c]
adjustMomentum glTarget e t r = zipWith (+) r ((e / 2) .* glTarget t)

adjustPosition :: Num c => c -> [c] -> [c] -> [c]
adjustPosition e r t = zipWith (+) t (e .* r)

acceptanceRatio :: Floating a => (t -> a) -> t -> t -> [a] -> [a] -> a
acceptanceRatio lTarget t0 t1 r0 r1 = auxilliaryTarget lTarget t1 r1
                                    / auxilliaryTarget lTarget t0 r0

auxilliaryTarget :: Floating a => (t -> a) -> t -> [a] -> a
auxilliaryTarget lTarget t r = exp (lTarget t - 0.5 * innerProduct r r)

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

