-- Various examples, using NUTS with dual-averaging.  Insert whatever trace 
-- (rosenbrockTrace, bnnTrace, etc.) you want into 'main' in order to spit out
-- some observations.
--
-- A convenient R script to display these traces:
--
-- require(ggplot2)
-- system('runhaskell Examples.hs > trace.csv')
-- d = read.csv('../tests/trace.dat', header = F)
-- ggplot(d, aes(V1, V2)) + geom_point(alpha = 0.05, col = 'darkblue')
--

module Numeric.MCMC.NUTS.Examples where

import Numeric.AD
import Numeric.MCMC.NUTS
import System.Random.MWC

logRosenbrock :: RealFloat a => [a] -> a
logRosenbrock [x0, x1] = negate (5 * (x1 - x0 ^ 2) ^ 2 + 0.05 * (1 - x0) ^ 2)

rosenbrockTrace :: IO [Parameters]
rosenbrockTrace = withSystemRandom . asGenST $ 
  nutsDualAveraging logRosenbrock (grad logRosenbrock) 10000 1000 [0.0, 0.0]

logHimmelblau :: RealFloat a => [a] -> a
logHimmelblau [x0, x1] = negate ((x0 ^ 2 + x1 - 11) ^ 2 + (x0 + x1 ^ 2 - 7) ^ 2)

himmelblauTrace :: IO [Parameters]
himmelblauTrace = withSystemRandom . asGenST $ 
  nutsDualAveraging logHimmelblau (grad logHimmelblau) 10000 1000 [0.0, 0.0]

logBnn :: RealFloat a => [a] -> a
logBnn [x0, x1] = -0.5 * (x0 ^ 2 * x1 ^ 2 + x0 ^ 2 + x1 ^ 2 - 8 * x0 - 8 * x1)

bnnTrace :: IO [Parameters]
bnnTrace = withSystemRandom . asGenST $ 
  nutsDualAveraging logBnn (grad logBnn) 10000 1000 [0.0, 0.0]

logBeale :: RealFloat a => [a] -> a
logBeale [x0, x1] 
  | and [x0 >= -4.5, x0 <= 4.5, x1 >= -4.5, x1 <= 4.5]
      = negate $  
          (1.5   - x0 + x0 * x1) ^ 2
        + (2.25  - x0 + x0 * x1 ^ 2) ^ 2
        + (2.625 - x0 + x0 * x1 ^ 3) ^ 2
  | otherwise = - (1 / 0)

bealeTrace :: IO [Parameters]
bealeTrace = withSystemRandom . asGenST $
  nutsDualAveraging logBeale (grad logBeale) 10000 1000 [0.0, 0.0]

printTrace :: Show a => [a] -> IO ()
printTrace = mapM_ (putStrLn . filter (`notElem` "[]") . show) 

main :: IO ()
main = bnnTrace >>= printTrace 

