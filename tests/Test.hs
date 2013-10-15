import Control.Lens
import Control.Monad
import Control.Monad.Primitive
import Data.Vector (singleton)
import Numeric.AD
import Numeric.MCMC.NUTS
import System.Random.MWC

import Debug.Trace

lTarget :: RealFloat a => [a] -> a
lTarget [x0, x1] = (-1) * (5 * (x1 - x0 ^ 2) ^ 2 + 0.05 * (1 - x0) ^ 2)

-- glTarget :: [Double] -> [Double]
-- glTarget = grad lTarget

glTarget [x, y] =
  let dx = 20 * x * (y - x ^ 2) + 0.1 * (1 - x)
      dy = -10 * (y - x ^ 2)
  in  [dx, dy]

t0   = [0.0, 0.0] :: [Double]
r0   = [0.0, 0.0] :: [Double]
logu = -0.12840   :: Double
v    = -1         :: Double
n    = 5          :: Int
m    = 1
madapt = 0
e    = 0.1        :: Double
eAvg = 1
h    = 0

t0da = [0.0, 0.0] :: [Double]
r0da = [0.0, 0.0] :: [Double]

runBuildTree :: PrimMonad m => Gen (PrimState m) -> m BuildTree
runBuildTree g = do
  liftM BuildTree $ buildTree lTarget glTarget g t0 r0 logu v n e

getThetas :: IO [Parameters]
getThetas = replicateM 1000 getTheta

getTheta :: IO Parameters
getTheta = do
  bt <- withSystemRandom . asGenIO $ \g ->
          buildTree lTarget glTarget g t0 r0 logu v n e  
  return $ bt^._5

getThetasDa :: IO [Parameters]
getThetasDa = replicateM 1000 getThetaDa

getThetaDa :: IO Parameters
getThetaDa = do
  bt <- withSystemRandom . asGenIO $ \g ->
          buildTreeDualAvg lTarget glTarget g t0 r0 logu v n e t0da r0da 
  return $ bt^._5

genMoveDa :: IO Parameters
genMoveDa = withSystemRandom . asGenIO $ \g -> do
    eps <- findReasonableEpsilon lTarget glTarget t0 g
    let daParams = basicDualAveragingParameters eps madapt
    blah <- nutsKernelDualAvg lTarget glTarget eps eAvg h m daParams t0 g
    return $  blah^._4

genMovesDa :: IO [Parameters]
genMovesDa = replicateM 1000 genMoveDa

genMove :: IO Parameters
genMove = withSystemRandom . asGenIO $ nutsKernel lTarget glTarget e t0 

genMoves :: IO [Parameters]
genMoves = replicateM 1000 genMove

main = do
  test <- withSystemRandom . asGenIO $ 
    nutsDualAveraging lTarget glTarget 5000 1000 t0
    -- nuts lTarget glTarget 5000 0.1 t0 
    -- genMovesDa
    -- genMoves
  mapM_ (putStrLn . filter (`notElem` "[]") . show) test
  -- mapM_ print test


