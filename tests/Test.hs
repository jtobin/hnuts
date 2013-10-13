import Control.Monad
import Control.Monad.Primitive
import Data.Vector (singleton)
import Numeric.AD
import Numeric.MCMC.NUTS
import System.Random.MWC

lTarget :: RealFloat a => [a] -> a
lTarget [x0, x1] = (-1) * (5 * (x1 - x0 ^ 2) ^ 2 + 0.05 * (1 - x0) ^ 2)

glTarget :: [Double] -> [Double]
glTarget = grad lTarget

-- glTarget [x, y] =
--   let dx = 20 * x * (y - x ^ 2) + 0.1 * (1 - x)
--       dy = -10 * (y - x ^ 2)
--   in  [dx, dy]

t0   = [0.0, 0.0] :: [Double]
r0   = [0.0, 0.0] :: [Double]
logu = -0.12840   :: Double
v    = -1         :: Double
n    = 5          :: Int
e    = 0.1        :: Double

runBuildTree :: PrimMonad m => Gen (PrimState m) -> m BuildTree
runBuildTree g = do
  liftM BuildTree $ buildTree lTarget glTarget g t0 r0 logu v n e

main = do
  test <- withSystemRandom . asGenIO $ nuts lTarget glTarget 20000 0.075 t0
  mapM_ print test



