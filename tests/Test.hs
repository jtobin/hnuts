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

t0 :: [Double]
t0 = [1.0, 1.0]

r0 :: [Double]
r0 = [0.0, 0.0]

logu = -0.12840 -- from octave
u    = exp logu
v    = -1 :: Double

n = 20   :: Int
e = 0.1 :: Double

runBuildTree :: PrimMonad m => Gen (PrimState m) -> m BuildTree
runBuildTree g = do
  liftM BuildTree $ buildTree lTarget glTarget g t0 r0 u v n e

main = do
  test <- create >>= nuts lTarget glTarget 1000 0.1 t0
  mapM_ print test



