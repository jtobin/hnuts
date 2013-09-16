import Numeric.AD
import Numeric.MCMC.NUTS
import System.Random.MWC

lTarget :: RealFloat a => [a] -> a
lTarget [x0, x1] = (-1) * (5 * (x1 - x0 ^ 2) ^ 2 + 0.05 * (1 - x0) ^ 2)

glTarget :: [Double] -> [Double]
glTarget = grad lTarget

inits :: [Double]
inits = [0.0, 0.0]

epochs :: Int
epochs = 100


