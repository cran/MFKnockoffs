library(testthat)
library(MFKnockoffs)

# Set a random seed for reproducibility when running R CMD CHECK.
# Note that this file is not run by devtools::test, so the tests will be truly
# random when run from RStudio.
set.seed(1234)

# Run the test suite.
test_check('MFKnockoffs')

# Reset the random seed.
rm(.Random.seed)