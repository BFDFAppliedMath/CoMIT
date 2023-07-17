# Tests---------------------------------------------------------
test_that("Forward Primer", {
  expect_equal(translateMutationForward("C->T"), "C/A")
  expect_equal(translateMutationForward("G->T"), "G/A")
})

test_that("Reverse Primer", {
  expect_equal(translateMutationReverse("C->T"), "G/T")
  expect_equal(translateMutationReverse("G->A"), "C/A")
  expect_equal(translateMutationReverse("A->G"), "T/G")
  expect_equal(translateMutationReverse("G->T"), "C/T")
  expect_equal(translateMutationReverse("T->C"), "A/C")
})
