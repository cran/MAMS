library(testthat)

toto <- mams()
toto5 <- mams(K = 5, nsim = 1000)

testthat::test_that("Simultaneous with user-define delta and delta0", {
    testthat::expect_no_warning(
    mams(method = "simultaneous", ushape = "obf",
                            nsim = 1000, J = 2,
                            delta = 0.8, delta0 = 0, sd = 1, 
                            p = NULL, p0 = NULL)
    )
})

testthat::test_that("Simultaneous with user-define p and p0", {
    testthat::expect_no_warning(
    mams(method = "simultaneous", ushape = "obf",
                            nsim = 1000, J = 2,
                            p = 0.73, p0 = 0.5)
    )
})

testthat::test_that("Simultaneous with user-define J=1 and p0=0.5", {
    testthat::expect_no_warning(
    mams(method = "simultaneous", ushape = "obf",
                            nsim = 1000, J = 1, r=1, r0=1,
                            p = 0.73, p0 = 0.5)
    )
})

testthat::test_that("Simultaneous with user-define J=1, K=5, delta, delta0", {
    testthat::expect_no_warning(
    mams(K=4, J=1, alpha=0.05, power=0.9, r=1, r0=1, p=NULL, 
                delta=0.545, delta0=0, sd=1, nsim = 1000)
    )
})

testthat::test_that("Testing mams switching from prob.scale to cohen.d", {
    
    testthat::expect_no_warning(
        toto5delta <- mams(obj=toto, K=5, delta=0.5, delta0=0, sd=1, 
                                nsim = 1000)
    )

    testthat::expect_no_warning(
        toto5p <- mams(obj=toto5delta, p=0.63, p0=0.5, nsim = 1000)
    )
})

test_that("Simulation with default parameters except K=5", {
    expect_silent(toto.sim <- mams.sim.simultaneous(obj = toto5, nsim=1000))
    expect_s3_class(toto.sim, "MAMS")
    expect_equal(toto.sim$K, 5)
})

test_that("Simultaneous simulation with specified `pv` parameters", {
    
    pv <- c(0.7, rep(0.5, 4))

    expect_silent(toto.sim.pv <- mams.sim.simultaneous(obj = toto5, 
    pv = pv, nsim = 1000))
    expect_equal(toto.sim.pv$input$par$pv, pv) 
})

test_that("Simultaneous simulation with specified `deltav` parameters", {
    deltav <- c(0.5, rep(0, 4))
    
    expect_silent(toto.sim.deltav <- mams.sim.simultaneous(obj = toto5, 
    deltav = deltav, nsim = 1000))
    expect_equal(toto.sim.deltav$input$par$deltav, deltav)
})

test_that("Simulation with incorrect number of deltav parameters", {
    deltav <- c(0.5, rep(0, 3)) 
    
    expect_error(
        toto.sim.deltav.wrong <- mams.sim.simultaneous(obj = toto5, 
                                                        deltav = deltav, 
                                                        nsim = 1000)
    )
})


test_that("mams error when provided object altered by mams_sim()", {
    toto.sim.deltav <- mams.sim.simultaneous(obj = toto5, 
                                                deltav = c(0.5, rep(0, 4)), 
                                                nsim = 1000)
    expect_error(mams(obj = toto.sim.deltav))
})

test_that("Example workflow without warnings", {
    expect_no_warning({
        toto.sim <- mams.sim.simultaneous(obj = toto5, nsim = 1000)

        toto.sim.pv <- mams.sim.simultaneous(obj = toto5, 
                                                    pv = c(0.7, rep(0.5, 4)), 
                                                    nsim = 1000)

        toto.sim.deltav <- mams.sim.simultaneous(obj = toto5, 
                                                    deltav = c(0.5, rep(0, 4)), 
                                                    nsim = 1000)
    })
})

test_that("mams() an error when nsim is less than 1000", {
    expect_error(
        mams(nsim = 999)
    )
})

test_that("mams handles additional parameters correctly", {
    expect_no_warning({
        toto5delta <- mams(obj = toto, K = 5, delta = 0.5, delta0 = 0, 
        sd = 1, nsim = 1000)
        toto5p <- mams(obj = toto5delta, p = 0.63, p0 = 0.5, nsim = 1000)
    })
    
    expect_equal(toto5delta$input$delta, 0.5)
    expect_equal(toto5delta$input$delta0, 0)
    expect_equal(toto5delta$input$sd, 1)
    expect_equal(toto5p$input$p, 0.63)
    expect_equal(toto5p$input$p0, 0.5)
})

test_that("Sequential simulations run without errors or warnings", {

    expect_no_warning(toto.sim <- mams.sim.simultaneous(obj = toto, 
                                                            nsim = 1000))

    expect_no_warning(mams.sim.simultaneous(obj = toto.sim, nsim = 1000))

    expect_no_warning(toto.sim3 <- mams.sim.simultaneous(obj = toto.sim, 
                                                    deltav = c(0.5, rep(0, 3)), 
                                                    nsim = 1000))
})

test_that("Parameter validation in mams_sim.simultaneous", {
    invalid_deltav <- "invalid_type"
    expect_error(
        mams.sim.simultaneous(obj = toto, deltav = invalid_deltav,
                                    nsim = 1000)
    )
})

test_that("Test summary function with complex pv or deltav", {
    
    deltav <- c(0.5, 0.4, 0, 0)
    toto.sim.deltav <- mams.sim.simultaneous(obj = toto, deltav = deltav, 
                                                    nsim = 1000)
    expect_no_warning(summary(toto.sim.deltav))

    pv <- c(0.73, 0.65, 0.5, 0.5)
    toto.sim.pv <- mams.sim.simultaneous(obj = toto, pv = pv, nsim = 1000)
    expect_no_warning(summary(toto.sim.pv))
})

test_that("Test automatic allocation for r and r0", {
expect_no_warning(toto1J  <- mams(toto, J = 1, nsim = 1000))
})

test_that("Default call test for automatically allocation of r and r0", {   
    expect_s3_class(toto, "MAMS")
    expect_equal(toto$J, 2)
    expect_equal(toto$input$r, 1:2)
    expect_equal(toto$input$r0, 1:2)
})

test_that("Passing obj as an unnamed argument without specifying J", {
    expect_no_warning(result_unnamed <- mams(toto, nsim = 1000))
    
    expect_s3_class(result_unnamed, "MAMS")
    expect_equal(result_unnamed$J, 2)
    expect_equal(result_unnamed$input$r, 1:2)
    expect_equal(result_unnamed$input$r0, 1:2)
})

test_that("Overriding J without specifying r and r0", {
    expect_no_warning(result_new_J <- mams(toto, J = 1, nsim = 1000))
    
    expect_s3_class(result_new_J, "MAMS")
    expect_equal(result_new_J$J, 1)
    expect_equal(result_new_J$input$r, 1)
    expect_equal(result_new_J$input$r0, 1)
})

test_that("Overriding J and providing custom r=1:2 and r0=2:3", {
expect_no_warning(result_custom_r <- mams(toto, J = 2, r = 1:2, r0 = 2:3))

    expect_equal(result_custom_r$J, 2)
    expect_equal(result_custom_r$input$r, 1:2)
    expect_equal(result_custom_r$input$r0, 2:3)
})

test_that("Check mams('method = simultaneous')
                            works with parameter H0=FALSE", {

    expect_no_warning(mams(H0 = FALSE))
})

test_that("Check summary works with parameter H0=FALSE", {
    expect_no_warning(summary(mams(H0 = FALSE)))
})


test_that("Check summary works with switching methods and parameters", {

    expect_no_warning({
    mams.fit1a = mams(method="sep",ushape="pocock", nsim=1000)
    mams.fit1b = mams(mams.fit1a, ushape="obf", nsim=1000)
    mams.fit1d = mams(mams.fit1b, J=3, K=3, delta = c(.25),sd=.5, delta0=0)
    print(mams.fit1d)
    summary(mams.fit1d)
    
    mams.fit0 = mams()
    mams.fit1a = mams(mams.fit0, method="sep",ushape="pocock", nsim=5000)
    mams.fit1b = mams.sim(mams.fit1a, deltav = c(.9,.5,.25,.1))
    summary(mams.fit1b)
    
    mams.fit1a = mams(mams.fit0, ushape="pocock", nsim=5000)
    mams.fit1b = mams.sim(mams.fit1a, deltav = c(.9,.5,.25,.1))
    summary(mams.fit1b)
    })

})

test_that("mams function calls with method='simultaneous' 
                                                without predefined delta0", {
    expect_no_error({
    mams.fit1a <- mams(method = "sep", ushape = "pocock", nsim = 1000)
    mams.fit1b <- mams(mams.fit1a, ushape = "obf", nsim = 1000)
        suppressWarnings({
            mams.fit1e <- mams(mams.fit1b, J = 2, K = 3, delta = c(0.25), 
                                sd = 0.5, nsim = 1000, method = "simultaneous")
        })
    })
})

test_that("Overriding J=2 without specifying r and r0", {
    expect_no_error({
    suppressWarnings({
mams.fit6e = mams(J=3, nsim=1000)
mams.fit6f = mams(mams.fit6e, delta0=0.075)
mams.fit6f$sim$H1$main$ess[1,4] == 30
mams.fit6g = mams(mams.fit6e, p0=0.521)
mams.fit6g$sim$H1$main$ess[1,4] == 30
    })
    })
})

test_that("mams checks for r0 and r vectors", {
    expect_error(
    mams(K = 3, J = 3, alpha = 0.05, power = 0.9, r = c(1, 2), 
    r0 = c(1, 2)),
    "Length of allocation ratios does not match number of stages."
    )

    expect_error(
    mams(K = 3, J = 2, alpha = 0.05, power = 0.9, r = c(1, 2), 
    r0 = c(3, 2)),
    "`r0` must be a monotonically increasing vector."
    )

    expect_error(
    mams(K = 3, J = 2, alpha = 0.05, power = 0.9, r = c(3, 2), 
    r0 = c(1, 2)),
    "`r` must be a monotonically increasing vector."
    )
    
    expect_error(
        mams(K = 3, J = 3, alpha = 0.05, power = 0.9, r = c(0, 2, 3), 
        r0 = c(1, 2, 3))
    )

    expect_error(
        mams(K = 3, J = 2, alpha = 0.05, power = 0.9, r = c(1, 2), 
        r0 = c(0, 2))
    )

    expect_error(
        mams(K = 3, J = 2, alpha = 0.05, power = 0.9, r = c(1, 2), 
        r0 = c(1.5, 2)),
        "First element of `r0` must be integers."
    )
})

test_that("Peter Greenstreet's example of allocation ratios", {
toto_1 <- mams(
    K = 3, J = 2,
    power = 0.8, r = c(3, 7), r0 = c(3, 7), ushape = "obf",
    lshape = "fixed", lfix = -9999, nstart = 1, nsim = 1000
)
expect_equal(toto_1$n, 9)
expect_equal(toto_1$rMat[1,1], 1)
expect_equal(toto_1$rMat[1,2], 2.333333, tolerance = 0.01)

})

### Extra tests for input parameters checks

test_that("sep method: allocation ratio vector length is checked", {
    expect_error(
    mams(K = 4, J = 2, alpha = 0.05, power = 0.9, r = c(1), r0 = c(1), 
            method = "sep"),
    "Length of allocation ratios does not match number of stages."
    )
    expect_error(
    mams(K = 4, J = 2, alpha = 0.05, power = 0.9, r = c(1, 2, 3), 
    r0 = c(1,2,3), 
    method = "sep"),
    "Length of allocation ratios does not match number of stages."
    )
})

test_that("sep method: allocation ratio vectors must be monotonically increasing
    and valid", {
# r0 not monotonically increasing
    expect_error(
    mams(K = 4, J = 2, alpha = 0.05, power = 0.9, r = c(1, 2), r0 = c(2, 1),
                method = "sep"),
    "`r0` must be a monotonically increasing vector."
    )
# r values must be >= 1
    expect_error(
    mams(K = 4, J = 2, alpha = 0.05, power = 0.9, r = c(0, 2), r0 = c(1, 2),
                method = "sep"),
    "`r` values must be >= 1."
    )
# r0[1] must be >= 1
    expect_error(
    mams(K = 4, J = 2, alpha = 0.05, power = 0.9, r = c(1, 2), r0 = c(0, 2),
                method = "sep"),
    "`r0\\[1\\]` must be >= 1."
    )
# All elements of r0 must be integers (via first element check message)
    expect_error(
    mams(K = 4, J = 2, alpha = 0.05, power = 0.9, r = c(1,2), 
                r0 = c(1.5, 2)),
    "First element of `r0` must be integers."
    )
})

test_that("sep method: effect size parameter check", {
# Both p and delta, sd provided => error
    expect_warning(expect_error(
    mams(K = 4, J = 2, alpha = 0.05, power = 0.9,
                r = 1:2, r0 = 1:2, 
                p = 0.75, delta = 0.5, sd = 1)
    ))

# p out-of-bound error: p < 0 or > 1
    expect_error(
    mams(K = 4, J = 2, alpha = 0.05, power = 0.9,
            r = 1:2, r0 = 1:2, 
            p = -0.1),
    "Treatment effect parameter 'p0' not within 0.5 and 1."
)

# sd must be positive when using delta scale.
    expect_warning(expect_error(
    mams(K = 4, J = 2, alpha = 0.05, power = 0.9,
            r = 1:2, r0 = 1:2, 
            delta = 0.5, sd = -1),
    "Standard deviation must be positive."
    ))
})

test_that("sep method: alpha must be between 0 and 1", {
    expect_error(
        mams(K = 4, J = 2, alpha = -0.05, power = 0.9, r = 1:2, r0 = 1:2),
        "Error rate or power not between 0 and 1."
    )
    expect_error(
    mams(K = 4, J = 2, alpha = 1.05, power = 0.9, r = 1:2, r0 = 1:2),
    "Error rate or power not between 0 and 1."
    )
})

test_that("K=3, J=1 yields correct power, type I error, and ESS", {
  design1a <- mams(
    K = 3, J = 1,
    alpha = 0.05, power = 0.8,
    delta = 2, delta0 = 0, sd = 10
  )

  # The power should be ≈ 0.8 
  eff_H1 <- as.numeric(design1a$sim$H1$main$efficacy[5, ])
  expect_equal(eff_H1, 0.8, tolerance = 0.02)

  # The p-value should be ≈ 0.05 
  eff_H0 <- as.numeric(design1a$sim$H0$main$efficacy[4, ])
  expect_equal(eff_H0, 0.05, tolerance = 0.01)

  # The ESS should be 422
  ess_H1 <- design1a$sim$H1$main$ess[1:4, 1]
  expect_true(all(ess_H1 == 422))
})

test_that("K=1, J=1 yields correct power, type I error, and ESS", {
  design1a_k1 <- mams(
    K = 1, J = 1,
    alpha = 0.05, power = 0.8,
    delta = 2, delta0 = 0, sd = 10
  )

  # The power should be ≈ 0.8 
  eff_H1_k1 <- as.numeric(design1a_k1$sim$H1$main$efficacy[3, ])
  expect_equal(eff_H1_k1, 0.8, tolerance = 0.02)

  # The p-value should be ≈ 0.05
  eff_H0_k1 <- as.numeric(design1a_k1$sim$H0$main$efficacy[2, ])
  expect_equal(eff_H0_k1, 0.05, tolerance = 0.01)

  # The ESS should be 310
  ess_H1_k1 <- design1a_k1$sim$H1$main$ess[1:2, 1]
  expect_true(all(ess_H1_k1 == 310))
})
