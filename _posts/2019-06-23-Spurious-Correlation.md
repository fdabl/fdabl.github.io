---
layout: post
title: "Spurious correlations and random walks"
date: 2019-06-29 11:00:00 +0100
categories: R
status: publish
published: true
# status: development
# published: false
---
 
The number of storks and the number of human babies delivered are positively correlated (Matthews, 2000). This is a classic example of a spurious correlation which most certainly has a causal explanation: a third variable, say economic development, is most likely to cause both an increase in storks and an increase in the number of human babies, hence the correlation.[^1] In this blog post, I discuss a more subtle case of spurious correlation, one that is not of causal but of statistical nature: *completely independent processes can be correlated substantially*.
 
## AR(1) processes and random walks
Moods, stockmarkets, the weather: everything changes, everything is in flux. The simplest model to describe change is an auto-regressive (AR) process of order one. Let $Y_t$ be a random variable where $t = [1, \ldots T]$ indexes discrete time. We write an AR(1) process as:
 
$$
Y_t = \phi \, Y_{t-1} + \epsilon_t \enspace ,
$$
 
where $\phi$ gives the correlation with the previous observation, and where $\epsilon_t \sim \mathcal{N}(0, \sigma^2)$. For $\phi = 1$ the process is called a *random walk*. We can simulate from these using the following code:
 

{% highlight r %}
simulate_ar <- function(n, phi, sigma = .1) {
  y <- rep(0, n)
  
  for (t in seq(2, n)) {
    y[t] <- phi*y[t-1] + rnorm(1, 0, sigma)
  }
  
  y
}
{% endhighlight %}
 
 
The following R code simulates data from three independent random walks and an AR(1) process with $\phi = 0.5$; the Figure below visualizes them.
 

{% highlight r %}
n <- 100
set.seed(1)
 
rw1 <- simulate_ar(n, phi = 1)
rw2 <- simulate_ar(n, phi = 1)
rw3 <- simulate_ar(n, phi = 1)
ar <- simulate_ar(n, phi = 0.5)
{% endhighlight %}
 
<img src="/assets/img/2019-06-23-Spurious-Correlation.Rmd/unnamed-chunk-3-1.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />
 
As we can see from the plot, the AR(1) process seems pretty well-behaved. This is in contrast to the three random walks: all of them have an initial upwards trend, after which the red line keeps on growing, while the blue line makes a downward jump. In contrast to AR(1) processes, random walks are *not stationary* since their variance is not constant across time. For some very good lecture notes on time-series analysis, see [here](https://www.economodel.com/time-series-analysis).
 
 
## Spurious correlations of random walks
If we look at the correlations of these three random walks across time points, we find that they are substantial:
 

{% highlight r %}
round(cor(cbind(red = rw1, green = rw2, blue = rw3)), 2)
{% endhighlight %}



{% highlight text %}
##         red green  blue
## red    1.00 -0.49 -0.29
## green -0.49  1.00  0.59
## blue  -0.29  0.59  1.00
{% endhighlight %}
 
I hope that this is at least a little bit of a shock. Upon reflection, however, it is clear that we are blundering: computing the correlation across time ignores the dependency between data points that is so typical of time-series data. To get more data about what is going on, we conduct a small simulation study.
 
In particular, we want to get an intuition of how this spurious correlation behaves with increasing sample sizes. We therefore simulate two independent random walks for sample sizes $n \in [50, 100, 200, 500, 1000, 2000]$ and compute their Pearson correlation, the test-statistic, and whether $p < \alpha$, where we set $\alpha$ to some an arbitrary value, say $\alpha = 0.05$. We repeated this 100 times and report the average of these quantities.
 

{% highlight r %}
library('dplyr')
set.seed(1)
 
times <- 100
ns <- c(50, 100, 200, 500, 1000, 2000)
 
comb <- expand.grid(times = seq(times), n = ns)
ncomb <- nrow(comb)
 
res <- matrix(NA, nrow = ncomb, ncol = 5)
colnames(res) <- c('ix', 'n', 'cor', 'tstat', 'pval')
 
for (i in seq(ncomb)) {
  n <- comb[i, 2]
  ix <- comb[i, 1]
  test <- cor.test(simulate_ar(n, phi = 1), simulate_ar(n, phi = 1))
  res[i, ] <- c(ix, n, test$estimate, test$statistic, test$p.value)
}
 
tab <- data.frame(res) %>% 
  group_by(n) %>% 
  summarize(
    avg_abs_corr = mean(abs(cor)),
    avg_abs_tstat = mean(abs(tstat)),
    percent_sig = mean(pval < .05)
  ) %>% data.frame
 
round(tab, 2)
{% endhighlight %}



{% highlight text %}
##      n avg_abs_corr avg_abs_tstat percent_sig
## 1   50         0.41          3.57        0.71
## 2  100         0.46          6.58        0.85
## 3  200         0.45          8.88        0.85
## 4  500         0.37         10.63        0.86
## 5 1000         0.41         17.05        0.88
## 6 2000         0.39         23.39        0.97
{% endhighlight %}
 
We observe that the average absolute correlation is very similar across $n$, but the test statistic grows with increased $n$, which naturally results in many more false rejections of the null hypothesis of no correlation between the two random walks.
 
To my knowledge, Granger and Newbold (1974) were the first to point out this puzzling fact.[^2] They regress one random walk onto the other instead of computing the Pearson correlation. (Note that the test statistic is the same). In a regression setting, we write:
 
$$
Y = \beta_0 + \beta_1 X + \epsilon \enspace ,
$$
 
where we assume that $\epsilon_i \sim \mathcal{N}(0, \sigma^2)$ (see also a [previous](https://fdabl.github.io/r/Curve-Fitting-Gaussian.html) blog post). This is evidently violated when performing linear regression on two random walks, as demonstrated by the residual plot below.
 
<img src="/assets/img/2019-06-23-Spurious-Correlation.Rmd/unnamed-chunk-6-1.png" title="plot of chunk unnamed-chunk-6" alt="plot of chunk unnamed-chunk-6" style="display: block; margin: auto;" />
 
Similar as above, we can have an AR(1) process on the residuals:
 
$$
\epsilon_t = \delta \epsilon_{t-1} + \eta_t \enspace ,
$$
 
and test whether $\delta = 0$. We can do so using the [Durbin-Watson test](https://en.wikipedia.org/wiki/Durbin%E2%80%93Watson_statistic), which yields:
 

{% highlight r %}
car::durbinWatsonTest(m)
{% endhighlight %}



{% highlight text %}
##  lag Autocorrelation D-W Statistic p-value
##    1       0.9357562    0.08623868       0
##  Alternative hypothesis: rho != 0
{% endhighlight %}
 
This indicates substantial autocorrelation, violating our modeling assumption of independent residuals. In the next section, we look at the deeper mathematical reasons for why we get such spurious correlation. In the Post Scriptum, we relax the constraint that $\phi = 1$ and look at how spurious correlation behaves for AR(1) processes.
 
 
<!-- In the next section, we will look more formally into the curious fact that two independent random walks are correlated. To understand why even with large $n$ the estimation goes awry, we have to make an excursion into asymptotia. -->
 
 
## Inconsistent estimation
The simulation results from the random walk simulations showed that the average (absolute) correlation stays roughly constant, while the test statistic increases with $n$. This indicates a problem with our estimator for the correlation. Because it is slightly easier to study, we focus on the regression parameter $\beta_1$ instead of the Pearson correlation. [Recall](https://fdabl.github.io/r/Curve-Fitting-Gaussian.html) that our regression estimate is
 
$$
\hat{\beta}_1 = \frac{\sum_{t=1}^N (x_t - \bar{x})(y_t - \bar{y})}{\sqrt{\sum_{t=1}^N (x_t - \bar{x})^2 \sum_{t=1}^N (y_t - \bar{y})^2}} \enspace ,
$$
 
where $\bar{x}$ and $\bar{y}$ are the empirical means of the realizations $x_t$ and $y_t$ of the AR(1) processes $X_t$ and $Y_t$, respectively. The test statistic associated with the null hypothesis $\beta_1 = 0$ is
 
$$
t_{\text{statistic}} := \frac{\hat{\beta_1} - 0}{se(\hat{\beta_1})} = \frac{\hat{\beta_1}}{\hat{\sigma} / \sqrt{\sum_{t=1}^N (x_t - \bar{x})^2}} \enspace ,
$$
 
where $\hat{\sigma}$ is the estimated standard deviation of the error. In simple linear regression, the test statistic follows a t-distribution with $n - 2$ degrees of freedom (it takes two parameters to fit a straight line). In the case of independent random walks, however, the test statistic does not have a limiting distribution; in fact, as $n \rightarrow \infty$, the distribution of $t_{\text{statistic}}$ diverges (Phillips, 1986).
 
To get an intuition for this, we plot the bootstrapped sampling distributions for $\beta_1$ and $t_{\text{statistic}}$, both for the case of regressing one independent AR(1) process onto another, and for random walk regression.
 

{% highlight r %}
regress_ar <- function(n, phi, sigma) {
  y <- simulate_ar(n, phi, sigma)
  x <- simulate_ar(n, phi, sigma)
  coef(summary(lm(y ~ x)))[2, c(1, 3)]
}
 
bootstrap_limit <- function(ns, phi, sigma, times = 200) {
  res <- matrix(NA, nrow = times * length(ns), ncol = 3)
  colnames(res) <- c('n', 'b1', 'tstat')
 
  ix <- 1
  for (n in ns) {
    for (i in seq(times)) {
      coefs <- regress_ar(n, phi, sigma)
      res[ix, ] <- c(n, coefs)
      ix <- ix + 1
    }
  }
 
  data.frame(res)
}
 
 
set.seed(1)
ns <- c(100, 200, 500, 1000, 2000)
res_ar <- bootstrap_limit(ns, .5, .1)
res_rw <- bootstrap_limit(ns,  1, .1)
{% endhighlight %}
 
The Figure below illustrates how things go wrong when regressing one independent random walk onto the other. In contrast to the estimate for the AR(1) regression, the estimate $\hat{\beta}_1$ does not decrease in the case of a random walk regression. Instead, it stays roughly within $[-0.75, 0.75]$ across all $n$. This shines further light on the initial simulation results that the average correlation stays roughly the same. Moreover, in contrast AR(1) regression for which the distribution of the test statistic does not change, the distribution of the test statistic for the random walk regression seems to diverge. This explains why we the proportion of false positives increases with $n$.
 
<img src="/assets/img/2019-06-23-Spurious-Correlation.Rmd/unnamed-chunk-9-1.png" title="plot of chunk unnamed-chunk-9" alt="plot of chunk unnamed-chunk-9" style="display: block; margin: auto;" />
 
Rigorous arguments of the above statements can be found in Phillips (1986) and Hamilton (1994, pp. 577).[^4] The explanations feature some nice asympotic arguments which I would love go into in detail; however, I'm currently in Santa Fe for a summer school that has a very tightly packed programme. On that note: it is [very, very cool](https://www.santafe.edu/engage/learn/schools/sfi-complex-systems-summer-school). You should definitely apply next year! In addition to the stimulating lectures, wonderful people, and exciting projects, the surroundings are stunning[^5].
 
<div style="text-align:center;">
  <img src="../assets/img/IAIA.jpeg" align="center" style="padding-top: 10px; padding-bottom: 10px;" width="720" height="620" />
</div>
 
 
<!-- ### Brownian Motion -->
<!-- The type of random walk we focused on in this blog post takes place in discrete, equidistant time steps.[^3] If we take the limit of $n \rightarrow \infty$, however, we move from a discrete time random walk to a continuous time Brownian motion. The gist of the argument is to make the difference $\Delta Y_t$ between time points $Y_{t+1}$ and $Y_t$ infinitesimally small. Recall that the Gaussian distribution is [closed under addition](https://fdabl.github.io/statistics/Two-Properties.html), and that -->
 
<!-- $$ -->
<!-- \begin{aligned} -->
<!-- Y_t &= \sum_{i=1}^t \eta_i \sim \mathcal{N}(0, t \cdot \sigma^2) \enspace \\[1em] -->
<!-- \Delta Y_t &= Y_{t+1} - Y_{t} = \sum_{i=1}^{t+1} \eta_i - \sum_{j=1}^t \eta_j = \eta_t \sim \mathcal{N}(0, \sigma^2) \enspace . -->
<!-- \end{aligned} -->
<!-- $$ -->
 
<!-- We may cut $\eta_t$ into $n$ pieces -->
 
<!-- $$ -->
<!-- \eta_t = \eta_{1t} + \eta_{2t} + \ldots + \eta_{nt} \enspace , -->
<!-- $$ -->
 
<!-- where $\eta_{it} \sim \mathcal{N}(0, \frac{1}{n})$. Therefore, as we increase $n$, the discrete-time process is defined at a finer and finer grid. For $n \rightarrow \infty$, this results into the continuous-time Brownian motion, which we denote as $W(t)$, where $W: t \in [0, 1] \rightarrow \mathbb{R}$. -->
 
 
<!-- ## Solutions -->
<!-- Hamilton (1994, p. 562) discusses three solutions. One of them is to *difference* the data before doing the regression, i.e., -->
 
<!-- $$ -->
<!-- \Delta Y_t = \beta_0 + \beta_1 \Delta X_t + \epsilon_t \enspace , -->
<!-- $$ -->
 
<!-- where $\Delta Y_t = Y_{t+1} - Y_t$. This does in fact work: -->
 
<!-- ```{r} -->
<!-- broom::tidy(lm(diff(rw1) ~ diff(rw2))) -->
<!-- ``` -->
 
<!-- ```{r, echo = FALSE} -->
<!-- n <- 1000 -->
<!-- dat <- matrix(0, nrow = n, ncol = 2) -->
<!-- B <- cbind( -->
<!--   c(.4, .2), -->
<!--   c(-.2, .4) -->
<!-- ) -->
 
<!-- for (i in seq(2, n)) { -->
<!--   z <- rnorm(1) -->
<!--   # dat[i, ] <- dat[i-1, ] %*% B + rnorm(2) -->
<!--   dat[i, ] <- c(.8, .4) * z + rnorm(2) -->
<!-- } -->
<!-- ``` -->
 
<!-- Why? Let $\eta_t$ and $\psi_t$ denote the errors of the two processes $Y$ and $X$, respectively, distributed according to zero-mean Gaussian with variances $\sigma_y$ and $\sigma_x$. We write -->
 
<!-- $$ -->
<!-- \Delta Y_t = \sum_{i=1}^{t+1} \eta_i - \sum_{i=1}^{t} \eta_i = \eta_{t+1} \sim \mathcal{N}(0, \sigma_y^2) \\[1em] -->
<!-- \Delta X_t = \sum_{i=1}^{t+1} \psi_i - \sum_{i=1}^{t} \eta_i = \psi_{t+1} \sim \mathcal{N}(0, \sigma_x^2) \enspace . -->
<!-- $$ -->
 
<!-- Now, since the respective differences are independent of each other, their correlation will be zero. -->
 
<!-- However, Hamilton notes that if the time-series are really stationary ($\vert \phi \lvert < 1$), then this can result in misspecified regression. Moreover, if $Y$ and $X$ are non-stationary but *cointegrated processes*, then this also will result in misspecification. -->
 
 
## Conclusion
"Correlation does not imply causation" is a common response to apparently spurious correlation. The idea is that we observe spurious associations because we do not have the full causal picture, as in the example of storks and human babies. In this blog post, we have seen that spurious correlation can be due to solely statistical reasons. In particular, we have seen that two independent random walks can be highly correlated. This can be diagnosed by looking at the residuals, which will *not* be independent and identically distributed, but will show a pronounced autocorrelation.
 
The mathematical explanation for the spurious correlation is not trivial. Using simulations, we found that the estimate of $\beta_1$ does not converge to the true value in the case of regressing one independent random walk onto another. Moreover, the test statistic diverges, meaning that with increasing sample size we are almost certain to reject the null hypothesis of no association. The spurious correlation occurs because our estimate is not consistent, which is a purely statistical explanation that does not invoke causal reasoning.
 
---
*I want to thank Toni Pichler and Andrea Bacilieri for helpful comments on this blog post.*
 
---
## Post Scriptum
<!-- ### Mean and variance of AR(1) and random walk -->
<!-- To better understand the differences between AR(1) processes and random walks, we look at their respective first two moments. We write out the process for some window of length $j$, and then recursively substitute: -->
 
<!-- $$ -->
<!-- \begin{aligned} -->
<!-- Y_t &= \phi \, Y_{t-1} + \epsilon_t \\[.5em] -->
<!--     &= \phi \, \left(\phi \, Y_{t-2} + \epsilon_{t-1}\right) + \epsilon_t \\[.5em] -->
<!--     &= \phi \, \left(\phi \, \left(\phi \, Y_{t-3} + \epsilon_{t-2}\right) + \epsilon_{t-1}\right) + \epsilon_t \\[.5em] -->
<!--     &= \vdots \\[.5em] -->
<!--     &= \phi^{j + 1} \, Y_{t - (j + 1)} + \sum_{i=t}^{t - (j + 1)} \phi^i \epsilon_{t-i} \\[.5em] -->
<!--     &= \sum_{i=0}^{t-1} \phi^i \epsilon_{t-i} \enspace , -->
<!-- \end{aligned} -->
<!-- $$ -->
 
<!-- where we assume that $Y_0 = 0$ is fixed. Let's compute the first two moments of this process. Exploiting linearity, we write: -->
 
 
<!-- $$ -->
<!-- \mathbb{E}[Y_t] = \mathbb{E}\left[\sum_{i=0}^{t-1} \phi^i \epsilon_{t-i}\right] = \sum_{i=0}^{t-1} \mathbb{E}\left[\phi^i \epsilon_{t-i}\right] = \sum_{i=0}^{t-1} \phi^i \mathbb{E}\left[\epsilon_{t-i}\right] = 0 \enspace . -->
<!-- $$ -->
 
<!-- This is also true for $\phi = 1$, i.e., a random walk. For the variance, we write: -->
 
<!-- $$ -->
<!-- \begin{aligned} -->
<!-- \text{Var}\left[Y_t\right] &= \mathbb{E}\left[\left(Y_t - \mathbb{E}[Y_t]\right)^2\right]  -->
<!-- = \mathbb{E}\left[Y_t^2\right] -->
<!-- = \mathbb{E}\left[\left(\sum_{i=0}^{t-1} \phi^i \epsilon_{t-i}\right)^2\right] \enspace , -->
<!-- \end{aligned} -->
<!-- $$ -->
 
<!-- where we split the quadratic into ["diagonal"](https://math.stackexchange.com/questions/125435/what-is-the-opposite-of-a-cross-term) terms and cross-terms, the latter of which have expectation zero by our assumption that the residuals are independent: -->
 
 
<!-- $$ -->
<!-- \begin{aligned} -->
<!-- \text{Var}\left[Y_t\right] &= \mathbb{E}\left[\sum_{i=0}^{t - 1} \left(\phi^i \epsilon_{t-i}\right)^2 + \sum_{i=0}^{t - 1} \sum_{j\neq i}^{t - 1} \left(\phi^i \epsilon_{t-i}\right) \left(\phi^j \epsilon_{t-j}\right)\right] \\[.5em] -->
<!-- &= \mathbb{E}\left[\sum_{i=0}^{t - 1} \left(\phi^i \epsilon_{t-i}\right)^2\right] \\[.5em] -->
<!-- &= \sum_{i=0}^{t - 1} \mathbb{E}\left[\left(\phi^i \epsilon_{t-i}\right)^2\right] \\[.5em] -->
<!-- &= \sum_{i=0}^{t - 1} \left(\phi^i\right)^2 \mathbb{E}\left[\epsilon_{t-i}^2\right] \\[.5em] -->
<!-- &= \sigma^2\sum_{i=0}^{t - 1} \left(\phi^2\right)^i \\[.5em] -->
<!-- &= \sigma^2 \frac{1}{1 - \phi^2} \enspace , -->
<!-- \end{aligned} -->
<!-- $$ -->
 
<!-- where the last line follows when $N \rightarrow \infty$ for $\vert\phi\vert < 1$ from a geometric series. For a random walk, however, this is not a geometric series anymore; it therefore does not converge, and the variance of a random walk does not exist. -->
 
### Spurious correlation of AR(1) processes
In the main text, we have looked at how the spurious correlation behaves for a random walk. Here, we study how the spurious correlation behaves as a function of $\phi \in [0, 1]$. We focus on sample sizes of $n = 200$, and adapt the simulation code from above.
 

{% highlight r %}
set.seed(1)
 
n <- 200
times <- 100
phis <- seq(0, 1, .02)
comb <- expand.grid(times = seq(times), n = n, phis)
ncomb <- nrow(comb)
 
res <- matrix(NA, nrow = ncomb, ncol = 6)
colnames(res) <- c('ix', 'n', 'phi', 'cor', 'tstat', 'pval')
 
for (i in seq(ncomb)) {
  ix <- comb[i, 1]
  n <- comb[i, 2]
  phi <- comb[i, 3]
  
  test <- cor.test(simulate_ar(n, phi = phi), simulate_ar(n, phi = phi))
  res[i, ] <- c(ix, n, phi, test$estimate, test$statistic, test$p.value)
}
 
dat <- data.frame(res) %>% 
  group_by(phi) %>% 
  summarize(
    avg_abs_corr = mean(abs(cor)),
    avg_abs_tstat = mean(abs(tstat)),
    percent_sig = mean(pval < .05)
  )
{% endhighlight %}
 
The Figure below shows that the issue of spurious correlation gets progressively worse as the AR(1) process approaches a random walk (i.e., $\phi = 1$). While this is true, the regression estimate remains consistent.
 
<img src="/assets/img/2019-06-23-Spurious-Correlation.Rmd/unnamed-chunk-11-1.png" title="plot of chunk unnamed-chunk-11" alt="plot of chunk unnamed-chunk-11" style="display: block; margin: auto;" />
 
 
## References
- Granger, C. W., & Newbold, P. ([1974](http://wolfweb.unr.edu/~zal/STAT758/Granger_Newbold_1974.pdf)). Spurious regressions in econometrics. *Journal of Econometrics, 2*(2), 111-120.
- Hamilton, J. D. ([1994](https://press.princeton.edu/titles/5386.html)). Time Series Analysis. P. Princeton, US: Princeton University Press.
- Kuiper, R. M., & Ryan, O. ([2018](https://www.tandfonline.com/doi/full/10.1080/10705511.2018.1431046)). Drawing conclusions from cross-lagged relationships: Re-considering the role of the time-interval. *Structural Equation Modeling: A Multidisciplinary Journal, 25*(5), 809-823.
- Phillips, P. C. ([1986](http://dido.econ.yale.edu/korora/phillips/pubs/art/a044.pdf)). Understanding spurious regressions in econometrics. *Journal of Econometrics, 33*(3), 311-340.
- Matthews, R. Storks deliver babies (p = 0.008) ([2000](https://onlinelibrary.wiley.com/doi/abs/10.1111/1467-9639.00013?casa_token=cWUllTD9P14AAAAA:PRERZz-uS2z9xX3DGt0-Qize94FuZuw-35s-2ECfUDY9Oi3J1m83cZh8EBHGlGh7fwQ2WHShOQuwB-YO)). *Teaching Statistics 22*(2), 36–38.
 
## Footnotes
[^1]: There are, of course, many [more](https://www.tylervigen.com/spurious-correlations).
[^2]: Thanks to Toni Pichler for drawing my attention to the fact that independent random walks are correlated, and Andrea Bacilieri for providing me with the classic references.
[^3]: For practical data analyses, this means that we need to sample data at equal time spaces. Moreover, the estimated coefficients are only valid with respect to that time inveral; for details, see Kuiper & Ryan (2018).
[^4]: Moreover, one way to avoid the spurious correlation is to *difference* the time-series. For other approaches, see Hamilton (1994, pp. 561).
[^5]: This awesome picture was made by Luther Seet.
