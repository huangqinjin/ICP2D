# Iterative Closest Point Algorithm for 2D Points

## Definition
$$\begin{align}
& p_i, i = 1,\cdots,n, \\
& q_i = sRp_i + t + N(0, \sigma_i^2I), i = 1,\cdots,n, \\
& \bar{p} = \frac{\sum_{i=1}^n\sigma_i^{-2}p_i}{\sum_{i=1}^n\sigma_i^{-2}} \\
& \bar{q} = \frac{\sum_{i=1}^n\sigma_i^{-2}q_i}{\sum_{i=1}^n\sigma_i^{-2}} \\
& \tilde{p}_i = p_i - \bar{p} \Longrightarrow \sum_{i=1}^n \sigma_i^{-2}\tilde{p}_i = 0 \\
& \tilde{q}_i = q_i - \bar{q} \Longrightarrow \sum_{i=1}^n \sigma_i^{-2}\tilde{q}_i = 0 \\
\end{align}$$


## [Objective](http://ceres-solver.org/nnls_covariance.html)
$$\begin{equation}
\min_{s, R, t} \frac12\sum_{i=1}^n\sigma_i^{-2}\|sRp_i + t - q_i\|^2
\end{equation}$$

Since
$$\begin{align*}
  & \|sRp_i + t - q_i\|^2 \\
= & \|(sR\tilde{p}_i - \tilde{q}_i) + (sR\bar{p} + t - \bar{q})\|^2 \\
= & \|sR\tilde{p}_i - \tilde{q}_i\|^2 + \|sR\bar{p} + t - \bar{q}\|^2 + 2(sR\bar{p} + t - \bar{q})^T(sR\tilde{p}_i - \tilde{q}_i),
\end{align*}$$

and
$$\begin{align*}
& \|sR\tilde{p}_i - \tilde{q}_i\|^2 = s^2\|\tilde{p}_i\|^2 + \|\tilde{q}_i\|^2 - 2s\tilde{q}_i^TR\tilde{p}_i, \\
& \sum_{i=1}^n\sigma_i^{-2}(sR\tilde{p}_i - \tilde{q}_i) = 0,
\end{align*}$$

the objective is equivalent to
$$\begin{align}
& \min (as^2 - 2rs) \\
& r = \sum_{i=1}^n\sigma_i^{-2}\tilde{q}_i^TR\tilde{p}_i \\
& a = \sum_{i=1}^n\sigma_i^{-2}\|\tilde{p}_i\|^2 \\
& t = \bar{q} - sR\bar{p}
\end{align}$$

Furthermore, we can prove that $ r_{\max} > 0 $, so the objective becomes
$$\begin{align}
& r_{\max} = \max_{R}\sum_{i=1}^n\sigma_i^{-2}\tilde{q}_i^TR\tilde{p}_i \\
& s = \frac{r_{\max}}{a} > 0
\end{align}$$


## Solution
By properties of trace, we have
$$\begin{equation}
r = \sum_{i=1}^n\mathrm{tr}(\sigma_i^{-2}\tilde{q}_i^TR\tilde{p}_i) = \mathrm{tr}\left(\sum_{i=1}^n\sigma_i^{-2}\tilde{p}_i\tilde{q}_i^TR\right)  = \mathrm{tr}(WR)
\end{equation}$$

where
$$\begin{equation}
W = \sum_{i=1}^n\sigma_i^{-2}\tilde{p}_i\tilde{q}_i^T.
\end{equation}$$

By SVD $ W = USV^T $, we have $ \mathrm{tr}(WR) = \mathrm{tr}(SV^TRU) $ and
$$\begin{equation}
\arg\max_{R}\mathrm{tr}(SV^TRU) = V\begin{pmatrix} I & 0 \\ 0 & \mathrm{sgn}\left|UV\right|\end{pmatrix}U^T.
\end{equation}$$


For 2D case i.e.
$$\begin{align}
W &= \begin{bmatrix}w_{11} & w_{12} \\ w_{21} & w_{22}\end{bmatrix}, \\
R &= \begin{bmatrix}\cos\theta & -\sin\theta \\ \sin\theta & \cos\theta\end{bmatrix}, \\
\end{align}$$

we have
$$\begin{gather}
\mathrm{tr}(WR) = (w_{11} + w_{22})\cos\theta + (w_{12} - w_{21})\sin\theta = k\cos(\theta - \varphi) \\
k = \sqrt{(w_{11} + w_{22})^2 + (w_{12} - w_{21})^2} \\
\varphi = \mathrm{arctan2}(w_{12} - w_{21}, w_{11} + w_{22}) \\
\end{gather}$$


## Iteration
Suppose we have
$$\begin{align}
& \bar{p}_n = \frac{\sum_{i=1}^n\sigma_i^{-2}p_i}{\sum_{i=1}^n\sigma_i^{-2}} \\
& \bar{q}_n = \frac{\sum_{i=1}^n\sigma_i^{-2}q_i}{\sum_{i=1}^n\sigma_i^{-2}} \\
& W_n = \sum_{i=1}^n\sigma_i^{-2}(p_i - \bar{p}_n)(q_i - \bar{q}_n)^T
\end{align}$$

Then
$$\begin{align}
\bar{p}_{n+1} = \frac{\sum_{i=1}^{n+1}\sigma_i^{-2}p_i}{\sum_{i=1}^{n+1}\sigma_i^{-2}} & = \bar{p}_n + \frac{\sigma_{n+1}^{-2}}{\sum_{i=1}^n\sigma_i^{-2}}(p_{n+1} - \bar{p}_{n+1}) \\
& = \bar{p}_n + \frac{\sigma_{n+1}^{-2}}{\sum_{i=1}^{n+1}\sigma_i^{-2}}(p_{n+1} - \bar{p}_n) \\
\bar{q}_{n+1} = \frac{\sum_{i=1}^{n+1}\sigma_i^{-2}q_i}{\sum_{i=1}^{n+1}\sigma_i^{-2}} & = \bar{q}_n + \frac{\sigma_{n+1}^{-2}}{\sum_{i=1}^n\sigma_i^{-2}}(q_{n+1} - \bar{q}_{n+1}) \\ 
& = \bar{q}_n + \frac{\sigma_{n+1}^{-2}}{\sum_{i=1}^{n+1}\sigma_i^{-2}}(q_{n+1} - \bar{q}_n) \\
\end{align}$$


Since
$$\begin{align*}
&(p_i - \bar{p}_{n+1})(q_i - \bar{q}_{n+1})^T \\
=& \left[p_i - \bar{p}_n - \frac{\sigma_{n+1}^{-2}(p_{n+1} - \bar{p}_{n+1})}{\sum_{i=1}^n\sigma_i^{-2}}\right]\left[q_i - \bar{q}_n - \frac{\sigma_{n+1}^{-2}(q_{n+1} - \bar{q}_{n+1})}{\sum_{i=1}^n\sigma_i^{-2}}\right]^T \\
=& (p_i - \bar{p}_n)(q_i - \bar{q}_n)^T + \frac{\sigma_{n+1}^{-4}}{(\sum_{i=1}^n\sigma_i^{-2})^2}(p_{n+1} - \bar{p}_{n+1})(q_{n+1} - \bar{q}_{n+1})^T \\
& - \frac{\sigma_{n+1}^{-2}(p_{n+1} - \bar{p}_{n+1})}{\sum_{i=1}^n\sigma_i^{-2}}(q_i - \bar{q}_n)^T - \frac{\sigma_{n+1}^{-2}(p_i - \bar{p}_n)}{\sum_{i=1}^n\sigma_i^{-2}}(q_{n+1} - \bar{q}_{n+1})^T \\
\end{align*}$$

and
$$
\sum_{i=1}^n\sigma_i^{-2}(p_i - \bar{p}_n) = \sum_{i=1}^n\sigma_i^{-2}(q_i - \bar{q}_n) = 0,
$$

we have
$$\begin{align}
& \sum_{i=1}^n\sigma_i^{-2}(p_i - \bar{p}_{n+1})(q_i - \bar{q}_{n+1})^T  \\
=& W_n + \frac{\sigma_{n+1}^{-4}}{\sum_{i=1}^n\sigma_i^{-2}}(p_{n+1} - \bar{p}_{n+1})(q_{n+1} - \bar{q}_{n+1})^T.
\end{align}$$

Thus
$$\begin{aligned}
 & W_{n+1} = \sum_{i=1}^{n+1}\sigma_i^{-2}(p_i - \bar{p}_{n+1})(q_i - \bar{q}_{n+1})^T  \\
=& \sum_{i=1}^n\sigma_i^{-2}(p_i - \bar{p}_{n+1})(q_i - \bar{q}_{n+1})^T + \sigma_{n+1}^{-2}(p_{n+1} - \bar{p}_{n+1})(q_{n+1} - \bar{q}_{n+1})^T  \\
=& W_n + \frac{\sum_{i=1}^{n+1}\sigma_i^{-2}}{\sum_{i=1}^n\sigma_i^{-2}}\sigma_{n+1}^{-2}(p_{n+1} - \bar{p}_{n+1})(q_{n+1} - \bar{q}_{n+1})^T \\
=& W_n + \frac{\sum_{i=1}^n\sigma_i^{-2}}{\sum_{i=1}^{n+1}\sigma_i^{-2}}\sigma_{n+1}^{-2}(p_{n+1} - \bar{p}_n)(q_{n+1} - \bar{q}_n)^T \\
\end{aligned}$$

