+++
date = '2025-03-14T00:38:48+08:00'
draft = false
title = 'StatPhy - Cluster Expansion'
math = true
+++

lecture notes from David Tong

<!--more-->
---
## 统计力学中的团簇展开 —— 图形化的计算方法

## 1 Introduction

**"Statistical mechanics is the art of turning the microscopic laws of physics into a description of Nature on a macroscopic scale." - David Tong**

### 1.1 为什么需要统计力学？

假设你破解了全部的物理学原理，知道了自然界的所有基本定律、基本粒子的特性以及它们之间的相互作用力，你该如何将这些知识转化成对你所处周围世界的理解？更具体地讲，假如有一个包含 $10^{23}$ 个粒子的盒子，已知它们的质量、电荷和相互作用，关于盒子本身，你能告诉我什么信息？

Newton 和 Schrödinger分别赐予了我们经典与量子情况下的粒子动力学方程，但是不论何种，对我们（包括人类与计算机）的有限算力来讲，联立 $10^{23}$ 个方程并尝试去解，是最坏也最不可能的策略。而且就算解出了运动轨迹/波函数，我们仍然无法截然宣称，这个盒子是热的还是冷的，是软的还是硬的？盒子会爆炸吗，还是会在下一秒坍缩成黑洞？**丰富的微观信息并没有给我们带来肯定的宏观答案**，因为真正的宏观系统不是由一个纯量子态描述的。它与环境接触，不断受到外部影响和冲击，每次系统经历一个小的扰动，就有可能转换到另一种状态。我们需要一个将微观、宏观联系起来的理论，一个放之四海而皆准的定则。它将告诉我们，统计不同于掰手指头数数，而是有趣的 LEGO 积木游戏，把零碎的片段，拼成完美的城堡。

### 1.2 统计力学从何而来，往何处去？

从1600年代到1900年代，科学家们不断发现着表现在不同物质上的统一的物理规律，比如我们熟知的 Boyle 定律（$p \propto \frac{1}{V}$）、 Charles 定律（$V \propto T$）和 Curie 定律（$\chi \propto \frac{1}{T}$）等。它将我们所感兴趣的宏观性质联系在一起，也隐含了物理学的一条重要信息——**普适性（Universality）**。不同物质的宏观性质有着相同规律，直觉上可以猜想它们来源于同一微观性质，因为那时一些人们认为，物质在微观上的组成是几乎一致的。经典热力学的成功没有掩盖 Ludwig Boltzmann 对微观、宏观联系的洞察力，他同时也是原子论的坚定捍卫者：“如果原子不存在，那么如何解释热现象？”基于原子论的观点， Boltzmann 发展了统计方法，奠定了现代统计力学的基础。

我国统计物理的先驱，王竹溪先生说:“统计物理学既是物理学，也是方法论。”从个例总结出一般规律、再去研究个例，是一切自然科学的方法基础。同时，现代统计物理不仅仅包含 Boltzmann 的统计力学，它并不偏倚还原论或是整体论，而是从**对称性（Symmetry）**和**尺度（Scaling）**出发，考察不同尺度的现象和本质，为研究相变和临界现象提供手段。例如，我们注意到，小的事物往往会影响大的事物，但很少影响非常大的事物；沿着关系链爬升，我们似乎会遗忘很小尺度上的事情。研究椋鸟成群飞行的动物学家无需理解 Higgs 玻色子的动力学，每个尺度各有分工，重整化群（Renormalisation Group）就是这样作为一个 formalism 的理论被提出。这不是本讲义讨论的主要内容，但会给我们带来思想上的启发。

## 2 统计物理的基本假设

### 2.1 基本概念和假设
统计物理通过研究大量微观粒子的统计行为，来解释宏观物理现象。为了描述这种联系，我们需要事先引入两个关键概念和两大基本假设：

- 微观态（Microstate）：微观态 $\omega_i\in\Omega=\{\omega_1,\omega_2,\dots\}$ 是对一个物理系统在**某一时刻**所有**微观自由度**（如粒子的位置、动量、量子态等）的**完整描述**。在经典力学中是相空间 $\{(\boldsymbol q,\boldsymbol p )\}$ 中的一个6N维向量 $(q_1,q_2,\dots,q_N;p_1,p_2,\dots,p_N)$ ；在量子力学中是系统的态矢量 $\vert \psi \rangle$ 或波函数 $\psi$。
- 系综（Ensemble）：系统在**一定宏观约束条件**下能够实现的**微观态**（也称**可访问态（accessible state）**）构成的一个**集合**。可以分为

- 微正则系综 $$\varepsilon_{micro}=\{\omega_i\in\Omega\vert E(\omega_i)=E, V(\omega_i)=V,N(\omega_i)=N\}$$
- 正则系综 $$\varepsilon_{canonical} =\{\omega_i\in\Omega\vert V(\omega_i)=V,N(\omega_i)=N\}$$
- 巨正则系综 $$\varepsilon_{grand} =\{\omega_i\in\Omega\vert V(\omega_i)=V\}$$

- 等概率（Equal-probability）假设：对处于平衡态（宏观量不再变化）的孤立系统，所有可访问的微观态的可能性都相同。
- 遍历（Ergodic）假设：在足够长的时间内，一个系统会以相同的机会访问其状态空间中的所有微观态。这意味着系统的时间平均与系综平均等价。

例如，一个具有2个可分辨粒子的系统，限定每个粒子的自旋只能向上或向下（1或-1），总共有$2^2$个可能的微观态，即(-1,-1),(-1,1),(1,-1),(1,1)，能量 $E=J\cdot S$ 分别为\{-2,0,0,2\}，假设访问每个微观态所用的时间为 $\Delta t$ ，能量的时间平均就是 $\langle E \rangle_{\rm{time}}=\frac{-2\cdot\Delta t+0+0+2\cdot\Delta t}{4\cdot\Delta t}=0$ ；而系综平均是 $\langle E \rangle_{\rm{ensemble}}=\frac{-2+0+0+2}{2^2}=0$ 。

这些假设奠定了从微观描述推导宏观性质的基础。

### 2.2 系统的可观测量

在一般的物理系统中，可观测量由算符 $\hat{\mathcal{O}}$ 表示，有 $\langle \hat{\mathcal{O}}\rangle = \sum_n p(n) \langle n \vert \hat{\mathcal{O}}\vert n \rangle $。其中 $p(n)$ 是第 $n$ 个微观态对应的概率，$\langle n \vert \hat{\mathcal{O}}\vert n \rangle$ 是在该状态下观测到的值。因此，要从宏观上理解 $\hat{\mathcal{O}}$ 的物理性质，我们必须搞清楚**微观态的分布 $p(n)$**。

### 2.3 正则系综的微观态分布

正则系综，即系统能够与一个大的热库（hot reservoir）达到平衡的情况，其粒子数N,体积V,温度T固定。给定 NVT 的值，我们可以唯一确定微观态的分布情况 $P_i=P(E_i)$。这里通过拉格朗日乘子法和最大熵原理推导它的形式，由约束条件：
- 归一化条件: $$\sum_i p_i = 1 $$
- 宏观能量不再变化\$$\sum_i p_i E_i = \langle E \rangle$$

以及 Gibbs 熵的定义 $S=-k_B\sum_ip_i\ln p_i$。构造拉格朗日函数：
$$\mathcal{L}= -\sum_i p_i \ln p_i -\alpha(\sum_i p_i-1)-\beta(\sum_i p_i E_i -\langle E\rangle)$$
要使得熵$\frac{S}{k_B}=\mathcal{L}$最大，则对每个 $p_i$ 求偏导且令其为0：
$$\frac{\partial \mathcal{L} }{\partial p_i} = -(1+\ln p_i)-\alpha-\beta E_i =0$$
解得$$p_i=\frac{e^{-\beta E_i}}{e^{\alpha+1}}\propto e^{-\beta E_i}$$
这就是 Boltzmann 分布，定义其中 $\beta=\frac{1}{k_B T}$。对概率进行归一化，我们得到
$$P_i = \frac{e^{-\beta E_i}}{\sum_i e^{-\beta E_i}} = \frac{e^{-\beta E_i}}{Z}$$
$Z$ 称为**配分函数（Partition function）**。当 $E_i=E\text{, }i=1,2,\dots,n$ 时，$P_i=\frac{1}{n}$，所有微观态的概率相同，退化为微正则系综的微观态分布。考虑粒子数可变的情况，则在 $E$ 上加一项 $-\mu N_i$ 表示将粒子放入系统后的能量改变量，得到巨正则系综的微观态分布 $P_i=\frac{e^{-\beta (E_i-\mu N_i)}}{\mathcal{Z}}$。

## 3 配分函数（Partition function）
配分函数不仅仅是一个归一化因子，而包含了物理系统的全部热力学信息。
$$\begin{align*}
   Z&=\sum_i e^{-\beta E_i}\\
   \mathcal{Z}&=\sum_{i} e^{\beta \mu N_i}e^{-\beta E_i}=\sum_{i}z^{N_i} e^{\beta E_i}
\end{align*}$$
对于已经写出配分函数的系统，可以用 **$\beta$ trick** 来获得一系列宏观可观测量。例如

$$\begin{align*}
   \langle E\rangle&=\sum_i P_i E_i=\frac{E_i e^{-\beta E_i}}{Z}= -\frac{d\ln Z}{d\beta}\\
    \langle N\rangle&=\sum_i P_i N_i=\frac{N_i e^{-\beta (E_i-\mu N_i)}}{\mathcal{Z}}= \frac{1}{\mu}\frac{d\ln \mathcal{Z}}{d\beta}
\end{align*}$$

### 3.1 From Quantum to Classical
以上介绍的是量子化的配分函数，系统只能取分立的能级，但我们可以相信的是：量子世界永不落后。我们一定可以从量子的情况推导经典系统中的连续配分函数
$$\begin{align}
    Z_1=\frac{1}{(2\pi\hbar)^3}\int d^3qd^3p e^{-\beta H(p,q)}
\end{align}$$
首先，系统的哈密顿量由 $\hat{H}=\frac{\hat{p}^2}{2m}+V(\hat{q})$ 给出，写出其配分函数
$$Z_1=\sum_n e^{-\beta E_n}=\sum_n\langle n\vert e^{-\beta \hat{H}} \vert n \rangle$$
插入两个恒等元 $\boldsymbol{1} = \int dq \vert q \rangle \langle q \vert$ 和 $\boldsymbol{1} = \int dp \vert p \rangle \langle p \vert$，得到
$$\begin{align*}
    Z_1 &=\sum_n \langle n \vert  \int dq \vert q \rangle \langle q \vert e^{-\beta \hat{H}}  \int dq' \vert q' \rangle \langle q' \vert n \rangle \\
    &= \int dq dq' \langle q \vert e^{-\beta \hat{H}} \vert q'\rangle \sum_n \langle q' \vert n \rangle \langle n\vert q \rangle
\end{align*}$$
其中 $\sum_n \vert n \rangle \langle n \vert = \boldsymbol{1}$也是恒等元，$\langle q'\vert q\rangle =\delta(q'-q)$ 是 $\delta$ 函数。原式可以化简成为

$$\begin{equation}
    Z_1=\int dq\langle q \vert e^{-\beta \hat{H}}\vert q \rangle
\end{equation}$$

到此为止我们完成了配分函数“对能量本征态求和”到“对位置本征态求积分”的转变。事实上，我们可以选择任何想要的本征态作为基底，也可以写成独立于基底的形式：$Z_1 = \rm{Tr }\text{ } e^{-\beta \hat{H}}$。

那么该如何向经典世界靠拢呢？令 $[\hat{q},\hat{p}]=i\hbar \rightarrow 0$，就得到指数项的分离形式 $$ e^{-\beta \hat{H}}=e^{-\beta \hat{p}^2/2m} e^{-\beta V(\hat{q})} +\mathcal{O}(\hbar)$$
再代入式(2)并忽略 $\mathcal{O}(\hbar)$ 项，注意这里符号上尖的细微差别

$$\begin{align*}
    Z_1&=\int dq\langle q \vert e^{-\beta \hat{p}^2/2m} e^{-\beta V(\hat{q})}\vert q \rangle\\
    &=\int dq e^{-\beta V(q)}\langle q\vert e^{-\beta \hat{p}^2/2m}\vert q \rangle\\
    &=\int dq dp dp' e^{-\beta V(q)}\langle q\vert p\rangle\langle p \vert e^{-\beta \hat{p}^2/2m}\vert p'\rangle \langle p' \vert q\rangle\\
    &=\frac{1}{2\pi \hbar}\int dq dp e^{-\beta H(p,q)}
\end{align*}$$

倒数第二步仍然是恒等元的插入，最后一步则利用了Fourier变换等式 $\langle q\vert p\rangle=\frac{1}{2\pi \hbar}e^{ipq/\hbar}$。

这是一维形式的结果，很容易推广至三维，也就是式 (1)。式中包含的常数大小并不重要，只要它的单位恰好能满足 $Z$ 是无量纲常数。而写成 $\hbar$ 可以提醒我们，这个式子是由量子情况过渡过来的。拥有了连续版本的配分函数，我们就可以愉快地处理接下来的相互作用气体（Interacting Gas）了。

## 4 团簇展开（Cluster Expansion）
### 4.1 真实气体的配分函数
对于理想气体（无相互作用），哈密顿量只包含动能项 $H=\sum_{i=1}^N \frac{p_i^2}{2m}$，代入配分函数，很容易积分得到
$$\begin{align*}
    Z_1&=\frac{1}{(2\pi\hbar)^3}\int d^3q d^3p e^{-\beta p^2/2m} \\
    &=V(\frac{mk_BT}{2\pi \hbar^2})^{3/2}\\
    Z_N&= \frac{Z_1^N}{N!} =\frac{V^N}{N!\lambda^{3N}}
\end{align*}$$
其中 $\lambda=(\frac{mk_BT}{2\pi \hbar^2})^{-1/2}$ 称为热波长。对于有相互作用 $U(r_{ij})$ 的真实气体，哈密顿量是
$$H=\sum_{i=1}^N \frac{p_i^2}{2m}+\sum_{i>j} U(r_{ij})$$
代入配分函数积分，可以写为
$$\begin{align*}
    Z_N &= \frac{1}{N!}\frac{1}{(2\pi \hbar)^{3N}} \int \prod_{i=1}^N d^3 p_i e^{-\beta \sum_j p_j^2/2m}\int\prod_{i=1}^{N} d^3 r_i e^{-\beta \sum_{j< k}{ U(r_{jk})}}\\
    &=\frac{1}{N!}\frac{1}{\lambda^{3N}}\int \prod_{i=1}^N d^3r_i e^{-\beta \sum_{j< k}{ U(r_{jk})}}
\end{align*}$$
相比理想气体即 $U(r_{ij})=0$ 的情况，我们无法直接从积分项中提取出有如体积 $V =\int d^3r$ 这样显而易见的信息。也许可以寄希望于 Taylor 展开：
$$\begin{align*}
    e^{-\beta \sum_{j< k}{ U(r_{jk})}}=1-\beta \sum_{j<k} U(r_{jk})+\frac{\beta^2}{2}\sum_{j<k,l<m} U(r_{jk})U(r_{lm})+\cdots
\end{align*}$$
在 $\beta U\ll 1$ 时，这自然是一个很好的展开，也对应于高温或近理想气体的行为。可是当我们考察偏离理想的情况，例如 $r\rightarrow 0\text{ 时，}U(r_{ij})\rightarrow \infty$，展开式将发散，高阶项将越来越大。在此引入 **Mayer f-function**：
$$f(r)=e^{-\beta U(r)}-1$$
可以发现 $r\rightarrow \infty \text{ 时，}f(r)\rightarrow 0\text{，}r \rightarrow 0\text{ 时，}f(r)\rightarrow -1$，$f_{ij}$是一个能使展开式收敛的好参量，我们计算时可以随时忽略高阶项。同时为了方便我们定义 $f_{ij}=f(r_{ij})$，从而配分函数可以写为
$$
\begin{equation}
\begin{split}
    \notag Z_N &=\frac{1}{N!}\frac{1}{\lambda^{3N}}\int \prod_{i=1}^N d^3r_i \prod_{j<k}(1+f_{jk})\\
    &=\frac{1}{N!}\frac{1}{\lambda^{3N}}\int \prod_{i=1}^N d^3r_i (1+\sum_{j<k}f_{jk} +\sum_{j<k,l<m} f_{jk}f_{lm}+\cdots) 
\end{split}
\end{equation}$$
可以看到，第一项 $\boldsymbol{1}$ 的积分后即 $\int \prod_{i=1}^N d^3r_i=V^{N}$，结果是理想气体配分函数 $Z_{ideal}=\frac{V^N}{N!\lambda^{3N}}$。形如第二项 $\boldsymbol{f_{jk}}$ 的积分是
$$\int \prod_{i=1}^N d^3r_i f_{12}=V^{N-2}\int d^3r_1 d^3 r_2 f(r_{12})=V^{N-1}\int d^3r f(r)$$
最后的等号作了 $r=r_1-r_2$ 的代换。总共有 $N(N-1)/2\approx \frac{N^2}{2}$个这样的项求和，结果是
$$Z^{(1)} = \frac{V^{N}}{N!\lambda^{3N}} \frac{N^2}{2V}\int d^3r f(r)=Z_{ideal}\frac{N^2}{2V}\int d^3r f(r)$$
同理可以继续往下算下去，不过计算高阶项的时候，有趣的事情发生了：我们可以用图形化的计算方法，来让计算变得更加简单，形式更加优美。
### 4.2 图形化的计算方法
我们注意到式(3)的表达式有一个特点，那就是形如  $f_{ij}f_{kl}f_{mn}\dots$ 的组合项，在整个求和中具有**唯一性**，且不管两两之间是否存在单个指标相同，整个项绝不会包含两个指标完全相同的因子，即**双指标的不重复性**。我们可以用**图**来定义这种规则：
- 在纸上画出 $N$ 个原子（请放心，我们不用真的画出 $10^{23}$ 数量级别的原子，在之后我们会将它们分解为一系列小的子图）。
- 用线段连接项 $f_{ij}f_{kl}f_{mn}\dots$ 中的指标对应的原子，即连接 $(i,j), (k,l),(m,n)\dots$等。

<!-- \begin{figure}[H]
    \centering
    \includegraphics[width=1.0\linewidth]{1.png}
    \caption{图中显示了 $N=4$ 的几个子图，分别代表不同的组合项}
    \label{wqwqqwqwq}
\end{figure} -->

这样的规则之下，我们已经得到的配分函数式(3)可以看作是对所有的子图求和然后积分，或是对所有的子图积分然后求和。我们把对子图 $G$ 的积分写作 $W[G]$ ，那么式(3)可以重新写成：
$$\begin{equation}
    Z_N =\frac{1}{N!\lambda^{3N}}\sum_{G} W[G]
\end{equation}$$

<!-- \begin{figure}[H]
    \centering
    \includegraphics[width=1\linewidth]{2.png}
    \caption{$N=5$ 的子图积分，可以分解成一个 3-cluster 和 2-cluster 的乘积}
    \label{fig:enter-label}
\end{figure} -->

绝大部分子图具有不连通的部分，如图(3)的 $f_{12}f_{34}$ 组合项，这样的项显然可以分解为更小的子图积分的乘积。我们把各个内部连通的部分叫做**团簇（clusters）**，具有 $l$ 个原子的团簇称为 $l-cluster$。$N$ 个原子的子图 $G$ 可以被分为 $m_l$ 个 $l-clusters$，即
$$\begin{equation}
    \sum_l^N m_l l=N
\end{equation}$$
有了 clusters 的存在，我们就有办法去处理更大规模的图。先考虑 3-cluster：
<!-- \begin{figure}[H]
    \centering
    \includegraphics[width=1.\linewidth]{3.png}
    \caption{3-cluster，总共包含四种不同的连通子图}
    \label{3cl}
\end{figure} -->
每个 3-cluster 会作为一类乘积因子 $U_3$ 出现在配分函数的各个组合项中，不管与它们组合的 $N-3$ 个原子构成的子图如何。并且由于组合项的唯一性，每个 cluster 只被计算一次：
<!-- \begin{figure}[H]
    \centering
    \includegraphics[width=0.8\linewidth]{4.png}
    \label{U3}
\end{figure} -->

一般地，对 $l-cluster$ 也可以写出乘积因子 $$U_l=\int \prod_{i=1}^l d^3r_i \sum_{G\in\{l-cluster\}} G$$
那么，整个配分函数就可以写成一系列 $U_l$ 的乘积，其中 $\{m_l\}$ 是满足约束关系式(5)的各种可能的组合。
$$\sum_G W[G]=N!\sum_{\{m_l\}} \prod_l \frac{U_l^{m_l}}{(l!)^{m_l}m_l!}$$

其中 $N! \prod_l \frac{1}{(l!)^{m_l}m_l!}$ 项表示将子图 $G$ 分割为一系列 $m_l$ 个 $l-clusters$ 的组合数。

<!-- \begin{figure}[H]
    \centering
    \includegraphics[width=1\linewidth]{5.png}
    \caption{可以验证，对于 $N=5, m_2=1,m_3=1$ 的情况，只包含 $U_3 U_2$ 项，并且系数应当是 $\frac{5!}{3!2!}=10$；从另一方面理解，我们只需要决定 5 个原子中哪两个位于右侧的 2-cluster，这样的组合数恰好为10}
    \label{fig:enter-labeqwqql}
\end{figure} -->

因此，最终的配分函数形式是
$$\begin{equation}
    Z_N=\frac{1}{\lambda^{3N}}\sum_{\{m_l\}}\prod_l  \frac{U_l^{m_l}}{(l!)^{m_l}m_l!}
\end{equation}$$

但是我们还不能满足的一点是，这样写出来的配分函数仍要考虑约束条件(5)下的各种可能性，这是很麻烦的，如果我们没有约束该多好！为了获得粒子数不受约束的情况，我们要从正则系综跳出，前往巨正则系综（可以证明，在热力学极限 $ N\rightarrow \infty,V\rightarrow \infty,N/V\rightarrow const.$ 的情况下，正则系综和巨正则系综等价）
$$\mathcal{Z}=\sum_{N}e^{\beta \mu N}Z_N$$
定义逸度（fugacity）$z=e^{\beta \mu}$。然后我们就可以进一步化简表达式(6)
$$\mathcal{Z}=\sum_{N}z^{N}Z_N=\sum_{m_l=0}^\infty \prod_{l=1}^\infty (\frac{z}{\lambda^3})^{m_l l} \frac{1}{m_l !} (\frac{U_l}{l!})^{m_l} =\prod_{l=1}^{\infty} \exp(\frac{U_l z^l}{\lambda^{3l}l!})$$

继续定义 $b_l=\frac{\lambda^3}{V}\frac{U_l}{l!\lambda^{3l}}$ ，则配分函数变为
$$\mathcal{Z}=\prod_{l=1}^{\infty} \exp(\frac{V}{\lambda^3}b_l z^l)=\exp(\frac{V}{\lambda^3}\sum_{l=1}^\infty b_l z^l)$$

从而，对所有子图可能的分割求和，被我们重新写成了对所有团簇的求和，这就是**团簇展开（cluster expansion）**。图形化的计算是统计力学和量子场论中的一种相当自然的 **trick**，两者都特别关注**连通图**，而在后者的领域中被称为**费曼图（Feymann diagram）**。

在 Cluster Expansion 中，只有连通簇的贡献才不会重复计算，直接影响宏观物理量；而在 Feynman 图中，连通图对应于物理上不可分割的过程。虽然它们应用的具体领域和形式略有不同，这两种图形方法都为解析求解和数值模拟提供了直观且强有力的工具，使得复杂系统的多体相互作用能够被系统地整理与求和，有助于我们识别、理解各阶扰动贡献。
