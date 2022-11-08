<!-- Use special library for math -->
<script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>

# Lecture 3


We say that $X$ is a Bernuoulli distributed if:
$$P(X=1)=p$$
$$P(X=0)=1-p:=q$$

<!-- An html time line with 2 parts and text in it -->
<div class="timeline" style="width: 100%; height: 100px; background-color: #ddd;">
  <div class="event" style="width: 50%; height: 100px; background-color: #4CAF50; float: left;">
    <div class="event-content" style="padding: 10px;">
      <h3>Event 1</h3>
      <p>Some text..</p>
    </div>
  </div>
  <div class="event" style="width: 50%; height: 100px; background-color: #2196F3; float: left;">
    <div class="event-content" style="padding: 10px;">
      <h3>Event 2</h3>
      <p>Some other text..</p>
    </div>
  </div>
</div>

Algotithm for Bernoulli distribution:
1. Generate a random number $U$ uniformly distributed on $[0,1]$
2. if u < p then return 1 else return 0

Consider now the eneral case of a random variable $X$ with a probability mass function $P(X=x_i)=p_i$ for $x=1,2,...,n$. We say that $X$ is a discrete random variable with a probability mass function $P(X=x)$ for $x=1,2,...,n$ if:
Divide the span $[0,1]$ into $n$ parts. Then the probability of $X$ being in the $i$th part is $p_i$.
Then if $U$ is generated to be in the $i$th part then $X=x_i$.

The cumulative distribution function of a discrete random variable $X$ is defined as:
$$F_X(t)=P(X\leq t)=\sum_{i:x_t\leq t}^n P(X=x_i)=\sum_{i:x_t\leq t}^n p_i$$
The sum over all $n$ terms is equal to 1.

<!-- ExampleBox green border-->
<div class="examplebox" style="border: 2px solid green; padding: 10px; margin: 10px;"> 

## Example
Let:
$$X = \begin{cases}3 \text{ with probability } 0.15\\ 5 \text{ with probability } 0.25\\ 9 \text{ with probability } 0.6\end{cases}$$
We have to calculate:
$$F_X(t)=p(X\leq t)$$
If this would have been continous we would have used an integral.
But since it is discrete we have to use a sum.
If $t<3$:
$$F_X(t)=P(X\leq t)=0$$
If $3\geq t < 5$:
$$F_X(t)=P(X\leq t)=0.15=0+0.15$$
If $5\geq t < 9$:
$$F_X(t)=P(X\leq t)=0.40=0.25+0.15$$
If $t\geq 9$:
$$F_X(t)=P(X\leq t)=1=0.6+0.25+0.15$$
Or:
$$F_X(t)=P(X\leq t)=\begin{cases}0 \text{ if } t<3\\ 0.15 \text{ if } 3\leq t<5\\ 0.40 \text{ if } 5\leq t<9\\ 1 \text{ if } t\geq 9\end{cases}$$

</div>

How to generate uniform discrete random numbers over 
<!-- Math bracket in markdown -->
$$I=\left\{0,1,2,...,n-1\right\}$$
The probabilty of $X$ being in the $i$th part is $p_i=\frac{1}{n}$.
Then if $U$ is generated to be in the $i$th part then $X=x_i$.
Let 
<!-- Floor operation markdown -->
$$X=\lfloor nU\rfloor$$
Where 
$$U\sim U(0,1)$$
Then $X$ is uniformly distributed over $I$.
What is the probabiltiy that $X=i$?
$$P(X=i)=P(\lfloor nU\rfloor=i)=P(U\leq \frac{i+1}{n})-P(U\leq \frac{i}{n})=\frac{1}{n}$$