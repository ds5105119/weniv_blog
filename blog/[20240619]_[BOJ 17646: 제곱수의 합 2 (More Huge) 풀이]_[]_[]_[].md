# BOJ 17646: 제곱수의 합 2 (More Huge) 풀이

<aside>
📒 알고리즘 분류: Miller–Rabin primality test(밀러–라빈 소수 판별법), 폴라드 로(Pollard's rho), 수학

</aside>

정수론 관련 문제이다. 생각보다는 간단한 방법으로 풀리는 문제이다. 훨씬 어려울줄 알고 논문을 읽었었는데, 배껴쓴것같은 논문도 찾게되었었다.

[Miller-Rabin 소수 판별법](https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test)과, [폴라드 로 알고리즘](https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm)을 사용하기 때문에 두 알고리즘은 위키피디아를 참조하면 된다.

### 두 개의 제곱수로 분리하는 알고리즘

---

페르마의 두 제곱수 정리(🔗[위키피디아 링크](https://en.wikipedia.org/wiki/Fermat%27s_theorem_on_sums_of_two_squares))는 홀수 소수 $p$에 대해 $p ≡ 1 \pmod 4$와 $p = x^2 + y^2$인 정수 $x, y$가 존재한다는 조건은 동치라는 것이다.

알고리즘은 2pass로 진행된다.

1. quadratic non-residue $c \mod p$를 찾는다.
2. $x = c^{(p-1)/4} \mod p$라 하자
3. $p, x$에 대해 gcd를 수행하고 두 나머지 $s, r$에 대해 $s, r < \sqrt p$이면 멈춘다. 이때 $r^2 + s^2 = p$를 만족한다.

$x$를 찾기 위해서 for 문을 사용하여 2부터 $\sqrt p$까지의 자연수 $c$에 대해 $c^{(p-1)/4} \mod p$를 계산하고 이것의 제곱을 $p$로 나눈 나머지가 -1인지 확인하면 된다.

$p ≡ 1 \pmod 4$에 대해, $p = a^2 + b^2 = (a + bi)(a - bi)$로 쓸 수 있다. 

다른 두 제곱수로 분해될 수 있는 소수 $q$에 대해, $q = (c + di)(c - di)$라고 한다면 $pq$는 다음과 같다.

$$
pq = ((ac - bd) + i(ad + bc))(ac-bd)-i(ad+bc))
$$

$$
\therefore pq = (ac-bd)^2 +(ad+bc)^2
$$

따라서, 정수 $n$의 소인수분해된 결과에 $p ≡ 3 \pmod 4$이고, $k$가 홀수인 $p^k$가 존재하지 않는다면, $n$은 두 제곱수로 분해될 수 있다.

```python
def div_prime(p):
  """
  :param p: 4k + 1 꼴로 나타내지는 소수와 2
  :return: p = a^2 + b^2 일 때, [a, b]
  """
  if p == 2:
      return 1, 1

  n, x, k, sq_value = p, 0, (p - 1) // 4, []

  for j in range(2, int(math.sqrt(p)) + 1):
      x = pow(j, k, p)
      if pow(x, 2, p) == p - 1:
          break

  while True:
      n, x = x, n % x
      if n * n < p:
          sq_value.append(n)
      if len(sq_value) == 2:
          break

  return sq_value
```

### 문제 솔루션

---

라그랑주의 네 제곱수 정리에 따르면 $n, k \in \mathbb{Z}^+ \cup \{0\}$에 대하여 $4^n(8k + 7)$ 은 세 제곱수로 나타낼 수 없다.

지금까지의 사실을 이용해 문제의 입력인 자연수 $n$에 대해

1. 한 개의 제곱수로 분리되는 경우 제곱근 함수로 체크
2. 두 개의 제곱수로 분리되는 경우 두 개의 제곱수로 분리
3. 자연수 $n$이 3개의 제곱수로 분리되는 경우 $n$에서 1, 4, 9, 16씩 제곱수를 뺀 값이 두 개의 제곱수로 분리되는지 확인
4. 자연수 $n$이 4개의 제곱수로 $4^k$를 n에서 빼고, 위 1, 2, 3 과정을 사용한다. 

### 전체 코드

---

```python
import sys
import math
import random
from collections import Counter
sys.setrecursionlimit(10 ** 6)

PRIME = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

def miller_rabin(n, p):
    d = n - 1
    while d % 2 == 0:
        if pow(p, d, n) == n - 1:
            return 1
        d >>= 1
    d = pow(p, d, n)
    return d == n - 1 or d == 1

def is_prime(n):
    """
    Miller Rabin Algorithm을 사용하여 소수를 판별한다.
    :param n: 자연수
    :return: 소수인 경우 True, 합성수인 경우 False를 반환
    """
    if n in PRIME:
        return True
    if n % 2 == 0 or n == 1:
        return False

    for a in PRIME:
        if not miller_rabin(n, a):
            return False
    return True

def pollard_rho(n):
    """
    :param n: 자연수
    :return: n의 소수인 약수 하나를 반환
    """
    if n & 1 == 0:
        return 2
    if is_prime(n):
        return n

    x = y = random.randint(2, n)
    c = random.randint(1, n - 2)
    d = 1

    while d == 1:
        x = (x ** 2 + c) % n
        y = (y ** 2 + c) % n
        y = (y ** 2 + c) % n
        d = math.gcd(x - y, n)
        if d == n:
            return pollard_rho(n)

    if is_prime(d):
        return d
    else:
        return pollard_rho(d)

def factorization(n):
    """
    :param n: 자연수
    :return: 소인수분해된 리스트를 반환
    """
    if n == 1:
        return [1]
    factor = []
    while n != 1:
        d = pollard_rho(n)
        factor.append(d)
        n //= d
    return factor

def div_prime(p):
    """
    :param p: 4k + 1 꼴로 나타내지는 소수와 2
    :return: p = a^2 + b^2 일 때, [a, b]
    """
    if p == 2:
        return 1, 1

    n, x, k, sq_value = p, 0, (p - 1) // 4, []

    for j in range(2, int(math.sqrt(p)) + 1):
        x = pow(j, k, p)
        if pow(x, 2, p) == p - 1:
            break

    while True:
        n, x = x, n % x
        if n * n < p:
            sq_value.append(n)
        if len(sq_value) == 2:
            break

    return sq_value

def sq2(n):
    """
    :param n: 2개로 분해될 수 있는 자연수
    :return: 2개로 분해된 값
    """
    factors = factorization(n)
    factors_count = Counter(factors)
    a, b = 1, 0

    for p, e in factors_count.items():
        if p % 4 == 3 and e % 2 == 1:
            return False
        elif (p % 4 == 1 or p == 2) and e % 2 == 1:
            n //= p
            c, d = div_prime(p)
            a, b = a * c - b * d, a * d + b * c

    if b != 0:
        return abs(a) * int(math.sqrt(n)), abs(b) * int(math.sqrt(n))
    else:
        return False

def solve(n):
    # 1개의 제곱수로 표현될 때
    if int(math.sqrt(n)) ** 2 == n:
        return [int(math.sqrt(n))]

    # 3개 이하의 제곱수로 표현될 때, n = 4^a * (8k + 7) 꼴이 아니다
    x, k = n, 0
    while x % 4 == 0:
        x //= 4
        k += 1

    if x % 8 != 7:
        # 2개의 제곱수로 표현될 때
        is_divided_by_2_sums = sq2(n)
        if is_divided_by_2_sums:
            return is_divided_by_2_sums

        # 같은 3개의 제곱수로 표현될 때
        if n % 3 == 0 and int(math.sqrt(n // 3)) ** 2 == n // 3:
            return [int(math.sqrt(n // 3)) for _ in range(3)]

        # 이외 3개의 제곱수로 표현될
        for i in range(1, int(math.sqrt(x))):
            p = sq2(x - (i ** 2))
            if p:
                return (2 ** k) * i, *[(2 ** k) * _ for _ in p]

    # 4개의 제곱수로 표현될 때
    else:
        return 2 ** k, *solve(n - 4 ** k)

def main():
    n = int(sys.stdin.readline())

    a = solve(n)
    a = [i for i in a if i]  # 0 제거
    a.sort()

    print(len(a))
    for i in a:
        print(i)

if __name__ == '__main__':
    main()

```

![Untitled](/img/BOJ 17646: 제곱수의 합 2 (More Huge) 풀이 img 0.png)
