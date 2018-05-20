#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

int64_t* factor(int64_t n, int64_t *count) {
	int64_t p, cnt = 0;
	int64_t *ret = (int64_t *)malloc(sizeof(int64_t)*(sqrt(n)+1));
	for (p = 2; p * p <= n; ++p) {
        while (n % p == 0) {
            ret[cnt++] = p;
            n /= p;
        }
	}
	if (n > 1) ret[cnt++] = n;
	if (count) *count = cnt;
	return ret;
}

int reciproc(int64_t q, int64_t p) {
    int64_t cnt, *m;
    int i, ans = 1;
    while (q > p)
        q -= p;
    m = factor(q, &cnt);
    if (cnt > 1)
        for (i = 0; i < cnt; ++i)
            ans *= reciproc(m[i], p);
    else if (q == 2)
        ans *= ((p * p - 1) / 8) & 1 ? -1 : 1;
    else if (q != 1) {
        ans *= ((p - 1) / 2 & 1) * ((q - 1) / 2 & 1) & 1 ? -1 : 1;
        ans *= reciproc(p, q);
    }
    free(m);
    return ans;
}

void extgcd(int64_t a, int64_t b, int64_t *d, int64_t *x, int64_t *y) {
    if (!b)
        *d = a, *x = 1, *y = 0;
    else {
        extgcd(b, a % b, d, y, x);
        *y -= *x * (a / b);
    }
}

int64_t inverse(int64_t a, int64_t p) {
    int64_t d, x, y;
    extgcd(a, p, &d, &x, &y);
    return d == 1 ? (x + p) % p : -1;
}

int legendre(int64_t a, int64_t p) {
    int64_t cnt, i, ret = 1;
    int64_t *m = factor(a, &cnt);
    for (i = 0; i < cnt; ++i)
        ret *= reciproc(m[i], p);
    free(m);
    return ret;
}

int64_t expmod(int64_t a, int64_t x, int64_t m) {
    int64_t sum = 1;
    while (x) {
        if (x & 1) {
            sum = (sum * a) % m;
            --x;
        }
        x /= 2;
        a = a * a % m;
    }
    return sum;
}

int64_t fastpow(int64_t a, int64_t n) {
    int64_t ans = 1, b = a;
    while (n) {
        if (n & 1)
            ans *= b;
        b *= b;
        n /= 2;
    }
    return ans;
}

int64_t solve(int64_t a, int64_t p, int64_t *result) {
	int64_t s = p - 1, t = 0, n, b, x, a1, c;
	int j;
	if (legendre(a, p) < 0)
		return 0;
	while (!(s & 1)) {
		s /= 2;
		t += 1;
	}
	c = t - 1;
	for (n = 1; legendre(n, p) == 1; ++n);
	b = expmod(n, s, p);
	a1 = inverse(a, p);
	x = expmod(a, (s + 1) / 2, p);
	printf("n = %lld, s = %lld, b = %lld\n", n, s, b);
    printf("a^-1 = %lld, b = %lld, x = %lld\n", a1, b, x);
	while (--t) {
        j = expmod(a1 * x * x, 1 << (t - 1), p) == 1 ? 0 : 1;
        x = (x * fastpow(b, j * (1 << (c - t)))) % p;
        printf("j = %d, x = %lld\n", j, x);
	}
	*result = x;
	return 1;
}

int main() {
    int64_t a = 315, p = 907, ans;
	if (solve(a, p, &ans))
        printf("Solution is: %lld, %lld.\n", ans, p - ans);
    else printf("No solution.\n");
	return 0;
}
