#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

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

int64_t gcd(int64_t a, int64_t b) {
    if (!a) return b;
    else return gcd(b % a, a);
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

int64_t gm_encrypt(char m, int64_t p, int64_t q) {
    int64_t n = p * q;
    int64_t x, y = n, c;
    for (x = 3; legendre(x, p) != -1 || legendre(x, q) != -1; ++x);
    srand((unsigned) time(NULL));
    while (gcd(y, n) != 1)
        y = rand() % 65535;
    c = y * y * (m ? x : 1);
    printf("x: %lld, y: %lld, c: %lld\n", x, y, c);
    return c;
}

char gm_decrypt(int64_t c, int64_t p, int64_t q) {
    int a, b;
    a = legendre(c, p);
    b = legendre(c, q);
    printf("legendre(c, p): %d, legendre(c, q): %d\n", a, b);
    if (a == 1 && b ==1)
        return 0;
    else return 1;
}


int main() {
    int64_t p = 613, q = 827;
    char d, m = 0;
    printf("message: %d\n", m);
    int64_t c = gm_encrypt(m, p, q);
    d = gm_decrypt(c, p, q);
    printf("decrypted: %d", d);
	return 0;
}
