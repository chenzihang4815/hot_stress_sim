/* Mixed integer/FP/branch/memory kernel. Not a bigger matmul:
   McPAT area is hypothesized to follow CPU config, not opcode mix. */
#include <stdio.h>
#include <stdint.h>

#define N 4096
#define ROUNDS 40

static uint32_t T[N];
static double F[N];

static uint32_t mix(uint32_t x, uint32_t k) {
    x ^= k;
    x *= 0x7feb352du;
    x ^= x >> 15;
    x *= 0x846ca68bu;
    x ^= x >> 16;
    return x;
}

int main(void) {
    uint32_t acc = 0x9e3779b9u;
    double facc = 0.0;
    int i, r, taken = 0;
    for (i = 0; i < N; i++) {
        T[i] = mix((uint32_t)i, 0x85ebca6bu);
        F[i] = (double)((i % 97) - 48) / 17.0;
    }
    for (r = 0; r < ROUNDS; r++) {
        for (i = 0; i < N; i++) {
            uint32_t x = T[(i * 17 + r) & (N - 1)];
            uint32_t y = mix(x, acc + (uint32_t)r);
            if ((y & 7u) == 0u) {
                acc ^= y * 3u;
                taken++;
            } else if (y & 1u) {
                acc += y;
            } else {
                acc = mix(acc, y);
            }
            T[i] = y;
            facc += F[i] * F[(i + r) & (N - 1)] - 0.015 * F[i];
            if (facc > 1e6) facc *= 0.5;
            if (facc < -1e6) facc *= 0.5;
        }
    }
    printf("acc=%u taken=%d facc=%.3f\n", acc, taken, facc);
    return 0;
}
