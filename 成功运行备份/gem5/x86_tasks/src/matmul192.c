#include <stdio.h>
#include <stdlib.h>
#define N 192
static double A[N][N], B[N][N], C[N][N];
int main(void) {
    int i, j, k;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++) {
            A[i][j] = (double)((i * 17 + j) % 13);
            B[i][j] = (double)((i + j * 9) % 11);
            C[i][j] = 0.0;
        }
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            for (k = 0; k < N; k++)
                C[i][j] += A[i][k] * B[k][j];
    printf("C00=%.1f C191=%.1f\n", C[0][0], C[N-1][N-1]);
    return 0;
}
