/*#include <stdio.h>
#include <malloc.h>
#include <xmmintrin.h>

#define N 4
#define M 10

inline float max(float a, float b) {
    return a > b ? a : b;
}

void multiply(float* a, float* b, float* R) {
    __m128 *xx = (__m128*)b;
    __m128* m128_result = (__m128*)R;
    __m128 mult, tmp;
    for(int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            mult = _mm_set1_ps(a[i * N + j]);
            for (int k = 0; k < N / 4; ++k) {
                tmp = _mm_mul_ps(mult, xx[N * j / 4 + k]);
                m128_result[N * i / 4 + k] = _mm_add_ps(m128_result[N * i / 4 + k], tmp);
            }
        }
    }
}

void sum(float* a, const float* b) {
    for(int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            a[i * N + j] += b[i * N + j];
}


int main() {
    float* matrix = malloc(sizeof(float) * N * N);
    float* I = malloc(sizeof(float) * N * N);
    float* B = malloc(sizeof(float) * N * N);
    float* R = malloc(sizeof(float) * N * N);
    float* deg = malloc(sizeof(float) * N * N);
    float* cur = malloc(sizeof(float) * N * N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i != j)
                matrix[i * N + j] = 5;
            else
                matrix[i * N + j] = 4;
        }
    }


    float ai = 0;
    float ab = 0;
    float cur_sum_i;
    float cur_sum_b;

    for (int i = 0; i < N; i++) {
        cur_sum_i = 0;
        cur_sum_b = 0;
        for (int j = 0; j < N; j++) {
            cur_sum_i += matrix[i * N + j];
            cur_sum_b += matrix[i * N + j];
        }
        ai = max(ai, cur_sum_i);
        ab = max(ab, cur_sum_b);
    }
    float p = ai * ab;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            B[i * N + j] = matrix[j * N + i] / p;
        }
    }

    multiply(B, matrix, R);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                R[i * N + j] = 1 - R[i * N + j];
                I[i * N + j] = 1;
            }
            else {
                R[i * N + j] *= -1;
                I[i * N + j] = 0;
            }
            deg[i * N + j] = R[i * N + j];
        }
    }


    for (int i = 1; i <= M; i++) {
        for (int h = 0; h < N; h++)
            for (int j = 0; j < N; j++)
                cur[h * N + j] = deg[h * N + j];
        sum(I, deg);
        multiply(R, cur, deg);
    }


    multiply(I, B, R);

    for (int h = 0; h < N; h++) {
        for (int j = 0; j < N; j++) {
            printf("%lf ", R[h * N + j]);
        }
        printf("\n");
    }

    free(matrix);
    free(B);
    free(R);
    free(deg);
    free(cur);
    free(I);
    return 0;
}*/

#include <malloc.h>
#include <xmmintrin.h>
#include <math.h>

#define N 1024
#define M 10

float max_(float a, float b) {
    return a > b ? a : b;
}

void multiply(float *a, float *b, float *R) {
    __m128 *xx = (__m128 *) b;
    __m128 *m128_result = (__m128 *) R;
    __m128 mult, tmp;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            R[i * N + j] = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            mult = _mm_set1_ps(a[i * N + j]);
            for (int k = 0; k < N / 4; ++k) {
                tmp = _mm_mul_ps(mult, xx[N * j / 4 + k]);
                m128_result[N * i / 4 + k] = _mm_add_ps(m128_result[N * i / 4 + k], tmp);
            }
        }
    }
}

/*void multiply(const float* a, const float* b, float* R) {
    for(int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i * N + j] = 0;
            for (int k = 0; k < N; k++)
                R[i * N + j] += a[i * N + k] * b[k * N + j];
        }
    }
}*/

void sum(float *a, const float *b) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            a[i * N + j] += b[i * N + j];
}


int main() {
    float *matrix = malloc(sizeof(float) * N * N);
    float *I = calloc(N * N, sizeof(float));
    float *B = calloc(N * N, sizeof(float));
    float *R = calloc(N * N, sizeof(float));
    float *deg = calloc(N * N, sizeof(float));
    float *cur = calloc(N * N, sizeof(float));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i != j)
                matrix[i * N + j] = 5.0f;
            else
                matrix[i * N + j] = 4.0f;
        }
    }

    float ai = 0;
    float ab = 0;
    float cur_sum_i;
    float cur_sum_b;

    for (int i = 0; i < N; i++) {
        cur_sum_i = 0;
        cur_sum_b = 0;
        for (int j = 0; j < N; j++) {
            cur_sum_i += fabsf(matrix[i * N + j]);
            cur_sum_b += fabsf(matrix[i * N + j]);
        }
        ai = max_(ai, cur_sum_i);
        ab = max_(ab, cur_sum_b);
    }
    float p = ai * ab;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            B[i * N + j] = matrix[j * N + i] / p;
        }
    }

    multiply(B, matrix, R);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                R[i * N + j] = 1 - R[i * N + j];
                I[i * N + j] = 1;
            } else {
                R[i * N + j] *= -1;
                I[i * N + j] = 0;
            }
            deg[i * N + j] = R[i * N + j];
        }
    }


    for (int i = 1; i <= M; i++) {
        for (int h = 0; h < N; h++)
            for (int j = 0; j < N; j++)
                cur[h * N + j] = deg[h * N + j];
        sum(I, deg);
        multiply(R, cur, deg);
    }

    multiply(I, B, R);

    free(matrix);
    free(B);
    free(R);
    free(deg);
    free(cur);
    free(I);
    return 0;
}
