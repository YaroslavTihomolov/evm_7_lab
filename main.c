#include <stdio.h>
#include <malloc.h>

#define N 3
#define M 10

inline float max(float a, float b) {
    return a > b ? a : b;
}

void multiply(float** a, float** b, float** R) {
    for(int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i][j] = 0;
            for (int k = 0; k < N; k++)
                R[i][j] += a[i][k] * b[k][j];
        }
    }
}

void sum(float** a, float** b) {
    for(int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            a[i][j] += b[i][j];
}


int main() {
    float** matrix = malloc(sizeof(float*) * N);
    float** I = malloc(sizeof(float*) * N);
    float** B = malloc(sizeof(float*) * N);
    float** R = malloc(sizeof(float*) * N);
    float** deg = malloc(sizeof(float*) * N);
    float** cur = malloc(sizeof(float*) * N);
    for (int i = 0; i < N; i++) {
        matrix[i] = malloc(sizeof(float) * N);
        I[i] = malloc(sizeof(float) * N);
        B[i] = malloc(sizeof(float) * N);
        R[i] = malloc(sizeof(float) * N);
        deg[i] = malloc(sizeof(float) * N);
        cur[i] = malloc(sizeof(float) * N);
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i != j)
                matrix[i][j] = 5;
            else
                matrix[i][j] = 4;
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
            cur_sum_i += matrix[i][j];
            cur_sum_b += matrix[j][i];
        }
        ai = max(ai, cur_sum_i);
        ab = max(ab, cur_sum_b);
    }
    float p = ai * ab;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            B[i][j] = matrix[j][i] / p;
        }
    }
    multiply(B, matrix, R);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                R[i][j] = 1 - R[i][j];
                I[i][j] = 1;
            }
            else {
                R[i][j] *= -1;
                I[i][j] = 0;
            }
            deg[i][j] = R[i][j];
        }
    }



    for (int i = 1; i <= M; i++) {
        for (int h = 0; h < N; h++)
            for (int j = 0; j < N; j++)
                cur[h][j] = deg[h][j];
        sum(I, deg);
        multiply(R, cur, deg);
    }

    multiply(I, B, R);

    for (int h = 0; h < N; h++) {
        for (int j = 0; j < N; j++) {
            printf("%f ", R[h][j]);
        }
        printf("\n");
    }


    for (int i = 0; i < N; i++) {
        free(matrix[i]);
        free(B[i]);
        free(R[i]);
    }

    free(matrix);
    free(B);
    free(R);
    return 0;
}
