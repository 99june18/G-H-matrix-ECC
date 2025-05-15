#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstring>

#define RUN_NUM        1000
#define CHIP_NUM       10
#define OECC_CW_LEN    136
#define OECC_REDUN_LEN 8

/* ────────── 전역 (136,128) Hamming H-행렬 ────────── */
unsigned int H_Matrix_SEC_Unbound[OECC_REDUN_LEN][OECC_CW_LEN];

void generator_oecc_H_matrix() {
    FILE *fp = fopen("H_matrix_SEC_Unbound.txt", "r");
    if (!fp) { perror("fopen"); exit(EXIT_FAILURE); }

    for (int r = 0; r < OECC_REDUN_LEN; ++r)
        for (int c = 0; c < OECC_CW_LEN; ++c)
            if (fscanf(fp, "%u", &H_Matrix_SEC_Unbound[r][c]) != 1) {
                fprintf(stderr, "Format error @ (%d,%d)\n", r, c);
                fclose(fp); exit(EXIT_FAILURE);
            }
    fclose(fp);
}

/* 1-bit 플립 */
void error_injection_SE(int chip, unsigned int cw[][OECC_CW_LEN]) {
    int pos = rand() % OECC_CW_LEN;
    cw[chip][pos] ^= 1;
}

/* 유틸 */
inline bool is_zero_syndrome(const unsigned int s[OECC_REDUN_LEN]) {
    for (int i = 0; i < OECC_REDUN_LEN; ++i) if (s[i]) return false;
    return true;
}
inline bool is_zero_codeword(const unsigned int cw[OECC_CW_LEN]) {
    for (int i = 0; i < OECC_CW_LEN; ++i) if (cw[i]) return false;
    return true;
}
void print_codeword(const unsigned int cw[OECC_CW_LEN]) {
    for (int c = 0; c < OECC_CW_LEN; ++c) printf("%u", cw[c]);
    putchar('\n');
}

/* 디코더: 정정 수행, S_after 반환 */
void decode_and_correct(int chip,
                        unsigned int cw[][OECC_CW_LEN],
                        unsigned int S_after[OECC_REDUN_LEN])
{
    unsigned int S[OECC_REDUN_LEN] = {0};

    /* S_before */
    for (int r = 0; r < OECC_REDUN_LEN; ++r)
        for (int c = 0; c < OECC_CW_LEN; ++c)
            S[r] ^= H_Matrix_SEC_Unbound[r][c] * cw[chip][c];

    printf("S_before = ");
    for (int r = OECC_REDUN_LEN - 1; r >= 0; --r) printf("%u", S[r]);
    putchar('\n');

    /* 1-bit 후보 탐색 */
    for (int err = 0; err < OECC_CW_LEN; ++err) {
        bool match = true;
        for (int r = 0; r < OECC_REDUN_LEN; ++r)
            if (S[r] != H_Matrix_SEC_Unbound[r][err]) { match = false; break; }
        if (match) { cw[chip][err] ^= 1; break; }
    }

    /* S_after 산출 */
    std::memset(S_after, 0, sizeof(unsigned int) * OECC_REDUN_LEN);
    for (int r = 0; r < OECC_REDUN_LEN; ++r)
        for (int c = 0; c < OECC_CW_LEN; ++c)
            S_after[r] ^= H_Matrix_SEC_Unbound[r][c] * cw[chip][c];

    printf("S_after  = ");
    for (int r = OECC_REDUN_LEN - 1; r >= 0; --r) printf("%u", S_after[r]);
    putchar('\n');
}

/* ─────────────────────────────── main ─────────────────────────────── */
int main() {
    generator_oecc_H_matrix();

    static unsigned int Chip_array[CHIP_NUM][OECC_CW_LEN];
    std::srand(static_cast<unsigned>(std::time(nullptr)));

    const int fault_chip = 0;
    int cnt_CE  = 0;   /* Correctable & fixed      */
    int cnt_SDC = 0;   /* Syndr=0 but data wrong   */
    int cnt_UE  = 0;   /* Uncorrectable (S≠0)      */

    for (int run = 0; run < RUN_NUM; ++run) {
        std::memset(Chip_array, 0, sizeof(Chip_array));   // golden codeword: all-zero

        /* ---- 오류 주입: 두 비트를 반드시 다르게 ---- */
        int p1 = rand() % OECC_CW_LEN;
        Chip_array[fault_chip][p1] ^= 1;
        /* ---- double bit error inject ---- */
        // int p2;
        // do { p2 = rand() % OECC_CW_LEN; } while (p2 == p1);
        // Chip_array[fault_chip][p2] ^= 1;

        printf("run %3d | Before codeword = ", run);
        print_codeword(Chip_array[fault_chip]);

        /* ---- 디코딩 & 정정 ---- */
        unsigned int S_after[OECC_REDUN_LEN];
        decode_and_correct(fault_chip, Chip_array, S_after);

        printf("run %3d | After  codeword = ", run);
        print_codeword(Chip_array[fault_chip]);

        bool syn0 = is_zero_syndrome(S_after);
        bool cw0  = is_zero_codeword(Chip_array[fault_chip]);

        if (syn0 && cw0) {            /* 정상 정정 → CE */
            ++cnt_CE;
            printf("==> CE  (correctable & clean)\n\n");
        } else if (syn0 && !cw0) {    /* 신드롬 0인데 데이터 오염 → SDC */
            ++cnt_SDC;
            printf("==> SDC (silent miscorrection)\n\n");
        } else {                      /* 신드롬 비0  → UE */
            ++cnt_UE;
            printf("==> UE  (uncorrectable)\n\n");
        }
    }

    /* ---- 통계 ---- */
    printf("=================================================\n");
    printf("Runs : %d\n", RUN_NUM);
    printf("CE   : %d\n", cnt_CE);
    printf("SDC  : %d\n", cnt_SDC);
    printf("UE   : %d\n", cnt_UE);
    printf("CE ratio  = %.2f%%\n",  100.0 * cnt_CE  / RUN_NUM);
    printf("SDC ratio = %.2f%%\n",  100.0 * cnt_SDC / RUN_NUM);
    printf("UE ratio  = %.2f%%\n",  100.0 * cnt_UE  / RUN_NUM);
    return 0;
}
