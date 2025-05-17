#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <string>

#define RUN_NUM        1000000        // 반복 횟수
#define CHIP_NUM       1
#define OECC_CW_LEN    288
#define OECC_REDUN_LEN 16
#define BOUND          16             // 1 바운드 = 16열

/* ────────── H-행렬 로드 (원본 그대로) ────────── */
unsigned int H_Matrix[OECC_REDUN_LEN][OECC_CW_LEN];

void load_H() {
    //FILE *fp = fopen("H_matrix_SEC_16bound.txt",  "r");
    //FILE *fp = fopen("H_matrix_SEC_16bound1.txt", "r");
    //FILE *fp = fopen("H_matrix_SEC_16bound2.txt", "r");
    //FILE *fp  = fopen("H_matrix_SEC_16bound3.txt", "r");
    FILE *fp = fopen("H_matrix_SEC_16bound4.txt", "r");     //지금까진 가장 좋은 성능임.
    //FILE *fp = fopen("H_matrix_SEC_16bound5.txt", "r");
    //FILE *fp = fopen("H_matrix_SEC_16bound6.txt", "r");
    if (!fp) { perror("fopen"); std::exit(EXIT_FAILURE); }

    for (int r = 0; r < OECC_REDUN_LEN; ++r)
        for (int c = 0; c < OECC_CW_LEN; ++c)
            if (fscanf(fp, "%u", &H_Matrix[r][c]) != 1) {
                fprintf(stderr,"Format error (%d,%d)\n", r, c);
                fclose(fp); std::exit(EXIT_FAILURE);
            }
    fclose(fp);
}

/* ────────── 오류 주입 함수들 ────────── */
inline void flip_1b(unsigned int cw[OECC_CW_LEN]){
    cw[ std::rand() % OECC_CW_LEN ] ^= 1;
}

inline void flip_2b(unsigned int cw[OECC_CW_LEN]){
    int p1 = std::rand() % OECC_CW_LEN;
    int p2; do{ p2 = std::rand() % OECC_CW_LEN; } while(p2==p1);
    cw[p1] ^= 1;  cw[p2] ^= 1;
}

inline void flip_3b(unsigned int cw[OECC_CW_LEN]){
    int p1 = std::rand() % OECC_CW_LEN;
    int p2; do{ p2 = std::rand() % OECC_CW_LEN; } while(p2==p1);
    int p3; do{ p3 = std::rand() % OECC_CW_LEN; } while(p3==p1 || p3==p2);
    cw[p1]^=1; cw[p2]^=1; cw[p3]^=1;
}

inline void flip_bounded_3b(unsigned int cw[OECC_CW_LEN], int start){
    int p1 = start + (std::rand() % BOUND);
    int p2; do{ p2 = start + (std::rand() % BOUND);} while(p2==p1);
    int p3; do{ p3 = start + (std::rand() % BOUND);} while(p3==p1||p3==p2);
    cw[p1]^=1; cw[p2]^=1; cw[p3]^=1;
}

/* 3-bit 오류 (각 비트가 **서로 다른 16-bit 바운드**에 위치)       */
/* start 바운드에서 1-bit, 나머지 두 비트는 다른 바운드에서 랜덤 선택 */
inline void flip_SE_SE_SE(unsigned int cw[OECC_CW_LEN], int start) {
    const int NUM_BLOCKS = OECC_CW_LEN / BOUND;          // 288 / 16 = 18

    /* ① 첫 비트: start 바운드 안에서 랜덤 */
    int p1 = start + (std::rand() % BOUND);
    int block1 = start / BOUND;

    /* ② 두 번째 비트: block1 이외의 바운드 선택 */
    int block2;
    do { block2 = std::rand() % NUM_BLOCKS; } while (block2 == block1);
    int p2 = block2 * BOUND + (std::rand() % BOUND);

    /* ③ 세 번째 비트: block1, block2 와 다른 바운드 선택 */
    int block3;
    do { block3 = std::rand() % NUM_BLOCKS; }
    while (block3 == block1 || block3 == block2);
    int p3 = block3 * BOUND + (std::rand() % BOUND);

    /* 비트 뒤집기 */
    cw[p1] ^= 1; cw[p2] ^= 1; cw[p3] ^= 1;
}


inline void flip_bounded_2b_2b(unsigned int cw[OECC_CW_LEN], int start){
    int p1 = start + (std::rand() % BOUND);
    int p2; do{ p2 = start + (std::rand() % BOUND);} while(p2==p1);
    int p3 = (start + 16) + (std::rand() % BOUND);
    int p4; do{ p4 = (start + 16) + (std::rand() % BOUND);} while(p3==p4);
    cw[p1]^=1; cw[p2]^=1; cw[p3]^=1; cw[p4]^=1;
}

/* ────────── 유틸 ────────── */
inline bool is_zero(const unsigned int *v,int n){
    for(int i=0;i<n;++i) if(v[i]) return false; return true;
}
std::string cw_str(const unsigned int cw[OECC_CW_LEN]){
    std::string s; s.reserve(OECC_CW_LEN);
    for(int i=0;i<OECC_CW_LEN;++i) s.push_back('0'+cw[i]);
    return s;
}

/* ────────── 바운드-내 신드롬 매칭 함수 ────────── */
bool is_in_bound(const unsigned int S[OECC_REDUN_LEN], int start){
    for(int col=start; col<start+BOUND; ++col){
        bool same=true;
        for(int r=0;r<OECC_REDUN_LEN;++r)
            if(S[r]!=H_Matrix[r][col]){ same=false; break; }
        if(same) return true;
    }
    return false;
}

/* ────────── 디코드 + 정정 ──────────
   반환값 : 0=DUE  1=CE  2=SDC
   S_before_out[16] 에 신드롬 저장                              */
int decode_correct(unsigned int cw[OECC_CW_LEN],
                   unsigned int S_before_out[OECC_REDUN_LEN])
{
    unsigned int S[OECC_REDUN_LEN]={0};

    /* S_before */
    for(int r=0;r<OECC_REDUN_LEN;++r)
        for(int c=0;c<OECC_CW_LEN;++c)
            S[r] ^= H_Matrix[r][c] * cw[c];

    std::memcpy(S_before_out, S, sizeof(unsigned int)*OECC_REDUN_LEN);

    /* 1-bit 정정 */
    for(int err=0; err<OECC_CW_LEN; ++err){
        bool match=true;
        for(int r=0;r<OECC_REDUN_LEN;++r)
            if(S[r]!=H_Matrix[r][err]){ match=false; break; }
        if(match){ cw[err]^=1; break; }
    }

    /* S_after */
    unsigned int Sa[OECC_REDUN_LEN]={0};
    for(int r=0;r<OECC_REDUN_LEN;++r)
        for(int c=0;c<OECC_CW_LEN;++c)
            Sa[r] ^= H_Matrix[r][c] * cw[c];

    bool syn0 = is_zero(Sa,OECC_REDUN_LEN);
    bool cw0  = is_zero(cw,OECC_CW_LEN);
    if(syn0 && cw0)       return 1;  // CE
    else if(syn0 && !cw0) return 2;  // SDC
    else                  return 0;  // DUE
}

/* ────────── main ────────── */
int main(){
    load_H();
    static unsigned int Chip[CHIP_NUM][OECC_CW_LEN];
    std::srand(static_cast<unsigned>(std::time(nullptr)));
    std::string last_sdc_before, last_sdc_after;   // ★ 추가

    int cnt_CE=0, cnt_SDC=0, cnt_DUE=0, cnt_in=0;
    const int start = 0;                 // 검사할 바운드 시작 열

    int cnt_in1=0, cnt_in2=0, cnt_in3=0, cnt_in4=0, cnt_in5=0, cnt_in6=0;

    for(int run=0; run<RUN_NUM; ++run){
        std::memset(Chip,0,sizeof(Chip));

        /* === 오류 시나리오 선택 === */
        // flip_1b(Chip[0]);
        // flip_2b(Chip[0]);
        //flip_3b(Chip[0]);
        //flip_bounded_3b(Chip[0], start);
        //flip_bounded_2b_2b(Chip[0], start);
        flip_SE_SE_SE(Chip[0], start);
        /* ========================== */

        std::string before = cw_str(Chip[0]);

        unsigned int S_before[OECC_REDUN_LEN];
        int status = decode_correct(Chip[0], S_before);

        bool inbound = is_in_bound(S_before, start);
        if(inbound) ++cnt_in;

        bool inbound1 = is_in_bound(S_before, 16);
        if(inbound1) ++cnt_in1;

        bool inbound2 = is_in_bound(S_before, 32);
        if(inbound2) ++cnt_in2;

        bool inbound3 = is_in_bound(S_before, 48);
        if(inbound3) ++cnt_in3;

        bool inbound4 = is_in_bound(S_before, 64);
        if(inbound4) ++cnt_in4;

        bool inbound5 = is_in_bound(S_before, 80);
        if(inbound5) ++cnt_in5;

        bool inbound6 = is_in_bound(S_before, 96);
        if(inbound6) ++cnt_in6;

        std::string after  = cw_str(Chip[0]);

        /* per-run 출력 */
        //printf("run %d | %s | in-bound:%s\n", run,
        //        status==1?"CE":status==2?"SDC":"DUE",
        //        inbound?"Y":"N");
        //printf("Before:%s\nAfter :%s\n\n", before.c_str(), after.c_str());

        /* 통계 누적 */
        if (status == 1) ++cnt_CE;
        else if (status == 2) {                 // SDC
            ++cnt_SDC;
            last_sdc_before = before;           // ★ 마지막 SDC 저장
            last_sdc_after  = after;
        }
        else ++cnt_DUE;
    }

    /* ---- 통계 출력 ---- */
    printf("\n================ summary ================\n");
    printf("Runs                       : %d\n", RUN_NUM);
    printf("CE   (Correct)             : %d\n", cnt_CE);
    printf("DUE  (Detect-only)         : %d\n", cnt_DUE);
    printf("SDC  (Miscorrect)          : %d\n", cnt_SDC);
    printf("CE ratio                   : %.4f%%\n",100.0*cnt_CE / RUN_NUM);
    printf("DUE ratio                  : %.4f%%\n",100.0*cnt_DUE/ RUN_NUM);
    printf("SDC ratio                  : %.4f%%\n",100.0*cnt_SDC/ RUN_NUM);
    printf("-----------------------------------------\n");
    printf("S_before matches bound(%d-%d): %d (%.4f%%)\n",
           start, start+BOUND-1, cnt_in, 100.0*cnt_in / RUN_NUM);
    printf("S_before matches bound(%d-%d): %d (%.4f%%)\n",
           16, 16+BOUND-1, cnt_in1, 100.0*cnt_in1 / RUN_NUM);
    printf("S_before matches bound(%d-%d): %d (%.4f%%)\n",
           32, 32+BOUND-1, cnt_in2, 100.0*cnt_in2 / RUN_NUM);
    printf("S_before matches bound(%d-%d): %d (%.4f%%)\n",
           48, 48+BOUND-1, cnt_in3, 100.0*cnt_in3 / RUN_NUM);
    printf("S_before matches bound(%d-%d): %d (%.4f%%)\n",
           64, 64+BOUND-1, cnt_in4, 100.0*cnt_in4 / RUN_NUM);
    printf("S_before matches bound(%d-%d): %d (%.4f%%)\n",
           80, 80+BOUND-1, cnt_in5, 100.0*cnt_in5 / RUN_NUM);
    printf("S_before matches bound(%d-%d): %d (%.4f%%)\n",
           96, 96+BOUND-1, cnt_in6, 100.0*cnt_in6 / RUN_NUM);

    if (cnt_SDC > 0) {                      // ★ 마지막 SDC 코드워드 출력
        printf("\nLast SDC case:\n");
        printf("Before: %s\n", last_sdc_before.c_str());
        printf("After : %s\n", last_sdc_after.c_str());
    }

    return 0;
}
