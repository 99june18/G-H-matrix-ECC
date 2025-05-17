#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <string>

#define RUN_NUM        1000000        // �ݺ� Ƚ��
#define CHIP_NUM       1
#define OECC_CW_LEN    288
#define OECC_REDUN_LEN 16
#define BOUND          16             // 1 �ٿ�� = 16��

/* �������������������� H-��� �ε� (���� �״��) �������������������� */
unsigned int H_Matrix[OECC_REDUN_LEN][OECC_CW_LEN];

void load_H() {
    //FILE *fp = fopen("H_matrix_SEC_16bound.txt",  "r");
    //FILE *fp = fopen("H_matrix_SEC_16bound1.txt", "r");
    //FILE *fp = fopen("H_matrix_SEC_16bound2.txt", "r");
    //FILE *fp  = fopen("H_matrix_SEC_16bound3.txt", "r");
    FILE *fp = fopen("H_matrix_SEC_16bound4.txt", "r");     //���ݱ��� ���� ���� ������.
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

/* �������������������� ���� ���� �Լ��� �������������������� */
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

/* 3-bit ���� (�� ��Ʈ�� **���� �ٸ� 16-bit �ٿ��**�� ��ġ)       */
/* start �ٿ�忡�� 1-bit, ������ �� ��Ʈ�� �ٸ� �ٿ�忡�� ���� ���� */
inline void flip_SE_SE_SE(unsigned int cw[OECC_CW_LEN], int start) {
    const int NUM_BLOCKS = OECC_CW_LEN / BOUND;          // 288 / 16 = 18

    /* �� ù ��Ʈ: start �ٿ�� �ȿ��� ���� */
    int p1 = start + (std::rand() % BOUND);
    int block1 = start / BOUND;

    /* �� �� ��° ��Ʈ: block1 �̿��� �ٿ�� ���� */
    int block2;
    do { block2 = std::rand() % NUM_BLOCKS; } while (block2 == block1);
    int p2 = block2 * BOUND + (std::rand() % BOUND);

    /* �� �� ��° ��Ʈ: block1, block2 �� �ٸ� �ٿ�� ���� */
    int block3;
    do { block3 = std::rand() % NUM_BLOCKS; }
    while (block3 == block1 || block3 == block2);
    int p3 = block3 * BOUND + (std::rand() % BOUND);

    /* ��Ʈ ������ */
    cw[p1] ^= 1; cw[p2] ^= 1; cw[p3] ^= 1;
}


inline void flip_bounded_2b_2b(unsigned int cw[OECC_CW_LEN], int start){
    int p1 = start + (std::rand() % BOUND);
    int p2; do{ p2 = start + (std::rand() % BOUND);} while(p2==p1);
    int p3 = (start + 16) + (std::rand() % BOUND);
    int p4; do{ p4 = (start + 16) + (std::rand() % BOUND);} while(p3==p4);
    cw[p1]^=1; cw[p2]^=1; cw[p3]^=1; cw[p4]^=1;
}

/* �������������������� ��ƿ �������������������� */
inline bool is_zero(const unsigned int *v,int n){
    for(int i=0;i<n;++i) if(v[i]) return false; return true;
}
std::string cw_str(const unsigned int cw[OECC_CW_LEN]){
    std::string s; s.reserve(OECC_CW_LEN);
    for(int i=0;i<OECC_CW_LEN;++i) s.push_back('0'+cw[i]);
    return s;
}

/* �������������������� �ٿ��-�� �ŵ�� ��Ī �Լ� �������������������� */
bool is_in_bound(const unsigned int S[OECC_REDUN_LEN], int start){
    for(int col=start; col<start+BOUND; ++col){
        bool same=true;
        for(int r=0;r<OECC_REDUN_LEN;++r)
            if(S[r]!=H_Matrix[r][col]){ same=false; break; }
        if(same) return true;
    }
    return false;
}

/* �������������������� ���ڵ� + ���� ��������������������
   ��ȯ�� : 0=DUE  1=CE  2=SDC
   S_before_out[16] �� �ŵ�� ����                              */
int decode_correct(unsigned int cw[OECC_CW_LEN],
                   unsigned int S_before_out[OECC_REDUN_LEN])
{
    unsigned int S[OECC_REDUN_LEN]={0};

    /* S_before */
    for(int r=0;r<OECC_REDUN_LEN;++r)
        for(int c=0;c<OECC_CW_LEN;++c)
            S[r] ^= H_Matrix[r][c] * cw[c];

    std::memcpy(S_before_out, S, sizeof(unsigned int)*OECC_REDUN_LEN);

    /* 1-bit ���� */
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

/* �������������������� main �������������������� */
int main(){
    load_H();
    static unsigned int Chip[CHIP_NUM][OECC_CW_LEN];
    std::srand(static_cast<unsigned>(std::time(nullptr)));
    std::string last_sdc_before, last_sdc_after;   // �� �߰�

    int cnt_CE=0, cnt_SDC=0, cnt_DUE=0, cnt_in=0;
    const int start = 0;                 // �˻��� �ٿ�� ���� ��

    int cnt_in1=0, cnt_in2=0, cnt_in3=0, cnt_in4=0, cnt_in5=0, cnt_in6=0;

    for(int run=0; run<RUN_NUM; ++run){
        std::memset(Chip,0,sizeof(Chip));

        /* === ���� �ó����� ���� === */
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

        /* per-run ��� */
        //printf("run %d | %s | in-bound:%s\n", run,
        //        status==1?"CE":status==2?"SDC":"DUE",
        //        inbound?"Y":"N");
        //printf("Before:%s\nAfter :%s\n\n", before.c_str(), after.c_str());

        /* ��� ���� */
        if (status == 1) ++cnt_CE;
        else if (status == 2) {                 // SDC
            ++cnt_SDC;
            last_sdc_before = before;           // �� ������ SDC ����
            last_sdc_after  = after;
        }
        else ++cnt_DUE;
    }

    /* ---- ��� ��� ---- */
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

    if (cnt_SDC > 0) {                      // �� ������ SDC �ڵ���� ���
        printf("\nLast SDC case:\n");
        printf("Before: %s\n", last_sdc_before.c_str());
        printf("After : %s\n", last_sdc_after.c_str());
    }

    return 0;
}
