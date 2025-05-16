#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <vector>
#include <string>

#define RUN_NUM        1000000
#define CHIP_NUM       1
#define OECC_CW_LEN    288
#define OECC_REDUN_LEN 16
#define BOUND          16          // 열 0~15 = 한 16-bit 바운드

/* ────────── H-행렬 로드 ────────── */
unsigned int H_Matrix[OECC_REDUN_LEN][OECC_CW_LEN];
void load_H() {
    //FILE *fp = fopen("H_matrix_SEC_16bound.txt", "r");
    // FILE *fp = fopen("H_matrix_SEC_16bound1.txt", "r");
    FILE *fp = fopen("H_matrix_SEC_16bound2.txt", "r");
    if (!fp) { perror("fopen"); exit(EXIT_FAILURE); }
    for (int r = 0; r < OECC_REDUN_LEN; ++r)
        for (int c = 0; c < OECC_CW_LEN; ++c)
            if (fscanf(fp, "%u", &H_Matrix[r][c]) != 1) {
                fprintf(stderr,"Format error @(%d,%d)\n",r,c);
                fclose(fp); exit(EXIT_FAILURE);
            }
    fclose(fp);
}

/* 유틸 */
inline void flip(int chip,unsigned int cw[][OECC_CW_LEN],int pos){ cw[chip][pos]^=1; }

inline bool is_zero(const unsigned int *v,int len){ for(int i=0;i<len;++i) if(v[i]) return false; return true; }

std::string cw_to_str(const unsigned int cw[OECC_CW_LEN]){
    std::string s; s.reserve(OECC_CW_LEN);
    for(int c=0;c<OECC_CW_LEN;++c) s.push_back('0'+cw[c]);
    return s;
}

/* 디코드+정정  (0=DUE 1=CE 2=SDC) */
int decode_correct(unsigned int cw[OECC_CW_LEN], bool *in_match){
    unsigned int S[OECC_REDUN_LEN]={0};
    for(int r=0;r<OECC_REDUN_LEN;++r)
        for(int c=0;c<OECC_CW_LEN;++c)
            S[r]^=H_Matrix[r][c]*cw[c];

    /* 바운드-내 매칭? */
    *in_match=false;
    for(int col=0;col<BOUND && !(*in_match);++col){
        bool same=true;
        for(int r=0;r<OECC_REDUN_LEN;++r)
            if(S[r]!=H_Matrix[r][col]){ same=false; break; }
        if(same) *in_match=true;
    }

    /* 1-bit 정정 */
    for(int err=0;err<OECC_CW_LEN;++err){
        bool same=true;
        for(int r=0;r<OECC_REDUN_LEN;++r)
            if(S[r]!=H_Matrix[r][err]){ same=false; break; }
        if(same){ cw[err]^=1; break; }
    }

    /* S_after */
    unsigned int Sa[OECC_REDUN_LEN]={0};
    for(int r=0;r<OECC_REDUN_LEN;++r)
        for(int c=0;c<OECC_CW_LEN;++c)
            Sa[r]^=H_Matrix[r][c]*cw[c];

    bool syn0=is_zero(Sa,OECC_REDUN_LEN);
    bool cw0 =is_zero(cw,OECC_CW_LEN);
    if(syn0 && cw0)      return 1;
    else if(syn0 && !cw0)return 2;
    else                 return 0;
}

/* ─────────────────────────── main ─────────────────────────── */
int main(){
    load_H();
    static unsigned int Chip[CHIP_NUM][OECC_CW_LEN];
    std::srand((unsigned)std::time(nullptr));

    int cnt_CE=0,cnt_SDC=0,cnt_DUE=0,cnt_in=0;
    std::vector<std::pair<std::string,std::string>> sdc_cases; // 저장

    for(int run=0;run<RUN_NUM;++run){
        std::memset(Chip,0,sizeof(Chip));

        /* Codeword 전체에서 random bit error 발생. */
        //int p1 = rand()%OECC_CW_LEN;
        //int p2; do{ p2=rand()%OECC_CW_LEN;} while(p1==p2);
        //int p3; do{ p3=rand()%OECC_CW_LEN;} while(p1==p2==p3);
        
        /* Bound 내에서의 error 발생, 첫 16b bound에 한정함. */
        int p1 = rand()%BOUND;
        int p2; do{ p2=rand()%BOUND;} while(p1==p2);
        int p3; do{ p3=rand()%BOUND;} while(p1==p2==p3);

        flip(0,Chip,p1); 
        flip(0,Chip,p2);
        flip(0,Chip,p3);

        std::string before = cw_to_str(Chip[0]);
        
        bool inbound;
        int status = decode_correct(Chip[0], &inbound);
        if(inbound) ++cnt_in;

        std::string after  = cw_to_str(Chip[0]);

        if(status==1) ++cnt_CE;
        else if(status==2){ ++cnt_SDC; sdc_cases.emplace_back(before,after); }
        else ++cnt_DUE;
    }

    /* 통계 출력 */
    printf("\n================ summary ================\n");
    printf("Runs                       : %d\n", RUN_NUM);
    printf("CE   (Correct)             : %d\n", cnt_CE);
    printf("DUE  (Detect-only)         : %d\n", cnt_DUE);
    printf("SDC  (Miscorrect)          : %d\n", cnt_SDC);
    printf("CE ratio                   : %.4f%%\n",100.0*cnt_CE /RUN_NUM);
    printf("DUE ratio                  : %.4f%%\n",100.0*cnt_DUE/RUN_NUM);
    printf("SDC ratio                  : %.4f%%\n",100.0*cnt_SDC/RUN_NUM);
    printf("=============== Bounded ? ==================\n");
    printf("S_before matches col0-15   : %d\n", cnt_in);
    printf("Matches ratio              : %.4f%%\n",100.0*cnt_in/RUN_NUM);

    /* SDC 사례 출력 */
    /* if(!sdc_cases.empty()){
        puts("\n--- SDC cases (Before > After) ---");
        for(const auto& p : sdc_cases){
            printf("Before: %s\nAfter : %s\n\n", p.first.c_str(), p.second.c_str());
        }
    }
    */
    
    return 0;
}
