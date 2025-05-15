# System ECC

Target scheme : RS16(18,16) + DEC <br>
Correct all of single, double, symbol error.

# build_SSC_DEC.py

Goal : Build an SSC-DEC H matrix for RS(18,16) over GF(2^16).

<코드 동작>
1. GF(16)의 2048개 primitive polynomial 중 하나를 선택
2. 선택한 primitive polynomial을 사용해 H matrix의 후보 column 2^16-1개를 생성
3. 생성된 column들 내에서 SSC-DEC 조건을 만족하는(Unity ECC paper의 Algorithm 1 참조) 18개의 column 조합을 탐색
   3-1. 성공 시 column의 조합 = 2 x 18 symbol 단위 matrix를 반환하고, binary로 변환한 32 x 288 bit matrix까지 반환한 후 탐색 종료
   3-2. 실패 시 fail 메시지를 띄우고 다음 primitive polynomial로 탐색 반복


Input  : GF_2^16__primitive_polynomial.txt

Outputs on success:
  H_indices.txt        chosen 18 column numbers
  H_symbol.txt         2 x 18 matrix, 16-bit hex values
  H_binary_32x288.txt  32 x 288 matrix, ASCII 0/1
