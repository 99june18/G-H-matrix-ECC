# On-die ECC

Target scheme : SEC or SEC-DED with bounded fault <br>
Correct all of single bit and if miscorrection happen, it happen bounded region of 16b width.<br><br>

# SEC_(136,128)
일반적인 on-die ECC인 (136,128) SEC의 H_matrix와 Fault simulation을 통해 ECC scheme의 신뢰성 측정 <br>

# SEC_DED(288,272)
궁국적으로 원하는 on-die ECC의 SEC-DED을 구현할 수 있는 16b bounded H matrix을 구현 <br>

# Bounded H matrix
row 0-4까지 5bit에서 1의 개수가 홀수 인 bit만 고르면 16개가 나오고, 이걸로 16column으로 나열하여, 16 bound 생성 <br>
이후 나머지 row 5-15까지는 16 bound 안에서 모든 column에서 동일한 숫자를 넣어서 bounded 성립하게 함 <br>
이때 1의 개수를 짝수개로 설정하여 column 전체로 봤을 때, 모든 column의 1의 개수가 홀수. <br>
그 말은 두개의 column을 xor하였을 때 1의 개수 짝수개. 즉 DED가 성립한다. <br>
H matrix의 종류가 다양해서 여러개 찾아서 최적의 H matrix 찾아낸다 <br>

# Implementation
1. git clone
2. make
3. ./Fault_sim_start