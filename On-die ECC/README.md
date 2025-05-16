# On-die ECC

Target scheme : SEC or SEC-DED with bounded fault <br>
Correct all of single bit and if miscorrection happen, it happen bounded region of 16b width.<br><br>

# SEC_(136,128)
일반적인 on-die ECC인 (136,128) SEC의 H_matrix와 Fault simulation을 통해 ECC scheme의 신뢰성 측정 <br>

# SEC_DED(288,272)
궁국적으로 원하는 on-die ECC의 SEC-DED을 구현할 수 있는 16b bounded H matrix을 구현 <br>

# Implementation
1. make
2. ./Fault_sim_start