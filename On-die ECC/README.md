# On-die ECC

Target scheme : SEC or SEC-DED with bounded fault <br>
Correct all of single bit and if miscorrection happen, it happen bounded region of 16b width.<br><br>

# SEC_(136,128)
�Ϲ����� on-die ECC�� (136,128) SEC�� H_matrix�� Fault simulation�� ���� ECC scheme�� �ŷڼ� ���� <br>

# SEC_DED(288,272)
�ñ������� ���ϴ� on-die ECC�� SEC-DED�� ������ �� �ִ� 16b bounded H matrix�� ���� <br>

# Bounded H matrix
row 0-4���� 5bit���� 1�� ������ Ȧ�� �� bit�� ���� 16���� ������, �̰ɷ� 16column���� �����Ͽ�, 16 bound ���� <br>
���� ������ row 5-15������ 16 bound �ȿ��� ��� column���� ������ ���ڸ� �־ bounded �����ϰ� �� <br>
�̶� 1�� ������ ¦������ �����Ͽ� column ��ü�� ���� ��, ��� column�� 1�� ������ Ȧ��. <br>
�� ���� �ΰ��� column�� xor�Ͽ��� �� 1�� ���� ¦����. �� DED�� �����Ѵ�. <br>
H matrix�� ������ �پ��ؼ� ������ ã�Ƽ� ������ H matrix ã�Ƴ��� <br>

# Implementation
1. git clone
2. make
3. ./Fault_sim_start