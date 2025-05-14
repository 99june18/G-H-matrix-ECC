# G-H-matrix-ECC

1. **Step 1** : Make H' matrix, which support RS16(18,16) + DEC (For system ECC) 
2. **Step 2** : Make H2 matrix, which support SEC or SECDED with bounded fault
3. **Step 3** : Make G2 matrix, with transpose the H2 matrix
4. **Step 4** : Organize G1 matrix, by using H'(G'), G2 matrix.