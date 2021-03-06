ghz=[{
"ghz":[11, 10],
"qasm":"""
OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
h q[11];
cx q[11],q[10];
measure q[11] -> c[0];
measure q[10] -> c[1];
""",
"status":[None, None],
"counts":{'0000000000000000': 501, '0000010000000000': 48, '0000100000000000': 44, '0000110000000000': 431}
},{
"ghz":[6, 11, 10],
"qasm":"""
OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
h q[6];
cx q[6],q[11];
cx q[11],q[10];
measure q[6] -> c[0];
measure q[11] -> c[1];
measure q[10] -> c[2];
""",
"status":[None, None],
"counts":{'0000000000000000': 359, '0000000001000000': 297, '0000010000000000': 26, '0000010001000000': 24, '0000100000000000': 23, '0000100001000000': 12, '0000110000000000': 141, '0000110001000000': 142}
},{
"ghz":[13, 12, 11, 10],
"qasm":"""
OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
h q[13];
h q[12];
h q[13];
cx q[12],q[13];
h q[12];
h q[13];
cx q[12],q[11];
cx q[11],q[10];
measure q[13] -> c[0];
measure q[12] -> c[1];
measure q[11] -> c[2];
measure q[10] -> c[3];
""",
"status":[None, None],
"counts":{'0000000000000000': 441, '0000010000000000': 11, '0000100000000000': 12, '0000110000000000': 7, '0001000000000000': 33, '0001010000000000': 8, '0001100000000000': 9, '0001110000000000': 55, '0010000000000000': 87, '0010010000000000': 3, '0010100000000000': 8, '0010110000000000': 32, '0011000000000000': 16, '0011010000000000': 24, '0011100000000000': 27, '0011110000000000': 251}
},{
"ghz":[4, 13, 12, 11, 10],
"qasm":"""
OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
h q[4];
h q[13];
h q[4];
cx q[13],q[4];
h q[13];
h q[4];
h q[12];
h q[13];
cx q[12],q[13];
h q[12];
h q[13];
cx q[12],q[11];
cx q[11],q[10];
measure q[4] -> c[0];
measure q[13] -> c[1];
measure q[12] -> c[2];
measure q[11] -> c[3];
measure q[10] -> c[4];
""",
"status":[None, None],
"counts":{'0000000000000000': 342, '0000000000010000': 41, '0000010000000000': 8, '0000010000010000': 2, '0000100000000000': 7, '0000100000010000': 1, '0000110000000000': 12, '0000110000010000': 7, '0001000000000000': 30, '0001000000010000': 3, '0001010000000000': 3, '0001010000010000': 2, '0001100000000000': 8, '0001100000010000': 2, '0001110000000000': 23, '0001110000010000': 32, '0010000000000000': 23, '0010000000010000': 33, '0010010000010000': 5, '0010100000000000': 1, '0010100000010000': 7, '0010110000000000': 6, '0010110000010000': 25, '0011000000000000': 5, '0011000000010000': 8, '0011010000000000': 2, '0011010000010000': 17, '0011100000000000': 4, '0011100000010000': 23, '0011110000000000': 49, '0011110000010000': 293}
},]