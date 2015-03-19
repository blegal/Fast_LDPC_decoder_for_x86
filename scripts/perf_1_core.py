#!/usr/bin/python
import subprocess
import os

C_LIST = []
for root, dirs, files in os.walk("../src/Constantes/"):
  for code in dirs:
    C_LIST.append( code )
  break
  
C_LIST=[
  '576x288', '1944x972', '2048x384', '2304x1152', '4000x2000'
]
GCC_LIST=['.icc']

for C in C_LIST:
	for COMPILER in GCC_LIST:
		subprocess.call('echo "EXECUTION OF LDPC DECODER FOR CODE [' + C + ']"',shell=True)
		subprocess.call('echo "xMS - 10 its - sse core"',shell=True)
		subprocess.call('./main.icc.' + C + ' -fixed -sse -OMS  1 -fer 10000000 -min 0.50 -max 0.51 -iter 10 -timer 30 -thread 1 | grep "(PERF1) Total Kernel throughput"',shell=True)
		subprocess.call('./main.icc.' + C + ' -fixed -sse -NMS 29 -fer 10000000 -min 0.50 -max 0.51 -iter 10 -timer 30 -thread 1 | grep "(PERF1) Total Kernel throughput"',shell=True)
		subprocess.call('echo "xMS - 10 its - avx cores"',shell=True)
		subprocess.call('./main.icc.' + C + ' -fixed -avx -OMS  1 -fer 10000000 -min 0.50 -max 0.51 -iter 10 -timer 30 -thread 1 | grep "(PERF1) Total Kernel throughput"',shell=True)
		subprocess.call('./main.icc.' + C + ' -fixed -avx -NMS 29 -fer 10000000 -min 0.50 -max 0.51 -iter 10 -timer 30 -thread 1 | grep "(PERF1) Total Kernel throughput"',shell=True)
		subprocess.call('echo "xMS - 20 its - sse core"',shell=True)
		subprocess.call('./main.icc.' + C + ' -fixed -sse -OMS  1 -fer 10000000 -min 0.50 -max 0.51 -iter 20 -timer 30 -thread 1 | grep "(PERF1) Total Kernel throughput"',shell=True)
		subprocess.call('./main.icc.' + C + ' -fixed -sse -NMS 29 -fer 10000000 -min 0.50 -max 0.51 -iter 20 -timer 30 -thread 1 | grep "(PERF1) Total Kernel throughput"',shell=True)
		subprocess.call('echo "xMS - 20 its - avx cores"',shell=True)
		subprocess.call('./main.icc.' + C + ' -fixed -avx -OMS  1 -fer 10000000 -min 0.50 -max 0.51 -iter 20 -timer 30 -thread 1 | grep "(PERF1) Total Kernel throughput"',shell=True)
		subprocess.call('./main.icc.' + C + ' -fixed -avx -NMS 29 -fer 10000000 -min 0.50 -max 0.51 -iter 20 -timer 30 -thread 1 | grep "(PERF1) Total Kernel throughput"',shell=True)
