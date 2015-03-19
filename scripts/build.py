#!/usr/bin/python
import subprocess
import os

C_LIST = []
for root, dirs, files in os.walk("../src/Constantes/"):
  for code in dirs:
    C_LIST.append( code )
  break
  
#C_LIST=[
#		'155x93',   '200x100',   '816x408', '1024x518', '1056x528', '1200x600',
#		'1248x624', '2640x1320', 
#		'4000x2000', '4896x2448', '8000x4000', '9972x4986', '20000x10000',
#	    '2388x597',   
# 802.11e
#	    '576x288', '960x480', '2304x1152',
# 802.11n
#		'1944x972', '1944x648', '1944x486', 
# 802.11an
#	    '2048x384',
# Codes DVB-NGH
#		'9216x4608', '9216x2304',
# Codes DVB-NGH
#		'16200x9000', '16200x7560', '16200x6480', '16200x5400', '16200x4320', '16200x2880',
# Codes DVB-S2
#		'64800x6480', '64800x7200', '64800x10800', '64800x16200', '64800x21600', '64800x32400',
#		'4000x2000']
GCC_LIST=['']

subprocess.call('rm main.icc main.icc.*',shell=True)
for C in C_LIST:
	subprocess.call('echo "CODE ' + C + '"  >> results_article',shell=True)
	subprocess.call('echo "#include \\"./' + C + '/constantes.h\\""     > ../src/Constantes/constantes.h',shell=True)
	subprocess.call('echo "#include \\"./' + C + '/constantes_sse.h\\"" > ../src/Constantes/constantes_sse.h',shell=True)
	for COMPILER in GCC_LIST:
		subprocess.call('echo "COMPILATION OF LDPC DECODER FOR CODE [' + C + ']"',shell=True)
		output = subprocess.check_output('make -f Makefile' + COMPILER + ' clean',shell=True)
		output = subprocess.check_output('make -f Makefile' + COMPILER + ' -j 8',shell=True)
		output = subprocess.check_output('mv main.icc main.icc.' + C + '',shell=True)
 		#check_call
