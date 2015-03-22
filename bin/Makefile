CC=icc

CFLAGS=-I../cpp_src -W -Wall -O3 -march=native -fast -ansi-alias \
			-fopenmp -finline -funroll-loops -no-prec-div -opt-prefetch -unroll-aggressive \
			-m64 -auto-ilp32 -I../src -xCORE-AVX2 -fma -mkl -I/opt/local/include/

LDFLAGS=-L/opt/local/lib -fopenmp -lboost_system-mt -lboost_timer-mt -mkl -lPcmMsr

EXEC=main.icc
#CFLAGS=-Wa,-q -I../cpp_src -W -Wall -ansi -pedantic -O3 -Wall -funroll-loops -ftree-vectorize \
#			-msse4a -march=native -mtune=native -ffast-math -fopenmp \
#			-fstrict-aliasing -fprefetch-loop-arrays -I../src

SRC= \
	../src/CBitGenerator/CBitGenerator.cpp \
	../src/CChanel/CChanel.cpp \
	../src/CChanel/CChanelAWGN_MKL.cpp \
	../src/CDecoder/template/CDecoder.cpp \
	../src/CDecoder/template/CDecoder_fixed.cpp \
	../src/CDecoder/template/CDecoder_fixed_AVX.cpp \
	../src/CDecoder/template/CDecoder_fixed_SSE.cpp \
	../src/CDecoder/OMS/CDecoder_OMS_fixed_SSE.cpp \
	../src/CDecoder/OMS/CDecoder_OMS_fixed_AVX.cpp \
 ../src/CDecoder/NMS/CDecoder_NMS_fixed_SSE.cpp \
 ../src/CDecoder/NMS/CDecoder_NMS_fixed_AVX.cpp \
 ../src/CEncoder/CFakeEncoder.cpp \
 ../src/CEncoder/Encoder.cpp \
 ../src/CEncoder/GenericEncoder.cpp \
 ../src/CErrorAnalyzer/CErrorAnalyzer.cpp \
 ../src/CFixPointConversion/CFastFixConversion.cpp \
 ../src/CFixPointConversion/CFixConversion.cpp \
 ../src/CTerminal/CTerminal.cpp \
 ../src/CTimer/CTimer.cpp \
 ../src/CTools/CTools.cpp \
 ../src/CTools/transpose_avx.cpp \
 ../src/CTrame/CTrame.cpp \
 ../src/main_p.cpp

OBJ= $(SRC:.cpp=.o)

all: $(EXEC)

main.icc: $(OBJ)
	@echo "[LINKING] $@"
	@$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	@echo "[C++] $<"
	@$(CC) $(CFLAGS) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
	find ../src/ -name "*.o" -exec rm {} \;
	find . -name "*.ic*" -exec rm {} \;

mrproper: clean
	rm $(EXEC)
