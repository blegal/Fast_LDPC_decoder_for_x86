# Fast_LDPC_decoder_for_x86

This is the source codes of the fast x86 LDPC decoder published (not yet) in
the IEEE TDPS journal, article.

In this git repository we published the source code of a LDPC decoder
implementation optimized for x86 target. This LDPC decoder implementation
efficiently takes advantage of the SIMD and SPMT programming model. The
approach used to achieve very high

Originally, this source code piece is part of a much larger project
that enables us to experiment LDPC decoding algorithm, data format, etc.
It explains why some piece of code are useless for you ;-)

In odrer to compile the LDPC decoder, you currently have to use Intel C++
compiler. Indeed, the ANWG channel model is implemented using the MKL library
function.

In order to compile the LDPC decoders, just open a terminal and go in
the "bin" directory.

> cd source_path
> cd bin

Then compile the source codes using "make"

> make

The output must look like this:

[C++] ../src/CBitGenerator/CBitGenerator.cpp
[C++] ../src/CChanel/CChanel.cpp
[C++] ../src/CChanel/CChanelAWGN_MKL.cpp
[C++] ../src/CDecoder/template/CDecoder.cpp
[C++] ../src/CDecoder/template/CDecoder_fixed.cpp
[C++] ../src/CDecoder/template/CDecoder_fixed_AVX.cpp
[C++] ../src/CDecoder/template/CDecoder_fixed_SSE.cpp
[C++] ../src/CDecoder/OMS/CDecoder_OMS_fixed_SSE.cpp
[C++] ../src/CDecoder/OMS/CDecoder_OMS_fixed_AVX.cpp
[C++] ../src/CDecoder/NMS/CDecoder_NMS_fixed_SSE.cpp
[C++] ../src/CDecoder/NMS/CDecoder_NMS_fixed_AVX.cpp
[C++] ../src/CEncoder/CFakeEncoder.cpp
[C++] ../src/CEncoder/Encoder.cpp
[C++] ../src/CEncoder/GenericEncoder.cpp
[C++] ../src/CErrorAnalyzer/CErrorAnalyzer.cpp
[C++] ../src/CFixPointConversion/CFastFixConversion.cpp
[C++] ../src/CFixPointConversion/CFixConversion.cpp
[C++] ../src/CTerminal/CTerminal.cpp
[C++] ../src/CTimer/CTimer.cpp
[C++] ../src/CTools/CTools.cpp
[C++] ../src/CTools/transpose_avx.cpp
[C++] ../src/CTrame/CTrame.cpp
[C++] ../src/main_p.cpp
[LINKING] main.icc

The compilation of the 576x288 LDPC decoder (default configuration) was successful,
the executable file is named "main.icc". To launch the LDPC decoder compiled
(576x288), just execute "main.icc" with some parameters:

OMS decoder (offset = 1/8)
> ./main.icc -fixed -avx -OMS 1 -min 0.5 -max 4.0 -iter 20
> ./main.icc -fixed -avx -OMS 1 -min 0.5 -max 4.0 -iter 20

NMS decoder (factor = 29/32)
> ./main.icc -fixed -avx -NMS 29 -min 0.5 -max 4.0 -iter 20
> ./main.icc -fixed -avx -OMS 1 -min 0.5 -max 4.0 -iter 20

To compile more LDPC decoders (code length, etc), execute
the build.py script from the bin directory:

> ../scripts/build.py
