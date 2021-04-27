# ccsx
efficient tool for generating  circular consensus sequences (ccs) from subreads. <br>
Ccsx works with both x86-64 CPUs and ARM CPUs supporting the NEON instruction sets(https://github.com/DLTcollab/sse2neon.git).  <br>
Efficiency is guaranteed with the use of bsalign(https://github.com/ruanjue/bsalign.git) <br>


# install
git clone https://github.com/110allan/ccsx.git <br>
#arm64
[git clone https://github.com/DLTcollab/sse2neon.git] <br>
cd ccsx <br>
make <br>
./ccsx<br>

# usage
Usage  : ccsx  [options] <INPUT> <OUTPUT> <br>
Generate circular consensus sequences (ccs) from subreads. <br>

Options:<br>
-h             Output this help <br>
-v             debug <br>
-m     <int>   Minimum length of subreads to use for generating CCS. [10] <br>
-M     <int>   Maximum length of subreads to use for generating CCS. [50000] <br>
-c     <int>   Minimum number of subreads required to generate CCS. [3] <br>
-j     <int>   Number of threads to use, 0 means autodetection. [2] <br>

Arguments:<br>
input          Input bam file.<br>
output         Output fasta file.<br>

