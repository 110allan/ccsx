# ccsx
efficient tool for generating  circular consensus sequences (ccs) from subreads. <br>
Ccsx works with both x86-64 CPUs and ARM CPUs supporting the NEON instruction sets(https://github.com/DLTcollab/sse2neon.git).  <br>
Efficiency is guaranteed with the use of bsalign(https://github.com/ruanjue/bsalign.git) <br>


# install
X86_64: <br>
git clone https://github.com/110allan/ccsx.git <br>
cd ccsx <br>
git clone  https://github.com/ruanjue/bsalign.git <br>
make <br>
./ccsx -h <br>

ARM64:<br>
git clone https://github.com/110allan/ccsx.git <br>
cd ccsx <br>
git clone https://github.com/DLTcollab/sse2neon.git <br>
git clone  https://github.com/ruanjue/bsalign.git <br>
make <br>
./ccsx -h <br>

# usage
Usage  : ccsx  [options] <INPUT> <OUTPUT> <br>
Generate circular consensus sequences (ccs) from subreads.<br>
<br>
Options:<br>
-h             Output this help <br>
-v             debug <br>
-m     <int>   Minimum total length of subreads in a hole to use for generating CCS. [5000] <br>
-M     <int>   Maximum total length of subreads in a hole to use for generating CCS. [500000] <br>
-c     <int>   Minimum number of subreads required to generate CCS. [3] <br>
-A             For fasta/fastq input,gzip allowed  <br>
-P             primitive bsalign,subread shred by default <br>
-X	   <str>   Exclude ZMWs from output file,a comma-separated list of ID <br>
-j     <int>   Number of threads to use. [2] <br>
<br>
Arguments:<br>
input          Input file.<br>
output         Output file.<br>

