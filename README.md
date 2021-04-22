# ccsx
efficient tool for generating  circular consensus sequences (ccs) from subreads.
Ccsx works with both x86-64 CPUs and ARM CPUs supporting the NEON instruction sets. 
Efficiency are ensured through the use of bsalign(https://github.com/ruanjue/bsalign.git)

# install
git clone https://github.com/110allan/ccsx.git
cd ccsx
make

# usage
Usage  : ccsx  [options] <INPUT> <OUTPUT>
Generate circular consensus sequences (ccs) from subreads.

Options:
-h             Output this help 
-v             debug 
-m     <int>   Minimum length of subreads to use for generating CCS. [10] 
-M     <int>   Maximum length of subreads to use for generating CCS. [50000] 
-c     <int>   Minimum number of subreads required to generate CCS. [3] 
-j     <int>   Number of threads to use, 0 means autodetection. [2] 

Arguments:
input          Input bam file.
output         Output fasta file.

