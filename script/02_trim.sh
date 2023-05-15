#!/usr/bin/bash
input_folder="/root/wangje/Project/吴霞/Data/01_Exrtract"
output_folder="/root/wangje/Project/吴霞/Data/02_trim/"
sample_prefix='BC230502'
core=6

tempfifo="my_temp_fifo"
mkfifo ${tempfifo}
rm -rf ${tempfifo}
exec ${core}<>${tempfifo}

for i in {1..96}
do
{
    trim_galore $input_folder/$sample_prefix-${i}_R2_barcode.fq.gz  -a AAAAAAAA --three_prime_clip_R1 1 -o $output_folder
}& 
done  >&${core} 
wait
exec ${core}>&-