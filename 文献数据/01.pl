use strict;
use warnings;
use IO::Compress::Gzip;

sub Extract_01 {
    my ($input_folder, $output_folder, $sample_id, $ligation_file, $RT_file) = @_;
    
    my $Read1 = "$input_folder/$sample_id" . "_1.fastq.gz";
    my $Read2 = "$input_folder/$sample_id" . "_2.fastq.gz";
    my $output_file = "$output_folder/$sample_id" . "_barcode.fastq.gz";
    my $f3 = new IO::Compress::Gzip $output_file or die "Cannot open $output_file for writing: $!\n";
    
    open my $ligation_fh, '<', $ligation_file or die "Cannot open $ligation_file: $!\n";
    my @ligation_list = <$ligation_fh>;
    chomp @ligation_list;
    close $ligation_fh;
    
    open my $RT_fh, '<', $RT_file or die "Cannot open $RT_file: $!\n";
    my @RT_list = <$RT_fh>;
    chomp @RT_list;
    close $RT_fh;
    
    my %names_dict;
    my @names;
    my @first_line;
    open my $read1_fh, '-|', "pyfastx.Fastq -i $Read1 --build_index=false" or die "Cannot open $Read1: $!\n";
    while (my ($name, $seq, $qual) = map { chomp; $_ } <$read1_fh>) {
        my $tmp_lig = substr($seq, 0, 10);
        if (grep { $_ eq $tmp_lig } @ligation_list) {
            my $tmp_RT = substr($seq, length($tmp_lig) + 14, 10);
            if (grep { $_ eq $tmp_RT } @RT_list) {
                my $umi = substr($seq, length($tmp_lig) + 6, 8);
                push @names, $name;
                push @first_line, "\@$name,$tmp_lig$tmp_RT,$umi";
            }
        }
    }
    close $read1_fh;
    
    @names_dict{@names} = @first_line;
    
    my $count = 0;
    open my $read2_fh, '-|', "pyfastx.Fastq -i $Read2 --build_index=false" or die "Cannot open $Read2: $!\n";
    while (my ($name, $seq, $qual) = map { chomp; $_ } <$read2_fh>) {
        if (exists $names_dict{$name}) {
            $f3->print($names_dict{$name}, "\n");
            $f3->print("$seq\n+\n$qual\n");
            $count++;
            my $percent = ($count / scalar(@names_dict)) * 100;
            print "\r" . "#" x int($percent) . ">>> $percent% ($count/" . scalar(@names_dict) . ")";
        }
    }
    close $read2_fh;
    
    $f3->close();
}

# 调用函数
Extract_01('input_folder', 'output_folder', 'sample_id', 'ligation_file.txt', 'RT_file.txt');
