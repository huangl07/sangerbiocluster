#!/usr/bin/perl
my $vep=0;
use Data::Dumper;
while($l=<STDIN>){
    chomp $l;
    if($l =~ /^#/){
        $l =~ s/^# //g;
        $l =~ s/\[\d+\]//g;
        @h=split(/\t/,$l);
		$output=0;

    } else{
        @a=split(/\t/,$l);
        @b=split(/\,/,$a[-1]);
		$type=scalar split(/\|/,$b[0],-1);
		#print $type,",";
	if($output == 0){
		if($type == 16) {
			print join("\t",@h[0..$#h-1],"allele\tfunctional_region\tPutative_impact\tgene_name\tgene_id\tFeature_type\tFeature_id\tTranscript_biotype\tExon_rank\tHGVS_c\tHGVS_p\tcDNA_position\tCDS_position\tProtein_position\tDistance_to_feature\n");
			$output=1;
			$vep=0;
		}else{
			print join("\t",@h[0..$#h-1],"allele\tfunctional_region\tPutative_impact\tgene_name\tgene_id\tFeature_type\tFeature_id\tTranscript_biotype\tExon_rank\tIntron_rank\tHGVS_c\tHGVS_p\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tDistance_to_feature\n");
			$output=1;
			$vep=1;
		}
	}
        @c=split(/\,/,$a[3]);
        foreach $b(@b){
            $d=(split(/\|/,$b,-1))[0];
            $b=~s/(\d+)\/(\d+)/($1\/$2)/g;
            @t=split(/\|/,$b,-1);
            for my $i (0 .. $#t) {
                if ($t[$i] eq "") {
                    $t[$i] = '--';
                }
            };
			$keyout=$d;
			if($vep == 1){
                if(length($d)>length($a[2])){
                    $keyout=~s/$a[2]//;
                }elsif(length($d) < length($a[2])){
                    $keyout="-";
  		        }
				print join("\t",@a[0..$#a-1],$keyout,@t[1..$#t-7]),"\n";
            } else {
            #print $d,",",join("\t",@t[0..$#t-1]),",";
            print join("\t",@a[0..$#a-1],$keyout,@t[1..$#t-1]),"\n";
            }
        }
    }
}
