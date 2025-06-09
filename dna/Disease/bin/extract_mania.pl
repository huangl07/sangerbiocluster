#!/usr/bin/perl
my $start=0;
while($l=<STDIN>){
    chomp $l;
    if($l =~ /^Gene 1\tGene 2/){
        $start = 1;
        print "Gene1\tGene2\tWeight\tType\tSource\n";
        next;
    };
    if($start == 0){
        next;
    }else{
        if( length($l) == 0 ){
            $start = 0;
            next;
        }else{
            print $l,"\n";
        };
    };
}
