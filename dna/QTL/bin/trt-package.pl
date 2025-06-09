#!/usr/bin/env perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
my @csv=glob("$fIn/*.csv");
my %trait;
my %YEAR;
my %LOC;
my %REP;
my %trtname;
my %sampleID;
foreach my $csv(@csv){
    `dos2unix $csv`;
    print "#######################";
    print $csv,"\n";
    print "########################";
    my $format=0;
    open In,$csv;
    my $nsamples=0;
    my @trait;
    my @year;
    my @locale;
    my @rep;
    while(<In>){
        chomp;
        s/\r//g;
        s/^\N{BOM}//g;
        next if ($_ eq ""||/^$/);
        $nsamples++;
        if($nsamples == 1){
            s/\,$//g;
            my (undef,@trt)=split(/\,/,$_);
            foreach my $trt(@trt){
                if($trt eq "year" || $trt eq "locale" || $trt eq "rep" || $trt eq "YEAR"||$trt eq "LOC"||$trt eq "REP"){
                    $format=1;
                }
                $trt=~s/\s+$//g;
                $trt=~s/^\s+//g;
                my $trtname=join("_",split(/\s+/,$trt));
                $trtname = "YEAR" if($trtname eq "year");
                $trtname = "LOC" if($trtname eq "loc");
                $trtname = "REP" if($trtname eq "rep");
                push @trait,$trtname;
            }
            next;
        }else{
            if($format == 1){
                s/\,$//g;
                my ($id,@value)=split(/\,/,$_);
                for(my $i=0;$i<@value;$i++){
                    $value[$i] = "NA" if($value[$i] eq "");
                }
                my $year="NA";
                my $local="NA";
                my $rep="NA";
                if($id eq ""){next;}
                for (my $i=0;$i<@value;$i++){
                    next if($i >= scalar @trait);
                    if($trait[$i] eq "YEAR"){$year=$value[$i];next;};
                    if($trait[$i] eq "LOC"){$local=$value[$i];next;}
                    if($trait[$i] eq "YEAR"){$rep=$value[$i];next;};
                    $trtname{$trait[$i]}=1;
                    $value[$i]=~s/\s+//g;
                    $value[$i]="NA" if($value[$i] eq "");
                    $year=~s/\.00//g;
                    $YEAR{$year}=1;
                    $LOC{$local}=1;
                    $REP{$rep}=1;
                    $trait{$id}{$year}{$local}{$rep}{$trait[$i]}=$value[$i];
                    $sampleID{$id}=1;
                }
            }else{
                s/\,$//g;
                my ($id,@value)=split(",",$_);
                for(my $i=0;$i<@value;$i++){
                    $value[$i] = "NA" if($value[$i] eq "");
                }
                if(!defined $id||$id eq ""){next;}
                if($id eq "year"){
                    @year=@value;
                    next;
                }
                if($id eq "locale"){
                    @locale=@value;
                    next;
                }
                if($id eq "rep"){
                    @rep=@value;
                    next;
                }
                for(my $i=0;$i<@value;$i++){
                    $year[$i]||="NA";
                    $locale[$i]||="NA";
                    $rep[$i]||="NA";
                    $value[$i]=~s/\s+//g;
                    $value[$i]="NA" if($value[$i] eq "");
                    $trtname{$trait[$i]}=1;
                    $year[$i]=~s/\.00//g;
                    $YEAR{$year[$i]}=1;
                    $LOC{$locale[$i]}=1;
                    $REP{$rep[$i]}=1;
                    $trait{$id}{$year[$i]}{$locale[$i]}{$rep[$i]}{$trait[$i]}=$value[$i];
                    $sampleID{$id}=1;
                }
            }
        }
    }
    close In;
}

open Out,">$fOut";
print Out join("\t","sampleID","YEAR","LOC","REP",sort keys %trtname),"\n";
foreach my $sampleID(sort keys %sampleID){
    #my @out;
    foreach my $year (keys %YEAR){
        foreach my $local(keys %LOC){
            foreach my $rep(keys %REP){
                my @out;
                foreach my $trt(sort keys %trtname){
                    $trait{$sampleID}{$year}{$local}{$rep}{$trt}||="NA";
                    push @out,$trait{$sampleID}{$year}{$local}{$rep}{$trt};
                }
                print Out join("\t",$sampleID,$year,$local,$rep,@out),"\n";
            }
        }
    }
}
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
  -h         Help

USAGE
        print $usage;
        exit;
}
