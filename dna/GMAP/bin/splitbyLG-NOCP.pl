#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$dOut,$fLG,$fKey,$type,$Pos,$chr);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"l:s"=>\$fLG,
				"d:s"=>\$dOut,
                "chr:s"=>\$chr,
				) or &USAGE;
&USAGE unless ($fIn and $dOut and $fLG);
mkdir $dOut if (!-d $dOut);
my $split="-";
if($chr){$split="-";}
$dOut=ABSOLUTE_DIR($dOut);
open In,$fIn;
my %indi;
my $head=`head -n 1 $fIn`;
chomp $head;
while(<In>){
	chomp;
	next if ($_ eq ""||/^$/ ||/MarkerID/);
	my ($id,@info)=split(/\s+/,$_);
	my $info=join("\t",@info);
	if($id =~ /X/){
		$id=~s/X/h/g;
	}
	my $scaid=(split(/$split/,$id))[0];
	push @{$indi{$scaid}},join("\t",$id,$info); 
}
close In;
open In,$fLG;
$/=">";
open L,">$dOut/lg.list";
while(<In>){
	chomp;
	next if ($_ eq ""||/^$/ ||/MarkerID/);
	my ($id,@info)=split(/\n/,$_);
	my $lgid=(split(/\s+/,$id))[0];
	my $l=join("\t",@info);
	my @scas=split(/\s+/,$l);
	print L "$lgid\t$dOut/$lgid.marker\n";
	open Out,">$dOut/$lgid.marker";
	print Out $head,"\n";
	foreach my $sca(@scas){
        if(!defined $indi{$sca}){
            print $sca,"\n";
            die;
        }
		print Out join("\n",@{$indi{$sca}}),"\n";
	}
	close Out;
}
close  In;
close L;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
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

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: huangl <long.huang\@majorbio.com> 

Options:
  -help			USAGE,
  -i	genotype file forced
  -l	linkage lg file
  -d	output dir
  -t	population type
  
   
USAGE
	print $usage;
	exit;
}

