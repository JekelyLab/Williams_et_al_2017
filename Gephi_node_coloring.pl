# This program colors the nodes of a gephi file based on single cell transcriptome data
#!/usr/bin/perl -w
#Gaspar Jekely Oct 2015
	
use strict;
use List::Util qw( min max );



	
my $usage = "Gephi_node_coloring.pl infile_gephi infile_gene_expression\n";
my $gephi = shift or die $usage;
my $expr = shift or die $usage;

open (GEPHI, "<$gephi") or die "Cannot open infile!\n";
open (EXPR, "<$expr") or die "Cannot open infile!\n";
open (GEPHI_OUT, ">$gephi"."_out.gexf") or die "Cannot open outfile!\n";

#open (my $in, '<', $expr) or die "Cannot open infile!\n";

my %hash={};
my @levels=();

foreach my $line (<EXPR>)
	{
	if ($line =~ /ERR/){
	chomp $line;
	my @fields = split (/\,/, $line);
	#print  $fields[0], "\n";	
	push (@levels, $fields[1]), 
	$hash{ $fields[0] } = $fields[1];
	#print $hash{ $fields[0] }; 
	}
	}
#print @levels;
my $max = max @levels;
#print "max= ", $max, "\n";

my @RGB=();
my %RGB={};
my $RGB1=255; my $RGB2=255; my $RGB3=255;

while (my $line = <GEPHI>)
	{
	if ($line =~ /node id/){
	print GEPHI_OUT $line;
	my @fields = split (/\"/, $line);
	$RGB1 = 255-$hash{ $fields[1] }/$max*255;
	}
	elsif ($line =~ /color/) {
	#print $line;
	print GEPHI_OUT "        <viz:color r=","\"", int($RGB1), "\" g=", "\"", int($RGB2), "\" b=\"", int($RGB3),"\"></viz:color>", "\n";
	#print "RGB: ", $RGB1, " ", $RGB2, " ", $RGB3, "\n";
	}
	else {print GEPHI_OUT $line;}
	}
	
