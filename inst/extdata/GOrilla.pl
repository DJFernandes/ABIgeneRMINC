#!/usr/bin/perl

# Adapted from script written by Aaron J. Mackey
#   https://gist.github.com/amackey/5286965

use strict;
use warnings;

use WWW::Mechanize;
use Getopt::Long;

# to install modules follow these commands
# perl -MCPAN -e shell
# install module [NAME] 

my $GOrillaURL = "http://cbl-gorilla.cs.technion.ac.il/";

my @organisms = qw(ARABIDOPSIS_THALIANA
		   SACCHAROMYCES_CEREVISIAE
		   CAENORHABDITIS_ELEGANS
		   DROSOPHILA_MELANOGASTER
		   DANIO_RERIO
		   HOMO_SAPIENS
		   MUS_MUSCULUS
		   RATTUS_NORVEGICUS
		 );
my %organisms; @organisms{@organisms} = (1) x @organisms;
my $organism = "MUS_MUSCULUS";

my @runmodes = qw(mhg hg);
my %runmodes; @runmodes{@runmodes} = (1) x @runmodes;
my $runmode = "hg";

my @ontologies = qw(proc func comp all);
my %ontologies; @ontologies{@ontologies} = (1) x @ontologies;
my $ontology = "proc";

my $pvalue = "0.001";
my $name = "";
my $email = "";
my $includedups = 0;
my $revigo = 1;
my $fast = 1;
# my ($targets, $background);
my $targets = "target.txt";
my $background = "background.txt";

my $result = GetOptions("organism=s" => \$organism,
		        "runmode=s" => \$runmode,
		        "targets=s" => \$targets,
		        "background=s" => \$background,
		        "ontology=s" => \$ontology,
			"pvalue=f" => \$pvalue,
			"name=s" => \$name,
			"email=s" => \$email,
			"includedups!" => \$includedups,
			"fast!" => \$fast,
		       );

die "No such organism $organism\n" unless $organisms{$organism};
die "No such runmode $runmode\n" unless $runmodes{$runmode};
die "No such ontology $ontology\n" unless $ontologies{$ontology};

die "Must supply both target and background files with runmode hg\n"
  unless ($runmode eq "mhg" || ($targets && $background));

die "Must supply target file with runmode mhg\n"
  unless ($runmode eq "hg" || $targets);

my $mech = WWW::Mechanize->new();

$mech->get($GOrillaURL);

$mech->form_name("gorilla");

$mech->select("species" => $organism);
$mech->set_fields("run_mode" => $runmode);
$mech->set_fields("target_file_name" => $targets);
if ($runmode eq "hg") {
#  $mech->set_file("background_file_name" => $background);
  $mech->set_fields("background_file_name" => $background);
}
$mech->set_fields("db" => $ontology);
$mech->select("pvalue_thresh" => $pvalue);
$mech->set_fields("analysis_name" => $name);
$mech->set_fields("user_email" => $email);
$mech->set_fields("output_excel" => 1);
$mech->set_fields("output_unresolved" => $includedups);
$mech->set_fields("output_revigo" => $revigo);
$mech->set_fields("fast_mode" => $fast);

$mech->click("run_gogo_button");

my $res = $mech->response();
my $base =  $res->base();
my ($id) = $base =~ m/id=(.*)/;

print "Results can be found at:
http://cbl-gorilla.cs.technion.ac.il/GOrilla/${id}/GOResults.html\n";

