#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $vcf_file;
my $output_file;
my $summary_file;
my $group_samples_file;
my $print_all_pos = 0;
my $help;

my %group_samples = ();
my %non_group_samples = ();
my %sample_min_depths = ();
my %sample_min_allele_pcts = ();
my %summary_counts = ();
my %sample_name_indexes = ();
my %sample_index_names = ();
my $gt_ad_index;

parse_args();
parse_group_samples();
parse_vcf();

if (defined($summary_file)) {
	print_summary();
}

exit(0);


sub parse_group_samples {
	my $group_samples_fh = open_fh($group_samples_file, 'r');
	
	while (my $line = <$group_samples_fh>) {
		if ($line =~ /^#/) {
			next();
		}

		chomp($line);
	
		my ($group, $sample, $min_depth, $min_allele_pct) = split(/\t/, $line);

		if (! defined($group) || ! defined($sample) || ! defined($min_depth) || ! defined($min_allele_pct)) {
			error("group samples file: $group_samples_file contains an invalid record\n\t$line");
		}

		if ($min_depth !~ /^[+-]?\d*\.?\d+$/) {
			error("group samples file: $group_samples_file min depth must be numeric\n\tmin depth: $min_depth\n\t$line");
		}

		if ($min_allele_pct !~ /^[+-]?\d*\.?\d+$/) {
			error("group samples file: $group_samples_file min allele pct must be numeric\n\tmin allele_pct: $min_allele_pct\n\t$line");
		}

		if ($min_depth < 0) {
			$min_depth = 0;
		}

		if ($min_allele_pct < 0) {
			$min_allele_pct = 0;
		}

		if ($min_allele_pct > 100) {
			$min_allele_pct = 100;
		}
	
		$group_samples{$group}{$sample}++;
		$sample_min_depths{$sample} = $min_depth;
		$sample_min_allele_pcts{$sample} = $min_allele_pct;
	}
	
	close($group_samples_fh);

	return(0);
}


sub parse_vcf {
	my $vcf_fh = open_fh($vcf_file, 'r');
	my $out_fh = open_fh($output_file, 'w');

	while (my $line = <$vcf_fh>) {
		chomp($line);

		if ($line =~ /^#/) {
			$summary_counts{'comments'}++;

			if ($line =~ /^#CHROM/) {
				my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split(/\t/, $line);

				foreach my $index (0..$#samples) {
					my $sample = $samples[$index];

					$sample_name_indexes{$sample} = $index;
					$sample_index_names{$index} = $sample;
				}

				check_missing_group_samples(\@samples);

				# print header
				print($out_fh "#chr\tpos");

				foreach my $group (sort keys %group_samples) {
					print($out_fh "\t${group}_allele");
				}

				print($out_fh "\n");
			}

			next();
		}

		$summary_counts{'vars'}++;

		my ($chr, $pos, $id, $ref_allele, $alt_alleles, $qual, $filter, $info, $format, @samples) = split(/\t/, $line);

		if ($info =~ /INDEL/) {
			$summary_counts{'indels'}++;

			next();
		}

		if (! defined($gt_ad_index)) {
			$gt_ad_index = check_vcf_format($format);
		}

		my @alleles = ();
	
		push(@alleles, $ref_allele);
	
		my @split_alt_alleles = split(/\,/, $alt_alleles);
	
		foreach my $split_alt_allele (@split_alt_alleles) {
			push(@alleles, $split_alt_allele);
		}


		my %group_sample_pri_alleles = ();
		my %group_cons_alleles = ();
		my %failed_check_samples = ();
		my $record_pass = 1;
	
		foreach my $group (sort keys %group_samples) {
			foreach my $sample (sort keys %{$group_samples{$group}}) {
				my $sample_index = $sample_name_indexes{$sample};
				my ($depth, $pri_allele_index, @allele_pcts) = proc_vcf_gt($samples[$sample_index]);
	
				if (! defined($pri_allele_index)) {
					$summary_counts{"undef primary allele: $sample"}++;
					$failed_check_samples{$sample}++;
					$record_pass = 0;
				}

				else {
					$group_sample_pri_alleles{$group}{$pri_allele_index}{$sample}++;

					if ($allele_pcts[$pri_allele_index] < $sample_min_allele_pcts{$sample}) {
						$summary_counts{"low primary allele%: $sample"}++;
						$failed_check_samples{$sample}++;
						$record_pass = 0;
					}
				}
	
				if ($depth < $sample_min_depths{$sample}) {
					$summary_counts{"low depth: $sample"}++;
					$failed_check_samples{$sample}++;
					$record_pass = 0;
				}
			}

			my $group_pri_allele_count = keys %{$group_sample_pri_alleles{$group}};

			if ($group_pri_allele_count > 1) {
				$summary_counts{"no group consensus allele: $group"}++;
				$record_pass = 0;
			}

			else {
				foreach my $pri_allele_index (keys %{$group_sample_pri_alleles{$group}}) {
					$group_cons_alleles{$group} = $pri_allele_index;
				}
			}
		}
	
		if ($record_pass == 0) {
			$summary_counts{'record fail'}++;

			next();
		}


		if ($print_all_pos) {
			$summary_counts{'record pass'}++;
		}

		else {
			my %cons_alleles = ();

			foreach my $group (keys %group_samples) {
				$cons_alleles{$group_cons_alleles{$group}}++;
			}

			my $group_count = keys %group_samples;
			my $cons_allele_count = keys %cons_alleles;

			if ($group_count == $cons_allele_count) {
				$summary_counts{'record pass'}++;
			}

			else {
				$summary_counts{'shared consensus alleles'}++;
				$summary_counts{'record fail'}++;

				next();
			}
		}


		print($out_fh "$chr\t$pos");

		foreach my $group (sort keys %group_samples) {
			print($out_fh "\t$alleles[$group_cons_alleles{$group}]");
		}

		print($out_fh "\n");
	}

	close($vcf_fh);
	close($out_fh);

	return(0);
}


sub check_missing_group_samples {
	my $samples_ref = shift();
	my %missing_samples = ();

	foreach my $group (keys %group_samples) {
		foreach my $sample (keys %{$group_samples{$group}}) {
			$missing_samples{$sample}++;
		}
	}

	foreach my $sample (@$samples_ref) {
		delete $missing_samples{$sample};
	}

	if (%missing_samples) {
		error(join(' ', 'the following samples:', join(', ', sort keys %missing_samples), "were not found in vcf file: $vcf_file"));
	}

	return(0);
}


sub check_vcf_format {
	my $format = shift();
	my @gt_tags = split(':', $format);

	foreach my $gt_index (0..$#gt_tags) {
		if ($gt_tags[$gt_index] eq 'AD') {
			return($gt_index);
		}
	}

	error("vcf FORMAT field does not include required AD tag");

	return(-1);
}


sub print_summary {
	my $summary_fh = open_fh($summary_file, 'w');

	foreach my $record (sort keys %summary_counts) {
		my $count = 0;

		if (exists($summary_counts{$record})) {
			$count = $summary_counts{$record};
		}

		print($summary_fh "$record: $count\n");
	}

	close($summary_fh);

	return(0);
}


sub proc_vcf_gt {
	my $gt_field = shift();
	my $depth;
	my $pri_allele_index = 0;
	my @ad_pct = ();
	my @split_gt_fields = split(/\:/, $gt_field);
	my $ads = $split_gt_fields[$gt_ad_index];
	my @ads = split(/\,/, $ads);

	foreach my $ad (@ads) {
		$depth += $ad;
	}

	foreach my $index (0..$#ads) {
		my $ad = $ads[$index];
		my $pct = 0;

		if ($depth > 0) {
			$pct = sprintf("%0.2f", $ad / $depth * 100);
		}

		push(@ad_pct, $pct);

		if ($pct > $ad_pct[$pri_allele_index]) {
			$pri_allele_index = $index;
		}
	}

	return($depth, $pri_allele_index, @ad_pct);
}


sub open_fh {
	my $file = shift();
	my $mode = shift();
	my $fh;

	if (! defined($file)) {
		error('open_fh file not defined');
	}

	if (! defined($mode)) {
		error("open_fh mode not defined for file: $file");
	}

	if ($mode ne 'r' && $mode ne 'w') {
		error("open_fh mode for file: $file must be either 'r' or 'w'");
	}

	if ($file =~ /\.gz$/) {
		if ($mode eq 'r') {
			open($fh, '-|', "gzip -dc $file") or error("can't read $file: $!");
		}

		elsif ($mode eq 'w') {
			open($fh, "| gzip -c > $file") or error("can't write $file: $!");
		}
	}

	elsif ($file =~ /\.bz2$/) {
		if ($mode eq 'r') {
			open($fh, '-|', "bzip2 -dc $file") or error("can't read $file: $!");
		}

		elsif ($mode eq 'w') {
			open($fh, "| bzip2 -c > $file") or error("can't write $file: $!");
		}
	}

	else {
		if ($mode eq 'r') {
			open($fh, '<', $file) or error("can't read $file: $!");
		}

		elsif ($mode eq 'w') {
			open($fh, '>', $file) or error("can't write $file: $!");
		}
	}

	return($fh);
}


sub error {
	my $msg = shift();

	if (defined($msg)) {
		print(STDERR "error: $msg\n");
	}

	exit(0);
}


sub arg_error {
	my $msg = shift();

	if (defined($msg)) {
		print(STDERR "error: $msg\n\n");
	}

	print("use option -h or --help to display help menu\n");

	exit(0);
}


sub parse_args {
	if ($#ARGV == -1) {
		pod2usage(-exitval => 0, -output => '-', -verbose => 1);
	}

	GetOptions ('v|vcf=s' => \$vcf_file,
				'g|group=s' => \$group_samples_file,
				'o|out=s' => \$output_file,
				's|sum=s' => \$summary_file,
				'a|all' => \$print_all_pos,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -output => '-', -verbose => 1);
	}

	if (! defined($vcf_file)) {
		arg_error('vcf file required');
	}

	if (! defined($group_samples_file)) {
		arg_error('group samples file required');
	}

	if (! defined($output_file)) {
		arg_error('output file required');
	}

	return(0);
}


__END__

=head1 NAME

abgd_consensus_alleles.pl

=head1 SYNOPSIS

abgd_consensus_alleles.pl [options]

=head1 DESCRIPTION

abgd_consensus_alleles.pl generates sample group consensus alleles from vcf variant files

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 -v --vcf     vcf file (required)
                FORMAT field must contain the AD tag

 -g --group   group samples file (required)
                tab-separated value file with the following fields
                #group  sample  min_depth   min_allele_pct

                example:
                a_cons  dur1    10    90
                a_cons  dur2    10    95
                a_cons  dur3    20    95
                b_cons  ipa1    10    95
                b_cons  ipa2    20    90
                ...

 -o --out     output file (required)

 -s --sum     summary file (default: disabled)

 -a --all     output all variant positions with valid consensus alleles, even if
                the consensus alleles are the same
                note: variant positions with non-consensus alleles are skipped
                default: only variant positions with differing alleles are output

 -h --help    display help menu

=cut
