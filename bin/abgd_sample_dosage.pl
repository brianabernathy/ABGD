#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $vcf_file;
my $group_cons_alleles_file;
my $output_file;
my $output_depths;
my $output_pcts;
my $norm_ref_sample;
my $sample_covs_file;
my $exclude_other = 0;
my $keep_norm_ref_no_alignment;
my $help;

my @group_names = ();
my %group_cons_alleles = ();
my %sample_covs = ();
my %sample_name_indexes = ();
my %sample_index_names = ();
my $gt_ad_index;
my $norm_ref_sample_index;

parse_args();
parse_group_cons_alleles();

if (defined($norm_ref_sample)) {
	parse_sample_covs();
}

parse_vcf();

exit(0);


sub parse_group_cons_alleles {
	my $group_cons_alleles_fh = open_fh($group_cons_alleles_file, 'r');
	
	while (my $line = <$group_cons_alleles_fh>) {
		chomp($line);
		my ($chr, $pos, @group_cons_alleles) = split(/\t/, $line);

		if ($line =~ /^#/) {
			foreach my $group_name (@group_cons_alleles) {
				$group_name =~ s/_allele//;

				push(@group_names, $group_name);
			}

			next();
		}

		foreach my $index (0..$#group_cons_alleles) {
			$group_cons_alleles{$chr}{$pos}{$index} = $group_cons_alleles[$index];
		}
	}
	
	close($group_cons_alleles_fh);

	if (! @group_names) {
		error('valid header not found in group consensus alleles file: $group_cons_alleles_file');
	}

	return(0);
}


sub parse_sample_covs {
	my $sample_covs_fh = open_fh($sample_covs_file, 'r');

	while (my $line = <$sample_covs_fh>) {
		if ($line =~ /^#/) {
			next();
		}

		chomp($line);
		my ($sample, $cov) = split(/\t/, $line);
		$sample_covs{$sample} = $cov;
	}

	if (! exists($sample_covs{$norm_ref_sample})) {
		error("norm ref sample: $norm_ref_sample not found in sample coverages file: $sample_covs_file");
	}

	return(0);
}


sub parse_vcf {
	my $vcf_fh = open_fh($vcf_file, 'r');
	my $out_fh = open_fh($output_file, 'w');

	while (my $line = <$vcf_fh>) {
		chomp($line);

		if ($line =~ /^#/) {
			if ($line =~ /^#CHROM/) {
				my ($chr, $pos, $id, $ref, $alts, $qual, $filter, $info, $format, @samples) = split(/\t/, $line);

				if (defined($norm_ref_sample)) {
					check_missing_group_samples(\@samples);
				}

				foreach my $index (0..$#samples) {
					my $sample = $samples[$index];

					if (defined($norm_ref_sample) && $sample eq $norm_ref_sample) {
						$norm_ref_sample_index = $index;
					}

					$sample_name_indexes{$sample} = $index;
					$sample_index_names{$index} = $sample;
				}

				if (defined($norm_ref_sample) && ! defined($norm_ref_sample_index)) {
					error("norm ref sample: $norm_ref_sample not found in vcf file: $vcf_file");
				}

				print_header($out_fh, $line);
			}

			next();
		}

		if (! defined($gt_ad_index)) {
			my ($chr, $pos, $id, $ref, $alts, $qual, $filter, $info, $format, @samples) = split(/\t/, $line);

			$gt_ad_index = check_vcf_format($format);
		}

		print_record($out_fh, $line);
	}

	close($vcf_fh);
	close($out_fh);

	return(0);
}


sub print_record {
	my $fh = shift();
	my $line = shift();

	my ($chr, $pos, $id, $ref, $alts, $qual, $filter, $info, $format, @sample_gts) = split(/\t/, $line);

	if (! exists($group_cons_alleles{$chr}{$pos})) {
		return(0);
	}

	if ($info =~ /INDEL/) {
		return(0);
	}

	my @group_cons_alleles = ();
	my @alleles = ();
	push(@alleles, $ref);
	my @split_alt_alleles = split(/\,/, $alts);

	foreach my $split_alt_allele (@split_alt_alleles) {
		push(@alleles, $split_alt_allele);
	}

	my ($norm_ref_depth, $norm_ref_ads_ref, $norm_ref_ad_pcts_ref);

	if (defined($norm_ref_sample)) {
		($norm_ref_depth, $norm_ref_ads_ref, $norm_ref_ad_pcts_ref) = proc_vcf_gt($sample_gts[$norm_ref_sample_index]);

		if (! defined($keep_norm_ref_no_alignment) && $norm_ref_depth == 0) {
			return(0);
		}
	}

	print($fh join("\t", $chr, $pos));

	foreach my $group_index (0..$#group_names) {
		my $group_cons_allele = '*';

		if (exists($group_cons_alleles{$chr}{$pos}{$group_index})) {
			$group_cons_allele = $group_cons_alleles{$chr}{$pos}{$group_index};
		}

		push(@group_cons_alleles, $group_cons_allele);

		print($fh "\t$group_cons_allele");
	}

	foreach my $sample_index (0..$#sample_gts) {
		my $sample_name = $sample_index_names{$sample_index};
		my $sample_gt = $sample_gts[$sample_index];
		my ($depth, $ads_ref, $ad_pcts_ref) = proc_vcf_gt($sample_gt);
		my %other_allele_indexes = ();

		foreach my $allele_index (0..$#alleles) {
			$other_allele_indexes{$allele_index}++;
		}

		print($fh "\t$depth");

		foreach my $group_index (0..$#group_names) {
			my $group_cons_allele = $group_cons_alleles[$group_index];
			my $group_cons_gt_allele_index;

			foreach my $allele_index (0..$#alleles) {
				my $allele = $alleles[$allele_index];

				if ($allele eq $group_cons_allele) {
					$group_cons_gt_allele_index = $allele_index;

					last();
				}
			}

			if (defined($group_cons_gt_allele_index)) {
				delete $other_allele_indexes{$group_cons_gt_allele_index};
			}

			my $ad = 0;

			if (defined($group_cons_gt_allele_index) && exists($$ads_ref[$group_cons_gt_allele_index])) {
				$ad = $$ads_ref[$group_cons_gt_allele_index];
			}

			if (defined($output_depths)) {
				print($fh "\t$ad");
			}

			if (defined($output_pcts)) {
				my $ad_pct = 0;

				if (defined($group_cons_gt_allele_index) && exists($$ad_pcts_ref[$group_cons_gt_allele_index])) {
					$ad_pct = $$ad_pcts_ref[$group_cons_gt_allele_index];
				}

				print($fh "\t$ad_pct");
			}


			if (defined($norm_ref_sample) && exists($sample_covs{$sample_name})) {
				my $sample_group_ad = $ad;
				my $norm_ref_group_ad = 0;
				my $norm_val = '0.00';

				if (defined($group_cons_gt_allele_index) && exists($$norm_ref_ads_ref[$group_cons_gt_allele_index])) {
					$norm_ref_group_ad = $$norm_ref_ads_ref[$group_cons_gt_allele_index];
				}

				my $sample_mean_cov = $sample_covs{$sample_name};
				my $norm_ref_mean_cov = $sample_covs{$norm_ref_sample};

				if ($sample_mean_cov > 0 && $norm_ref_mean_cov > 0 && $norm_ref_group_ad > 0) {
					$norm_val = sprintf("%0.2f", (($sample_group_ad / $sample_mean_cov) / ($norm_ref_group_ad / $norm_ref_mean_cov)));
				}

				print($fh "\t$norm_val");
			}
		}

		if (! $exclude_other) {
			my $other_ad = 0;
			my $other_ad_pct = 0;
			my $norm_ref_other_ad = 0;

			foreach my $other_allele_index (keys %other_allele_indexes) {
				if (exists($$ads_ref[$other_allele_index])) {
					$other_ad += $$ads_ref[$other_allele_index];
				}

				if (exists($$ad_pcts_ref[$other_allele_index])) {
					$other_ad_pct += $$ad_pcts_ref[$other_allele_index];
				}

				if (exists($$norm_ref_ads_ref[$other_allele_index])) {
					$norm_ref_other_ad += $$norm_ref_ads_ref[$other_allele_index];
				}
			}

			if (defined($output_depths)) {
				print($fh "\t$other_ad");
			}

			if (defined($output_pcts)) {
				print($fh "\t$other_ad_pct");
			}

			if (defined($norm_ref_sample) && exists($sample_covs{$sample_name})) {
				my $other_norm_val = '0.00';
				my $sample_mean_cov = $sample_covs{$sample_name};
				my $norm_ref_mean_cov = $sample_covs{$norm_ref_sample};

				if ($sample_mean_cov > 0 && $norm_ref_mean_cov > 0 && $norm_ref_other_ad > 0) {
					$other_norm_val = sprintf("%0.2f", (($other_ad / $sample_mean_cov) / ($norm_ref_other_ad / $norm_ref_mean_cov)));
				}

				print($fh "\t$other_norm_val");
			}
		}
	}
				
	print($fh "\n");

	return(0);
}


sub print_header {
	my $fh = shift();
	my $line = shift();

	my ($chr, $pos, $id, $ref, $alts, $qual, $filter, $info, $format, @samples) = split(/\t/, $line);

	print($fh join("\t", '#chr', 'pos'));

	foreach my $group (@group_names) {
		print($fh "\t${group}_allele");
	}

	foreach my $sample (@samples) {
		print($fh "\t${sample}_total_depth");

		foreach my $group (@group_names, 'other') {
			if ($group eq 'other' && $exclude_other) {
				next();
			}

			if (defined($output_depths)) {
				print($fh "\t${sample}_${group}_depth");
			}

			if (defined($output_pcts)) {
				print($fh "\t${sample}_${group}_pct");
			}

			if (defined($norm_ref_sample) && exists($sample_covs{$sample})) {
				print($fh "\t${sample}_${group}_norm");
			}
		}
	}
				
	print($fh "\n");

	return(0);
}


sub check_missing_group_samples {
	my $samples_ref = shift();
	my %missing_cov_samples = ();
	my %missing_vcf_samples = ();

	foreach my $cov_sample (keys %sample_covs) {
		$missing_vcf_samples{$cov_sample}++;
	}

	foreach my $vcf_sample (@$samples_ref) {
		delete $missing_vcf_samples{$vcf_sample};
		$missing_cov_samples{$vcf_sample}++;
	}

	foreach my $cov_sample (keys %sample_covs) {
		delete $missing_cov_samples{$cov_sample};
	}

	if (%missing_cov_samples) {
		warning(join(' ', 'the following samples:', join(', ', sort keys %missing_cov_samples), "were not found in the sample coverages file: $sample_covs_file, but were found in the vcf file: $vcf_file, these samples will not be normalized"));
	}

	if (%missing_vcf_samples) {
		warning(join(' ', 'the following samples:', join(', ', sort keys %missing_vcf_samples), "were not found in the vcf file: $vcf_file, but were found in the sample coverages file: $sample_covs_file, these samples will not be processed"));
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


sub proc_vcf_gt {
	my $gt_field = shift();
	my $depth;
	my @ad_pcts = ();
	my @split_gt_fields = split(/\:/, $gt_field);
	my $ads = $split_gt_fields[$gt_ad_index];

	$ads =~ s/\./0/g;

	my @ads = split(/\,/, $ads);

	foreach my $ad (@ads) {
		$depth += $ad;
	}

	foreach my $ad (@ads) {
		my $pct = '0.00';

		if ($depth > 0) {
			$pct = sprintf("%0.2f", $ad / $depth * 100);
		}
		push(@ad_pcts, $pct);
	}

	return($depth, \@ads, \@ad_pcts);
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


sub warning {
	my $msg = shift();

	if (defined($msg)) {
		print(STDERR "warning: $msg\n");
	}

	return(0);
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
				'a|alleles=s' => \$group_cons_alleles_file,
				'o|out=s' => \$output_file,
				'd|depth' => \$output_depths,
				'p|pct' => \$output_pcts,
				'n|norm=s' => \$norm_ref_sample,
				'c|cov=s' => \$sample_covs_file,
				'e|excl' => \$exclude_other,
				'k|keep' => \$keep_norm_ref_no_alignment,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -output => '-', -verbose => 1);
	}

	if (! defined($vcf_file)) {
		arg_error('vcf file required');
	}

	if (! defined($group_cons_alleles_file)) {
		arg_error('group consensus alleles file required');
	}

	if (! defined($output_file)) {
		arg_error('output file required');
	}

	if (defined($norm_ref_sample) && ! defined($sample_covs_file)) {
		arg_error('sample genome coverage file required for normalization');
	}

	if (defined($sample_covs_file) && ! defined($norm_ref_sample)) {
		arg_error('normalization reference sample required for normalization');
	}

	if (! defined($output_depths) && ! defined($output_pcts) && ! defined($norm_ref_sample)) {
		arg_error('one or more of -d/--depth, -p/--pct, or -n/--norm options required');
	}

	return(0);
}


__END__

=head1 NAME

abgd_sample_dosage.pl

=head1 SYNOPSIS

abgd_sample_dosage.pl [options]

=head1 DESCRIPTION

abgd_sample_dosage.pl generates sample group dosage percentages and normalized values

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 -v --vcf       samples vcf file (required)
                  FORMAT field must contain the AD tag

 -a --alleles   group consensus alleles file (required)
                  generated by abgd_consensus_alleles.pl

                  tab-separated value file with the following fields
                  #chr	pos	group1_cons_allele	group2_cons_allele	groupX_cons_allele

 -o --out       output file (required)

 one or more of -d/--depth, -p/--pct, or -n/--norm options required

 -d --depth     output sample allele depths

 -p --pct       output sample allele percentages

 -n --norm      normalization reference sample (output sample norm values)

 -c --cov       sample genome coverage file (required for normalization)
                 only samples included in the file will be normalized

                  tab-separated value file with the following fields
                  #sample	coverage

                  example:
                  norm.ref	100.5
                  sample1	30.5
                  sample2	30.2
                  sample3	30.2

 -e --excl      exclude 'other' group stats for samples

 -k --keep      output variant sites where the norm ref sample has no alignments
                  default: skip these records

 -h --help      display help menu

=cut
