package InterProScan_SDK::InterProScan_SDKImpl;
use strict;
use Bio::KBase::Exceptions;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org 
our $VERSION = "0.1.0";

=head1 NAME

InterProScan_SDK

=head1 DESCRIPTION

A KBase module: InterProScan_SDK

=cut

#BEGIN_HEADER
use Bio::KBase::AuthToken;
use Bio::KBase::workspace::Client;
use GenomeAnnotationAPI::GenomeAnnotationAPIClient;
use Config::IniFiles;
use Data::Dumper;
use File::Path;
use DateTime;
use JSON::XS;
use Cwd;
our $currentcontext;

#Initialization function for call
sub util_initialize_call {
	my ($self,$params,$ctx) = @_;
	$currentcontext = $ctx;
	return $params;
}

sub util_currentuser {
	return $currentcontext->user_id();
}

sub util_token {
	return $currentcontext->token();
}

sub util_provenance {
	return $currentcontext->provenance();
}

sub util_error {
	my ($self,$message) = @_;
	Carp::confess($message);
}

sub util_to_fasta {
	my ($self,$seqName, $seq, $len) = @_;
	# default to 80 characters of sequence per line
	$len = 80 unless $len;
	my $formatted_seq = ">$seqName\n";
	while (my $chunk = substr($seq, 0, $len, "")) {
		$formatted_seq .= "$chunk\n";
	}
	return $formatted_seq;
}

sub util_scratchdir {
	my ($self) = @_;
	return $self->{scratch};
}

sub util_validate_args {
	my ($self,$args,$mandatoryArguments,$optionalArguments,$substitutions) = @_;
	print "Retrieving input parameters.\n";
	if (!defined($args)) {
	    $args = {};
	}
	if (ref($args) ne "HASH") {
		$self->util_error("Arguments not hash");	
	}
	if (defined($substitutions) && ref($substitutions) eq "HASH") {
		foreach my $original (keys(%{$substitutions})) {
			$args->{$original} = $args->{$substitutions->{$original}};
		}
	}
	if (defined($mandatoryArguments)) {
		for (my $i=0; $i < @{$mandatoryArguments}; $i++) {
			if (!defined($args->{$mandatoryArguments->[$i]})) {
				push(@{$args->{_error}},$mandatoryArguments->[$i]);
			}
		}
	}
	$self->util_error("Mandatory arguments ".join("; ",@{$args->{_error}})." missing.") if (defined($args->{_error}));
	if (defined($optionalArguments)) {
		foreach my $argument (keys(%{$optionalArguments})) {
			if (!defined($args->{$argument})) {
				$args->{$argument} = $optionalArguments->{$argument};
			}
		}	
	}
	return $args;
}

sub util_configure_ws_id {
	my ($self,$ws,$id) = @_;
	my $input = {};
 	if ($ws =~ m/^\d+$/) {
 		$input->{wsid} = $ws;
	} else {
		$input->{workspace} = $ws;
	}
	if ($id =~ m/^\d+$/) {
		$input->{objid} = $id;
	} else {
		$input->{name} = $id;
	}
	return $input;
}

sub util_runexecutable {
	my ($self,$Command) = @_;
	my $OutputArray;
	push(@{$OutputArray},`$Command`);
	return $OutputArray;
}

sub util_from_json {
	my ($self,$data) = @_;
    if (!defined($data)) {
    	die "Data undefined!";
    }
    return decode_json $data;
}

sub util_get_genome {
	my ($self,$ref) = @_;
	my $output = $self->util_ga_client()->get_genome_v1({
		genomes => [{
			"ref" => $ref
		}],
		ignore_errors => 1,
		no_data => 0,
		no_metadata => 1
	});
	return $output->{genomes}->[0]->{data};
}

sub util_ga_client {
	my ($self,$input) = @_;
	if (!defined($self->{_gaclient})) {
		$self->{_gaclient} = new GenomeAnnotationAPI::GenomeAnnotationAPIClient($ENV{ SDK_CALLBACK_URL });
	}
	return $self->{_gaclient};
}

sub func_annotate_genome_with_interpro_pipeline {
	my ($self,$params) = @_;
    $params = $self->util_validate_args($params,["workspace","genome_id","genome_output_id"],{
    	genome_workspace => $params->{workspace},
    });
    my $annofunc = "Annotate Genome with InterPro Pipeline";
  	my $timestamp = DateTime->now()->datetime();
    #Step 1: Get genome from workspace
    my $wsClient = Bio::KBase::workspace::Client->new($self->{'workspace-url'},token=>$self->util_token());
    my $genome = $self->util_get_genome($params->{genome_workspace}."/".$params->{genome_id});
    #Step 2: Print protein FASTA file
    File::Path::mkpath $self->util_scratchdir();
    my $filename = $self->util_scratchdir()."/protein.fa";
    open ( my $fh, ">", $filename) || $self->util_error("Failure to open file: $filename, $!");
    my $genehash = {};
    foreach my $gene (@{$genome->{features}}) {
    	if (defined($gene->{protein_translation}) && length($gene->{protein_translation}) > 0) {
    		$genehash->{$gene->{id}} = $gene;
    		print $fh $self->util_to_fasta($gene->{id}, $gene->{protein_translation});
    	}
    }
    close($fh);
    #Step 3: Run interpro
    my $orig_cwd = cwd;
    chdir $self->util_scratchdir();
    system("/data/interproscan/interproscan.sh");
    system("/data/interproscan/interproscan.sh -i ".$self->util_scratchdir()."/protein.fa -f tsv -o ".$self->util_scratchdir()."/protein.tsv --disable-precalc -iprscan -iprlookup -hm");
    chdir $orig_cwd;
    #Step 4: Parsing interpro results
    $filename = $self->util_scratchdir()."/protein.tsv";
    my $numftr = 0;
    my $numdomains = 0;
    my $domainhash = {};
    my $ftrhash = {};
    open ( my $fh, "<", $filename) || $self->util_error("Failure to open file: $filename, $!");
    while (my $line = <$fh>) {
    	chomp($line);
    	my $array = [split(/\t/,$line)];
    	if (@{$array} < 13 || !defined($genehash->{$array->[0]})) {
    		next;
    	}
    	if (!defined($ftrhash->{$array->[0]})) {
    		$ftrhash->{$array->[0]} = 1;
    		$numftr++;
    	}
    	if (!defined($domainhash->{$array->[11]})) {
    		$domainhash->{$array->[11]} = 1;
    		$numdomains++;
    	}
    	my $ftr = $genehash->{$array->[0]};
    	if (!defined($ftr->{ontology_terms}->{InterPro}->{$array->[11]})) {
			$ftr->{ontology_terms}->{InterPro}->{$array->[11]} = {
				 evidence => [],
				 id => $array->[11],
				 term_name => $array->[12],
				 ontology_ref => "7537/36/2",#TODO: Need to make an interpro ontology and then set this ref to that
				 term_lineage => [],
			};
		}
		my $found = 0;
		for (my $k=0; $k < @{$ftr->{ontology_terms}->{InterPro}->{$array->[11]}->{evidence}}; $k++) {
			if ($ftr->{ontology_terms}->{InterPro}->{$array->[11]}->{evidence}->[$k]->{method} eq $annofunc) {
				$ftr->{ontology_terms}->{InterPro}->{$array->[11]}->{evidence}->[$k]->{timestamp} = $timestamp;
				$ftr->{ontology_terms}->{InterPro}->{$array->[11]}->{evidence}->[$k]->{method_version} = $VERSION;
				$ftr->{ontology_terms}->{InterPro}->{$array->[11]}->{evidence}->[$k]->{alignment_evidence} = [[$array->[6]+0,$array->[7]+0,abs($array->[7]-$array->[6]),$array->[8]+0]];
				$found = 1;
				last;
			}
		}
		if ($found == 0) {
			push(@{$ftr->{ontology_terms}->{InterPro}->{$array->[11]}->{evidence}},{
				method => $annofunc,
				method_version => $VERSION,
				timestamp => $timestamp,
				alignment_evidence => [[$array->[6]+0,$array->[7]+0,abs($array->[7]-$array->[6]),$array->[8]+0]]
			});
		}
    }
    close($fh);
    #Step 5: Saving the genome and report
    my $gaoutput = $self->util_ga_client()->save_one_genome_v1({
		workspace => $params->{workspace},
        name => $params->{genome_output_id},
        data => $genome,
        provenance => $self->util_provenance(),
        hidden => 0
	});
    my $reportObj = {
		'objects_created' => [],
		'text_message' => $numftr." annotated with ".$numdomains." distinct interpro domains by interpro scan!"
	};
    my $info = $wsClient->save_objects({
    	workspace => $params->{workspace},
    	objects => [{
    		type => "KBaseReport.Report",
    		data => $reportObj,
    		name => $params->{genome_output_id}.".annotate_genome_with_interpro_pipeline.report",
    		hidden => 1,
    		provenance => $self->util_provenance(),
    		meta => {}
    	}]
    });
   	return {
		report_name => $params->{genome_output_id}.".annotate_genome_with_interpro_pipeline.report",
		report_ref => $params->{workspace}.'/'.$params->{genome_output_id}.".annotate_genome_with_interpro_pipeline.report"
	};
}
#END_HEADER

sub new
{
    my($class, @args) = @_;
    my $self = {
    };
    bless $self, $class;
    #BEGIN_CONSTRUCTOR
    my $config_file = $ENV{ KB_DEPLOYMENT_CONFIG };
    my $cfg = Config::IniFiles->new(-file=>$config_file);
    $self->{'kbase-endpoint'} = $cfg->val('InterProScan_SDK','kbase-endpoint');
    $self->{'workspace-url'} = $cfg->val('InterProScan_SDK','workspace-url');
    $self->{'job-service-url'} = $cfg->val('InterProScan_SDK','job-service-url');
    $self->{'shock-url'} = $cfg->val('InterProScan_SDK','shock-url');
    $self->{'handle-service-url'} = $cfg->val('InterProScan_SDK','handle-service-url');
    $self->{'scratch'} = $cfg->val('InterProScan_SDK','scratch');
    if (!defined($self->{'workspace-url'})) {
		die "no workspace-url defined";
	}
    #END_CONSTRUCTOR

    if ($self->can('_init_instance'))
    {
	$self->_init_instance();
    }
    return $self;
}

=head1 METHODS



=head2 annotate_genome_with_interpro_pipeline

  $output = $obj->annotate_genome_with_interpro_pipeline($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is an InterProScan_SDK.AnnotateGenomeWithInterproPipelineParams
$output is an InterProScan_SDK.StandardFunctionOutput
AnnotateGenomeWithInterproPipelineParams is a reference to a hash where the following keys are defined:
	workspace has a value which is an InterProScan_SDK.workspace_name
	genome_workspace has a value which is an InterProScan_SDK.workspace_name
	genome_id has a value which is an InterProScan_SDK.genome_id
	genome_output_id has a value which is an InterProScan_SDK.genome_id
workspace_name is a string
genome_id is a string
StandardFunctionOutput is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is an InterProScan_SDK.Report_ref
Report_ref is a string

</pre>

=end html

=begin text

$params is an InterProScan_SDK.AnnotateGenomeWithInterproPipelineParams
$output is an InterProScan_SDK.StandardFunctionOutput
AnnotateGenomeWithInterproPipelineParams is a reference to a hash where the following keys are defined:
	workspace has a value which is an InterProScan_SDK.workspace_name
	genome_workspace has a value which is an InterProScan_SDK.workspace_name
	genome_id has a value which is an InterProScan_SDK.genome_id
	genome_output_id has a value which is an InterProScan_SDK.genome_id
workspace_name is a string
genome_id is a string
StandardFunctionOutput is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is an InterProScan_SDK.Report_ref
Report_ref is a string


=end text



=item Description



=back

=cut

sub annotate_genome_with_interpro_pipeline
{
    my $self = shift;
    my($params) = @_;

    my @_bad_arguments;
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_genome_with_interpro_pipeline:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotate_genome_with_interpro_pipeline');
    }

    my $ctx = $InterProScan_SDK::InterProScan_SDKServer::CallContext;
    my($output);
    #BEGIN annotate_genome_with_interpro_pipeline
    $self->util_initialize_call($params,$ctx);
	$output = $self->func_annotate_genome_with_interpro_pipeline($params);
    #END annotate_genome_with_interpro_pipeline
    my @_bad_returns;
    (ref($output) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"output\" (value was \"$output\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_genome_with_interpro_pipeline:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotate_genome_with_interpro_pipeline');
    }
    return($output);
}




=head2 version 

  $return = $obj->version()

=over 4

=item Parameter and return types

=begin html

<pre>
$return is a string
</pre>

=end html

=begin text

$return is a string

=end text

=item Description

Return the module version. This is a Semantic Versioning number.

=back

=cut

sub version {
    return $VERSION;
}

=head1 TYPES



=head2 bool

=over 4



=item Description

A binary boolean


=item Definition

=begin html

<pre>
an int
</pre>

=end html

=begin text

an int

=end text

=back



=head2 genome_id

=over 4



=item Description

A string representing a Genome id.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 workspace_name

=over 4



=item Description

A string representing a workspace name.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 Report_ref

=over 4



=item Description

The workspace ID for a Report object
@id ws KBaseReport.Report


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 StandardFunctionOutput

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is an InterProScan_SDK.Report_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is an InterProScan_SDK.Report_ref


=end text

=back



=head2 AnnotateGenomeWithInterproPipelineParams

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace has a value which is an InterProScan_SDK.workspace_name
genome_workspace has a value which is an InterProScan_SDK.workspace_name
genome_id has a value which is an InterProScan_SDK.genome_id
genome_output_id has a value which is an InterProScan_SDK.genome_id

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace has a value which is an InterProScan_SDK.workspace_name
genome_workspace has a value which is an InterProScan_SDK.workspace_name
genome_id has a value which is an InterProScan_SDK.genome_id
genome_output_id has a value which is an InterProScan_SDK.genome_id


=end text

=back



=cut

1;
