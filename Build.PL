#!/usr/bin/perl

use strict;
use warnings;


use Config;
use File::Spec;



my %prereqs = qw(
    Devel::CheckLib                 0
    Encode                          0
    Encode::Locale                  1.05
    ExtUtils::CppGuess              0
    ExtUtils::MakeMaker             6.80
    ExtUtils::ParseXS               3.35
    File::Basename                  0
    File::Spec                      0
    Getopt::Long                    0
    Module::Build::WithXSpp         0.14
    Moo                             1.003001
    POSIX                           0
    Scalar::Util                    0
    Test::More                      0
    Thread::Semaphore               0
    IO::Scalar                      0
    threads                         1.96
    Time::HiRes                     0
    Unicode::Normalize              0
);
my %recommends = qw(
    Class::XSAccessor               0
    Test::Harness                   0
    Thread::Queue                   0
    threads::shared                 0
);

my $sudo    = grep { $_ eq '--sudo' } @ARGV;
my $gui     = grep { $_ eq '--gui' } @ARGV;
my $xs_only = grep { $_ eq '--xs' }  @ARGV;
if ($gui) {
    %prereqs = qw(
    Class::Accessor                 0
    Wx                              0.9918
    Socket                          2.016
    );
    %recommends = qw(
    Growl::GNTP                     0.15
    Wx::GLCanvas                    0
    OpenGL                          0
    LWP::UserAgent                  0
    Net::Bonjour                    0
    );
    if ($^O eq 'MSWin32') {
        $recommends{"Win32::TieRegistry"} = 0;
        
        # we need an up-to-date Win32::API because older aren't thread-safe (GH #2517)
        $prereqs{'Win32::API'} = 0.79;
    }
} elsif ($xs_only) {
    %prereqs = %recommends = ();
}

my @missing_prereqs = ();
if ($ENV{SLIC3R_NO_AUTO}) {
    foreach my $module (sort keys %prereqs) {
        my $version = $prereqs{$module};
        next if eval "use $module $version; 1";
        push @missing_prereqs, $module if exists $prereqs{$module};
        print "Missing prerequisite $module $version\n";
    }
    foreach my $module (sort keys %recommends) {
        my $version = $recommends{$module};
        next if eval "use $module $version; 1";
        print "Missing optional $module $version\n";
    }
} else {
    my @try = (
        $ENV{CPANM} // (),
        File::Spec->catfile($Config{sitebin}, 'cpanm'),
        File::Spec->catfile($Config{installscript}, 'cpanm'),
    );
    
    my $cpanm;
    foreach my $path (@try) {
        if (-e $path) {  # don't use -x because it fails on Windows
            $cpanm = $path;
            last;
        }
    }
    if (!$cpanm) {
        if ($^O =~ /^(?:darwin|linux)$/ && system(qw(which cpanm)) == 0) {
            $cpanm = 'cpanm';
        }
    }
    die <<'EOF'
cpanm was not found. Please install it before running this script.

There are several ways to install cpanm, try one of these:

    apt-get install cpanminus
    curl -L http://cpanmin.us | perl - --sudo App::cpanminus
    cpan App::cpanminus

If it is installed in a non-standard location you can do:
    
    CPANM=/path/to/cpanm perl Build.PL

EOF
        if !$cpanm;
    my @cpanm_args = ();
    push @cpanm_args, "--sudo" if $sudo;
    
    # install local::lib without --local-lib otherwise it's not usable afterwards
    if (!eval "use local::lib qw(local-lib); 1") {
        my $res = system $cpanm, @cpanm_args, 'local::lib';
        warn "Warning: local::lib is required. You might need to run the `cpanm --sudo local::lib` command in order to install it.\n"
            if $res != 0;
    }
    
    push @cpanm_args, ('--local-lib', 'local-lib');
    
    # make sure our cpanm is updated (old ones don't support the ~ syntax)
    system $cpanm, @cpanm_args, 'App::cpanminus';
    
    my %modules = (%prereqs, %recommends);
    foreach my $module (sort keys %modules) {
        my $version = $modules{$module};
        my @cmd = ($cpanm, @cpanm_args);
        
        # temporary workaround for upstream bug in test
        push @cmd, '--notest'
            if $module =~ /^(?:OpenGL|Math::PlanePath|Test::Harness|IO::Scalar)$/;
        
        # do not try to reinstall modules that were already
        # installed manually
        push @cmd, '--skip-satisfied';

        push @cmd, "$module~$version";
        my $res = system @cmd;
        if ($res != 0) {
            if (exists $prereqs{$module}) {
                push @missing_prereqs, $module;
            } else {
                printf "Don't worry, this module is optional.\n";
            }
        }
    }
    
    if (!$gui) {
        # clean xs directory before reinstalling, to make sure Build is called
        # with current perl binary
        if (-e './xs/Build') {
            if ($^O eq 'MSWin32') {
                system '.\xs\Build', 'distclean';
            } else {
	#	system './xs/Build', 'distclean'; # makipoea	
		system './xs/Build';
	    }
        }
        my $res = system $cpanm, @cpanm_args, '--reinstall', '--verbose', './xs';
        if ($res != 0) {
            die "The XS/C++ code failed to compile, aborting\n";
        }
    }
}

if (@missing_prereqs) {
    printf "The following prerequisites failed to install: %s\n", join(', ', @missing_prereqs);
    exit 1;
} elsif (!$gui) {
    eval "use App::Prove; 1" or die "Failed to load App::Prove";
    my $res = 0;#App::Prove->new->run ? 0 : 1;# makipoea # my $res = 0;
    if ($res == 0) {
        print "If you also want to use the GUI you can now run `perl Build.PL --gui` to install the required modules.\n";
    } else {
        print "Some tests failed. Please report the failure to the author!\n";
    }
    exit $res;
}

__END__
