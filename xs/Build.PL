#!/usr/bin/perl -w

use strict;
use warnings;

use Config;
use Devel::CheckLib;
use ExtUtils::CppGuess;
use Module::Build::WithXSpp;

my $cpp_guess = ExtUtils::CppGuess->new;
my $mswin = $^O eq 'MSWin32';
my $linux = $^O eq 'linux';



# prevent an annoying concatenation warning by Devel::CheckLib
$ENV{LD_RUN_PATH} //= "";

# HAS_BOOL         : stops Perl/lib/CORE/handy.h from doing "#  define bool char" for MSVC
# NOGDI            : prevents inclusion of wingdi.h which defines functions Polygon() and Polyline() in global namespace
# BOOST_ASIO_DISABLE_KQUEUE : prevents a Boost ASIO bug on OS X: https://svn.boost.org/trac/boost/ticket/5339
my @cflags = qw(-DHAS_BOOL -DNOGDI -DSLIC3RXS -DBOOST_ASIO_DISABLE_KQUEUE -Dexprtk_disable_rtl_io_file  -Dexprtk_disable_rtl_vecops -Dexprtk_disable_string_capabilities -Dexprtk_disable_enhanced_features); # makipoea j'ai retiré le flag : -Dexprtk_disable_return_statement
push @cflags, "-DSLIC3R_BUILD_COMMIT=$ENV{SLIC3R_GIT_VERSION}" if defined $ENV{SLIC3R_GIT_VERSION};

push @cflags, "-I/usr/include" if -d "/usr/include"; #makipoea


# std=c++11 Enforce usage of C++11 (required now). Minimum compiler supported: gcc 4.9, clang 3.3, MSVC 14.0
if ($cpp_guess->is_gcc)	{
	# GCC is pedantic with c++11 std, so use -std=gnu++11 to be able to use M_PI
	push @cflags, qw(-std=gnu++11);
} else {
	push @cflags, qw(-std=c++11);
}

my @ldflags = ();

#makipoea
my @LIBS = $cpp_guess->is_msvc
    ? qw(-LIBPATH:src/libslic3r)
    : qw(-Lsrc/libslic3r);
#makipoea
push @LIBS, '-lglpk';# makipoea

if ($linux && (defined $ENV{SLIC3R_STATIC} && $ENV{SLIC3R_STATIC})) {
    push @ldflags, qw(-static-libgcc -static-libstdc++);
    if ($ENV{TRAVIS}) {
        # On the build server, link to the actual static libraries to make sure we get them in the list.
    	push @ldflags, qw(/usr/lib/gcc/x86_64-linux-gnu/4.9/libstdc++.a /usr/lib/gcc/x86_64-linux-gnu/4.9/libgcc.a);
    }
    # ExtUtils::CppGuess has a hard-coded -lstdc++, so we filter it out
    {
        no strict 'refs';
        no warnings 'redefine';
        my $func = "ExtUtils::CppGuess::_get_lflags";
        my $orig = *$func{CODE};
        *{$func} = sub {
            my $lflags = $orig->(@_);
            $lflags =~ s/\s*-lstdc\+\+//;
            return $lflags;
        };
    }
}

if ($^O eq 'darwin') {
    push @cflags, qw(-stdlib=libc++);
    push @ldflags, qw(-framework IOKit -framework CoreFoundation -lc++);
    
    # Due to a bug/misconfiguration/stupidity, boost 1.52 and libc++ don't like each
    # other much: a compilation error "Constexpr function never produces a constant
    # expression" pops up when trying to compile anything that uses
    # boost/chrono/duration.hpp (namely boost/thread for us). This is a workaround
    # that prevents this from happening, not needed with newer Boost versions.
    # See here for more details: https://svn.boost.org/trac/boost/ticket/7671
    push @cflags, qw(-DBOOST_THREAD_DONT_USE_CHRONO -DBOOST_NO_CXX11_RVALUE_REFERENCES -DBOOST_THREAD_USES_MOVE);
    
    # ExtUtils::CppGuess has a hard-coded -lstdc++, so we filter it out
    {
        no strict 'refs';
        no warnings 'redefine';
        my $func = "ExtUtils::CppGuess::_get_lflags";
        my $orig = *$func{CODE};
        *{$func} = sub {
            my $lflags = $orig->(@_);
            $lflags =~ s/\s*-lstdc\+\+//;
            return $lflags;
        };
    }
}
if ($mswin) {
    push @cflags, qw(-DNOMINMAX);
}

my @INC  = qw(-Isrc/libslic3r);
#@LIBS = $cpp_guess->is_msvc ? qw(-LIBPATH:src/libslic3r) : qw(-Lsrc/libslic3r); #makipoea

# makipoea la ligne du haut commencait pas in my @LIBS.....

# search for Boost in a number of places
my @boost_include = ();
if (defined $ENV{BOOST_INCLUDEDIR}) {
    push @boost_include, $ENV{BOOST_INCLUDEDIR}
} elsif (defined $ENV{BOOST_DIR}) {
    # User might have provided a path relative to the Build.PL in the main directory
    foreach my $subdir ($ENV{BOOST_DIR}, "../$ENV{BOOST_DIR}", "$ENV{BOOST_DIR}/include", "../$ENV{BOOST_DIR}/include") {
        next if $subdir =~ m{^\.\.//};
        printf "Checking %s: ", $subdir;
        if (-d $subdir) {
            push @boost_include, $subdir;
            print "found\n";
        } else {
            print "not found\n";
        }
    }
    die "Invalid BOOST_DIR: could not find Boost headers\n" if !@boost_include;
} else {
    # Boost library was not defined by the environment.
    # Try to guess at some default paths.
    if ($mswin) {
        for my $path (glob('C:\dev\boost*\include'), glob ('C:\boost*\include')) {
            push @boost_include, $path;
        }
        if (! @boost_include) {
            # No boost\include. Try to include the boost root.
            for my $path (glob('C:\dev\boost*'), glob ('C:\boost*')) {
                push @boost_include, $path;
            }
        }
    } else {
        push @boost_include, grep { -d $_ }
            qw(/opt/local/include /usr/local/include /opt/include /usr/include);
    }
}

#makipoea
check_lib(
    lib     => 'glpk',
    header  => 'glpk.h',
) or die "GLPK library not found, please install it.\n";
#makipoea

my @boost_libs = ();
if (defined $ENV{BOOST_LIBRARYPATH}) {
    push @boost_libs, $ENV{BOOST_LIBRARYPATH}
} elsif (defined $ENV{BOOST_DIR}) {
    # User might have provided a path relative to the Build.PL in the main directory
    foreach my $subdir ("$ENV{BOOST_DIR}/stage/lib", "../$ENV{BOOST_DIR}/stage/lib", "$ENV{BOOST_DIR}/lib", "../$ENV{BOOST_DIR}/lib") {
        next if $subdir =~ m{^\.\.//};
        printf "Checking %s: ", $subdir;
        if (-d $subdir) {
            push @boost_libs, $subdir;
            print "found\n";
        } else {
            print "not found\n";
        }
    }
    die "Invalid BOOST_DIR: could not find Boost libs\n" if !@boost_libs;
} else {
    # Boost library was not defined by the environment.
    # Try to guess at some default paths.
    if ($mswin) {
        for my $path (
            glob('C:\dev\boost*\lib'),       glob ('C:\boost*\lib'), 
            glob('C:\dev\boost*\stage\lib'), glob ('C:\boost*\stage\lib')) {
            push @boost_libs, $path;
        }
    } else {
        push @boost_libs, grep { -d $_ }
            qw(/opt/local/lib /usr/local/lib /opt/lib /usr/lib /lib);
    }
}
# In order to generate the -l switches we need to know how Boost libraries are named
my $have_boost = 0;
my $have_boost_optional = 0;
my @boost_libraries = qw(system thread filesystem);  # we need these
my @boost_optional_libraries = qw(nowide);  # we need these, but if they aren't present we can deal
# check without explicit lib path (works on Linux)
if (! $mswin) {
    $have_boost = 1
        if check_lib(
            lib     => [ map "boost_${_}", @boost_libraries ],
        );
    $have_boost_optional = 1
        if check_lib(
            lib     => [ map "boost_${_}", @boost_optional_libraries ],
        );
}

if (!$ENV{SLIC3R_STATIC} && $have_boost) {
    # The boost library was detected by check_lib on Linux.
    push @LIBS, map "-lboost_${_}", @boost_libraries;
    if (!$ENV{SLIC3R_STATIC} && $have_boost_optional) {
        push @LIBS, map "-lboost_${_}", @boost_optional_libraries;
    }
} else {
    # Either static linking, or check_lib could not be used to find the boost libraries.
    my $lib_prefix = 'libboost_';
    my $lib_ext = $Config{lib_ext};
    PATH: foreach my $path (@boost_libs) {
        # Try to find the boost system library.
        my @files = glob "$path/${lib_prefix}system*$lib_ext";
        next if !@files;
  
        if ($files[0] =~ /${lib_prefix}system([^.]*)$lib_ext$/) {
            # Suffix contains the version number, the build type etc.
            my $suffix = $1;
            # Verify existence of all required boost libraries at $path.
            for my $lib (map "${lib_prefix}${_}${suffix}${lib_ext}", @boost_libraries) {
                # If the library file does not exist, try next library path.
                -f "$path/$lib" or next PATH;
            }
            if (! $cpp_guess->is_msvc) {
                # Test the correctness of boost libraries by linking them to a minimal C program.
                check_lib(
                    lib     => [ map "boost_${_}${suffix}", @boost_libraries ],
                    INC     => join(' ', map "-I$_", @INC,  @boost_include),
                    LIBS    => "-L$path",
                ) or next;
            }
            push @INC, (map " -I$_", @boost_include);  # TODO: only use the one related to the chosen lib path
            if ($ENV{SLIC3R_STATIC} || $cpp_guess->is_msvc) {
                push @LIBS, map "${path}/${lib_prefix}$_${suffix}${lib_ext}", @boost_libraries;
            } else {
                push @LIBS, " -L$path", (map " -lboost_$_$suffix", @boost_libraries);
            }
            $have_boost = 1;
            last;
        }
    }
}
die <<'EOF' if !$have_boost;
Slic3r requires the Boost libraries. Please make sure they are installed.

If they are installed, this script should be able to locate them in several
standard locations. If this is not the case, you might want to supply their 
path through the BOOST_INCLUDEDIR and BOOST_LIBRARYPATH environment variables:

    BOOST_INCLUDEDIR=/usr/local/include BOOST_LIBRARYPATH=/usr/lib perl Build.PL

If you just compiled Boost in its source directory without installing it in the
system you can just provide the BOOST_DIR variable pointing to that directory.

EOF

if ($ENV{SLIC3R_DEBUG}) {
    # only on newer GCCs: -ftemplate-backtrace-limit=0
    push @cflags, '-DSLIC3R_DEBUG';
    push @cflags, $cpp_guess->is_msvc ? '-Gd' : '-g';
} else {
    # Disable asserts in the release builds.
    push @cflags, '-DNDEBUG', '-O';
}
if ($cpp_guess->is_gcc) {
    # our templated XS bindings cause undefined-var-template warnings
    push @cflags, qw(-Wno-undefined-var-template);
}

my $build = Module::Build::WithXSpp->new(
    module_name     => 'Slic3r::XS',
    dist_abstract   => 'XS code for Slic3r',
    build_requires => {qw(
        ExtUtils::ParseXS           3.35
        ExtUtils::Typemaps          1.00
        ExtUtils::Typemaps::Default 1.05
        ExtUtils::XSpp              0.18
        Module::Build               0.3601
        Test::More                  0
    )},
    configure_requires => {qw(
        ExtUtils::CppGuess          0.07
        Module::Build               0.38
        Module::Build::WithXSpp     0.13
    )},
    extra_compiler_flags => [ @INC, @cflags ],
    extra_linker_flags => [ @LIBS, @ldflags ],
    
    # Provides extra C typemaps that are auto-merged
    extra_typemap_modules => {
        'ExtUtils::Typemaps::Basic' => '1.05',
    },
    
    # for MSVC builds
    early_includes => [qw(
        cstring
        cstdlib
        ostream
        sstream
        libslic3r/GCodeSender.hpp
    )]
);

$build->create_build_script;

__END__
