#!/usr/bin/perl
#
# Configuration script for RIP4 code
# 
# Be sure to run as ./configure (to avoid getting a system configure command by mistake)
#

$sw_perl_path = perl ;
$sw_netcdf_path = "" ;
$sw_netcdff_lib = "";
$sw_ncarg_lib = "";
$sw_ldflags=""; 
$sw_compileflags=""; 
$sw_os = "ARCH" ;           # ARCH will match any
$sw_mach = "ARCH" ;         # ARCH will match any
$sw_ncarg_version="6"; 
$sw_compiler="unknown"; 

while ( substr( $ARGV[0], 0, 1 ) eq "-" )
 {
  if ( substr( $ARGV[0], 1, 5 ) eq "perl=" )
  {
    $sw_perl_path = substr( $ARGV[0], 6 ) ;
  }
  if ( substr( $ARGV[0], 1, 7 ) eq "netcdf=" )
  {
    $sw_netcdf_path = substr( $ARGV[0], 8 ) ;
  }
  if ( substr( $ARGV[0], 1, 8 ) eq "netcdff=" )
  {
    $sw_netcdff_lib = substr( $ARGV[0], 9 ) ;
  }
  if ( substr( $ARGV[0], 1, 3 ) eq "os=" )
  {
    $sw_os = substr( $ARGV[0], 4 ) ;
  }
  if ( substr( $ARGV[0], 1, 5 ) eq "mach=" )
  {
    $sw_mach = substr( $ARGV[0], 6 ) ;
  }
  if ( substr( $ARGV[0], 1, 5 ) eq "vrsn=" )
  {
    $sw_ncarg_version = substr( $ARGV[0], 6 ) ;
  }
  if ( substr( $ARGV[0], 1, 6 ) eq "cmplr=" )
  {
    $sw_compiler = substr( $ARGV[0], 7 ) ;
  }
  if ( substr( $ARGV[0], 1, 8 ) eq "ldflags=" )
  {
    $sw_ldflags = substr( $ARGV[0], 9 ) ;
# multiple options separated by spaces are passed in from sh script
# separated by ! instead. Replace with spaces here.
    $sw_ldflags =~ s/!/ /g ;
  }
  if ( substr( $ARGV[0], 1, 13 ) eq "compileflags=" )
  {
    $sw_compileflags = substr( $ARGV[0], 14 ) ;
    $sw_compileflags =~ s/!/ /g ;
  }
  shift @ARGV ;
 }

# The required ncarg libraries vary according to version number
printf " sw_ncarg_version = %s \n",$sw_ncarg_version;
if ( $sw_ncarg_version == 6 ) {
  $sw_ncarg_lib = "-lcairo -lfontconfig -lpixman-1 -lfreetype -lexpat -lpthread -lbz2 -lXrender";
}
elsif ( $sw_ncarg_version == 5 ) {
  $sw_ncarg_lib = "-lXpm";
}
# If we can guess the ncarg compiler, add the libraries to the list
  if ( $sw_compiler == "gfortran" ) {
    $sw_ncarg_lib = $sw_ncarg_lib . " -lgfortran -lgcc" ;
  }
  elsif ( $sw_compiler == "pgf90" ) {
    $sw_ncarg_lib = $sw_ncarg_lib . " -lpgftnrtl -lpgc" ;
  }
printf " library = %s \n",$sw_ncarg_lib;

# parse the configure.rip file

$validresponse = 0 ;

# Display the choices to the user and get selection
until ( $validresponse ) {
  printf "------------------------------------------------------------------------\n" ;
  printf "Please select from among the following supported platforms.\n\n" ;

  $opt = 1 ;
  open CONFIGURE_DEFAULTS, "< ./arch/configure.defaults" 
      or die "Cannot open ./arch/configure.defaults for reading" ;
  while ( <CONFIGURE_DEFAULTS> )
  {
    if ( substr( $_, 0, 5 ) eq "#ARCH" && ( index( $_, $sw_os ) >= 0 ) && ( index( $_, $sw_mach ) >= 0 ) )
    {
      $optstr[$opt] = substr($_,6) ;
      $optstr[$opt] =~ s/^[ 	]*// ;
      if ( substr( $optstr[$opt], 0,4 ) ne "NULL" )
      {
        printf "  %2d.  %s",$opt,$optstr[$opt] ;
        $opt++ ;
      }
    }
  }
  close CONFIGURE_DEFAULTS ;

  $opt -- ;

  printf "\nEnter selection [%d-%d] : ",1,$opt ;
  $response = <STDIN> ;

  if ( $response == -1 ) { exit ; }

  if ( $response >= 1 && $response <= $opt ) 
  { $validresponse = 1 ; }
  else
  { printf("\nInvalid response (%d)\n",$response);}
}
printf "------------------------------------------------------------------------\n" ;

$optchoice = $response ;

open CONFIGURE_DEFAULTS, "< ./arch/configure.defaults" 
      or die "Cannot open ./arch/configure.defaults for reading" ;
$latchon = 0 ;
while ( <CONFIGURE_DEFAULTS> )
{
  if ( substr( $_, 0, 5 ) eq "#ARCH" && $latchon == 1 )
  {
    $latchon = 0 ;
  }
  if ( $latchon == 1 )
  {
    $_ =~ s/CONFIGURE_PERL_PATH/$sw_perl_path/g ;
    $_ =~ s/CONFIGURE_NETCDF_PATH/$sw_netcdf_path/g ;
    $_ =~ s/CONFIGURE_LDFLAGS/$sw_ldflags/g ;
    $_ =~ s/CONFIGURE_COMPILEFLAGS/$sw_compileflags/g ;
    if ( $sw_netcdf_path ) 
      { $_ =~ s/CONFIGURE_WRFIO_NF/wrfio_nf/g ;
	$_ =~ s:CONFIGURE_NETCDF_FLAG:-DNETCDF: ;
	$_ =~ s:CONFIGURE_NETCDF_LIB_PATH:-L../external/io_netcdf -lwrfio_nf -L$sw_netcdf_path/lib -lnetcdf: ;
	 }
    else                   
      { $_ =~ s/CONFIGURE_WRFIO_NF//g ;
	$_ =~ s:CONFIGURE_NETCDF_FLAG::g ;
	$_ =~ s:CONFIGURE_NETCDF_LIB_PATH::g ;
	 }


    @machopts = ( @machopts, $_ ) ;
  }
  if ( substr( $_, 0, 5 ) eq "#ARCH" && $latchon == 0 )
  {
    $x=substr($_,6) ;
    $x=~s/^[     ]*// ;
    if ( $x eq $optstr[$optchoice] )
    {
      $latchon = 1 ;
    }
  }
}
close CONFIGURE_DEFAULTS ;

#printf "------------------------------------------------------------------------\n" ;
#foreach $f ( @machopts )
#{
#  if ( substr( $f , 0 , 8 ) eq "external" ) { last ; }
#  print $f ;
#}
#printf "------------------------------------------------------------------------\n" ;
#printf "\nYou have chosen: %s",$optstr[$optchoice] ;
#printf "Listed above are the default options for this platform.\n" ;
#printf "Settings are written to the file configure.rip here in the top-level\n" ;
#printf "directory.  If you wish to change settings, please edit that file.\n" ;
#printf "If you wish to change the default options, edit the file:\n\n" ;
#printf "     arch/configure.defaults\n" ;
#printf "\n" ;

open CONFIGURE_WRF, "> configure.rip" or die "cannot append configure.rip" ;
open ARCH_PREAMBLE, "< arch/preamble" or die "cannot open arch/preamble" ;
my @preamble;
# apply substitutions to the preamble...
while ( <ARCH_PREAMBLE> )
  {
  $_ =~ s:CONFIGURE_NETCDFF_LIB:$sw_netcdff_lib:g;
  $_ =~ s:CONFIGURE_NCARG_LIB:$sw_ncarg_lib:g;
  @preamble = ( @preamble, $_ ) ;
  }
close ARCH_PREAMBLE ;
print CONFIGURE_WRF @preamble  ;
close ARCH_PREAMBLE ;
printf CONFIGURE_WRF "# Settings for %s", $optstr[$optchoice] ;
print CONFIGURE_WRF @machopts  ;
open ARCH_POSTAMBLE, "< arch/postamble" or die "cannot open arch/postamble" ;
while ( <ARCH_POSTAMBLE> ) { print CONFIGURE_WRF } ;
close ARCH_POSTAMBLE ;
close CONFIGURE_WRF ;

printf "Configuration successful. To build RIP4, type: compile \n" ;
printf "------------------------------------------------------------------------\n" ;


