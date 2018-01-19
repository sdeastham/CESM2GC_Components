#!/bin/sh -f

echo 

if [ $# -ne 1 ]; then
    echo "Invoke archive_baseline.sh -help for usage." 
    exit 1
fi

if [ $1 == "-help" ]; then
cat << EOF1
NAME

	archive_baseline.sh - archive pretag baselines to set locations on
                              hobart and yellowstone.


SYNOPSIS

	archive_baseline.sh TAGNAME
	  [-help]


ENVIROMENT VARIABLES

	CESM_TESTDIR - Directory that contains the CESM finished results you wish to archive.
	CAM_TESTDIR - Directory that contains the CAM finished results you wish to archive.
	CAM_FC      - Compiler used, only used on hobart (PGI,NAG), where the compiler
                      name is appended to the archive directory.


BASELINE ARCHIVED LOCATION

	hobart:     /fs/cgd/csm/models/atm/cam/pretag_bl/TAGNAME_pgi
	              /fs/cgd/csm/models/atm/cam/pretag_bl/TAGNAME_nag
        yellowstone:  /glade/p/cesm/cseg/models/atm/cam/pretag_bl/TAGNAME



HOW TO USE ARCHIVE BASELINES

	Set BL_TESTDIR to the archived baseline you wish to load.


WORK FLOW

	This is an example for hobart.

	Modify your sandbox with the changes you want.
        setenv CAM_FC PGI
        setenv CAM_TESTDIR /scratch/cluster/fischer/cam5_2_06
        Run the cam test suite.
        Make your trunk tag
	archive_baseline.sh cam5_2_06
        
	Create a new sandbox.
        setenv CAM_FC PGI
	setenv CAM_TESTDIR /scratch/cluster/fischer/cam5_2_07
        setenv BL_TESTDIR /fs/cgd/csm/models/atm/cam/pretag_bl/cam5_2_06_pgi
        Run the cam test suite.
        Make your trunk tag
        archive_baseline.sh cam5_2_07	
   

WARNING
	
	System changes can cause answer changes. So you may need to create new baselines
        if you are getting unexpected baseline failures.

EOF1
exit
fi

if [ -z "$CAM_TESTDIR" ]; then
  echo "ERROR: please set CAM_TESTDIR"
  echo
  exit 1
fi


hostname=`hostname`
case $hostname in

  ho*)
    echo "server: hobart"
    if [ -z "$CAM_FC" ]; then
      CAM_FC="PGI"
    fi
    test_file_list="tests_pretag_hobart_${CAM_FC,,}"
    baselinedir="/fs/cgd/csm/models/atm/cam/pretag_bl/$1_${CAM_FC,,}"
  ;;
  ys*)
    echo "server: yellowstone"
    if [ -z "$CAM_FC" ]; then
      CAM_FC="INTEL"
    fi
    if [ -z "$CESM_TESTDIR" ]; then
      echo "ERROR: please set CESM_TESTDIR (test-root in the create_test)"
      echo
      exit 1
    fi
    test_file_list="tests_pretag_yellowstone"
    baselinedir="/glade/p/cesm/cseg/models/atm/cam/pretag_bl/$1"
  ;;
  * ) echo "ERROR: machine $hostname not currently supported"; exit 1 ;;
esac

echo " "
echo "CESM Archiving to /glade/p/cesmdata/cseg/ccsm_baselines/$1"
echo " "
case $hostname in
  ys*)
    ../../../../cime/scripts/Tools/bless_test_results  -t '' -c '' --hist-only  -r $CESM_TESTDIR -b $1 -f -s
esac

echo
echo "Archiving to ${baselinedir}"
echo
if [ -d ${baselinedir} ]; then
   echo "ERROR: Baseline $baselinedir already exists."
   exit 1
fi

mkdir $baselinedir

if [ ! -d ${baselinedir} ]; then
   echo "ERROR: Failed to make ${baselinedir}"
   exit 1
fi

echo "Archiving the following directories."
test_list=""
while read input_line; do
    test_list="${input_line} "
  for test_id in ${test_list}; do
      master_line=`grep $test_id input_tests_master`
       str1=${master_line%% *}
       temp=${master_line#$str1 }
       str2=${temp%% *}

       temp=${temp#$str2 }
       str3=${temp%% *}
       temp=${temp#$str3 }
       str4=${temp%% *}
       temp=${temp#$str4 }
       str5=${temp%% *}

       temp=${str2%%.*}
       scr1=${temp#"TBL"}
       scr1=${temp#"TBL"}



       if grep -c TBL ${str2} > /dev/null; then
         case="TSM${scr1}.$str3.$str4.$str5"
         ls -ld ${CAM_TESTDIR}/${case}
         cp -rp ${CAM_TESTDIR}/${case} ${baselinedir}/${case}
         chmod -R a+r ${baselinedir}
         chmod -R g+w ${baselinedir}
       fi

  done

done < ${test_file_list}