#!/usr/bin/sh

# Get the OOPMH-LIB root directory from a makefile
# OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
EXPECTED_NUM_TESTS=1

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation


# Validation for bubble unsteady
#-----------------------------------------
mkdir bubble_unsteady
echo "Running bubble unsteady validation "
../bubble_unsteady -n 4 -r 0.5 -c 0.01 -q 1 -o bubble_unsteady/ > \
    OUTPUT_bubble_unsteady
echo "done"
echo " " >> validation.log
echo "Bubble unsteady validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat bubble_unsteady/* > bubble_unsteady.dat

diff=$(zcmp ../validata/bubble_unsteady.dat.gz bubble_unsteady.dat )
if [ $? != 0 ]; then
    echo "[ERROR] Compare failed to run." >> validation.log
elif [ $diff ]; then
    echo "[FAILED]" >> validation.log
else
    echo "[OK]" >> validation.log
fi
#-----------------------------------------

cd ..

#######################################################################

#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
ACTUAL_NUM_TESTS_PASSED=$(grep OK Validation/validation.log | wc -l)
ACTUAL_NUM_TESTS_FAILED=$(grep FAILED Validation/validation.log | wc -l)
ACTUAL_NUM_TESTS_ERROR=$(grep ERROR Validation/validation.log | wc -l)
if [ $ACTUAL_NUM_TESTS_ERROR > 0 ]; then
    echo "Error in testing scripts. Check validation.log."
    return 1
fi
if [ $ACTUAL_NUM_TESTS_PASSED > $EXPECTED_NUM_TESTS ]; then
if [ $ACTUAL_NUM_TESTS_FAILED > 0 ]; then
    echo "Failed test. Check validation.log."
    return 1
fi
if [ $ACTUAL_NUM_TESTS_PASSED > $EXPECTED_NUM_TESTS ]; then
    echo "Passed more tests than expected! Check EXPECTED_NUM_TESTS and validation.log."
    return 2
fi
if [ $ACTUAL_NUM_TESTS_PASSED = $EXPECTED_NUM_TESTS ]; then
    echo "Passed all $(ACTUAL_NUM_TESTS_PASSED) tests."
    return 0
fi
if [ $ACTUAL_NUM_TESTS_PASSED < $EXPECTED_NUM_TESTS ]; then
    echo "Passed less tests than expected! Check EXPECTED_NUM_TESTS and validation.log."
    return 1
fi

# Never get here
exit 10
