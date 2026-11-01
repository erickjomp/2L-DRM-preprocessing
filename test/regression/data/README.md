# Regression test data — do not modify

Input files of the small cropped AMZ2 CPM 4 km case (16 x 14 cells over the
Ecuadorian Andes, 60 levels, hourly from 2018-05-31 00:00 to 2018-06-02 00:00),
used by the regression tests (`test/regression`). Only the variables read by
prepro2LDRM-slabvars are kept. `p_lowest_level.dat.gz` holds the pressure of the
lowest model level of the 48 files used, for the NaN checks.

Any change to these files changes the results, so the regression tests will fail
against the stored references in `test/regression/reference/`. To regenerate them,
use `test/regression/make_test_data.py` and update the references with
`--update-reference` in the same commit.
