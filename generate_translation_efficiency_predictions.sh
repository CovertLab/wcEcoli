# Script to run fitting used to generate translation efficiency predictions.
# Requires setup outlined in requirements.txt and in README.md.

make clean compile

## Set B0 - RNA mass per cell
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-mass-fraction-rna \
--save-cell-specs B00

## Set B1 - mRNA mass fraction
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-mass-fraction-mrna \
--save-cell-specs B01

## Set B2 - RNA sequencing, Covert 2004
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-rna-seq=Covert \
--save-cell-specs B02

## Set B3 - RNA sequencing, Cho 2011
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-rna-seq=Cho \
--save-cell-specs B03

## Set B4 - RNA sequencing, Dong 2011
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-rna-seq=Dong \
--save-cell-specs B04

## Set B5 - mRNA half lives, with kas
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-rna-half-life=with_kas \
--save-cell-specs B05

## Set B6 - mRNA half lives, without kas
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-rna-half-life=without_kas \
--save-cell-specs B06

## Set B7 - r-protein half lives
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-r-protein-degradation \
--save-cell-specs B07

## Set B8 - protein mass per cell
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-mass-fraction-protein \
--save-cell-specs B08

## Set B9 - transcription elongation rate
# TODO: add after --flat-elongation is separated

## Set B10 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-rnap-activity=0.1 \
--save-cell-specs B10

## Set B11 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-rnap-activity=0.2 \
--save-cell-specs B11

## Set B12 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-rnap-activity=0.3 \
--save-cell-specs B12

## Set B13 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-rnap-activity=0.4 \
--save-cell-specs B13

## Set B14 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-rnap-activity=0.5 \
--save-cell-specs B14

## Set B15 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-rnap-activity=0.6 \
--save-cell-specs B15

## Set B16 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-rnap-activity=0.7 \
--save-cell-specs B16

## Set B17 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-rnap-activity=0.8 \
--save-cell-specs B17

## Set B18 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-rnap-activity=0.9 \
--save-cell-specs B18

## Set B19 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-rnap-activity=1.0 \
--save-cell-specs B19

## Set B20 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--disable-rnap-fraction-increase --alternate-rnap-activity=0.1 \
--save-cell-specs B20

## Set B21 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--disable-rnap-fraction-increase --alternate-rnap-activity=0.2 \
--save-cell-specs B21

## Set B22 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--disable-rnap-fraction-increase --alternate-rnap-activity=0.3 \
--save-cell-specs B22

## Set B23 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--disable-rnap-fraction-increase --alternate-rnap-activity=0.4 \
--save-cell-specs B23

## Set B24 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--disable-rnap-fraction-increase --alternate-rnap-activity=0.5 \
--save-cell-specs B24

## Set B25 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--disable-rnap-fraction-increase --alternate-rnap-activity=0.6 \
--save-cell-specs B25

## Set B26 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--disable-rnap-fraction-increase --alternate-rnap-activity=0.7 \
--save-cell-specs B26

## Set B27 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--disable-rnap-fraction-increase --alternate-rnap-activity=0.8 \
--save-cell-specs B27

## Set B28 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--disable-rnap-fraction-increase --alternate-rnap-activity=0.9 \
--save-cell-specs B28

## Set B29 - RNA polymerase activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--disable-rnap-fraction-increase --alternate-rnap-activity=1.0 \
--save-cell-specs B29

## Set B30 - translation elongation rate
# TODO: add after --flat-elongation is separated

## Set B31 - ribosome activity
python runscripts/manual/runFitter.py \
--debug \
--disable-ribosome-fitting --disable-rnapoly-fitting --flat-elongation \
--alternate-ribosome-activity \
--save-cell-specs B31
