# fagales_traits

## Usage
`mkdir output
cat discard_terms.csv unit_columns.csv > all_fields_dropped.csv
python3 ./trait_processing/traiter_process_terms.py ./trait_processing/Fagales_2023-01-26/Fagales_2023-01-26.csv ./trait_processing/Manual_trait_extraction_newheader/Fagales_fill-in_controlled_fields.csv ./trait_processing/output/out ./trait_processing/output/codeguide.csv ./trait_processing/output/distancematrix ./trait_processing/all_fields_dropped.csv ./trait_processing/concatenate_terms.csv ./trait_processing/range_terms_quantitative.csv ./trait_processing/trait_processing/range_terms_count.csv 0.95`

Note it is assumed that the base repository directory is the working directory. If you don't want to use the second manual dataset, you can just define an empty data frame. Note the first two CSVs are trait data and the rest are controlled fields and terms, which are documented in these Excel files: `./trait_processing/Trait_field_reconciliation.xlsx` and `./trait_processing/Trait_term_reconciliation.xlsx`. 


Description of folders as follows:

## trait_processing folder
Contains the main analysis script and field/term vocabularies, as well as the following subfolders:

### Fagales_2023-01-26
Automated extraction results (HTML is interactive output; csv is flattened output).

### Manual_trait_extraction
Manual extraction results ("controlled fields" files have been reformatted to follow automated extraction").

### Manual_trait_extraction_newheader
Same as above but matched headers to automated extraction.

### output
Directory containing script output.

## trait_ordination
Contains ordinations demonstrated in the manuscript. Code is at `pcoa_traitdata.r`, the CSVs are local copies of the output directory described above, and PDFs show output graphs.

## manual_scoring
Contains scored errors in `out_droppedmissing_Betulaceae` and result graphs in PDFs.



