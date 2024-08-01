

import ipyrad.analysis as ipa
import pandas as pd

converter = ipa.vcf_to_hdf5(
    name="Ecoli_LD5K",
    data="vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_subsample.vcf",
    ld_block_size=5000,
)

# run the converter
converter.run()